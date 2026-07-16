#' Resolve Parallel Configuration
#'
#' Validates and normalizes the user-facing parallel knobs and returns a list
#' describing whether/when parallel dispatch should happen. The decision is
#' centralized here so callers cannot accidentally nest two parallel layers.
#'
#' @param use_parallel Logical master switch.
#' @param n_cores Integer number of workers, or `NULL` for auto.
#' @param parallel_strategy Character: `"grid"`, `"restart"`, or `"auto"`.
#' @param requested Character: which parallel layer is asking
#'   (`"grid"` or `"restart"`). Ignored unless `parallel_strategy != "auto"`.
#' @param seed Integer seed or `NULL`.
#' @param n_tasks Integer: number of independent tasks the active parallel
#'   layer will dispatch (e.g. `length(k_grid)` or `nstart`). Used only when
#'   `n_cores = NULL` to avoid spawning more workers than tasks.
#' @return A list with `active` (logical), `n_cores` (integer), `strategy`
#'   (resolved character), and `seed`.
#' @noRd
mv_parallel_config <- function(use_parallel, n_cores, parallel_strategy = "auto",
                               requested = NULL, seed = NULL, n_tasks = NULL) {
  stopifnot(
    is.logical(use_parallel), length(use_parallel) == 1, !is.na(use_parallel)
  )

  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
      stop("'n_cores' must be a positive integer (or NULL for automatic).")
    }
    n_cores <- as.integer(n_cores)
  }

  stopifnot(is.null(seed) || (is.numeric(seed) && length(seed) == 1))
  if (!is.null(seed)) seed <- as.integer(seed)

  if (!is.null(n_tasks)) {
    stopifnot(is.numeric(n_tasks), length(n_tasks) == 1, n_tasks >= 1)
    n_tasks <- as.integer(n_tasks)
  }

  strategy <- match.arg(parallel_strategy, c("auto", "grid", "restart"))

  # future may not be installed; if use_parallel is requested but unavailable
  # we signal a clear error rather than silently falling back.
  if (use_parallel && !requireNamespace("future", quietly = TRUE)) {
    stop(
      "'use_parallel = TRUE' requires the 'future' package. ",
      "Install it with install.packages(\"future\") or set use_parallel = FALSE."
    )
  }

  active <- FALSE
  if (use_parallel) {
    if (is.null(requested)) {
      active <- TRUE
    } else if (strategy == "auto") {
      active <- TRUE
    } else {
      active <- (strategy == requested)
    }
  }

  if (active && is.null(n_cores)) {
    n_cores <- mv_optimal_n_cores(n_tasks = n_tasks)
  }

  list(active = active, n_cores = n_cores %||% 1L, strategy = strategy, seed = seed)
}

#' Choose an Optimal Number of Parallel Workers
#'
#' Picks a sensible default worker count for a clustering job, balancing
#' available hardware against the actual amount of parallelizable work.
#'
#' The heuristic:
#' - Start from `future::availableCores()` (honors `mc.cores`,
#'   `OMP_NUM_THREADS`, cgroups, etc.).
#' - Never exceed `future::availableWorkers()`.
#' - Cap at `max_workers` so we don't spawn hundreds of PSOCK processes on
#'   large machines where oversubscription would hurt more than help.
#' - Scale to the parallelizable unit count (`n_tasks`): there is no point in
#'   more workers than tasks (it only adds startup/serialization overhead), so
#'   the result is `min(cores, workers, max_workers, n_tasks)` (at least 1).
#'
#' @param n_tasks Integer: number of independent tasks the parallel layer will
#'   dispatch (e.g. `length(k_grid)` or `nstart`). `NULL` means "unknown yet",
#'   in which case the cap is ignored.
#' @param max_workers Integer: hard upper bound on workers (default 8). Tune
#'   down on memory-constrained hosts; tune up for very large grids.
#' @return A positive integer number of workers to use.
#' @noRd
mv_optimal_n_cores <- function(n_tasks = NULL, max_workers = 8L) {
  if (!requireNamespace("future", quietly = TRUE)) {
    return(1L)
  }
  cores <- tryCatch(
    as.integer(future::availableCores()),
    error = function(e) 1L
  )
  workers <- tryCatch(
    # availableWorkers() returns a CHARACTER vector of hostnames; the worker
    # count is its LENGTH, not an integer coercion of the names.
    length(future::availableWorkers()),
    error = function(e) cores
  )
  n <- min(cores, workers, max_workers)
  if (!is.null(n_tasks)) {
    n <- min(n, max(1L, as.integer(n_tasks)))
  }
  max(1L, n)
}

#' Run a Future-Lapply Over a List
#'
#' The single parallel dispatch primitive used throughout Ampharos. It:
#' 1. Picks a Windows-safe plan (multisession => PSOCK workers).
#' 2. Seeds each task deterministically from `(seed, task_index)` inside the
#'    worker (see `mv_task_seed()`); `future.seed` is `FALSE` so future performs
#'    no additional RNG handling and results depend only on our seeding.
#' 3. Restores the previous plan afterwards (no leaked workers).
#'
#' The `FUN` argument must be self-contained: it should depend only on the
#' objects passed via `...` plus the Ampharos internals, which are made
#' available inside each worker through `mv_parallel_worker_setup()`. The small
#' `mv_*` helper closures used by the dispatcher are bundled automatically by
#' future's default global detection; callers should keep any large objects
#' passed via `...` minimal to limit serialization cost.
#'
#' @param X A vector or list to iterate over.
#' @param FUN Function applied to each element of `X`.
#' @param config List from `mv_parallel_config()`.
#' @param ... Extra arguments forwarded to `FUN` (must be serializable).
#' @return A list of length `length(X)` with the results.
#' @noRd
mv_future_lapply <- function(X, FUN, config, ...) {
  stopifnot(is.list(config), is.logical(config$active))

  # Shared per-task wrapper. Crucially, the *worker environment* is set up and
  # the RNG is seeded deterministically from (seed, task index) BEFORE FUN runs.
  # Because the seed depends only on (seed, idx) and not on execution order, the
  # sequential path (below) produces byte-identical results to the parallel path:
  # the same seed always yields the same fit, regardless of how tasks are
  # scheduled across workers. A single parallel layer is enforced elsewhere, so
  # workers themselves never spawn sub-workers.
  #
  # The caller's RNG state (kind + .Random.seed) is saved and restored so the
  # sequential path leaves no global RNG side-effects behind.
  run_task <- function(x, idx, ...) {
    # Worker setup MUST run first: under future.globals = FALSE the mv_*
    # helpers are only resolvable on the worker after mv_parallel_worker_setup
    # attaches the private environment. Restoring RNG state after setup is fine.
    mv_parallel_worker_setup()
    rng_save <- mv_rng_state_save()
    on.exit(mv_rng_state_restore(rng_save), add = TRUE)
    if (!is.null(config$seed)) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(mv_task_seed(config$seed, idx))
    }
    FUN(x, ...)
  }

  if (!config$active) {
    return(lapply(seq_along(X), function(i) run_task(X[[i]], i, ...)))
  }

  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop(
      "'use_parallel = TRUE' requires the 'future' and 'future.apply' packages. ",
      "Install them with install.packages(c(\"future\", \"future.apply\")) ",
      "or set use_parallel = FALSE."
    )
  }

  n_workers <- config$n_cores
  n_tasks <- length(X)

  # Avoid spawning more workers than tasks (pure overhead otherwise).
  n_workers <- min(n_workers, max(1L, n_tasks))

  # Temporarily switch to a Windows-safe multisession plan. Restore on exit.
  # NOTE: future::plan() is global; concurrent plan changes in the same
  # session are not supported and would race.
  prev_plan <- future::plan(
    future::multisession,
    workers = n_workers,
    .init = FALSE
  )
  on.exit({
    ok <- tryCatch(
      future::plan(prev_plan, .init = FALSE),
      error = function(e) FALSE
    )
    if (isFALSE(ok)) {
      warning("Failed to restore previous future::plan(); default plan may have changed.")
    }
  }, add = TRUE)

  # future.seed = FALSE: we rely entirely on our own per-task seeding
  # (mv_task_seed), so future does not touch the RNG. We keep future's default
  # global detection (future.globals) so the small mv_* helper closures used by
  # run_task are bundled automatically; only large objects passed via ... should
  # be kept minimal by callers to limit serialization cost.
  future.apply::future_lapply(
    seq_along(X),
    function(i, ...) run_task(X[[i]], i, ...),
    ...,
    future.seed = FALSE,
    future.scheduling = FALSE
  )
}

#' Save/Restore RNG State
#'
#' Captures the current RNG kind and `.Random.seed` so a code block can mutate
#' the RNG and then restore exactly what the caller had. Returns a list that is
#' fed back to `mv_rng_state_restore()`.
#'
#' @return A list with elements `kind` and `seed` (or `NULL` seed when the RNG
#'   had not yet been initialized).
#' @noRd
mv_rng_state_save <- function() {
  kind <- RNGkind()
  seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  list(kind = kind, seed = seed)
}

#' @param state List produced by `mv_rng_state_save()`.
#' @rdname mv_rng_state_save
#' @noRd
mv_rng_state_restore <- function(state) {
  if (is.null(state)) return(invisible(NULL))
  RNGkind(state$kind[1], state$kind[2], state$kind[3])
  if (is.null(state$seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", state$seed, envir = .GlobalEnv)
  }
  invisible(NULL)
}

#' Derive a Deterministic Per-Task Seed
#'
#' Maps a user seed and a 1-based task index to a positive integer seed that is
#' stable across R sessions and independent of execution order. Used so the
#' sequential and parallel code paths draw identical random numbers for a given
#' task index.
#'
#' The mix is `((base * 1000003) + idx * 2 - 1) mod .Machine$integer.max`. This
#' is a simple, fast, order-independent hash: collisions are only possible if
#' two distinct `(base, idx)` pairs map to the same value modulo
#' `.Machine$integer.max` (about 2.1e9), which requires either ~2.1e9 tasks or
#' adversarial base/idx choices. For the typical scale of this package (grids
#' of tens of `k` values, hundreds of `nstart` restarts) collisions are
#' effectively impossible. If far larger task counts are ever needed, switch to
#' a 64-bit hash (e.g. `digest::digest`).
#'
#' @param base_seed Integer user seed.
#' @param idx Positive integer task index (1-based).
#' @return A positive integer seed.
#' @noRd
mv_task_seed <- function(base_seed, idx) {
  stopifnot(is.numeric(base_seed), is.finite(base_seed),
            is.numeric(idx), idx >= 1)
  # Mix the two integers; .Machine$integer.max keeps it representable.
  cap <- .Machine$integer.max
  derived <- (as.integer(base_seed) * 1000003L + as.integer(idx) * 2L - 1L) %% cap
  if (derived <= 0L) derived <- derived + cap
  derived
}

#' Set Up Worker Environment
#'
#' Makes the Ampharos internals available inside a parallel worker. Because
#' `mv_future_lapply()` runs with `future.globals = FALSE`, the dispatched task
#' closure relies on these helpers being resolvable on the worker's search
#' path; this function guarantees that. If the `Ampharos` package is already
#' *attached* on the worker (so its symbols resolve as bare names), nothing is
#' sourced and the function returns immediately. Otherwise it sources the
#' `R/` directory into a private environment (never `.GlobalEnv`) located via
#' the `AMPHAROS_ROOT` environment variable (set by the test harness and by
#' `mv_with_parallel`) or by walking up from the worker's working directory.
#'
#' The private environment is attached just above `.GlobalEnv` on the worker so
#' the `mv_*` functions are visible to the dispatched task without polluting or
#' conflicting with the worker's global environment or other packages.
#'
#' @noRd
mv_parallel_worker_setup <- function() {
  # If the package is attached, its exported (and, for our own internal calls,
  # internal) symbols already resolve as bare names on the worker; skip the
  # fragile re-sourcing entirely. We test for attachment, not just namespace
  # load, because future does not attach packages on workers.
  if ("package:Ampharos" %in% search()) {
    return(invisible(NULL))
  }

  root <- Sys.getenv("AMPHAROS_ROOT", "")
  if (root == "" || !dir.exists(file.path(root, "R"))) {
    root <- getwd()
    while (!dir.exists(file.path(root, "R")) && dirname(root) != root) {
      root <- dirname(root)
    }
  }

  r_dir <- file.path(root, "R")
  if (!dir.exists(r_dir)) {
    return(invisible(NULL))
  }

  source_order <- c(
    "matrix_utils.R",
    "Matrix_Base.R",
    "matrix_init_whiten.R",
    "matrix_init_kmeans.R",
    "matrix_init_dbscan.R",
    "matrix_init_emrefine.R",
    "Matrix_Noise_BR.R",
    "KS_Score.R",
    "Matrix_Noise.R",
    "matrix_parallel.R"
  )

  # Source into a private, detached environment to avoid clobbering the
  # worker's .GlobalEnv or other attached packages. The parent is baseenv()
  # (not emptyenv()) so base primitives/operators (`<-`, `{`, `if`, `for`, ...)
  # resolve during sys.source; the ampharos helpers are all defined inside this
  # env. Note we also source matrix_parallel.R so mv_task_seed / mv_rng_state_*
  # are available to the run_task closure.
  private_env <- new.env(parent = baseenv())
  for (f in source_order) {
    fp <- file.path(r_dir, f)
    if (file.exists(fp)) {
      sys.source(fp, envir = private_env)
    }
  }
  # Attach just above .GlobalEnv so dispatched tasks can find the helpers.
  do.call(attach, list(what = private_env, name = "ampharos_worker_env",
                       pos = 2L, warn.conflicts = FALSE))
  invisible(NULL)
}

#' Null-Coalescing Helper
#'
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
