# Benchmark helpers for HC noise clustering.
#
# Usage:
#   devtools::load_all()
#   source("scripts/benchmark_estimate_k_hc.R")
#
# The functions below assume the Ampharos objects are already loaded into the
# current session by `devtools::load_all()`.

mv_hc_benchmark_simulate <- function(n_per_group,
                                     contamination = 0,
                                     g = length(n_per_group),
                                     component_means = NULL,
                                     row_sd = 0.35,
                                     col_sd = 0.35,
                                     contamination_type = c("column_replace", "full_uniform"),
                                     contamination_range = c(-15, 15),
                                     seed = NULL) {
  contamination_type <- match.arg(contamination_type)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(n_per_group) == 1L) {
    n_per_group <- rep.int(as.integer(n_per_group), g)
  }

  if (!is.numeric(n_per_group) || any(n_per_group < 1)) {
    stop("'n_per_group' must be a positive integer vector or scalar.")
  }

  n_per_group <- as.integer(n_per_group)
  g <- as.integer(g)

  if (is.null(component_means)) {
    component_means <- mv_hc_benchmark_default_means(g = g)
  }

  if (length(component_means) != g) {
    stop("'component_means' must contain exactly one mean matrix per group.")
  }

  row_count <- nrow(component_means[[1]])
  col_count <- ncol(component_means[[1]])

  if (any(vapply(component_means, nrow, integer(1)) != row_count) ||
      any(vapply(component_means, ncol, integer(1)) != col_count)) {
    stop("All matrices in 'component_means' must have the same dimensions.")
  }

  clean_samples <- vector("list", sum(n_per_group))
  truth_labels <- integer(sum(n_per_group))
  sample_index <- 1L

  for (component_id in seq_len(g)) {
    for (replicate_id in seq_len(n_per_group[component_id])) {
      clean_samples[[sample_index]] <- mv_hc_benchmark_simulate_one(
        mean_matrix = component_means[[component_id]],
        row_sd = row_sd,
        col_sd = col_sd
      )
      truth_labels[sample_index] <- component_id
      sample_index <- sample_index + 1L
    }
  }

  contamination_count <- as.integer(round(sum(n_per_group) * contamination))
  contamination_count <- max(0L, contamination_count)

  noise_samples <- list()
  if (contamination_count > 0L) {
    for (noise_index in seq_len(contamination_count)) {
      noise_samples[[noise_index]] <- mv_hc_benchmark_simulate_noise(
        row_count = row_count,
        col_count = col_count,
        contamination_type = contamination_type,
        contamination_range = contamination_range
      )
    }
  }

  x_list <- c(clean_samples, noise_samples)
  truth <- c(truth_labels, rep.int(0L, contamination_count))

  list(
    x_list = x_list,
    truth = truth,
    contamination = contamination,
    contamination_count = contamination_count,
    contamination_type = contamination_type,
    g = g,
    n_per_group = n_per_group,
    component_means = component_means
  )
}

mv_hc_benchmark_default_means <- function(g, row_count = 2L, col_count = 3L, spread = 2) {
  g <- as.integer(g)
  offsets <- seq_len(g) - (g + 1) / 2

  lapply(offsets, function(offset) {
    matrix(
      offset * spread,
      nrow = row_count,
      ncol = col_count
    )
  })
}

mv_hc_benchmark_simulate_one <- function(mean_matrix, row_sd = 0.35, col_sd = 0.35) {
  row_count <- nrow(mean_matrix)
  col_count <- ncol(mean_matrix)

  row_scale <- diag(row_sd, row_count)
  col_scale <- diag(col_sd, col_count)

  mean_matrix + row_scale %*% matrix(rnorm(row_count * col_count), row_count, col_count) %*% col_scale
}

mv_hc_benchmark_simulate_noise <- function(row_count,
                                           col_count,
                                           contamination_type = c("column_replace", "full_uniform"),
                                           contamination_range = c(-15, 15)) {
  contamination_type <- match.arg(contamination_type)

  if (contamination_type == "full_uniform") {
    return(matrix(
      runif(row_count * col_count, min = contamination_range[1], max = contamination_range[2]),
      nrow = row_count,
      ncol = col_count
    ))
  }

  noisy_matrix <- matrix(rnorm(row_count * col_count), nrow = row_count, ncol = col_count)
  noisy_column <- sample.int(col_count, 1L)
  noisy_matrix[, noisy_column] <- runif(row_count, min = contamination_range[1], max = contamination_range[2])
  noisy_matrix
}

mv_hc_benchmark_generate_permutations <- function(values) {
  values <- as.integer(values)

  if (length(values) <= 1L) {
    return(list(values))
  }

  permutations <- list()
  position <- 1L

  for (index in seq_along(values)) {
    remainder <- values[-index]
    tails <- mv_hc_benchmark_generate_permutations(remainder)

    for (tail in tails) {
      permutations[[position]] <- c(values[index], tail)
      position <- position + 1L
    }
  }

  permutations
}

mv_hc_benchmark_align_labels <- function(predicted, truth, g = max(truth, predicted, na.rm = TRUE)) {
  predicted <- as.integer(predicted)
  truth <- as.integer(truth)
  g <- as.integer(g)

  if (g < 1L) {
    stop("'g' must be at least 1.")
  }

  if (g > 8L) {
    stop("Exact label alignment is only implemented for g <= 8.")
  }

  if (length(predicted) != length(truth)) {
    stop("'predicted' and 'truth' must have the same length.")
  }

  candidate_mappings <- mv_hc_benchmark_generate_permutations(seq_len(g))
  best_accuracy <- -Inf
  best_labels <- predicted
  best_mapping <- seq_len(g)

  for (mapping in candidate_mappings) {
    relabeled <- predicted

    for (component_id in seq_len(g)) {
      relabeled[predicted == component_id] <- mapping[component_id]
    }

    current_accuracy <- mean(relabeled == truth)

    if (is.finite(current_accuracy) && current_accuracy > best_accuracy) {
      best_accuracy <- current_accuracy
      best_labels <- relabeled
      best_mapping <- mapping
    }
  }

  list(
    labels = best_labels,
    mapping = best_mapping,
    accuracy = best_accuracy
  )
}

mv_hc_benchmark_score <- function(predicted, truth, g = max(truth, predicted, na.rm = TRUE)) {
  aligned <- mv_hc_benchmark_align_labels(predicted, truth, g = g)

  overall_accuracy <- mean(aligned$labels == truth)

  non_noise_mask <- truth != 0L
  noise_mask <- truth == 0L
  predicted_noise_mask <- aligned$labels == 0L

  non_noise_accuracy <- if (any(non_noise_mask)) {
    mean(aligned$labels[non_noise_mask] == truth[non_noise_mask])
  } else {
    NA_real_
  }

  noise_accuracy <- if (any(noise_mask)) {
    mean(predicted_noise_mask[noise_mask])
  } else {
    NA_real_
  }

  noise_precision <- if (any(predicted_noise_mask)) {
    mean(truth[predicted_noise_mask] == 0L)
  } else {
    NA_real_
  }

  noise_recall <- noise_accuracy

  noise_f1 <- if (is.na(noise_precision) || is.na(noise_recall) ||
      (noise_precision + noise_recall) == 0) {
    NA_real_
  } else {
    2 * noise_precision * noise_recall / (noise_precision + noise_recall)
  }

  balanced_accuracy <- mean(c(non_noise_accuracy, noise_accuracy), na.rm = TRUE)

  list(
    overall_accuracy = overall_accuracy,
    balanced_accuracy = balanced_accuracy,
    non_noise_accuracy = non_noise_accuracy,
    noise_accuracy = noise_accuracy,
    noise_precision = noise_precision,
    noise_recall = noise_recall,
    noise_f1 = noise_f1,
    aligned_labels = aligned$labels,
    label_mapping = aligned$mapping
  )
}

mv_hc_benchmark_grid_for_dims <- function(row_count, col_count, n_points = 60L) {
  dimension <- row_count * col_count

  if (!is.finite(dimension) || dimension <= 0) {
    return(10^seq(-16, -1, length.out = n_points))
  }

  center_log10 <- -0.75 * dimension
  half_width <- max(8, ceiling(dimension))
  lower_log10 <- max(log10(.Machine$double.xmin), center_log10 - half_width)
  upper_log10 <- center_log10 + half_width

  grid <- 10^seq(lower_log10, upper_log10, length.out = n_points)
  grid <- grid[is.finite(grid) & grid > 0]

  if (length(grid) < 2L) {
    grid <- 10^seq(-16, -1, length.out = n_points)
  }

  sort(unique(grid))
}

mv_hc_benchmark_default_k_grid <- function(x_list, n_points = 60L) {
  mv_hc_benchmark_grid_for_dims(
    row_count = nrow(x_list[[1]]),
    col_count = ncol(x_list[[1]]),
    n_points = n_points
  )
}

mv_hc_benchmark_get_k_grid <- function(x_list, k_grid = NULL) {
  if (!is.null(k_grid)) {
    k_grid <- sort(unique(as.numeric(k_grid)))
    k_grid <- k_grid[is.finite(k_grid) & k_grid > 0]

    if (length(k_grid) >= 2L) {
      return(k_grid)
    }
  }

  mv_hc_benchmark_default_k_grid(x_list)
}

mv_hc_benchmark_collapse_ranges <- function(k_values = NULL, indices, k_grid = NULL) {
  if (is.null(k_values) && !is.null(k_grid)) {
    k_values <- k_grid
  }

  if (is.null(k_values)) {
    stop("Either 'k_values' or 'k_grid' must be provided.")
  }

  if (!length(indices)) {
    return(data.frame(
      k_min = numeric(0),
      k_max = numeric(0),
      n_values = integer(0)
    ))
  }

  indices <- sort(unique(as.integer(indices)))
  run_breaks <- c(0L, which(diff(indices) != 1L), length(indices))
  range_count <- length(run_breaks) - 1L

  data.frame(
    k_min = vapply(seq_len(range_count), function(range_index) {
      start_index <- run_breaks[range_index] + 1L
      k_values[indices[start_index]]
    }, numeric(1)),
    k_max = vapply(seq_len(range_count), function(range_index) {
      end_index <- run_breaks[range_index + 1L]
      k_values[indices[end_index]]
    }, numeric(1)),
    n_values = vapply(seq_len(range_count), function(range_index) {
      run_breaks[range_index + 1L] - run_breaks[range_index]
    }, integer(1))
  )
}

mv_hc_benchmark_estimate_k <- function(x_list,
                                       truth,
                                       g,
                                       initializations = c("kmeans", "emrefine", "dbscan"),
                                       contamination = NA_real_,
                                       max_iter = 100,
                                       tol = 1e-06,
                                       nstart = 100,
                                       adaptive_grid = TRUE,
                                       k_grid = NULL,
                                       noise_pi_init = 0.05,
                                       use_parallel = FALSE,
                                       n_cores = NULL,
                                       verbose = FALSE) {
  k_grid <- mv_hc_benchmark_get_k_grid(x_list, k_grid)

  results <- vector("list", length(initializations))

  for (initialization_index in seq_along(initializations)) {
    initialization <- initializations[initialization_index]
    fit_error <- NULL
    fit_time <- system.time({
      fit <- tryCatch(
        mv_noise_fit(
          x_list = x_list,
          g = g,
          noise_type = "hc",
          max_iter = max_iter,
          tol = tol,
          nstart = nstart,
          estimate_k = TRUE,
          k_grid = k_grid,
          adaptive_grid = adaptive_grid,
          noise_pi_init = noise_pi_init,
          init = initialization,
          verbose = verbose,
          use_parallel = use_parallel,
          n_cores = n_cores
        ),
        error = function(e) {
          fit_error <<- conditionMessage(e)
          NULL
        }
      )
    })

    if (is.null(fit)) {
      results[[initialization_index]] <- data.frame(
        initialization = initialization,
        contamination = contamination,
        selected_k = NA_real_,
        overall_accuracy = NA_real_,
        balanced_accuracy = NA_real_,
        non_noise_accuracy = NA_real_,
        noise_accuracy = NA_real_,
        noise_precision = NA_real_,
        noise_recall = NA_real_,
        noise_f1 = NA_real_,
        logLik = NA_real_,
        iterations = NA_integer_,
        converged = FALSE,
        elapsed_sec = unname(fit_time["elapsed"]),
        status = "failed",
        error_message = fit_error,
        k_grid_n = length(k_grid),
        k_grid_min = min(k_grid),
        k_grid_max = max(k_grid),
        stringsAsFactors = FALSE
      )
      next
    }

    score <- mv_hc_benchmark_score(fit$cluster, truth, g = g)
    selected_k <- fit$k_selection$selected_k

    results[[initialization_index]] <- data.frame(
      initialization = initialization,
      contamination = contamination,
      selected_k = selected_k,
      overall_accuracy = score$overall_accuracy,
      balanced_accuracy = score$balanced_accuracy,
      non_noise_accuracy = score$non_noise_accuracy,
      noise_accuracy = score$noise_accuracy,
      noise_precision = score$noise_precision,
      noise_recall = score$noise_recall,
      noise_f1 = score$noise_f1,
      logLik = tail(fit$logLik, 1L),
      iterations = fit$iterations,
      converged = fit$converged,
      elapsed_sec = unname(fit_time["elapsed"]),
      status = "ok",
      error_message = NA_character_,
      k_grid_n = length(k_grid),
      k_grid_min = min(k_grid),
      k_grid_max = max(k_grid),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

mv_hc_benchmark_exhaustive_k <- function(x_list,
                                         truth,
                                         g,
                                         initializations = c("kmeans", "emrefine", "dbscan"),
                                         k_grid = NULL,
                                         select_by = c("balanced_accuracy", "overall_accuracy", "noise_f1"),
                                         max_iter = 100,
                                         tol = 1e-06,
                                         nstart = 100,
                                         noise_pi_init = 0.05,
                                         use_parallel = FALSE,
                                         n_cores = NULL,
                                         verbose = FALSE) {
  select_by <- match.arg(select_by)

  k_grid <- mv_hc_benchmark_get_k_grid(x_list, k_grid)

  all_results <- list()

  for (initialization in initializations) {
    per_k_results <- vector("list", length(k_grid))

    for (k_index in seq_along(k_grid)) {
      candidate_k <- k_grid[k_index]
      fit_error <- NULL
      fit_time <- system.time({
        fit <- tryCatch(
          mv_noise_fit(
            x_list = x_list,
            g = g,
            noise_type = "hc",
            max_iter = max_iter,
            tol = tol,
            nstart = nstart,
            estimate_k = FALSE,
            noise_k = candidate_k,
            noise_pi_init = noise_pi_init,
            init = initialization,
            verbose = verbose,
            use_parallel = use_parallel,
            n_cores = n_cores
          ),
          error = function(e) {
            fit_error <<- conditionMessage(e)
            NULL
          }
        )
      })

      if (is.null(fit)) {
        per_k_results[[k_index]] <- data.frame(
          initialization = initialization,
          k = candidate_k,
          overall_accuracy = NA_real_,
          balanced_accuracy = NA_real_,
          non_noise_accuracy = NA_real_,
          noise_accuracy = NA_real_,
          noise_precision = NA_real_,
          noise_recall = NA_real_,
          noise_f1 = NA_real_,
          logLik = NA_real_,
          iterations = NA_integer_,
          converged = FALSE,
          elapsed_sec = unname(fit_time["elapsed"]),
          status = "failed",
          error_message = fit_error,
          stringsAsFactors = FALSE
        )
        next
      }

      score <- mv_hc_benchmark_score(fit$cluster, truth, g = g)

      per_k_results[[k_index]] <- data.frame(
        initialization = initialization,
        k = candidate_k,
        overall_accuracy = score$overall_accuracy,
        balanced_accuracy = score$balanced_accuracy,
        non_noise_accuracy = score$non_noise_accuracy,
        noise_accuracy = score$noise_accuracy,
        noise_precision = score$noise_precision,
        noise_recall = score$noise_recall,
        noise_f1 = score$noise_f1,
        logLik = tail(fit$logLik, 1L),
        iterations = fit$iterations,
        converged = fit$converged,
        elapsed_sec = unname(fit_time["elapsed"]),
        status = "ok",
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    per_k_table <- do.call(rbind, per_k_results)
    score_values <- per_k_table[[select_by]]
    finite_scores <- is.finite(score_values)

    if (!any(finite_scores)) {
      best_value <- NA_real_
      best_indices <- integer(0)
      best_ranges <- data.frame(
        k_min = numeric(0),
        k_max = numeric(0),
        n_values = integer(0)
      )
    } else {
      best_value <- max(score_values[finite_scores])
      best_indices <- which(finite_scores & abs(score_values - best_value) <= sqrt(.Machine$double.eps))
      best_ranges <- mv_hc_benchmark_collapse_ranges(k_values = per_k_table$k, indices = best_indices)
    }

    all_results[[initialization]] <- list(
      results = per_k_table,
      select_by = select_by,
      best_value = best_value,
      best_indices = best_indices,
      best_k_ranges = best_ranges
    )
  }

  all_results
}

mv_hc_benchmark_run <- function(n_per_group,
                                contamination_levels = c(0, 0.05, 0.1, 0.2),
                                initializations = c("kmeans", "emrefine", "dbscan"),
                                g = length(n_per_group),
                                replicates = 10,
                                component_means = NULL,
                                row_sd = 0.35,
                                col_sd = 0.35,
                                contamination_type = c("column_replace", "full_uniform"),
                                contamination_range = c(-15, 15),
                                max_iter = 100,
                                tol = 1e-06,
                                nstart = 100,
                                adaptive_grid = TRUE,
                                k_grid = NULL,
                                noise_pi_init = 0.05,
                                use_parallel = FALSE,
                                n_cores = NULL,
                                seed = NULL,
                                verbose = FALSE) {
  contamination_type <- match.arg(contamination_type)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  benchmark_means <- if (is.null(component_means)) {
    mv_hc_benchmark_default_means(g = g)
  } else {
    component_means
  }

  if (is.null(k_grid)) {
    k_grid <- mv_hc_benchmark_grid_for_dims(
      row_count = nrow(benchmark_means[[1]]),
      col_count = ncol(benchmark_means[[1]]),
      n_points = 60L
    )
  } else {
    k_grid <- mv_hc_benchmark_get_k_grid(benchmark_means, k_grid)
  }

  auto_selection_results <- list()
  exhaustive_results <- list()

  run_index <- 1L

  for (contamination_level in contamination_levels) {
    for (replicate_index in seq_len(replicates)) {
      sim_seed <- if (is.null(seed)) NULL else seed + run_index
      simulated <- mv_hc_benchmark_simulate(
        n_per_group = n_per_group,
        contamination = contamination_level,
        g = g,
        component_means = component_means,
        row_sd = row_sd,
        col_sd = col_sd,
        contamination_type = contamination_type,
        contamination_range = contamination_range,
        seed = sim_seed
      )

      auto_selection_results[[run_index]] <- mv_hc_benchmark_estimate_k(
        x_list = simulated$x_list,
        truth = simulated$truth,
        g = g,
        initializations = initializations,
        contamination = contamination_level,
        max_iter = max_iter,
        tol = tol,
        nstart = nstart,
        adaptive_grid = adaptive_grid,
        k_grid = k_grid,
        noise_pi_init = noise_pi_init,
        use_parallel = use_parallel,
        n_cores = n_cores,
        verbose = verbose
      )

      exhaustive_results[[run_index]] <- mv_hc_benchmark_exhaustive_k(
        x_list = simulated$x_list,
        truth = simulated$truth,
        g = g,
        initializations = initializations,
        k_grid = k_grid,
        max_iter = max_iter,
        tol = tol,
        nstart = nstart,
        noise_pi_init = noise_pi_init,
        use_parallel = use_parallel,
        n_cores = n_cores,
        verbose = verbose
      )

      run_index <- run_index + 1L
    }
  }

  list(
    auto_selection = do.call(rbind, auto_selection_results),
    exhaustive = exhaustive_results
  )
}

mv_hc_benchmark_save_results <- function(results,
                                        output_dir,
                                        prefix = "hc_benchmark") {
  if (missing(results) || is.null(results)) {
    stop("'results' must be a benchmark result object.")
  }

  if (missing(output_dir) || !nzchar(output_dir)) {
    stop("'output_dir' must be a non-empty path.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  auto_path <- file.path(output_dir, paste0(prefix, "_auto_selection.csv"))
  exhaustive_path <- file.path(output_dir, paste0(prefix, "_exhaustive.csv"))
  ranges_path <- file.path(output_dir, paste0(prefix, "_best_k_ranges.csv"))

  if (is.data.frame(results)) {
    utils::write.csv(results, auto_path, row.names = FALSE)
    return(invisible(list(auto_selection = auto_path)))
  }

  if (!is.list(results) || is.null(results$auto_selection) || is.null(results$exhaustive)) {
    stop("'results' must be a data frame or a benchmark run object returned by mv_hc_benchmark_run().")
  }

  utils::write.csv(results$auto_selection, auto_path, row.names = FALSE)

  exhaustive_rows <- list()
  range_rows <- list()
  row_index <- 1L

  for (run_id in seq_along(results$exhaustive)) {
    run_result <- results$exhaustive[[run_id]]

    if (!is.list(run_result) || is.null(run_result)) {
      next
    }

    for (initialization in names(run_result)) {
      init_result <- run_result[[initialization]]

      if (is.null(init_result$results)) {
        next
      }

      init_table <- init_result$results
      init_table$run_id <- run_id
      init_table$initialization <- initialization
      init_table$select_by <- init_result$select_by
      init_table$best_value <- init_result$best_value
      exhaustive_rows[[row_index]] <- init_table
      row_index <- row_index + 1L

      if (!is.null(init_result$best_k_ranges) && nrow(init_result$best_k_ranges) > 0) {
        best_ranges <- init_result$best_k_ranges
        best_ranges$run_id <- run_id
        best_ranges$initialization <- initialization
        best_ranges$select_by <- init_result$select_by
        best_ranges$best_value <- init_result$best_value
        range_rows[[length(range_rows) + 1L]] <- best_ranges
      }
    }
  }

  if (length(exhaustive_rows) > 0) {
    exhaustive_table <- do.call(rbind, exhaustive_rows)
    utils::write.csv(exhaustive_table, exhaustive_path, row.names = FALSE)
  }

  if (length(range_rows) > 0) {
    ranges_table <- do.call(rbind, range_rows)
    utils::write.csv(ranges_table, ranges_path, row.names = FALSE)
  }

  invisible(list(
    auto_selection = auto_path,
    exhaustive = if (length(exhaustive_rows) > 0) exhaustive_path else NULL,
    best_k_ranges = if (length(range_rows) > 0) ranges_path else NULL
  ))
}
