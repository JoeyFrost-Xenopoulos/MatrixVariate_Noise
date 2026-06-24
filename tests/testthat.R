library(testthat)

# Source all R files in the package
r_dir <- file.path(dirname(dirname(sys.frame(1)$ofile)), "R")
for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

test_check("Ampharos")
