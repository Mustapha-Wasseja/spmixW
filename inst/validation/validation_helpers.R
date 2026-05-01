# validation_helpers.R — shared helpers for all validation scripts

library(spmixW)

# Check if a value is within tolerance of truth
check_param <- function(name, estimate, truth, tol, matlab_est = NA) {
  pass <- abs(estimate - truth) <= tol
  status <- if (pass) "PASS" else "FAIL"
  cat(sprintf("  %-20s | True=%-8.4f | R=%-8.4f | MATLAB=%-8s | tol=%-5.2f | %s\n",
              name, truth, estimate,
              if (is.na(matlab_est)) "N/A" else sprintf("%.4f", matlab_est),
              tol, status))
  pass
}

# Make a k-nearest-neighbour W from random coordinates
make_knn_W <- function(N, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  coords <- cbind(runif(N), runif(N))
  W <- as.matrix(normw(make_knw(coords, k, row_normalise = FALSE)))
  W
}

# Generate panel fixed effects matching MATLAB pattern
make_fe <- function(N, Time) {
  sfe <- (1:N) / N         # region FE
  tfe <- (1:Time) / Time   # time FE
  rep(sfe, times = Time) + rep(tfe, each = N)
}

# Print header
print_header <- function(title) {
  cat("\n", strrep("=", 70), "\n")
  cat(" ", title, "\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("  %-20s | %-8s | %-8s | %-8s | %-5s | %s\n",
              "Parameter", "True", "R Est", "MATLAB", "Tol", "Status"))
  cat(strrep("-", 70), "\n")
}
