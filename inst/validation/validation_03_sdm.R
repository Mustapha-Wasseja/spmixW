# Validation 3: SDM Panel (Chapter 2)
source("validation_helpers.R")

results <- list()

set.seed(10203040)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
rho_true <- 0.7
beta_true <- c(1, 1)
theta_true <- c(-1, -1)

W <- make_knn_W(N, k = 6, seed = 10203040)
Wbig <- kronecker(diag(Time), W)
X <- matrix(rnorm(nobs * k), ncol = k)
WX <- Wbig %*% X
fe <- make_fe(N, Time)
evec <- rnorm(nobs)

A <- diag(nobs) - rho_true * Wbig
y <- as.numeric(solve(A) %*% (X %*% beta_true + WX %*% theta_true + fe + evec))

print_header("Validation 3: SDM Panel — Both FE")
res <- sdm_panel(y, X, W, N, Time, ndraw = 5000, nomit = 1500,
                 prior = list(model = 3, rval = 0))

# Full beta vector: [beta1, beta2, theta1, theta2]
r <- c(
  check_param("beta_1", res$beta[1], 1.0, 0.15),
  check_param("beta_2", res$beta[2], 1.0, 0.15),
  check_param("theta_1", res$beta[3], -1.0, 0.15),
  check_param("theta_2", res$beta[4], -1.0, 0.15),
  check_param("rho", res$rho, 0.7, 0.08)
)
results <- c(results, r)

# Effects identity
cat("\n  Effects identity (direct + indirect = total):\n")
for (j in seq_len(res$p)) {
  dd <- res$direct_draws[, j] + res$indirect_draws[, j] - res$total_draws[, j]
  ok <- max(abs(dd)) < 1e-10
  results <- c(results, ok)
  cat(sprintf("    x%d: %s (max err = %.2e)\n", j, if (ok) "PASS" else "FAIL", max(abs(dd))))
}

# Indirect should be negative (theta dominates)
cat("  Indirect effects should be negative (theta=-1 dominates):\n")
for (j in seq_len(res$p)) {
  ie <- res$indirect[j, 1]
  ok <- ie < 0
  results <- c(results, ok)
  cat(sprintf("    x%d indirect = %.4f: %s\n", j, ie, if (ok) "PASS" else "FAIL"))
}

cat(sprintf("\nValidation 3 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
