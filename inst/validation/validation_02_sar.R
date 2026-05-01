# Validation 2: SAR Panel (Chapter 1, sar_panel_gd.m)
source("validation_helpers.R")

results <- list()

set.seed(10203040)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
rho_true <- 0.6
beta_true <- c(1, 1)
sige_true <- 1

W <- make_knn_W(N, k = 5, seed = 10203040)
Wbig <- kronecker(diag(Time), W)
X <- matrix(rnorm(nobs * k), ncol = k)
fe <- make_fe(N, Time)
evec <- rnorm(nobs, sd = sqrt(sige_true))

A <- diag(nobs) - rho_true * Wbig
y <- as.numeric(solve(A) %*% (X %*% beta_true + fe + evec))

print_header("Validation 2: SAR Panel — Homoscedastic, Both FE")
res <- sar_panel(y, X, W, N, Time, ndraw = 5000, nomit = 1500,
                 prior = list(model = 3, rval = 0))

r <- c(
  check_param("beta_1", res$beta[1], 1.0, 0.10, 0.9845),
  check_param("beta_2", res$beta[2], 1.0, 0.10, 1.0100),
  check_param("rho", res$rho, 0.6, 0.05, 0.5990),
  check_param("sigma^2", res$sige, 1.0, 0.25, 0.8694)
)
results <- c(results, r)

# Effects checks
cat("\n  Effects checks:\n")
for (j in seq_len(res$p)) {
  di <- res$direct[j, 1]
  ii <- res$indirect[j, 1]
  ti <- res$total[j, 1]
  cat(sprintf("    x%d: direct=%.4f, indirect=%.4f, total=%.4f\n", j, di, ii, ti))

  # direct + indirect = total (algebraic identity)
  dd <- res$direct_draws[, j] + res$indirect_draws[, j] - res$total_draws[, j]
  r <- c(r, max(abs(dd)) < 1e-10)
  cat(sprintf("    direct + indirect = total: %s (max err = %.2e)\n",
              if (max(abs(dd)) < 1e-10) "PASS" else "FAIL", max(abs(dd))))

  # direct > raw beta (spatial multiplier)
  r <- c(r, di > res$beta[j])
  cat(sprintf("    direct (%.3f) > beta (%.3f): %s\n",
              di, res$beta[j], if (di > res$beta[j]) "PASS" else "FAIL"))
}
results <- c(results, r)

# Total/Direct ratio check (~1/(1-rho) = 2.5)
ratio <- mean(res$total[, 1]) / mean(res$direct[, 1])
cat(sprintf("  Total/Direct ratio: %.2f (expected ~%.2f)\n", ratio, 1/(1-rho_true)))
r <- c(abs(ratio - 1/(1-rho_true)) < 0.5)
results <- c(results, r)

cat(sprintf("\nValidation 2 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
