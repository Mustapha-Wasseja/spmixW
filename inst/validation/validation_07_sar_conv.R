# Validation 7: SAR Convex Combination (Chapter 3, sar_conv_panel_gd.m)
source("validation_helpers.R")

results <- list()

# Scaled-down version: N=200 instead of N=3111
set.seed(10203444)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
rho_true <- 0.6
beta_true <- c(1, 1)
gamma_true <- c(0.3, 0.7, 0.0)

# Three distinct W matrices
W1 <- make_knn_W(N, k = 5, seed = 10001)   # "contiguity-like"
W2 <- make_knn_W(N, k = 6, seed = 20002)   # "6-nn"
W3 <- make_knn_W(N, k = 8, seed = 30003)   # "distance-like"

Wbigs <- lapply(list(W1, W2, W3), function(W) kronecker(diag(Time), W))
Wc <- gamma_true[1] * Wbigs[[1]] + gamma_true[2] * Wbigs[[2]] + gamma_true[3] * Wbigs[[3]]

X <- matrix(rnorm(nobs * k), ncol = k)
fe <- make_fe(N, Time)
evec <- rnorm(nobs)

A <- diag(nobs) - rho_true * as.matrix(Wc)
y <- as.numeric(solve(A) %*% (X %*% beta_true + fe + evec))

# ---- Run with M=3 ----
print_header("Validation 7a: SAR Conv Panel — M=3, true gamma=(0.3, 0.7, 0.0)")
res3 <- sar_conv_panel(y, X, list(W1, W2, W3), N, Time,
                       ndraw = 20000, nomit = 5000,
                       prior = list(model = 3))

r <- c(
  check_param("rho", res3$rho, 0.6, 0.15),
  check_param("beta_1", res3$beta[1], 1.0, 0.15),
  check_param("beta_2", res3$beta[2], 1.0, 0.15)
)
results <- c(results, r)

# Gamma ordering: gamma2 > gamma1 > gamma3
cat(sprintf("  gamma = (%.3f, %.3f, %.3f)\n", res3$gamma[1], res3$gamma[2], res3$gamma[3]))
r_g <- c(
  res3$gamma[2] > res3$gamma[1],  # gamma2 should be largest
  res3$gamma[1] > res3$gamma[3]   # gamma3 should be smallest
)
cat(sprintf("  gamma2 > gamma1: %s\n", if (r_g[1]) "PASS" else "FAIL"))
cat(sprintf("  gamma1 > gamma3: %s\n", if (r_g[2]) "PASS" else "FAIL"))
results <- c(results, r_g)

# Acceptance rates
cat(sprintf("  rho acceptance rate: %.3f\n", res3$rho_acc_rate))
cat(sprintf("  gamma acceptance rate: %.3f\n", res3$gamma_acc_rate))
r_acc <- c(
  res3$rho_acc_rate > 0.15 && res3$rho_acc_rate < 0.70,
  res3$gamma_acc_rate > 0.05
)
cat(sprintf("  rho acc in [0.15, 0.70]: %s\n", if (r_acc[1]) "PASS" else "FAIL"))
cat(sprintf("  gamma acc > 0.05: %s\n", if (r_acc[2]) "PASS" else "FAIL"))
results <- c(results, r_acc)

# ---- Run with M=2 (only true matrices) ----
print_header("Validation 7b: SAR Conv Panel — M=2, true gamma=(0.3, 0.7)")
res2 <- sar_conv_panel(y, X, list(W1, W2), N, Time,
                       ndraw = 15000, nomit = 5000,
                       prior = list(model = 3))

r <- c(
  check_param("rho (M=2)", res2$rho, 0.6, 0.15),
  check_param("gamma_1 (M=2)", res2$gamma[1], 0.3, 0.15),
  check_param("gamma_2 (M=2)", res2$gamma[2], 0.7, 0.15)
)
results <- c(results, r)

cat(sprintf("\nValidation 7 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
