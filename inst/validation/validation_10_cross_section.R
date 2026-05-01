# Validation 10: Cross-Sectional Mode (T=1, Chapter 7)
source("validation_helpers.R")

results <- list()

# ---- SAR cross-section ----
set.seed(10203444)
N <- 300; Time <- 1; k <- 2
nobs <- N
rho_true <- 0.6
beta_true <- c(1, 1)
intercept_true <- 2

W <- make_knn_W(N, k = 5, seed = 10203444)
X <- cbind(1, matrix(rnorm(N * k), ncol = k))  # intercept in first column
evec <- rnorm(N)

A <- diag(N) - rho_true * W
y <- as.numeric(solve(A) %*% (X %*% c(intercept_true, beta_true) + evec))

print_header("Validation 10a: SAR Cross-Section (T=1, model=0)")
res_sar <- sar_panel(y, X, W, N, Time, ndraw = 5000, nomit = 1500,
                     prior = list(model = 0, rval = 0))

r <- c(
  check_param("intercept", res_sar$beta[1], 2.0, 0.5),
  check_param("beta_1", res_sar$beta[2], 1.0, 0.15),
  check_param("beta_2", res_sar$beta[3], 1.0, 0.15),
  check_param("rho", res_sar$rho, 0.6, 0.08)
)
results <- c(results, r)

# ---- SAR Conv cross-section ----
print_header("Validation 10b: SAR Conv Cross-Section (T=1, M=2)")
W1 <- make_knn_W(N, k = 5, seed = 10001)
W2 <- make_knn_W(N, k = 8, seed = 20002)
gamma_true <- c(0.7, 0.3)
Wc <- gamma_true[1] * W1 + gamma_true[2] * W2

X2 <- matrix(rnorm(N * k), ncol = k)
A2 <- diag(N) - rho_true * Wc
y2 <- as.numeric(solve(A2) %*% (X2 %*% beta_true + evec))

res_conv <- sar_conv_panel(y2, X2, list(W1, W2), N, Time,
                           ndraw = 15000, nomit = 5000,
                           prior = list(model = 0))

r <- c(
  check_param("rho (conv)", res_conv$rho, 0.6, 0.2),
  check_param("beta_1 (conv)", res_conv$beta[1], 1.0, 0.2),
  check_param("beta_2 (conv)", res_conv$beta[2], 1.0, 0.2)
)
results <- c(results, r)

# Gamma ordering should be roughly correct
cat(sprintf("  gamma = (%.3f, %.3f) — true=(0.7, 0.3)\n", res_conv$gamma[1], res_conv$gamma[2]))
r <- res_conv$gamma[1] > res_conv$gamma[2]
cat(sprintf("  gamma1 > gamma2: %s\n", if (r) "PASS" else "FAIL"))
results <- c(results, r)

cat(sprintf("\nValidation 10 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
