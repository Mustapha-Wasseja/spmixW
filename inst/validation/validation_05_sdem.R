# Validation 5: SDEM Panel (Chapter 2, sdem_panel_gd.m)
source("validation_helpers.R")

results <- list()

set.seed(10203040)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
lambda_true <- 0.7
beta_true <- c(1, 1)
theta_true <- c(-1, -1)
sige_true <- 0.1

W <- make_knn_W(N, k = 6, seed = 10203040)
Wbig <- kronecker(diag(Time), W)
X <- matrix(rnorm(nobs * k), ncol = k)
WX <- Wbig %*% X
fe <- make_fe(N, Time)
evec <- rnorm(nobs, sd = sqrt(sige_true))

B <- diag(nobs) - lambda_true * Wbig
u <- as.numeric(solve(B) %*% evec)
y <- as.numeric(X %*% beta_true + WX %*% theta_true + fe + u)

print_header("Validation 5: SDEM Panel — Both FE")
res <- sdem_panel(y, X, W, N, Time, ndraw = 5000, nomit = 1500,
                  prior = list(model = 3, rval = 0))

r <- c(
  check_param("beta_1", res$beta[1], 1.0, 0.10, 0.9885),
  check_param("beta_2", res$beta[2], 1.0, 0.10, 1.0094),
  check_param("theta_1", res$beta[3], -1.0, 0.15, -1.0459),
  check_param("theta_2", res$beta[4], -1.0, 0.15, -0.9455),
  check_param("lambda", res$rho, 0.7, 0.08, 0.6880)
)
results <- c(results, r)

# SDEM effects: direct = beta, indirect = theta
cat("\n  SDEM effects check (direct ~ beta, indirect ~ theta):\n")
for (j in seq_len(res$p)) {
  de <- res$direct[j, 1]
  ie <- res$indirect[j, 1]
  r2 <- c(abs(de - beta_true[j]) < 0.15, abs(ie - theta_true[j]) < 0.2)
  results <- c(results, r2)
  cat(sprintf("    x%d: direct=%.4f (true=%.1f) %s, indirect=%.4f (true=%.1f) %s\n",
              j, de, beta_true[j], if (r2[1]) "PASS" else "FAIL",
              ie, theta_true[j], if (r2[2]) "PASS" else "FAIL"))
}

cat(sprintf("\nValidation 5 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
