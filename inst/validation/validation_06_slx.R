# Validation 6: SLX Panel (Chapter 2)
source("validation_helpers.R")

results <- list()

set.seed(10203040)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
beta_true <- c(1, 1)
theta_true <- c(-1, -1)

W <- make_knn_W(N, k = 6, seed = 10203040)
Wbig <- kronecker(diag(Time), W)
X <- matrix(rnorm(nobs * k), ncol = k)
WX <- Wbig %*% X
fe <- make_fe(N, Time)
evec <- rnorm(nobs)

y <- as.numeric(X %*% beta_true + WX %*% theta_true + fe + evec)

print_header("Validation 6: SLX Panel — Both FE")
res <- slx_panel(y, X, W, N, Time, ndraw = 5000, nomit = 1500,
                 prior = list(model = 3, rval = 0))

r <- c(
  check_param("beta_1", res$beta[1], 1.0, 0.10),
  check_param("beta_2", res$beta[2], 1.0, 0.10)
)
results <- c(results, r)

# SLX: direct = beta, indirect = theta (exact, no spatial multiplier)
cat("\n  SLX effects (exact: direct=beta, indirect=theta):\n")
for (j in seq_len(res$p)) {
  de <- res$direct[j, 1]
  ie <- res$indirect[j, 1]
  r2 <- c(abs(de - beta_true[j]) < 0.10, abs(ie - theta_true[j]) < 0.15)
  results <- c(results, r2)
  cat(sprintf("    x%d: direct=%.4f (true=%.1f) %s, indirect=%.4f (true=%.1f) %s\n",
              j, de, beta_true[j], if (r2[1]) "PASS" else "FAIL",
              ie, theta_true[j], if (r2[2]) "PASS" else "FAIL"))
}

cat(sprintf("\nValidation 6 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
