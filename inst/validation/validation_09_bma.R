# Validation 9: SAR Convex BMA (Chapter 6)
source("validation_helpers.R")

results <- list()

# M=2, DGP uses only W1 (gamma = (1, 0))
set.seed(221010)
N <- 100; Time <- 5; k <- 2
nobs <- N * Time
rho_true <- 0.5
beta_true <- c(1, -0.5)

W1 <- make_knn_W(N, k = 5, seed = 11111)
W2 <- make_knn_W(N, k = 8, seed = 22222)

W1big <- kronecker(diag(Time), W1)
X <- matrix(rnorm(nobs * k), ncol = k)
fe <- make_fe(N, Time)
evec <- rnorm(nobs)

A <- diag(nobs) - rho_true * W1big
y <- as.numeric(solve(A) %*% (X %*% beta_true + fe + evec))

print_header("Validation 9: SAR Convex BMA — M=2, true DGP uses W1 only")
suppressMessages(
  res <- sar_conv_bma(y, X, list(W1, W2), N, Time,
                      ndraw = 10000, nomit = 3000,
                      prior = list(model = 3))
)

# Number of models = 2^2 - 1 = 3
r <- res$nmodels == 3
cat(sprintf("  Number of models evaluated: %d (expected 3): %s\n",
            res$nmodels, if (r) "PASS" else "FAIL"))
results <- c(results, r)

# Model probs
cat("  Model probabilities:\n")
for (i in seq_len(res$nmodels)) {
  subset_str <- paste(which(res$subsets[i, ] == 1), collapse = ",")
  cat(sprintf("    Model %d (W={%s}): prob=%.4f, logm=%.2f\n",
              i, subset_str, res$probs[i], res$logm[i]))
}

# Model with W1 alone should have highest or near-highest prob
# (W1 is the true generator)
# Model 1 = {W1}, Model 2 = {W2}, Model 3 = {W1,W2}
# W2-only model should have lowest probability
r <- res$probs[2] < max(res$probs[1], res$probs[3])
cat(sprintf("  W2-only prob < best model prob: %s\n", if (r) "PASS" else "FAIL"))
results <- c(results, r)

# BMA beta recovery
r <- c(
  check_param("BMA beta_1", res$beta[1], 1.0, 0.25),
  check_param("BMA beta_2", res$beta[2], -0.5, 0.25)
)
results <- c(results, r)

# Probs sum to 1
r <- abs(sum(res$probs) - 1) < 1e-10
cat(sprintf("  Probs sum to 1: %s\n", if (r) "PASS" else "FAIL"))
results <- c(results, r)

cat(sprintf("\nValidation 9 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
