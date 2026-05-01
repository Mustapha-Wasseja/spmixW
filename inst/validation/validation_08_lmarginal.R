# Validation 8: Model Comparison (Chapter 5)
source("validation_helpers.R")

results <- list()

# Generate data from SDM model
set.seed(221010)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
rho_true <- 0.4
beta_true <- c(1, 2)
theta_true <- c(1.5, -0.5)

W <- make_knn_W(N, k = 5, seed = 221010)
Wbig <- kronecker(diag(Time), W)
X <- matrix(rnorm(nobs * k), ncol = k)
WX <- Wbig %*% X
fe <- make_fe(N, Time)
evec <- rnorm(nobs)

A <- diag(nobs) - rho_true * Wbig
y <- as.numeric(solve(A) %*% (X %*% beta_true + WX %*% theta_true + fe + evec))

# Demean for model comparison
dm <- demean_panel(y, X, N, Time, model = 3)

print_header("Validation 8: Model Comparison — SDM DGP")
lm_res <- lmarginal_panel(dm$ywith, dm$xwith, W, N, Time)

cat(sprintf("  Log-marginals: SLX=%.2f, SDM=%.2f, SDEM=%.2f\n",
            lm_res$logm_slx, lm_res$logm_sdm, lm_res$logm_sdem))
cat(sprintf("  Probabilities: SLX=%.4f, SDM=%.4f, SDEM=%.4f\n",
            lm_res$probs[1], lm_res$probs[2], lm_res$probs[3]))

# SDM should win (data was generated from SDM)
r <- c(
  lm_res$logm_sdm > lm_res$logm_slx,
  lm_res$probs[2] > 0.3  # SDM prob should be substantial
)
cat(sprintf("  SDM logm > SLX logm: %s\n", if (r[1]) "PASS" else "FAIL"))
cat(sprintf("  SDM prob > 0.3: %s\n", if (r[2]) "PASS" else "FAIL"))
results <- c(results, r)

# Probs sum to 1
r <- abs(sum(lm_res$probs) - 1) < 1e-10
cat(sprintf("  Probs sum to 1: %s\n", if (r) "PASS" else "FAIL"))
results <- c(results, r)

cat(sprintf("\nValidation 8 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
