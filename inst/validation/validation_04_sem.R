# Validation 4: SEM Panel (Chapter 2)
source("validation_helpers.R")

results <- list()

set.seed(10203040)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
lambda_true <- 0.7
beta_true <- c(1, 1)

W <- make_knn_W(N, k = 6, seed = 10203040)
Wbig <- kronecker(diag(Time), W)
X <- matrix(rnorm(nobs * k), ncol = k)
fe <- make_fe(N, Time)
evec <- rnorm(nobs)

B <- diag(nobs) - lambda_true * Wbig
u <- as.numeric(solve(B) %*% evec)
y <- as.numeric(X %*% beta_true + fe + u)

print_header("Validation 4: SEM Panel — Both FE")
res <- sem_panel(y, X, W, N, Time, ndraw = 5000, nomit = 1500,
                 prior = list(model = 3, rval = 0))

r <- c(
  check_param("beta_1", res$beta[1], 1.0, 0.10),
  check_param("beta_2", res$beta[2], 1.0, 0.10),
  check_param("lambda", res$rho, 0.7, 0.08)
)
results <- c(results, r)

cat(sprintf("\nValidation 4 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
