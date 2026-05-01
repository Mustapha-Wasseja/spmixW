# Validation 1: OLS Panel (Chapter 1, ols_panel_gd.m)
source("validation_helpers.R")

results <- list()

# ---- DGP matching MATLAB rng(10203040) ----
set.seed(10203040)
N <- 200; Time <- 10; k <- 2
nobs <- N * Time
beta_true <- c(1, 1)
sige_true <- 5

X <- matrix(rnorm(nobs * k), ncol = k)
fe <- make_fe(N, Time)
evec <- rnorm(nobs, sd = sqrt(sige_true))
y <- as.numeric(X %*% beta_true + fe + evec)

# ---- Run (a): Homoscedastic, both FE ----
print_header("Validation 1a: OLS Panel — Homoscedastic, Both FE")
res_a <- ols_panel(y, X, N, Time, ndraw = 2500, nomit = 500,
                   prior = list(model = 3, rval = 0))

r <- c(
  check_param("beta_1", res_a$beta[1], 1.0, 0.15, 0.9654),
  check_param("beta_2", res_a$beta[2], 1.0, 0.15, 1.0214),
  check_param("sigma^2", res_a$sige, 5.0, 1.5, 4.3532)
)
results <- c(results, r)

# ---- Run (b): Heteroscedastic + outliers ----
print_header("Validation 1b: OLS Panel — Heteroscedastic, Outliers in T=5,6")
y_out <- y
# Add outliers to periods 5 and 6
idx5 <- ((5-1)*N + 1):(5*N)
idx6 <- ((6-1)*N + 1):(6*N)
y_out[idx5] <- y_out[idx5] + 10
y_out[idx6] <- y_out[idx6] + 10

res_b <- ols_panel(y_out, X, N, Time, ndraw = 2500, nomit = 500,
                   prior = list(model = 1, rval = 5))

r <- c(
  check_param("beta_1", res_b$beta[1], 1.0, 0.25, 0.9079),
  check_param("beta_2", res_b$beta[2], 1.0, 0.15, 1.0228)
)
# Check that vi for outlier periods are larger
vi_mean <- res_b$vmean
vi_t5 <- mean(vi_mean[idx5])
vi_other <- mean(vi_mean[-c(idx5, idx6)])
cat(sprintf("  vi(t=5) / vi(other) = %.2f (should be > 1.5)\n", vi_t5 / vi_other))
r <- c(r, vi_t5 / vi_other > 1.5)
results <- c(results, r)

# ---- Run (c): Informative prior pulling beta toward 0.5 ----
print_header("Validation 1c: OLS Panel — Tight Prior beta=(0.5,0.5)")
res_c <- ols_panel(y, X, N, Time, ndraw = 2500, nomit = 500,
                   prior = list(model = 1, rval = 0,
                                beta_mean = c(0.5, 0.5),
                                beta_var = diag(2) * 0.001))

# Beta should be pulled toward 0.5 (between 0.5 and 1.0)
r <- c(
  check_param("beta_1 (biased)", res_c$beta[1], 0.5, 0.3, 0.6285),
  check_param("beta_2 (biased)", res_c$beta[2], 0.5, 0.3, 0.6496)
)
cat(sprintf("  beta_1 < 1.0 (pulled toward prior): %s\n",
            if (res_c$beta[1] < 1.0) "PASS" else "FAIL"))
r <- c(r, res_c$beta[1] < 1.0)
results <- c(results, r)

cat(sprintf("\nValidation 1 Summary: %d/%d passed\n", sum(unlist(results)), length(unlist(results))))
