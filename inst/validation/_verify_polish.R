devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

set.seed(42)
coords <- cbind(runif(80), runif(80))
W1 <- make_knw(coords, k = 4)
W2 <- make_knw(coords, k = 8)

# Generate data cleanly
panel <- simulate_panel(
  N = 80, T = 10,
  W = list(W1, W2),
  gamma = c(0.7, 0.3),
  rho = 0.4,
  beta = c(1, 1),
  theta = c(-0.5, 0.3),
  seed = 42
)

# Estimate SDM convex with named W matrices
res <- spmodel(y ~ x1 + x2,
               data = panel,
               W = list(geography = W1, trade = W2),
               model = "sdm",
               id = "region", time = "year",
               effects = "twoway",
               ndraw = 8000, nomit = 2000, thin = 2)

# CHECK 1-4: print
print(res)

# CHECK 5: all accessors work without warnings
cat("\n--- Accessor checks ---\n")
cat("rho mean:", mean(res$rho_draws), "\n")
cat("rho scalar:", res$rho, "\n")
cat("beta:", res$beta, "\n")
cat("gamma:", res$gamma, "\n")
cat("sigma2:", res$sigma2, "\n")
cat("beta_draws dims:", dim(as.matrix(res$beta_draws)), "\n")
cat("sigma_draws length:", length(res$sigma_draws), "\n")
cat("gamma_draws dims:", dim(res$gamma_draws), "\n")
cat("rho_draws length:", length(res$rho_draws), "\n")
cat("rho_draws is numeric?", is.numeric(res$rho_draws), "\n")

# CHECK 6: tidy() returns proper variable names
cat("\n--- Tidy output ---\n")
print(tidy.spmixW(res))

# CHECK 7: simulate_panel returns clean data frame
cat("\n--- Simulated data preview ---\n")
print(head(panel))
cat("Columns:", names(panel), "\n")
cat("Dimensions:", nrow(panel), "x", ncol(panel), "\n")

cat("\n=== All checks passed! ===\n")
