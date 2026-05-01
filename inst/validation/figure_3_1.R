# Reproduce Figure 3.1 from panelg_manual.pdf Chapter 3
# SAR convex combination with 3 W matrices on US counties data

devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

plot_dir <- "C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW/inst/validation/plots"

# ============================================================
# 1. Read the US counties shapefile
# ============================================================
cat("Reading US counties shapefile...\n")
shp_path <- "C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/toolbox_panelg/demo_data/uscounties_projected.shp"

if (!requireNamespace("sf", quietly = TRUE)) {
  stop("Package 'sf' is required. Install with: install.packages('sf')")
}

counties <- sf::st_read(shp_path, quiet = TRUE)
coords_full <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(counties)))
N_full <- nrow(coords_full)
cat(sprintf("Full shapefile: N = %d counties\n", N_full))

# ============================================================
# 2. Subsample if N is too large for memory/speed
# ============================================================
N_target <- 500
if (N_full > N_target) {
  cat(sprintf("Subsampling %d counties for computational feasibility...\n", N_target))
  set.seed(10203444)
  idx <- sort(sample(N_full, N_target))
  latt <- coords_full[idx, 2]  # Y = latitude
  long <- coords_full[idx, 1]  # X = longitude
  N <- N_target
} else {
  latt <- coords_full[, 2]
  long <- coords_full[, 1]
  N <- N_full
}
cat(sprintf("Using N = %d counties\n", N))

# ============================================================
# 3. Build 3 weight matrices (matching MATLAB DGP)
# ============================================================
cat("Building weight matrices...\n")
coords <- cbind(long, latt)

# W1: Delaunay contiguity approximated by ~12 nearest neighbours
# (Delaunay on projected coords typically has 5-7 neighbours on average;
#  we use a kNN approximation since R's Delaunay is expensive for large N)
W_cont <- as.matrix(make_knw(coords, k = 12, row_normalise = TRUE))

# W2: 6 nearest neighbours
W_6nn <- as.matrix(make_knw(coords, k = 6, row_normalise = TRUE))

# W3: Inverse distance with 6nn cutoff
dmat <- as.matrix(dist(coords))
diag(dmat) <- 1  # avoid division by zero
W3_raw <- (1 / dmat) * (as.matrix(W_6nn) > 0)  # zero beyond 6nn
diag(W3_raw) <- 0
W_dist <- as.matrix(normw(W3_raw))

cat(sprintf("W dimensions: %d x %d\n", N, N))

# ============================================================
# 4. DGP matching MATLAB sar_conv_panel_gd.m
# ============================================================
set.seed(10203444)
Time <- 20
k <- 2
nobs <- N * Time
rho_true <- 0.6
beta_true <- c(1, 1)
sige_true <- 1
gamma_true <- c(0.3, 0.7, 0.0)

X <- matrix(rnorm(nobs * k), ncol = k)
evec <- rnorm(nobs) * sqrt(sige_true)

# Fixed effects
SFE <- rep((1:N) / N, times = Time)
TFE <- rep((1:Time) / Time, each = N)

# Build Wc and generate y
Wc <- gamma_true[1] * W_cont + gamma_true[2] * W_6nn + gamma_true[3] * W_dist
Wbig <- kronecker(diag(Time), Wc)

cat("Generating y from SAR DGP: y = (I - rho*Wc)^{-1} (X*beta + SFE + TFE + e)...\n")
A <- Matrix::Diagonal(nobs) - rho_true * Matrix::Matrix(Wbig, sparse = TRUE)
y <- as.numeric(Matrix::solve(A, X %*% beta_true + SFE + TFE + evec))

cat(sprintf("\nDGP: N=%d, T=%d, NT=%d\n", N, Time, nobs))
cat(sprintf("True: rho=%.2f, beta=(%.1f, %.1f), gamma=(%.1f, %.1f, %.1f)\n",
            rho_true, beta_true[1], beta_true[2],
            gamma_true[1], gamma_true[2], gamma_true[3]))

# ============================================================
# 5. Estimate with sar_conv_panel
# ============================================================
cat("\nEstimating SAR convex combination panel...\n")
cat("ndraw=20000, nomit=5000, thin=4\n\n")

t0 <- proc.time()
result <- sar_conv_panel(y, X, list(W_cont, W_6nn, W_dist), N, Time,
                         ndraw = 20000, nomit = 5000,
                         prior = list(model = 3, thin = 4))
elapsed <- (proc.time() - t0)[3]
cat(sprintf("Estimation time: %.1f seconds\n\n", elapsed))

# ============================================================
# 6. Print results and compare to manual
# ============================================================
cat("============================================================\n")
cat("ESTIMATION RESULTS vs MANUAL (Table 3.1)\n")
cat("============================================================\n\n")

# Manual values (from N=3111, ndraw=50000)
manual_beta1 <- 1.001; manual_beta2 <- 1.003
manual_rho <- 0.611
manual_g1 <- 0.324; manual_g2 <- 0.658; manual_g3 <- 0.018

cat(sprintf("%-15s  %-8s  %-10s  %-10s\n", "Parameter", "True", "R (N=500)", "MATLAB (N=3111)"))
cat(strrep("-", 55), "\n")
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10.4f\n", "beta_1", 1.0, result$beta[1], manual_beta1))
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10.4f\n", "beta_2", 1.0, result$beta[2], manual_beta2))
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10.4f\n", "rho", 0.6, result$rho, manual_rho))
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10.4f\n", "gamma_1", 0.3, result$gamma[1], manual_g1))
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10.4f\n", "gamma_2", 0.7, result$gamma[2], manual_g2))
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10.4f\n", "gamma_3", 0.0, result$gamma[3], manual_g3))
cat(sprintf("%-15s  %-8.3f  %-10.4f  %-10s\n", "sigma^2", 1.0, result$sige, "0.953"))

cat(sprintf("\nMH rho acceptance:   %.1f%%\n", result$rho_acc_rate * 100))
cat(sprintf("MH gamma acceptance: %.1f%%\n", result$gamma_acc_rate * 100))

# Credible intervals for gamma
gsave <- as.matrix(result$gdraw)
bsave <- as.matrix(result$bdraw)
psave <- as.numeric(result$pdraw)
ssave <- as.numeric(result$sdraw)

cat("\nGamma credible intervals:\n")
for (m in 1:3) {
  gl <- quantile(gsave[, m], 0.025)
  gu <- quantile(gsave[, m], 0.975)
  cat(sprintf("  gamma_%d: [%.4f, %.4f]  %s\n", m, gl, gu,
              if (gl > 0.01) "significant" else "includes ~0"))
}

# ============================================================
# 7. Figure 3.1: MCMC monitoring draws (last 1000 draws)
# ============================================================
matlab_blue   <- "#0072BD"
matlab_orange <- "#D95319"

n_draws <- nrow(gsave)
last_n <- min(1000, n_draws)
draw_range <- (n_draws - last_n + 1):n_draws
iter_x <- seq_len(last_n)

png(file.path(plot_dir, "figure_3_1.png"), width = 800, height = 600, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

# Top-left: gammas
plot(iter_x, gsave[draw_range, 1], type = "l", col = matlab_blue,
     ylim = c(0, 1), xlab = "draws", ylab = expression(gamma),
     main = "gammas")
lines(iter_x, gsave[draw_range, 2], col = matlab_orange)
lines(iter_x, gsave[draw_range, 3], col = "gray50")
abline(h = gamma_true[1], lty = 2, col = matlab_blue, lwd = 0.8)
abline(h = gamma_true[2], lty = 2, col = matlab_orange, lwd = 0.8)
legend("right",
       c(expression(gamma[1]*" (W"[cont]*")"),
         expression(gamma[2]*" (W"[6*nn]*")"),
         expression(gamma[3]*" (W"[dist]*")")),
       col = c(matlab_blue, matlab_orange, "gray50"),
       lty = 1, lwd = 1.5, cex = 0.65, bty = "n")

# Top-right: rho
plot(iter_x, psave[draw_range], type = "l", col = matlab_blue,
     xlab = "draws", ylab = expression(rho), main = "rho")
abline(h = rho_true, lty = 2, col = "red", lwd = 1)

# Bottom-left: beta
plot(iter_x, bsave[draw_range, 1], type = "l", col = matlab_blue,
     ylim = range(c(bsave[draw_range, 1], bsave[draw_range, 2])),
     xlab = "draws", ylab = expression(beta), main = "beta")
lines(iter_x, bsave[draw_range, 2], col = matlab_orange)
abline(h = beta_true[1], lty = 2, col = matlab_blue, lwd = 0.8)
abline(h = beta_true[2], lty = 2, col = matlab_orange, lwd = 0.8)
legend("topright", c(expression(beta[1]), expression(beta[2])),
       col = c(matlab_blue, matlab_orange), lty = 1, lwd = 1.5,
       cex = 0.7, bty = "n")

# Bottom-right: sigma
plot(iter_x, ssave[draw_range], type = "l", col = matlab_blue,
     xlab = "draws", ylab = expression(sigma^2), main = "sigma")
abline(h = sige_true, lty = 2, col = "red", lwd = 1)

mtext("Figure 3.1: MCMC monitoring draws (last 1000 iterations)",
      outer = TRUE, cex = 1, font = 2)

dev.off()
cat("\nSaved: figure_3_1.png\n")

# Effects
if (!is.null(result$direct)) {
  cat("\nEffects decomposition:\n")
  cat(sprintf("  %-5s  %-10s  %-10s  %-10s\n", "Var", "Direct", "Indirect", "Total"))
  for (j in seq_len(result$p)) {
    cat(sprintf("  x%-4d  %-10.4f  %-10.4f  %-10.4f\n",
                j, result$direct[j, 1], result$indirect[j, 1], result$total[j, 1]))
  }
}

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
