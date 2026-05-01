# Taylor order sweep: test orders 4-8 on the same N=500 DGP
# Also run exact log-det as gold standard benchmark

devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

# ============================================================
# DGP (identical to rho_diagnostic.R Test B)
# ============================================================
shp_path <- "C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/toolbox_panelg/demo_data/uscounties_projected.shp"
counties <- sf::st_read(shp_path, quiet = TRUE)
coords_full <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(counties)))

set.seed(10203444)
N <- 500
idx <- sort(sample(nrow(coords_full), N))
latt <- coords_full[idx, 2]
long <- coords_full[idx, 1]
coords <- cbind(long, latt)

W_cont <- as.matrix(make_knw(coords, k = 12, row_normalise = TRUE))
W_6nn  <- as.matrix(make_knw(coords, k = 6, row_normalise = TRUE))
dmat <- as.matrix(dist(coords)); diag(dmat) <- 1
W3_raw <- (1 / dmat) * (W_6nn > 0); diag(W3_raw) <- 0
W_dist <- as.matrix(normw(W3_raw))

set.seed(10203444)
Time <- 20; k <- 2; nobs <- N * Time
rho_true <- 0.6; beta_true <- c(1, 1)
gamma_true <- c(0.3, 0.7, 0.0)

X <- matrix(rnorm(nobs * k), ncol = k)
evec <- rnorm(nobs)
SFE <- rep((1:N) / N, times = Time)
TFE <- rep((1:Time) / Time, each = N)

Wc <- gamma_true[1] * W_cont + gamma_true[2] * W_6nn + gamma_true[3] * W_dist
Wbig <- Matrix::Matrix(kronecker(diag(Time), Wc), sparse = TRUE)
A <- Matrix::Diagonal(nobs) - rho_true * Wbig
y <- as.numeric(Matrix::solve(A, X %*% beta_true + SFE + TFE + evec))

Wlist <- list(W_cont, W_6nn, W_dist)

# ============================================================
# First: verify the Taylor approximation quality at rho=0.6
# Compare Taylor log-det vs exact eigenvalue log-det
# ============================================================
cat("============================================================\n")
cat("Taylor approximation accuracy at rho=0.6, gamma=(0.3, 0.7, 0)\n")
cat("============================================================\n\n")

# Exact log-det via eigenvalues of Wc
eig_vals <- Re(eigen(Wc, only.values = TRUE)$values)
exact_lndet <- sum(log(1 - 0.6 * eig_vals))
cat(sprintf("Exact log|I - 0.6*Wc| = %.6f\n", exact_lndet))

for (ord in 4:8) {
  traces_ord <- log_det_taylor(Wlist, max_order = ord)
  taylor_val <- eval_taylor_lndet(traces_ord, rho = 0.6, gamma = c(0.3, 0.7, 0))
  err <- taylor_val - exact_lndet
  cat(sprintf("Order %d Taylor:          %.6f  (error = %+.4f, rel = %.2f%%)\n",
              ord, taylor_val, err, 100 * abs(err / exact_lndet)))
}

# ============================================================
# Sweep: sar_conv_panel at Taylor orders 4-8
# ============================================================
cat("\n============================================================\n")
cat("Taylor order sweep: sar_conv_panel estimates\n")
cat("============================================================\n\n")

orders <- 4:8
results_table <- data.frame(
  Order = integer(),
  rho = numeric(),
  bias = numeric(),
  gamma1 = numeric(),
  gamma2 = numeric(),
  gamma3 = numeric(),
  time_sec = numeric()
)

for (ord in orders) {
  cat(sprintf("Running order %d...\n", ord))
  t0 <- proc.time()
  res <- sar_conv_panel(y, X, Wlist, N, Time,
                        ndraw = 15000, nomit = 5000,
                        prior = list(model = 3, thin = 4, taylor_order = ord))
  elapsed <- (proc.time() - t0)[3]

  results_table <- rbind(results_table, data.frame(
    Order = ord,
    rho = res$rho,
    bias = res$rho - rho_true,
    gamma1 = res$gamma[1],
    gamma2 = res$gamma[2],
    gamma3 = res$gamma[3],
    time_sec = elapsed
  ))
  cat(sprintf("  rho=%.4f  bias=%+.4f  gamma=(%.3f, %.3f, %.3f)  time=%.0fs\n",
              res$rho, res$rho - rho_true,
              res$gamma[1], res$gamma[2], res$gamma[3], elapsed))
}

cat("\n============================================================\n")
cat("SWEEP RESULTS\n")
cat("============================================================\n")
cat(sprintf("%-7s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n",
            "Order", "rho", "bias", "gamma1", "gamma2", "gamma3", "time(s)"))
cat(strrep("-", 62), "\n")
for (i in seq_len(nrow(results_table))) {
  r <- results_table[i, ]
  cat(sprintf("%-7d  %-8.4f  %+.4f  %-8.4f  %-8.4f  %-8.4f  %-8.0f\n",
              r$Order, r$rho, r$bias, r$gamma1, r$gamma2, r$gamma3, r$time_sec))
}
cat(strrep("-", 62), "\n")
cat(sprintf("%-7s  %-8.4f  %+.4f\n", "True", rho_true, 0.0))
cat(sprintf("%-7s  %-8.4f  %+.4f\n", "SAR*", 0.5997, -0.0003))
cat("\n* SAR with known true Wc (from rho_diagnostic Test A)\n")
