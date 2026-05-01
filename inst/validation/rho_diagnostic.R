# Diagnostic: Is the rho bias from W_cont approximation or convex estimation?

devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

# ============================================================
# 1. Read shapefile and subsample (same as figure_3_1.R)
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

# ============================================================
# 2. Build W matrices — kNN-12 proxy for contiguity
# ============================================================
W_cont_knn <- as.matrix(make_knw(coords, k = 12, row_normalise = TRUE))
W_6nn <- as.matrix(make_knw(coords, k = 6, row_normalise = TRUE))

# Inverse distance with 6nn cutoff
dmat <- as.matrix(dist(coords)); diag(dmat) <- 1
W3_raw <- (1 / dmat) * (W_6nn > 0); diag(W3_raw) <- 0
W_dist <- as.matrix(normw(W3_raw))

# ============================================================
# 3. Build Delaunay contiguity via spdep::tri2nb
# ============================================================
cat("Building Delaunay contiguity via spdep::tri2nb...\n")
nb_del <- spdep::tri2nb(coords)
listw_del <- spdep::nb2listw(nb_del, style = "W")
W_cont_del <- as.matrix(listw_del$weights |>
  (\(w) {
    n <- length(w)
    m <- matrix(0, n, n)
    for (i in seq_len(n)) {
      nbs <- nb_del[[i]]
      m[i, nbs] <- w[[i]]
    }
    m
  })())

cat(sprintf("Delaunay: mean neighbours = %.1f\n", mean(spdep::card(nb_del))))
cat(sprintf("kNN-12:   exactly 12 neighbours per unit\n\n"))

# ============================================================
# 4. Generate data using SAME DGP for all tests
#    Use the kNN-12 W_cont (same as figure_3_1.R)
# ============================================================
set.seed(10203444)
Time <- 20; k <- 2; nobs <- N * Time
rho_true <- 0.6; beta_true <- c(1, 1)
gamma_true <- c(0.3, 0.7, 0.0)

X <- matrix(rnorm(nobs * k), ncol = k)
evec <- rnorm(nobs)
SFE <- rep((1:N) / N, times = Time)
TFE <- rep((1:Time) / Time, each = N)

Wc_knn <- gamma_true[1] * W_cont_knn + gamma_true[2] * W_6nn + gamma_true[3] * W_dist
Wbig_knn <- Matrix::Matrix(kronecker(diag(Time), Wc_knn), sparse = TRUE)
A_knn <- Matrix::Diagonal(nobs) - rho_true * Wbig_knn
y_knn <- as.numeric(Matrix::solve(A_knn, X %*% beta_true + SFE + TFE + evec))

# ============================================================
# TEST A: Standard SAR with known Wc (kNN proxy)
# ============================================================
cat("TEST A: sar_panel() with known Wc (kNN-12 based), same y\n")
t0 <- proc.time()
res_A <- sar_panel(y_knn, X, Wc_knn, N, Time, ndraw = 10000, nomit = 3000,
                   prior = list(model = 3, rval = 0))
cat(sprintf("  rho = %.4f  (true = 0.6)  time = %.0fs\n\n",
            res_A$rho, (proc.time() - t0)[3]))

# ============================================================
# TEST B: Convex combination with kNN-12 proxy (same as figure_3_1.R)
# ============================================================
cat("TEST B: sar_conv_panel() with kNN-12 proxy for Wcont, same y\n")
t0 <- proc.time()
res_B <- sar_conv_panel(y_knn, X, list(W_cont_knn, W_6nn, W_dist), N, Time,
                        ndraw = 20000, nomit = 5000,
                        prior = list(model = 3, thin = 4))
cat(sprintf("  rho = %.4f  gamma = (%.3f, %.3f, %.3f)  time = %.0fs\n\n",
            res_B$rho, res_B$gamma[1], res_B$gamma[2], res_B$gamma[3],
            (proc.time() - t0)[3]))

# ============================================================
# TEST C: Generate NEW data using Delaunay Wcont, then estimate convex
# ============================================================
cat("TEST C: sar_conv_panel() with Delaunay Wcont — new y from Delaunay DGP\n")
Wc_del <- gamma_true[1] * W_cont_del + gamma_true[2] * W_6nn + gamma_true[3] * W_dist
Wbig_del <- Matrix::Matrix(kronecker(diag(Time), Wc_del), sparse = TRUE)
A_del <- Matrix::Diagonal(nobs) - rho_true * Wbig_del

set.seed(10203444)
X_c <- matrix(rnorm(nobs * k), ncol = k)
evec_c <- rnorm(nobs)
y_del <- as.numeric(Matrix::solve(A_del, X_c %*% beta_true + SFE + TFE + evec_c))

t0 <- proc.time()
res_C <- sar_conv_panel(y_del, X_c, list(W_cont_del, W_6nn, W_dist), N, Time,
                        ndraw = 20000, nomit = 5000,
                        prior = list(model = 3, thin = 4))
cat(sprintf("  rho = %.4f  gamma = (%.3f, %.3f, %.3f)  time = %.0fs\n\n",
            res_C$rho, res_C$gamma[1], res_C$gamma[2], res_C$gamma[3],
            (proc.time() - t0)[3]))

# ============================================================
# TEST D: Standard SAR with known Wc (Delaunay based)
# ============================================================
cat("TEST D: sar_panel() with known Wc (Delaunay based), Delaunay y\n")
t0 <- proc.time()
res_D <- sar_panel(y_del, X_c, Wc_del, N, Time, ndraw = 10000, nomit = 3000,
                   prior = list(model = 3, rval = 0))
cat(sprintf("  rho = %.4f  (true = 0.6)  time = %.0fs\n\n",
            res_D$rho, (proc.time() - t0)[3]))

# ============================================================
# Summary
# ============================================================
cat("============================================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("============================================================\n")
cat(sprintf("%-50s  rho = %.4f\n", "A: SAR with known Wc (kNN-12 proxy)", res_A$rho))
cat(sprintf("%-50s  rho = %.4f\n", "B: Conv with kNN-12 Wcont (same y as A)", res_B$rho))
cat(sprintf("%-50s  rho = %.4f\n", "C: Conv with Delaunay Wcont (Delaunay y)", res_C$rho))
cat(sprintf("%-50s  rho = %.4f\n", "D: SAR with known Wc (Delaunay, same y as C)", res_D$rho))
cat(strrep("-", 65), "\n")
cat("True rho = 0.6000\n\n")

if (abs(res_A$rho - 0.6) < 0.05) {
  cat("DIAGNOSIS: Test A close to truth => kNN-12 proxy is fine.\n")
  cat("  Bias in Test B is from convex combination estimation.\n")
} else {
  cat("DIAGNOSIS: Test A also biased => problem is the kNN-12 proxy,\n")
  cat("  not the convex combination estimator.\n")
}

if (abs(res_C$rho - 0.6) < abs(res_B$rho - 0.6)) {
  cat("  Delaunay (Test C) improves rho estimate over kNN proxy (Test B).\n")
} else {
  cat("  Delaunay (Test C) does NOT improve over kNN proxy (Test B).\n")
}
