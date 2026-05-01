devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

# Same DGP as sweep
shp_path <- "C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/toolbox_panelg/demo_data/uscounties_projected.shp"
counties <- sf::st_read(shp_path, quiet = TRUE)
coords_full <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(counties)))
set.seed(10203444)
N <- 500
idx <- sort(sample(nrow(coords_full), N))
coords <- cbind(coords_full[idx, 1], coords_full[idx, 2])

W_cont <- as.matrix(make_knw(coords, k = 12, row_normalise = TRUE))
W_6nn  <- as.matrix(make_knw(coords, k = 6, row_normalise = TRUE))
dmat <- as.matrix(dist(coords)); diag(dmat) <- 1
W3_raw <- (1 / dmat) * (W_6nn > 0); diag(W3_raw) <- 0
W_dist <- as.matrix(normw(W3_raw))

set.seed(10203444)
Time <- 20; nobs <- N * Time
X <- matrix(rnorm(nobs * 2), ncol = 2)
evec <- rnorm(nobs)
SFE <- rep((1:N) / N, times = Time)
TFE <- rep((1:Time) / Time, each = N)
Wc <- 0.3 * W_cont + 0.7 * W_6nn
Wbig <- Matrix::Matrix(kronecker(diag(Time), Wc), sparse = TRUE)
A <- Matrix::Diagonal(nobs) - 0.6 * Wbig
y <- as.numeric(Matrix::solve(A, X %*% c(1, 1) + SFE + TFE + evec))

Wlist <- list(W_cont, W_6nn, W_dist)

for (ord in c(4, 6, 8)) {
  cat(sprintf("Order %d: ", ord))
  t0 <- proc.time()
  res <- sar_conv_panel(y, X, Wlist, N, Time,
                        ndraw = 15000, nomit = 5000,
                        prior = list(model = 3, thin = 4, taylor_order = ord))
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("rho=%.4f (bias=%+.4f)  gamma=(%.3f, %.3f, %.3f)  time=%.0fs\n",
              res$rho, res$rho - 0.6,
              res$gamma[1], res$gamma[2], res$gamma[3], elapsed))
}
