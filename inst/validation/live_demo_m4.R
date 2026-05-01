# Live demo: SAR convex combination with M=4 weight matrices
# Produces 6 diagnostic PNG plots and full console output

devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

plot_dir <- "C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW/inst/validation/plots"

# Colour palette for 4 weight matrices
cols4 <- c("#2166AC", "#B2182B", "#4DAF4A", "#984EA3")
W_names <- c("W1: 4-nearest neighbours",
             "W2: 6-nearest neighbours",
             "W3: 10-nearest neighbours",
             "W4: inverse-distance (6nn)")

# ============================================================
# DGP
# ============================================================
set.seed(42)
N <- 150
Time <- 15
k <- 2
nobs <- N * Time

latt <- runif(N)
long <- runif(N)

W1 <- as.matrix(make_knw(cbind(latt, long), k = 4))
W2 <- as.matrix(make_knw(cbind(latt, long), k = 6))
W3 <- as.matrix(make_knw(cbind(latt, long), k = 10))

# W4: inverse-distance with 6nn cutoff
dmat <- as.matrix(dist(cbind(latt, long))) + diag(N)
W4_raw <- (1 / dmat) * (W2 > 0)
W4 <- as.matrix(normw(W4_raw))

# True parameters
rho_true   <- 0.5
beta_true  <- c(1.5, -0.8)
gamma_true <- c(0.4, 0.35, 0.25, 0.0)

# Build true Wc
Wc <- gamma_true[1] * W1 + gamma_true[2] * W2 +
      gamma_true[3] * W3 + gamma_true[4] * W4
Wbig <- kronecker(diag(Time), Wc)

X <- matrix(rnorm(nobs * k), ncol = k)
SFE <- rep((1:N) / N, times = Time)
TFE <- rep((1:Time) / Time, each = N)
evec <- rnorm(nobs)

A <- diag(nobs) - rho_true * Wbig
y <- as.numeric(solve(A) %*% (X %*% beta_true + SFE + TFE + evec))

cat("============================================================\n")
cat("DGP: SAR convex combination panel\n")
cat(sprintf("N=%d, T=%d, NT=%d, k=%d, M=4\n", N, Time, nobs, k))
cat(sprintf("True rho = %.2f\n", rho_true))
cat(sprintf("True beta = (%.2f, %.2f)\n", beta_true[1], beta_true[2]))
cat(sprintf("True gamma = (%.2f, %.2f, %.2f, %.2f)\n",
            gamma_true[1], gamma_true[2], gamma_true[3], gamma_true[4]))
cat("============================================================\n\n")

# ============================================================
# Estimation
# ============================================================
cat("Running SAR convex combination panel estimation...\n")
cat("ndraw=15000, nomit=5000, thin=2\n\n")

t0 <- proc.time()
result <- sar_conv_panel(y, X, list(W1, W2, W3, W4), N, Time,
                         ndraw = 15000, nomit = 5000,
                         prior = list(model = 3, thin = 2))
elapsed <- (proc.time() - t0)[3]
cat(sprintf("Estimation time: %.1f seconds\n\n", elapsed))

# ============================================================
# Console output
# ============================================================
cat("============================================================\n")
cat("ESTIMATION RESULTS\n")
cat("============================================================\n\n")

# Coefficients
cat("Posterior coefficient estimates:\n")
bsave <- as.matrix(result$bdraw)
for (j in seq_along(result$beta)) {
  bm <- result$beta[j]
  bs <- sd(bsave[, j])
  bl <- quantile(bsave[, j], 0.025)
  bu <- quantile(bsave[, j], 0.975)
  cat(sprintf("  beta_%d:  Mean=%.4f  SD=%.4f  [%.4f, %.4f]\n", j, bm, bs, bl, bu))
}

cat(sprintf("\nrho:       Mean=%.4f  SD=%.4f  [%.4f, %.4f]\n",
            result$rho, sd(as.numeric(result$pdraw)),
            quantile(as.numeric(result$pdraw), 0.025),
            quantile(as.numeric(result$pdraw), 0.975)))
cat(sprintf("sigma^2:   Mean=%.4f\n", result$sige))

# Gamma
cat("\nConvex combination weights (gamma):\n")
gsave <- as.matrix(result$gdraw)
for (m in 1:4) {
  gm <- result$gamma[m]
  gs <- sd(gsave[, m])
  gl <- quantile(gsave[, m], 0.025)
  gu <- quantile(gsave[, m], 0.975)
  t_val <- gm / gs
  cat(sprintf("  gamma_%d: Mean=%.4f  SD=%.4f  t=%.2f  [%.4f, %.4f]  (true=%.2f)  %s\n",
              m, gm, gs, t_val, gl, gu, gamma_true[m], W_names[m]))
}

# Acceptance rates
cat(sprintf("\nMH acceptance rates:\n"))
cat(sprintf("  rho:   %.1f%%\n", result$rho_acc_rate * 100))
cat(sprintf("  gamma: %.1f%%\n", result$gamma_acc_rate * 100))

# Effects
cat("\nEffects decomposition:\n")
cat(sprintf("  %-5s  %-10s  %-10s  %-10s\n", "Var", "Direct", "Indirect", "Total"))
for (j in seq_len(result$p)) {
  cat(sprintf("  x%-4d  %-10.4f  %-10.4f  %-10.4f\n",
              j, result$direct[j, 1], result$indirect[j, 1], result$total[j, 1]))
}

# ============================================================
# PLOT 1: Gamma posteriors (2x2)
# ============================================================
png(file.path(plot_dir, "gamma_posteriors.png"), width = 800, height = 600, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (m in 1:4) {
  d <- density(gsave[, m], from = 0, to = 1)
  plot(d, main = W_names[m], xlab = bquote(gamma[.(m)]),
       ylab = "Density", col = cols4[m], lwd = 2.5, xlim = c(0, 1))
  polygon(d, col = adjustcolor(cols4[m], alpha.f = 0.15), border = NA)
  abline(v = gamma_true[m], col = "red", lty = 2, lwd = 2)
  abline(v = result$gamma[m], col = "gray30", lty = 3, lwd = 1.5)
  legend("topright",
         legend = c(sprintf("True = %.2f", gamma_true[m]),
                    sprintf("Est  = %.3f", result$gamma[m])),
         col = c("red", "gray30"), lty = c(2, 3), lwd = c(2, 1.5),
         cex = 0.7, bty = "n")
}
mtext("Posterior Densities of Convex Weights", outer = TRUE, cex = 1.1, font = 2)
dev.off()
cat("\nSaved: gamma_posteriors.png\n")

# ============================================================
# PLOT 2: Rho diagnostics (1x2)
# ============================================================
png(file.path(plot_dir, "rho_diagnostics.png"), width = 800, height = 400, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

rho_draws <- as.numeric(result$pdraw)
plot(rho_draws, type = "l", col = adjustcolor("gray40", 0.5),
     main = expression("Trace plot: " * rho),
     xlab = "Iteration (post burn-in)", ylab = expression(rho))
abline(h = rho_true, col = "red", lwd = 2, lty = 2)
abline(h = result$rho, col = cols4[1], lwd = 1.5, lty = 3)
legend("topright", c(sprintf("True = %.2f", rho_true), sprintf("Mean = %.3f", result$rho)),
       col = c("red", cols4[1]), lty = c(2, 3), lwd = c(2, 1.5), cex = 0.7, bty = "n")

d <- density(rho_draws)
plot(d, main = expression("Posterior density: " * rho),
     xlab = expression(rho), ylab = "Density", col = cols4[1], lwd = 2.5)
polygon(d, col = adjustcolor(cols4[1], alpha.f = 0.15), border = NA)
abline(v = rho_true, col = "red", lty = 2, lwd = 2)
abline(v = result$rho, col = "gray30", lty = 3, lwd = 1.5)
legend("topright", c(sprintf("True = %.2f", rho_true), sprintf("Mean = %.3f", result$rho)),
       col = c("red", "gray30"), lty = c(2, 3), lwd = c(2, 1.5), cex = 0.7, bty = "n")

dev.off()
cat("Saved: rho_diagnostics.png\n")

# ============================================================
# PLOT 3: Beta posteriors (1x2)
# ============================================================
png(file.path(plot_dir, "beta_posteriors.png"), width = 800, height = 400, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

for (j in 1:2) {
  d <- density(bsave[, j])
  plot(d, main = bquote("Posterior density: " * beta[.(j)]),
       xlab = bquote(beta[.(j)]), ylab = "Density",
       col = cols4[j], lwd = 2.5)
  polygon(d, col = adjustcolor(cols4[j], alpha.f = 0.15), border = NA)
  abline(v = beta_true[j], col = "red", lty = 2, lwd = 2)
  abline(v = result$beta[j], col = "gray30", lty = 3, lwd = 1.5)
  legend("topright",
         c(sprintf("True = %.2f", beta_true[j]),
           sprintf("Mean = %.4f", result$beta[j])),
         col = c("red", "gray30"), lty = c(2, 3), lwd = c(2, 1.5),
         cex = 0.7, bty = "n")
}
dev.off()
cat("Saved: beta_posteriors.png\n")

# ============================================================
# PLOT 4: Effects summary with error bars
# ============================================================
png(file.path(plot_dir, "effects_summary.png"), width = 800, height = 500, res = 120)
par(mar = c(5, 4, 3, 1))

p <- result$p
eff_means <- rbind(result$direct[, 1], result$indirect[, 1], result$total[, 1])

# 90% credible intervals
eff_lo <- eff_hi <- matrix(0, 3, p)
for (j in seq_len(p)) {
  eff_lo[1, j] <- quantile(result$direct_draws[, j], 0.05)
  eff_hi[1, j] <- quantile(result$direct_draws[, j], 0.95)
  eff_lo[2, j] <- quantile(result$indirect_draws[, j], 0.05)
  eff_hi[2, j] <- quantile(result$indirect_draws[, j], 0.95)
  eff_lo[3, j] <- quantile(result$total_draws[, j], 0.05)
  eff_hi[3, j] <- quantile(result$total_draws[, j], 0.95)
}

eff_cols <- c(cols4[1], cols4[3], "gray40")
bp <- barplot(eff_means, beside = TRUE,
              names.arg = paste0("x", seq_len(p)),
              col = eff_cols,
              main = "Effects Decomposition with 90% Credible Intervals",
              ylab = "Effect", ylim = c(min(eff_lo) * 1.15, max(eff_hi) * 1.15))

# Add error bars
for (j in seq_len(p)) {
  for (e in 1:3) {
    x_pos <- bp[e, j]
    arrows(x_pos, eff_lo[e, j], x_pos, eff_hi[e, j],
           angle = 90, code = 3, length = 0.04, lwd = 1.5)
  }
}

legend("topright", c("Direct", "Indirect", "Total"),
       fill = eff_cols, cex = 0.8, bty = "n")
dev.off()
cat("Saved: effects_summary.png\n")

# ============================================================
# PLOT 5: Gamma traces (all overlaid)
# ============================================================
png(file.path(plot_dir, "gamma_trace.png"), width = 800, height = 400, res = 120)
par(mar = c(4, 4, 3, 1))

n_draws <- nrow(gsave)
plot(gsave[, 1], type = "l", col = adjustcolor(cols4[1], 0.5),
     ylim = c(0, 1), main = "MCMC Traces: Convex Weights",
     xlab = "Iteration (post burn-in)", ylab = expression(gamma))
for (m in 2:4) {
  lines(gsave[, m], col = adjustcolor(cols4[m], 0.5))
}
# Add horizontal lines for true values
for (m in 1:4) {
  abline(h = gamma_true[m], col = cols4[m], lty = 2, lwd = 1.5)
}
legend("right", paste0("W", 1:4, " (true=", gamma_true, ")"),
       col = cols4, lty = 1, lwd = 2, cex = 0.65, bty = "n")

dev.off()
cat("Saved: gamma_trace.png\n")

# ============================================================
# PLOT 6: Acceptance rates summary
# ============================================================
png(file.path(plot_dir, "acceptance_rates.png"), width = 600, height = 400, res = 120)
par(mar = c(2, 2, 3, 2))

plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
title(main = "Metropolis-Hastings Acceptance Rates", cex.main = 1.2)

text(5, 8.5, sprintf("rho acceptance rate:    %.1f%%", result$rho_acc_rate * 100),
     cex = 1.3, font = 2, col = cols4[1])
text(5, 7, sprintf("gamma acceptance rate:  %.1f%%", result$gamma_acc_rate * 100),
     cex = 1.3, font = 2, col = cols4[3])

# Guideline
text(5, 5, "Target range: 20-60% for good mixing", cex = 1, col = "gray50")

# Status indicators
rho_ok <- result$rho_acc_rate > 0.15 && result$rho_acc_rate < 0.70
gam_ok <- result$gamma_acc_rate > 0.05 && result$gamma_acc_rate < 0.70
text(5, 3.5, sprintf("rho:   %s", if (rho_ok) "OK" else "NEEDS TUNING"),
     cex = 1.1, col = if (rho_ok) "#4DAF4A" else "#B2182B", font = 2)
text(5, 2.5, sprintf("gamma: %s", if (gam_ok) "OK" else "NEEDS TUNING"),
     cex = 1.1, col = if (gam_ok) "#4DAF4A" else "#B2182B", font = 2)

# Estimation metadata
text(5, 0.8, sprintf("N=%d  T=%d  M=%d  ndraw=%d  nomit=%d  thin=%d  time=%.0fs",
                      N, Time, 4, 15000, 5000, 2, elapsed),
     cex = 0.7, col = "gray60")

dev.off()
cat("Saved: acceptance_rates.png\n")

cat("\n============================================================\n")
cat("All 6 plots saved to:\n")
cat(plot_dir, "\n")
cat("============================================================\n")
