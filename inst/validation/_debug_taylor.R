devtools::load_all("C:/Users/musta/Dropbox/Econometrics Papers/EAC Paper/spatial_econo_paper/functions/R_package/spmixW")

# The key question: in sar_conv_panel, the panel scaling is done as:
#   traces_panel$T2 <- Time * traces$T2
#   traces_panel$T3 <- Time * traces$T3
#   traces_panel$T4 <- Time * traces$T4
# But the new eval_taylor_lndet uses traces$traces[["2"]], NOT traces$T2
# So the panel scaling is NOT being applied to the actual computation!

# Check: are traces$T2 and traces$traces[["2"]] the same object?
N <- 20
W1 <- matrix(0, N, N); W2 <- matrix(0, N, N)
for (i in 1:N) { W1[i, (i %% N) + 1] <- 0.5; W2[i, ((i-2) %% N) + 1] <- 0.5 }

tr <- log_det_taylor(list(W1, W2), max_order = 4)

cat("T2 field:", head(tr$T2), "\n")
cat("traces[['2']] field:", head(tr$traces[["2"]]), "\n")
cat("Identical:", identical(tr$T2, tr$traces[["2"]]), "\n")

# Now simulate what sar_conv_panel does:
traces_panel <- tr
traces_panel$T2 <- 20 * tr$T2   # this modifies the backward-compat field
traces_panel$T3 <- 20 * tr$T3
traces_panel$T4 <- 20 * tr$T4

# But does eval_taylor_lndet use traces$traces or traces$T2?
# It uses traces$traces[["2"]] — which was NOT scaled!
cat("\nAfter panel scaling:\n")
cat("traces_panel$T2:", head(traces_panel$T2), "\n")
cat("traces_panel$traces[['2']]:", head(traces_panel$traces[["2"]]), "\n")
cat("Panel scaling applied to traces list:", identical(traces_panel$T2, traces_panel$traces[["2"]]), "\n")

cat("\n*** BUG CONFIRMED: panel scaling modifies $T2 but eval_taylor_lndet reads $traces[['2']] ***\n")
cat("The fix: sar_conv_panel must scale traces$traces, not traces$T2/T3/T4\n")
