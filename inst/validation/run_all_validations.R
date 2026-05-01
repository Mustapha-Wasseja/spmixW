# run_all_validations.R — Master script to run all validation checks
# Usage: Rscript run_all_validations.R

cat("=================================================================\n")
cat("spmixW Package Validation Suite\n")
cat("Reproducing key examples from LeSage panelg_manual.pdf\n")
cat("=================================================================\n\n")

scripts <- c(
  "validation_01_ols.R",
  "validation_02_sar.R",
  "validation_03_sdm.R",
  "validation_04_sem.R",
  "validation_05_sdem.R",
  "validation_06_slx.R",
  "validation_07_sar_conv.R",
  "validation_08_lmarginal.R",
  "validation_09_bma.R",
  "validation_10_cross_section.R"
)

all_results <- list()
for (script in scripts) {
  cat(sprintf("\n>>> Running %s ...\n", script))
  tryCatch({
    source(script, local = TRUE)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  })
}

cat("\n=================================================================\n")
cat("VALIDATION SUITE COMPLETE\n")
cat("=================================================================\n")
