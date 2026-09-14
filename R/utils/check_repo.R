# =============================================================================
# Repository reproducibility checks
# IOTC BSH/SMA Catch Reconstruction
# =============================================================================

required_dirs <- c(
  "R",
  "R/processing",
  "R/reconstruction",
  "R/diagnostics",
  "R/plotting",
  "docs",
  "metadata",
  "outputs"
)

required_files <- c(
  "README.md",
  ".gitignore",
  "docs/Reconstructing_hist_BSH_SMA_catches_IOTC_updated.qmd",
  "docs/refs.bib"
)

cat("\n--- Directory check ---\n")

for (x in required_dirs) {
  status <- if (dir.exists(x)) "OK" else "MISSING"
  cat(sprintf("%-45s %s\n", x, status))
}

cat("\n--- File check ---\n")

for (x in required_files) {
  status <- if (file.exists(x)) "OK" else "MISSING"
  cat(sprintf("%-70s %s\n", x, status))
}


# -----------------------------------------------------------------------------
# Key analytical outputs
# -----------------------------------------------------------------------------

key_outputs <- c(
  "outputs/final_BSH_SMA_annual_reconstructed_catch.csv",
  "outputs/sma/sma_scenario_base.rds"
)

cat("\n--- Key analytical outputs ---\n")

for (x in key_outputs) {
  status <- if (file.exists(x)) "OK" else "MISSING"
  cat(sprintf("%-70s %s\n", x, status))
}


# -----------------------------------------------------------------------------
# Look for common temporary files that should not be committed
# -----------------------------------------------------------------------------

all_files <- list.files(
  ".",
  recursive = TRUE,
  all.files = TRUE,
  full.names = TRUE
)

temp_files <- all_files[
  grepl(
    "\\.vdoc|\\.Rhistory$|\\.RData$|\\.knit\\.md$|\\.utf8\\.md$|~$",
    all_files,
    ignore.case = TRUE
  )
]

cat("\n--- Temporary files found ---\n")

if (length(temp_files) == 0) {
  cat("None\n")
} else {
  cat(paste(temp_files, collapse = "\n"))
  cat("\n")
}


# -----------------------------------------------------------------------------
# Session information
# -----------------------------------------------------------------------------

cat("\n--- R session ---\n")
print(sessionInfo())