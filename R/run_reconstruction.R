# =============================================================================
# IOTC BSH / SMA historical catch reconstruction
#
# Master production workflow
#
# Runs the minimum analytical sequence required to regenerate the
# final BSH and SMA reconstructed annual catch series.
#
# Diagnostic, exploratory, and manuscript plotting scripts are NOT
# included here unless they are required by the reconstruction.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Import and standardize IOTC data
# -----------------------------------------------------------------------------

source("R/processing/Initial_data_import.R")
source("R/processing/harmonize_iotc_codes.R")
source("R/processing/build_reconstruction_strata.R")


# =============================================================================
# Runtime options
# =============================================================================

run_cross_validation <- FALSE


# -----------------------------------------------------------------------------
# 2. Blue shark reconstruction
# -----------------------------------------------------------------------------

source("R/reconstruction/bsh_shark_nei_disaggregation.R")
source("R/reconstruction/bsh_nei_scenarios.R")
source("R/reconstruction/bsh_nei_historical_sensitivity.R")
source("R/reconstruction/bsh_ratio_reconstruction.R")


  cat(
    "\n\n=========================================================\n",
    " Running R/reconstruction/bsh_gam_species_composition.R",
    "\n=========================================================\n"
  )

source("R/reconstruction/bsh_gam_species_composition.R")
source("R/reconstruction/reconstruct_bsh_scenarios.R")

cat ("\n\n=========================================================\n",
"\nFinished  Blue shark reconstruction \n",
"\n\n=========================================================\n" )
# -----------------------------------------------------------------------------
# 3. Shortfin mako reconstruction
# -----------------------------------------------------------------------------

source("R/reconstruction/sma_reporting_states.R")
source("R/reconstruction/sma_hscr_species_composition.R")
source("R/reconstruction/sma_ratio_nonreporting.R")
source("R/reconstruction/sma_gam_species_composition.R")
source("R/reconstruction/reconstruct_sma_scenarios.R")

cat ("\n\n=========================================================\n",
   "\nFinished  mako shark reconstruction \n",
   "\n\n=========================================================\n")
# -----------------------------------------------------------------------------
# 4. Confirm final output exists
# -----------------------------------------------------------------------------

final_file <- "outputs/final_BSH_SMA_annual_reconstructed_catch.csv"

if (!file.exists(final_file)) {
  stop(
    "Pipeline completed, but final output was not created:\n",
    final_file
  )
}

message(
  "\nReconstruction complete.\n",
  "Final output: ", final_file
)