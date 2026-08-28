# =============================================================================
# bsh_gam_species_composition.R
#
# GAM-based reconstruction of blue shark species composition
#
# PURPOSE
# -------
# Estimate:
#
#   P(BSH | identified shark catch)
#
# for Fleet x Fishery x Year strata with incomplete species reporting.
#
# This provides a statistical alternative to the
# Hierarchical species-composition reconstruction (HSCR).
#
# The model is fitted only to strata containing species-specific shark
# information. It is then used to predict the BSH proportion in strata
# containing generic shark catch.
#
# IMPORTANT
# ---------
# This script addresses UNDER-REPORTING TO SPECIES.
#
# It does NOT estimate completely unreported shark catch.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(purrr)
library(tibble)


# =============================================================================
# 1. Species definitions
# =============================================================================

identified_shark_codes <- c(
  "ALV",  # Thresher
  "BSH",  # Blue shark
  "BTH",  # Bigeye thresher
  "FAL",  # Silky shark
  "LMA",  # Longfin mako
  "MAK",  # Mako sharks
  "OCS",  # Oceanic whitetip shark
  "POR",  # Porbeagle
  "PTH",  # Pelagic thresher
  "RHN",  # Whale shark
  "SMA",  # Shortfin mako
  "SPK",  # Great hammerhead
  "SPL",  # Scalloped hammerhead
  "SPN",  # Hammerhead sharks nei
  "SPY",  # Hammerhead sharks, etc. nei
  "SPZ",  # Smooth hammerhead
  "THR"   # Thresher sharks nei
)

generic_shark_codes <- c(
  "SKH"
)

# Other generic categories retained for later sensitivity testing
excluded_generic_shark_codes <- c(
  "MSK",
  "AG20",
  "SHXX",
  "SHRK",
  "AG22",
  "SKX",
  "SKHS"
)


# =============================================================================
# 2. Reporting periods
# =============================================================================

period_levels <- c(
  "1950_1963",
  "1964_1980",
  "1981_1994",
  "1995_2004",
  "2005_2014",
  "2015_2025"
)

add_reporting_period <- function(dat) {

  dat |>
    mutate(
      period = case_when(
        YEAR <= 1963 ~ "1950_1963",
        YEAR <= 1980 ~ "1964_1980",
        YEAR <= 1994 ~ "1981_1994",
        YEAR <= 2004 ~ "1995_2004",
        TRUE         ~ "2015_2025"
      ),

      period = factor(
        period,
        levels = period_levels,
        ordered = TRUE
      )
    )
}


# =============================================================================
# 3. Build Fleet x Fishery x Year shark-composition dataset
# =============================================================================

gam_shark_data <- rc_dat |>
  add_reporting_period() |>
  group_by(
    YEAR,
    period,

    FLEET_CODE,
    FLEET,

    FISHERY_GROUP_CODE,
    FISHERY_GROUP,

    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(

    # -------------------------------------------------------------
    # Species-specific sharks
    # -------------------------------------------------------------

    bsh_mt = sum(
      CATCH[SPECIES_CODE == "BSH"],
      na.rm = TRUE
    ),

    sma_mt = sum(
      CATCH[SPECIES_CODE == "SMA"],
      na.rm = TRUE
    ),

    identified_shark_mt = sum(
      CATCH[SPECIES_CODE %in% identified_shark_codes],
      na.rm = TRUE
    ),

    generic_shark_mt = sum(
      CATCH[SPECIES_CODE %in% generic_shark_codes],
      na.rm = TRUE
    ),

    # -------------------------------------------------------------
    # Target species
    # -------------------------------------------------------------

    yft_mt = sum(
      CATCH[SPECIES_CODE == "YFT"],
      na.rm = TRUE
    ),

    bet_mt = sum(
      CATCH[SPECIES_CODE == "BET"],
      na.rm = TRUE
    ),

    skj_mt = sum(
      CATCH[SPECIES_CODE == "SKJ"],
      na.rm = TRUE
    ),

    alb_mt = sum(
      CATCH[SPECIES_CODE == "ALB"],
      na.rm = TRUE
    ),

    swo_mt = sum(
      CATCH[SPECIES_CODE == "SWO"],
      na.rm = TRUE
    ),

    sbf_mt = sum(
      CATCH[SPECIES_CODE == "SBF"],
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# =============================================================================
# 4. Derive explanatory variables
# =============================================================================

gam_shark_data <- gam_shark_data |>
  mutate(

    target_primary_mt =
      yft_mt +
      bet_mt +
      skj_mt +
      alb_mt +
      swo_mt,

    target_plus_sbf_mt =
      target_primary_mt +
      sbf_mt,

    tropical_tuna_mt =
      yft_mt +
      bet_mt +
      skj_mt,

    temperate_target_mt =
      alb_mt +
      swo_mt,

    # -------------------------------------------------------------
    # Target composition
    # -------------------------------------------------------------

    prop_tropical =
      if_else(
        target_plus_sbf_mt > 0,
        tropical_tuna_mt / target_plus_sbf_mt,
        0
      ),

    prop_alb =
      if_else(
        target_plus_sbf_mt > 0,
        alb_mt / target_plus_sbf_mt,
        0
      ),

    prop_swo =
      if_else(
        target_plus_sbf_mt > 0,
        swo_mt / target_plus_sbf_mt,
        0
      ),

    prop_sbf =
      if_else(
        target_plus_sbf_mt > 0,
        sbf_mt / target_plus_sbf_mt,
        0
      ),

    # -------------------------------------------------------------
    # Response
    # -------------------------------------------------------------

    p_bsh =
      if_else(
        identified_shark_mt > 0,
        bsh_mt / identified_shark_mt,
        NA_real_
      )
  )