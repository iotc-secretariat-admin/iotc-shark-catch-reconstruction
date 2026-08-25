# =============================================================================
# bsh_shark_nei_disaggregation.R
#
# Disaggregate generic shark retained catch to blue shark (BSH).
#
# PURPOSE
# -------
# This script addresses under-reporting to species:
#
#   Shark catch was reported
#       ↓
#   Catch was reported as a generic shark category
#       ↓
#   Estimate the proportion attributable to BSH
#
# This is distinct from complete non-reporting, which is handled by the
# ratio/GAM reconstruction.
#
# IMPORTANT
# ---------
# The script does NOT add estimates for strata with no shark reporting.
# It only partitions generic shark catch that already exists in RC.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# -----------------------------------------------------------------------------
# 1. Generic shark categories
#
# These are the categories currently being treated as potentially containing
# blue shark.
# -----------------------------------------------------------------------------

generic_shark_codes <- c(
  "SKH",    # Sharks various nei
  "SHRK",   # Sharks various nei
  "AG22",   # Sharks various nei
  "SKX",    # Sharks, rays, skates, etc. nei
  "SKHS"    # Sharks finned
)


# Not currently included.
# Retained here explicitly for later sensitivity analyses.
excluded_generic_shark_codes <- c(
  "MSK",    # Sharks mackerel and porbeagles nei
  "AG20",   # Sharks mackerel, porbeagles nei
  "SHXX"    # Sharks nei other than oceanic whitetip shark and blue shark
)


# =============================================================================
# 2. First make an inventory of ALL shark-related species codes
#
# Do this before estimating species composition.
#
# We do not want to assume that BSH and SMA are the only identified sharks.
# =============================================================================

shark_species_inventory <- rc_dat |>
  filter(
    str_detect(
      coalesce(SPECIES_CATEGORY, ""),
      regex("shark", ignore_case = TRUE)
    ) |
      str_detect(
        coalesce(SPECIES, ""),
        regex("shark", ignore_case = TRUE)
      ) |
      SPECIES_CODE %in% c(
        "BSH",
        "SMA",
        generic_shark_codes,
        excluded_generic_shark_codes
      )
  ) |>
  distinct(
    SPECIES_CATEGORY_CODE,
    SPECIES_CATEGORY,
    SPECIES_CODE,
    SPECIES,
    SPECIES_SCIENTIFIC
  ) |>
  arrange(
    SPECIES_CATEGORY_CODE,
    SPECIES_CODE
  )

print(shark_species_inventory, n = Inf)

# =============================================================================
# 3. Define reconstruction periods
#
# These periods are provisional and represent reporting regimes rather than
# biological periods.
# =============================================================================

rc_shark_dat <- rc_dat |>
  mutate(
    period = case_when(
      YEAR <= 1963 ~ "1950_1963",
      YEAR <= 1980 ~ "1964_1980",
      YEAR <= 1994 ~ "1981_1994",
      YEAR <= 2004 ~ "1995_2004",
      YEAR <= 2014 ~ "2005_2014",
      TRUE         ~ "2015_2025"
    )
  )


# =============================================================================
# 4. Identified shark species
#
# IMPORTANT:
# Replace this vector after reviewing shark_species_inventory.
#
# Include species-specific shark codes.
# Do NOT include generic shark categories.
# =============================================================================

identified_shark_codes <- c(
  "BSH",
  "SMA"

  # Add other species-specific shark codes after reviewing:
  # shark_species_inventory
)


# =============================================================================
# 5. Species-specific shark composition by Fleet x Fishery x Year
# =============================================================================

shark_composition_strata <- rc_shark_dat |>
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

    # Species-specific BSH
    bsh_mt = sum(
      CATCH[SPECIES_CODE == "BSH"],
      na.rm = TRUE
    ),

    # Species-specific SMA
    sma_mt = sum(
      CATCH[SPECIES_CODE == "SMA"],
      na.rm = TRUE
    ),

    # All identified shark species
    identified_shark_mt = sum(
      CATCH[SPECIES_CODE %in% identified_shark_codes],
      na.rm = TRUE
    ),

    # Generic sharks potentially requiring disaggregation
    shark_nei_mt = sum(
      CATCH[SPECIES_CODE %in% generic_shark_codes],
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>
  mutate(

    # BSH fraction among identified shark catch
    p_bsh = if_else(
      identified_shark_mt > 0,
      bsh_mt / identified_shark_mt,
      NA_real_
    ),

    p_sma = if_else(
      identified_shark_mt > 0,
      sma_mt / identified_shark_mt,
      NA_real_
    )
  )

# =============================================================================
# 6. Species-composition donors
# =============================================================================

bsh_comp_donors <- shark_composition_strata |>
  filter(
    identified_shark_mt > 0,
    !is.na(p_bsh)
  )



# Build the hierarchical BSH composition estimates
# =============================================================================
# 7. Level 1:
#    Fleet x Fishery x Period
# =============================================================================

comp_L1 <- bsh_comp_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_CODE,
    period
  ) |>
  summarise(

    donor_n = n(),

    n_donor_years =
      n_distinct(YEAR),

    identified_shark_mt =
      sum(identified_shark_mt),

    bsh_mt =
      sum(bsh_mt),

    # pooled composition
    p_bsh =
      bsh_mt / identified_shark_mt,

    .groups = "drop"
  ) |>
  filter(
    donor_n >= 2,
    identified_shark_mt > 0
  ) |>
  rename(
    p_bsh_L1 = p_bsh,
    n_L1 = donor_n
  )