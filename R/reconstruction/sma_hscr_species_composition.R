#sma_hscr_species_composition.R

# =============================================================================
# sma_hscr_species_composition.R
#
# Hierarchical Species-Composition Reconstruction (HSCR) for shortfin mako.
#
# Two separate taxonomic allocation problems:
#
#   A. MAK -> SMA
#      MAK = Isurus spp.
#
#   B. SKH -> SMA
#      SKH = Various sharks NEI
#
# These allocations are estimated independently.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)


# =============================================================================
# 1. Definitions
# =============================================================================

identified_shark_codes <- c(
  "ALV",
  "BSH",
  "BTH",
  "FAL",
  "LMA",
  "OCS",
  "POR",
  "PTH",
  "RHN",
  "SMA",
  "SPK",
  "SPL",
  "SPN",
  "SPY",
  "SPZ",
  "THR"
)


# =============================================================================
# 2. Build donor composition dataset
# =============================================================================

sma_comp_strata <- rc_dat |>
  mutate(
    period = case_when(
      YEAR <= 1963 ~ "1950_1963",
      YEAR <= 1980 ~ "1964_1980",
      YEAR <= 1994 ~ "1981_1994",
      YEAR <= 2004 ~ "1995_2004",
      YEAR <= 2014 ~ "2005_2014",
      TRUE         ~ "2015_2025"
    )
  ) |>
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

    sma_mt =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    lma_mt =
      sum(
        CATCH[SPECIES_CODE == "LMA"],
        na.rm = TRUE
      ),

    mak_mt =
      sum(
        CATCH[SPECIES_CODE == "MAK"],
        na.rm = TRUE
      ),

    generic_shark_mt =
      sum(
        CATCH[SPECIES_CODE == "SKH"],
        na.rm = TRUE
      ),

    identified_shark_mt =
      sum(
        CATCH[
          SPECIES_CODE %in%
            identified_shark_codes
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(

    identified_mako_mt =
      sma_mt +
        lma_mt,

    p_sma_given_mako =
      if_else(
        identified_mako_mt > 0,
        sma_mt /
          identified_mako_mt,
        NA_real_
      ),

    p_sma_given_shark =
      if_else(
        identified_shark_mt > 0,
        sma_mt /
          identified_shark_mt,
        NA_real_
      )
  )


# =============================================================================
# PART A
# MAK -> SMA
# =============================================================================


# =============================================================================
# 3. MAK donor level 1:
# same fleet x fishery x period
# =============================================================================

mak_L1 <- sma_comp_strata |>
  filter(
    identified_mako_mt > 0
  ) |>
  group_by(
    FLEET_CODE,
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L1 =
      n(),

    donor_first_year_L1 =
      min(YEAR),

    donor_last_year_L1 =
      max(YEAR),

    donor_sma_mt_L1 =
      sum(sma_mt),

    donor_lma_mt_L1 =
      sum(lma_mt),

    donor_identified_mako_mt_L1 =
      donor_sma_mt_L1 +
        donor_lma_mt_L1,

    p_sma_mak_L1 =
      donor_sma_mt_L1 /
        donor_identified_mako_mt_L1,

    .groups = "drop"
  )


# =============================================================================
# 4. MAK donor level 2:
# same fleet x fishery group x period
# =============================================================================

mak_L2 <- sma_comp_strata |>
  filter(
    identified_mako_mt > 0
  ) |>
  group_by(
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L2 =
      n(),

    donor_first_year_L2 =
      min(YEAR),

    donor_last_year_L2 =
      max(YEAR),

    donor_sma_mt_L2 =
      sum(sma_mt),

    donor_lma_mt_L2 =
      sum(lma_mt),

    p_sma_mak_L2 =
      donor_sma_mt_L2 /
      (donor_sma_mt_L2 +
         donor_lma_mt_L2),

    .groups = "drop"
  )


# =============================================================================
# 5. MAK donor level 3:
# fishery x period
# =============================================================================

mak_L3 <- sma_comp_strata |>
  filter(
    identified_mako_mt > 0
  ) |>
  group_by(
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L3 =
      n(),

    n_fleets_L3 =
      n_distinct(
        FLEET_CODE
      ),

    donor_first_year_L3 =
      min(YEAR),

    donor_last_year_L3 =
      max(YEAR),

    donor_sma_mt_L3 =
      sum(sma_mt),

    donor_lma_mt_L3 =
      sum(lma_mt),

    p_sma_mak_L3 =
      donor_sma_mt_L3 /
      (donor_sma_mt_L3 +
         donor_lma_mt_L3),

    .groups = "drop"
  )


# =============================================================================
# 6. MAK donor level 4:
# fishery group x period
# =============================================================================

mak_L4 <- sma_comp_strata |>
  filter(
    identified_mako_mt > 0
  ) |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L4 =
      n(),

    n_fleets_L4 =
      n_distinct(
        FLEET_CODE
      ),

    donor_first_year_L4 =
      min(YEAR),

    donor_last_year_L4 =
      max(YEAR),

    donor_sma_mt_L4 =
      sum(sma_mt),

    donor_lma_mt_L4 =
      sum(lma_mt),

    p_sma_mak_L4 =
      donor_sma_mt_L4 /
      (donor_sma_mt_L4 +
         donor_lma_mt_L4),

    .groups = "drop"
  )


# =============================================================================
# 7. MAK allocation candidates
# =============================================================================

sma_mak_candidates <- sma_comp_strata |>
  filter(
    mak_mt > 0
  ) |>

  left_join(
    mak_L1,
    by = c(
      "FLEET_CODE",
      "FISHERY_CODE",
      "period"
    )
  ) |>

  left_join(
    mak_L2,
    by = c(
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>

  left_join(
    mak_L3,
    by = c(
      "FISHERY_CODE",
      "period"
    )
  ) |>

  left_join(
    mak_L4,
    by = c(
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>

  mutate(

    mak_source = case_when(

      !is.na(p_sma_mak_L1) ~
        "fleet_fishery_period",

      !is.na(p_sma_mak_L2) ~
        "fleet_group_period",

      !is.na(p_sma_mak_L3) &
        n_fleets_L3 >= 2 ~
        "multi_fleet_fishery_period",

      !is.na(p_sma_mak_L4) &
        n_fleets_L4 >= 2 ~
        "multi_fleet_group_period",

      TRUE ~
        "no_supported_donor"
    ),

    p_sma_from_mak =
      case_when(

        mak_source ==
          "fleet_fishery_period" ~
          p_sma_mak_L1,

        mak_source ==
          "fleet_group_period" ~
          p_sma_mak_L2,

        mak_source ==
          "multi_fleet_fishery_period" ~
          p_sma_mak_L3,

        mak_source ==
          "multi_fleet_group_period" ~
          p_sma_mak_L4,

        TRUE ~
          NA_real_
      ),

    sma_from_mak =
      mak_mt *
        p_sma_from_mak
  )

sma_mak_candidates <- sma_mak_candidates |>
  mutate(

    mak_allocated =
      !is.na(p_sma_from_mak),

    allocated_mak_mt =
      if_else(
        mak_allocated,
        mak_mt,
        0
      ),

    unallocated_mak_mt =
      if_else(
        mak_allocated,
        0,
        mak_mt
      )
  )

glimpse(sma_mak_candidates)
sma_mak_summary <- sma_mak_candidates |>
  summarise(

    mak_mt =
      sum(mak_mt, na.rm = TRUE),

    allocated_mak_mt =
      sum(allocated_mak_mt, na.rm = TRUE),

    unallocated_mak_mt =
      sum(unallocated_mak_mt, na.rm = TRUE),

    estimated_sma_from_mak =
      sum(sma_from_mak, na.rm = TRUE),

    pct_mak_supported =
      100 * allocated_mak_mt / mak_mt
  )
# =============================================================================
# PART B
# SKH -> SMA
# =============================================================================


# =============================================================================
# 8. SKH donor L1
# =============================================================================

skh_L1 <- sma_comp_strata |>
  filter(
    identified_shark_mt > 0
  ) |>
  group_by(
    FLEET_CODE,
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L1 =
      n(),

    donor_first_year_L1 =
      min(YEAR),

    donor_last_year_L1 =
      max(YEAR),

    donor_sma_mt_L1 =
      sum(sma_mt),

    donor_identified_shark_mt_L1 =
      sum(
        identified_shark_mt
      ),

    p_sma_skh_L1 =
      donor_sma_mt_L1 /
        donor_identified_shark_mt_L1,

    .groups = "drop"
  )


# =============================================================================
# 9. SKH donor L2
# =============================================================================

skh_L2 <- sma_comp_strata |>
  filter(
    identified_shark_mt > 0
  ) |>
  group_by(
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L2 =
      n(),

    donor_sma_mt_L2 =
      sum(sma_mt),

    donor_identified_shark_mt_L2 =
      sum(
        identified_shark_mt
      ),

    p_sma_skh_L2 =
      donor_sma_mt_L2 /
        donor_identified_shark_mt_L2,

    .groups = "drop"
  )


# =============================================================================
# 10. SKH donor L3
# =============================================================================

skh_L3 <- sma_comp_strata |>
  filter(
    identified_shark_mt > 0
  ) |>
  group_by(
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L3 =
      n(),

    n_fleets_L3 =
      n_distinct(
        FLEET_CODE
      ),

    donor_sma_mt_L3 =
      sum(sma_mt),

    donor_identified_shark_mt_L3 =
      sum(
        identified_shark_mt
      ),

    p_sma_skh_L3 =
      donor_sma_mt_L3 /
        donor_identified_shark_mt_L3,

    .groups = "drop"
  )


# =============================================================================
# 11. SKH candidates / N1 and N2
# =============================================================================

sma_skh_candidates <- sma_comp_strata |>
  filter(
    generic_shark_mt > 0
  ) |>

  left_join(
    skh_L1,
    by = c(
      "FLEET_CODE",
      "FISHERY_CODE",
      "period"
    )
  ) |>

  left_join(
    skh_L2,
    by = c(
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>

  left_join(
    skh_L3,
    by = c(
      "FISHERY_CODE",
      "period"
    )
  ) |>

  mutate(

    # -------------------------------------------------------------
    # N1: stricter empirical HSCR
    # -------------------------------------------------------------

    sma_N1_source = case_when(

      !is.na(p_sma_skh_L1) ~
        "same_fleet_same_fishery",

      !is.na(p_sma_skh_L2) ~
        "same_fleet_same_group",

      !is.na(p_sma_skh_L3) &
        n_fleets_L3 >= 2 ~
        "multi_fleet_same_fishery",

      TRUE ~
        "insufficient_support"
    ),

    p_sma_N1 = case_when(

      !is.na(p_sma_skh_L1) ~
        p_sma_skh_L1,

      !is.na(p_sma_skh_L2) ~
        p_sma_skh_L2,

      !is.na(p_sma_skh_L3) &
        n_fleets_L3 >= 2 ~
        p_sma_skh_L3,

      TRUE ~
        NA_real_
    ),

    sma_from_skh_N1 =
      generic_shark_mt *
        p_sma_N1,


    # -------------------------------------------------------------
    # N2: broader hierarchical borrowing
    # -------------------------------------------------------------

    sma_N2_source = case_when(

      !is.na(p_sma_skh_L1) ~
        "same_fleet_same_fishery",

      !is.na(p_sma_skh_L2) ~
        "same_fleet_same_group",

      !is.na(p_sma_skh_L3) ~
        "same_fishery",

      TRUE ~
        "no_donor"
    ),

    p_sma_N2 =
      coalesce(
        p_sma_skh_L1,
        p_sma_skh_L2,
        p_sma_skh_L3
      ),

    sma_from_skh_N2 =
      generic_shark_mt *
        p_sma_N2
  )


# =============================================================================
# 12. Summaries
# =============================================================================

sma_mak_candidates <- sma_mak_candidates |>
  mutate(

    mak_allocated =
      !is.na(p_sma_from_mak),

    allocated_mak_mt =
      if_else(
        mak_allocated,
        mak_mt,
        0
      ),

    unallocated_mak_mt =
      if_else(
        mak_allocated,
        0,
        mak_mt
      )
  )

sma_mak_summary <- sma_mak_candidates |>
  summarise(

    mak_mt =
      sum(mak_mt, na.rm = TRUE),

    allocated_mak_mt =
      sum(allocated_mak_mt, na.rm = TRUE),

    unallocated_mak_mt =
      sum(unallocated_mak_mt, na.rm = TRUE),

    estimated_sma_from_mak =
      sum(sma_from_mak, na.rm = TRUE),

    pct_mak_supported =
      100 * allocated_mak_mt / mak_mt
  )


sma_skh_candidates <- sma_skh_candidates |>
  mutate(

    skh_supported_N1 =
      !is.na(p_sma_N1),

    skh_supported_N2 =
      !is.na(p_sma_N2),

    allocated_skh_N1_mt =
      if_else(
        skh_supported_N1,
        generic_shark_mt,
        0
      ),

    allocated_skh_N2_mt =
      if_else(
        skh_supported_N2,
        generic_shark_mt,
        0
      ),

    unallocated_skh_N1_mt =
      if_else(
        skh_supported_N1,
        0,
        generic_shark_mt
      ),

    unallocated_skh_N2_mt =
      if_else(
        skh_supported_N2,
        0,
        generic_shark_mt
      )
  )

sma_skh_summary <- sma_skh_candidates |>
  summarise(

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    allocated_skh_N1_mt =
      sum(allocated_skh_N1_mt, na.rm = TRUE),

    allocated_skh_N2_mt =
      sum(allocated_skh_N2_mt, na.rm = TRUE),

    sma_N1_mt =
      sum(sma_from_skh_N1, na.rm = TRUE),

    sma_N2_mt =
      sum(sma_from_skh_N2, na.rm = TRUE),

    unallocated_N1_mt =
      sum(unallocated_skh_N1_mt, na.rm = TRUE),

    unallocated_N2_mt =
      sum(unallocated_skh_N2_mt, na.rm = TRUE),

    pct_skh_supported_N1 =
      100 * allocated_skh_N1_mt /
        generic_shark_mt,

    pct_skh_supported_N2 =
      100 * allocated_skh_N2_mt /
        generic_shark_mt
  )


print(sma_mak_summary)
print(sma_skh_summary)


# =============================================================================
# 13. Save
# =============================================================================

saveRDS(
  sma_mak_candidates,
  "outputs/sma/sma_mak_candidates.rds"
)

saveRDS(
  sma_skh_candidates,
  "outputs/sma/sma_skh_candidates.rds"
)

write_csv(
  sma_mak_summary,
  "outputs/sma/sma_mak_summary.csv"
)

write_csv(
  sma_skh_summary,
  "outputs/sma/sma_skh_summary.csv"
)