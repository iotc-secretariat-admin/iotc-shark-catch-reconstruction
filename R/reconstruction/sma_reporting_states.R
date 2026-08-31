# sma_reporting_states.R

# =============================================================================
# sma_reporting_states.R
#
# Construct Fleet x Fishery x Year reporting states for shortfin mako.
#
# Reporting problems:
#
#   1. SMA reported directly
#   2. MAK reported (Isurus spp.; requires mako-species allocation)
#   3. SKH reported (generic sharks; requires shark-composition allocation)
#   4. No shark reporting despite target fishing
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

target_species_codes <- c(
  "YFT",
  "BET",
  "SKJ",
  "ALB",
  "SWO",
  "SBF"
)

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
        YEAR <= 2014 ~ "2005_2014",
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
# 2. Aggregate RC to Fleet x Fishery x Year
# =============================================================================

sma_reporting_states <- rc_dat |>
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

    rc_total =
      sum(CATCH, na.rm = TRUE),

    rc_target =
      sum(
        CATCH[
          SPECIES_CODE %in%
            target_species_codes
        ],
        na.rm = TRUE
      ),

    rc_yft =
      sum(
        CATCH[SPECIES_CODE == "YFT"],
        na.rm = TRUE
      ),

    rc_bet =
      sum(
        CATCH[SPECIES_CODE == "BET"],
        na.rm = TRUE
      ),

    rc_skj =
      sum(
        CATCH[SPECIES_CODE == "SKJ"],
        na.rm = TRUE
      ),

    rc_alb =
      sum(
        CATCH[SPECIES_CODE == "ALB"],
        na.rm = TRUE
      ),

    rc_swo =
      sum(
        CATCH[SPECIES_CODE == "SWO"],
        na.rm = TRUE
      ),

    rc_sbf =
      sum(
        CATCH[SPECIES_CODE == "SBF"],
        na.rm = TRUE
      ),

    rc_sma =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    rc_lma =
      sum(
        CATCH[SPECIES_CODE == "LMA"],
        na.rm = TRUE
      ),

    rc_mak =
      sum(
        CATCH[SPECIES_CODE == "MAK"],
        na.rm = TRUE
      ),

    rc_shark_nei =
      sum(
        CATCH[SPECIES_CODE == "SKH"],
        na.rm = TRUE
      ),

    rc_identified_sharks =
      sum(
        CATCH[
          SPECIES_CODE %in%
            identified_shark_codes
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  )


# =============================================================================
# 3. Add CA SMA information
# =============================================================================

ca_sma_fy <- ca_dat |>
  filter(
    SPECIES_CODE == "SMA"
  ) |>
  group_by(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    ca_sma_mt =
      sum(
        CATCH[
          CATCH_UNIT_CODE == "MT"
        ],
        na.rm = TRUE
      ),

    ca_sma_no =
      sum(
        CATCH[
          CATCH_UNIT_CODE == "NO"
        ],
        na.rm = TRUE
      ),

    ca_sma_records =
      n(),

    .groups = "drop"
  )


sma_reporting_states <- sma_reporting_states |>
  left_join(
    ca_sma_fy,
    by = c(
      "YEAR",
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "FISHERY_CODE"
    )
  )


# =============================================================================
# 4. Add effort information
# =============================================================================

ef_fy <- ef_dat |>
  group_by(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    effort_units =
      paste(
        sort(
          unique(
            EFFORT_UNIT_CODE
          )
        ),
        collapse = ";"
      ),

    n_effort_units =
      n_distinct(
        EFFORT_UNIT_CODE
      ),

    effort_records =
      n(),

    has_hooks =
      any(
        EFFORT_UNIT_CODE ==
          "HOOKS"
      ),

    has_trips =
      any(
        EFFORT_UNIT_CODE ==
          "TRIPS"
      ),

    has_fdays =
      any(
        EFFORT_UNIT_CODE ==
          "FDAYS"
      ),

    .groups = "drop"
  )


sma_reporting_states <- sma_reporting_states |>
  left_join(
    ef_fy,
    by = c(
      "YEAR",
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "FISHERY_CODE"
    )
  )


# =============================================================================
# 5. Reporting flags
# =============================================================================

sma_reporting_states <- sma_reporting_states |>
  mutate(

    has_target =
      rc_target > 0,

    has_rc_sma =
      rc_sma > 0,

    has_rc_lma =
      rc_lma > 0,

    has_mak =
      rc_mak > 0,

    has_shark_nei =
      rc_shark_nei > 0,

    has_identified_sharks =
      rc_identified_sharks > 0,

    has_ca_sma =
      coalesce(
        ca_sma_mt > 0 |
          ca_sma_no > 0,
        FALSE
      ),

    has_effort =
      coalesce(
        effort_records > 0,
        FALSE
      )
  )


# =============================================================================
# 6. Primary SMA reporting state
# =============================================================================

sma_reporting_states <- sma_reporting_states |>
  mutate(

    sma_reporting_state = case_when(

      # -----------------------------------------------------------
      # Direct SMA reporting
      # -----------------------------------------------------------

      has_rc_sma ~
        "S1_SMA_reported_RC",

      # -----------------------------------------------------------
      # Genus-level mako reporting
      # -----------------------------------------------------------

      !has_rc_sma &
        has_mak ~
        "S2_MAK_reported_no_SMA",

      # -----------------------------------------------------------
      # Generic sharks reported
      # -----------------------------------------------------------

      !has_rc_sma &
        !has_mak &
        has_shark_nei ~
        "S3_generic_shark_reported",

      # -----------------------------------------------------------
      # Target fishing but no shark reporting
      # -----------------------------------------------------------

      !has_rc_sma &
        !has_mak &
        !has_shark_nei &
        has_target ~
        "S4_target_fishing_no_shark_reporting",

      # -----------------------------------------------------------
      # CA evidence despite no RC SMA
      # -----------------------------------------------------------

      !has_rc_sma &
        has_ca_sma ~
        "S5_CA_SMA_no_RC_SMA",

      TRUE ~
        "S6_other_no_SMA_information"
    )
  )



sma_reporting_states <- sma_reporting_states |>
  mutate(

    requires_mak_allocation =
      rc_mak > 0,

    requires_skh_allocation =
      rc_shark_nei > 0,

    requires_nonreport_reconstruction =
      rc_target > 0 &
      rc_sma == 0 &
      rc_mak == 0 &
      rc_shark_nei == 0
  )
# =============================================================================
# 7. Reconstruction problem
# =============================================================================

sma_reporting_states <- sma_reporting_states |>
  mutate(

    sma_problem = case_when(

      sma_reporting_state ==
        "S1_SMA_reported_RC" ~
        "reported",

      sma_reporting_state ==
        "S2_MAK_reported_no_SMA" ~
        "mako_species_allocation",

      sma_reporting_state ==
        "S3_generic_shark_reported" ~
        "generic_shark_allocation",

      sma_reporting_state ==
        "S4_target_fishing_no_shark_reporting" ~
        "potential_non_reporting",

      sma_reporting_state ==
        "S5_CA_SMA_no_RC_SMA" ~
        "RC_CA_reporting_conflict",

      TRUE ~
        "insufficient_information"
    )
  )


# =============================================================================
# 8. Summary
# =============================================================================

sma_state_summary <- sma_reporting_states |>
  group_by(
    sma_reporting_state
  ) |>
  summarise(

    n_strata =
      n(),

    total_target_mt =
      sum(
        rc_target,
        na.rm = TRUE
      ),

    total_sma_mt =
      sum(
        rc_sma,
        na.rm = TRUE
      ),

    total_mak_mt =
      sum(
        rc_mak,
        na.rm = TRUE
      ),

    total_shark_nei_mt =
      sum(
        rc_shark_nei,
        na.rm = TRUE
      ),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    .groups = "drop"
  )

print(
  sma_state_summary,
  n = Inf
)


# =============================================================================
# 9. Early longline diagnostics
# =============================================================================

sma_early_ll_summary <- sma_reporting_states |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    YEAR <= 1980
  ) |>
  group_by(
    sma_reporting_state
  ) |>
  summarise(

    n_strata =
      n(),

    target_mt =
      sum(
        rc_target,
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        rc_sma,
        na.rm = TRUE
      ),

    mak_mt =
      sum(
        rc_mak,
        na.rm = TRUE
      ),

    shark_nei_mt =
      sum(
        rc_shark_nei,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  arrange(
    desc(target_mt)
  )

print(
  sma_early_ll_summary,
  n = Inf
)


# =============================================================================
# 10. Fleet-level early longline summary
# =============================================================================

sma_early_ll_fleets <- sma_reporting_states |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    YEAR <= 1980
  ) |>
  group_by(
    FLEET_CODE,
    FLEET
  ) |>
  summarise(

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    target_mt =
      sum(rc_target),

    sma_mt =
      sum(rc_sma),

    mak_mt =
      sum(rc_mak),

    shark_nei_mt =
      sum(rc_shark_nei),

    n_years =
      n_distinct(YEAR),

    n_years_sma =
      n_distinct(
        YEAR[rc_sma > 0]
      ),

    .groups = "drop"
  ) |>
  arrange(
    desc(target_mt)
  )

print(
  sma_early_ll_fleets,
  n = Inf
)


# =============================================================================
# 11. Save
# =============================================================================

dir.create(
  "outputs/sma",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  sma_reporting_states,
  "outputs/sma/sma_reporting_states.rds"
)

write_csv(
  sma_state_summary,
  "outputs/sma/sma_reporting_state_summary.csv"
)

write_csv(
  sma_early_ll_summary,
  "outputs/sma/sma_early_ll_summary.csv"
)

write_csv(
  sma_early_ll_fleets,
  "outputs/sma/sma_early_ll_fleets.csv"
)