# =============================================================================
# sma_hscr_nonreporting.R
#
# HSCR reconstruction of completely unreported SMA catches.
#
# Applies only to strata with target fishing and no:
#
#   SMA
#   MAK
#   SKH
#
# Ratio modeled as:
#
#   SMA / target catch
#
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)


# =============================================================================
# 1. Donors
# =============================================================================

sma_ratio_donors <- sma_reporting_states |>
  filter(
    rc_sma > 0,
    rc_target > 0
  ) |>
  mutate(
    sma_target_ratio =
      rc_sma /
        rc_target
  )


# =============================================================================
# 2. Candidates
# =============================================================================

sma_ratio_candidates <- sma_reporting_states |>
  filter(
    sma_problem ==
      "potential_non_reporting",
    rc_target > 0
  )

# fix the donor class. “all longline” is still too broad a fallback donor pool. LLO LLD LLF are all wildly different
add_longline_class <- function(dat) {

  dat |>
    mutate(

      longline_class = case_when(

        FISHERY_GROUP_CODE != "LL" ~
          NA_character_,

        FISHERY_CODE == "LLD" ~
          "LL_deep_freezing",

        FISHERY_CODE == "LLF" ~
          "LL_fresh",

        FISHERY_CODE == "LLO" ~
          "LL_other",

        # If SLL occurs as a fishery code in the object
        FISHERY_CODE == "SLL" ~
          "LL_shark_target",

        TRUE ~
          "LL_other_unspecified"
      )
    )
}

sma_ratio_donors <-
  add_longline_class(
    sma_ratio_donors
  )

sma_ratio_candidates <-
  add_longline_class(
    sma_ratio_candidates
  )



# =============================================================================
# 3. L1: fleet x fishery x period
# =============================================================================

ratio_L1 <- sma_ratio_donors |>
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

    donor_target_mt_L1 =
      sum(rc_target),

    donor_sma_mt_L1 =
      sum(rc_sma),

    ratio_L1 =
      donor_sma_mt_L1 /
        donor_target_mt_L1,

    .groups = "drop"
  )


# =============================================================================
# 4. L2: fleet x same longline class x period
# =============================================================================

ratio_L2 <- sma_ratio_donors |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    !is.na(longline_class)
  ) |>
  group_by(
    FLEET_CODE,
    longline_class,
    period
  ) |>
  summarise(

    n_L2 =
      n(),

    donor_first_year_L2 =
      min(YEAR),

    donor_last_year_L2 =
      max(YEAR),

    donor_target_mt_L2 =
      sum(
        rc_target,
        na.rm = TRUE
      ),

    donor_sma_mt_L2 =
      sum(
        rc_sma,
        na.rm = TRUE
      ),

    ratio_L2 =
      donor_sma_mt_L2 /
        donor_target_mt_L2,

    .groups = "drop"
  )


# =============================================================================
# 5. L3: fishery x period
# =============================================================================

ratio_L3 <- sma_ratio_donors |>
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

    donor_target_mt_L3 =
      sum(rc_target),

    donor_sma_mt_L3 =
      sum(rc_sma),

    ratio_L3 =
      donor_sma_mt_L3 /
        donor_target_mt_L3,

    .groups = "drop"
  )


# =============================================================================
# 6. L4: all longline x period
# =============================================================================

ratio_L4 <- sma_ratio_donors |>
  filter(
    FISHERY_GROUP_CODE ==
      "LL"
  ) |>
  group_by(
    longline_class,
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

    donor_target_mt_L4 =
      sum(rc_target),

    donor_sma_mt_L4 =
      sum(rc_sma),

    ratio_L4 =
      donor_sma_mt_L4 /
        donor_target_mt_L4,

    .groups = "drop"
  )


# =============================================================================
# 7. Attach donor hierarchy and select reconstruction ratio
#
# IMPORTANT
# ---------
# Longline-wide donor levels (L2 and L4) may only be used for candidate
# strata that are themselves longline fisheries.
#
# This prevents longline donor ratios from being assigned to GN, LI, PS,
# or other fisheries simply because the fleet and period match.
#
# No coalesce() is used here because donor eligibility depends on the
# candidate fishery and must be evaluated explicitly.
# =============================================================================

sma_ratio_candidates <- sma_ratio_candidates |>

  # ---------------------------------------------------------------------------
  # L1: same fleet x same fishery x period
  # ---------------------------------------------------------------------------

  left_join(
    ratio_L1,
    by = c(
      "FLEET_CODE",
      "FISHERY_CODE",
      "period"
    )
  ) |>

  # ---------------------------------------------------------------------------
  # L2: same fleet x broader longline x period
  #
  # NOTE:
  # These donor values are joined to all rows but are only eligible below
  # when the candidate itself belongs to the same longline class i.e. LLO or LLF FISHERY_GROUP_CODE == "LL".
  # ---------------------------------------------------------------------------

left_join(
  ratio_L2,
  by = c(
    "FLEET_CODE",
    "longline_class",
    "period"
  )
)|>

  # ---------------------------------------------------------------------------
  # L3: same fishery x period, potentially borrowing across fleets
  # ---------------------------------------------------------------------------

  left_join(
    ratio_L3,
    by = c(
      "FISHERY_CODE",
      "period"
    )
  ) |>

  # ---------------------------------------------------------------------------
  # L4: broader longline x period, potentially borrowing across fleets
  #
  # Again, this level is only eligible for LL candidate strata.
  # ---------------------------------------------------------------------------

  left_join(
    ratio_L4,
    by = c(
          "longline_class",
          "period")
  ) |>

  mutate(

    # -------------------------------------------------------------------------
    # Selected donor source
    # -------------------------------------------------------------------------

    ratio_source = case_when(

      !is.na(ratio_L1) ~
        "fleet_fishery_period",

      FISHERY_GROUP_CODE == "LL" &
        !is.na(ratio_L2) ~
        "fleet_longline_period",

      !is.na(ratio_L3) ~
        "fishery_period",

      FISHERY_GROUP_CODE == "LL" &
        !is.na(ratio_L4) ~
        "longline_period",

      TRUE ~
        "no_donor"
    ),


    # -------------------------------------------------------------------------
    # Selected ratio
    #
    # Explicitly tied to ratio_source rather than simply taking the first
    # non-missing value.
    # -------------------------------------------------------------------------

    selected_ratio = case_when(

      ratio_source ==
        "fleet_fishery_period" ~
        ratio_L1,

      ratio_source ==
        "fleet_longline_period" ~
        ratio_L2,

      ratio_source ==
        "fishery_period" ~
        ratio_L3,

      ratio_source ==
        "longline_period" ~
        ratio_L4,

      TRUE ~
        NA_real_
    ),


    # -------------------------------------------------------------------------
    # Evidential support
    # -------------------------------------------------------------------------

    ratio_support = case_when(

      ratio_source ==
        "fleet_fishery_period" ~
        "same_fleet_same_fishery",

      ratio_source ==
        "fleet_longline_period" ~
        "same_fleet_broader_longline",

      ratio_source ==
        "fishery_period" &
        !is.na(n_fleets_L3) &
        n_fleets_L3 == 1 ~
        "single_fleet_proxy",

      ratio_source ==
        "fishery_period" &
        !is.na(n_fleets_L3) &
        n_fleets_L3 >= 2 ~
        "multi_fleet_same_fishery",

      ratio_source ==
        "longline_period" &
        !is.na(n_fleets_L4) &
        n_fleets_L4 == 1 ~
        "single_fleet_proxy",

      ratio_source ==
        "longline_period" &
        !is.na(n_fleets_L4) &
        n_fleets_L4 >= 2 ~
        "multi_fleet_longline",

      ratio_source ==
        "no_donor" ~
        "none",

      TRUE ~
        "other"
    ),


    # -------------------------------------------------------------------------
    # First year represented by selected donor
    # -------------------------------------------------------------------------

    selected_donor_first_year = case_when(

      ratio_source ==
        "fleet_fishery_period" ~
        donor_first_year_L1,

      ratio_source ==
        "fleet_longline_period" ~
        donor_first_year_L2,

      ratio_source ==
        "fishery_period" ~
        donor_first_year_L3,

      ratio_source ==
        "longline_period" ~
        donor_first_year_L4,

      TRUE ~
        NA_real_
    ),


    # -------------------------------------------------------------------------
    # Last year represented by selected donor
    # -------------------------------------------------------------------------

    selected_donor_last_year = case_when(

      ratio_source ==
        "fleet_fishery_period" ~
        donor_last_year_L1,

      ratio_source ==
        "fleet_longline_period" ~
        donor_last_year_L2,

      ratio_source ==
        "fishery_period" ~
        donor_last_year_L3,

      ratio_source ==
        "longline_period" ~
        donor_last_year_L4,

      TRUE ~
        NA_real_
    ),


    # -------------------------------------------------------------------------
    # Temporal support
    # -------------------------------------------------------------------------

    temporal_support = case_when(

      ratio_source ==
        "no_donor" ~
        "no_donor",

      is.na(selected_ratio) ~
        "no_donor",

      YEAR <
        selected_donor_first_year ~
        "backcast",

      YEAR >
        selected_donor_last_year ~
        "forward_extrapolation",

      TRUE ~
        "within_donor_period"
    ),


    # -------------------------------------------------------------------------
    # Complete non-reporting reconstruction
    # -------------------------------------------------------------------------

    sma_nonreport_estimate = case_when(

      !is.na(selected_ratio) &
        rc_target > 0 ~
        rc_target *
          selected_ratio,

      TRUE ~
        NA_real_
    )
  )

# =============================================================================
# 8. Summary & diagnostics 
# =============================================================================


# =============================================================================
# Inspect extreme SMA donor ratios
# =============================================================================

sma_ratio_donors |>
  arrange(
    desc(sma_target_ratio)
  ) |>
  select(
    YEAR,
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    FISHERY,
    rc_target,
    rc_sma,
    sma_target_ratio
  ) |>
  print(n = 100)

sma_ratio_donors |>
  group_by(
    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(

    n = n(),

    target_mt =
      sum(rc_target),

    sma_mt =
      sum(rc_sma),

    pooled_ratio =
      sma_mt / target_mt,

    median_annual_ratio =
      median(
        sma_target_ratio,
        na.rm = TRUE
      ),

    max_annual_ratio =
      max(
        sma_target_ratio,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  arrange(
    desc(pooled_ratio)
  ) |>
  print(n = Inf)
# 
# check the pooled ratios across the types of LL for the different periods. ? any big changes? 
sma_ratio_donors |>
  filter(
    FISHERY_GROUP_CODE == "LL"
  ) |>
  group_by(
    longline_class,
    period
  ) |>
  summarise(

    n_strata =
      n(),

    n_fleets =
      n_distinct(FLEET_CODE),

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

    pooled_ratio =
      sma_mt /
        target_mt,

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    .groups = "drop"
  ) |>
  print(n = Inf)

# =======================================================
#  Make a summary 
# ====================================================== 


sma_nonreport_summary <- sma_ratio_candidates |>
  group_by(
    ratio_support,
    temporal_support
  ) |>
  summarise(

    n_strata =
      n(),

    target_mt =
      sum(rc_target),

    reconstructed_sma_mt =
      sum(
        sma_nonreport_estimate,
        na.rm = TRUE
      ),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    .groups = "drop"
  ) |>
  arrange(
    desc(
      reconstructed_sma_mt
    )
  )

print(
  sma_nonreport_summary,
  n = Inf
)


# =============================================================================
# 9. Save
# =============================================================================

saveRDS(
  sma_ratio_candidates,
  "outputs/sma/sma_nonreport_candidates.rds"
)

write_csv(
  sma_nonreport_summary,
  "outputs/sma/sma_nonreport_summary.csv"
)



# check and re-check. 
sma_nonreport_summary

sma_reconstruction_summary # this comes from R\reconstruction\reconstruct_sma_scenarios.R

sma_scenario_base |>
  count(invalid_nonreport_overlap)

sma_ratio_candidates |>
  count(
    ratio_source,
    ratio_support,
    temporal_support,
    sort = TRUE
  )

