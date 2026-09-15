# =============================================================================
# bsh_ratio_reconstruction.R
# Blue shark ratio-based reconstruction
#
# Estimate unreported BSH catch using BSH : target-catch ratios from
# comparable reporting strata.
#
# Input:
#   reconstruction_strata
#
# Output:
#   bsh_ratio_results
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)

# -----------------------------------------------------------------------------
# 1. Define analytical periods
# -----------------------------------------------------------------------------

bsh_ratio_dat <- reconstruction_strata |>
  mutate(
    period = case_when(
      YEAR <= 1969 ~ "1950_1969",
      YEAR <= 1980 ~ "1970_1980",
      YEAR <= 1994 ~ "1981_1994",
      YEAR <= 2004 ~ "1995_2004",
      YEAR <= 2014 ~ "2005_2014",
      TRUE         ~ "2015_2025"
    )
  )


# -----------------------------------------------------------------------------
# 2. Define donor strata
#
# A donor must:
#   - report positive BSH catch
#   - contain positive target catch
#
# For the initial analysis, restrict to longline.
# -----------------------------------------------------------------------------

bsh_donors <- bsh_ratio_dat |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    rc_bsh > 0,
    rc_target > 0
  ) |>
  mutate(
    ratio_bsh_target = rc_bsh / rc_target
  )


# -----------------------------------------------------------------------------
# 3. Inspect ratio distribution
# -----------------------------------------------------------------------------

bsh_ratio_summary <- bsh_donors |>
  group_by(period) |>
  summarise(
    n_strata = n(),
    n_fleets = n_distinct(FLEET_CODE),

    median_ratio = median(
      ratio_bsh_target,
      na.rm = TRUE
    ),

    mean_ratio = mean(
      ratio_bsh_target,
      na.rm = TRUE
    ),

    trimmed_mean_ratio = mean(
      ratio_bsh_target,
      trim = 0.10,
      na.rm = TRUE
    ),

    q25_ratio = quantile(
      ratio_bsh_target,
      0.25,
      na.rm = TRUE
    ),

    q75_ratio = quantile(
      ratio_bsh_target,
      0.75,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

print(bsh_ratio_summary)


bsh_ratio_by_fleet <- bsh_donors |>
  group_by(
    period,
    FLEET_CODE,
    FLEET,
    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(
    n = n(),

    target_mt = sum(
      rc_target,
      na.rm = TRUE
    ),

    bsh_mt = sum(
      rc_bsh,
      na.rm = TRUE
    ),

    pooled_ratio =
      bsh_mt / target_mt,

    median_stratum_ratio =
      median(
        ratio_bsh_target,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  arrange(
    period,
    desc(target_mt)
  )

print(bsh_ratio_by_fleet, n = Inf)

# -----------------------------------------------------------------------------
# 4. Level 1: Fleet x Fishery x Period
# -----------------------------------------------------------------------------

# ratio_L1 <- bsh_donors |>
#   group_by(
#     FLEET_CODE,
#     FISHERY_CODE,
#     period
#   ) |>
#   summarise(
#     donor_n = n(),
#     donor_target_mt = sum(rc_target),
#     donor_bsh_mt = sum(rc_bsh),

#     donor_ratio =
#       donor_bsh_mt / donor_target_mt,

#     .groups = "drop"
#   ) |>
#   filter(
#     donor_n >= 2,
#     donor_target_mt > 0
#   ) |>
#   rename(
#     ratio_L1 = donor_ratio,
#     n_L1 = donor_n
#   )

ratio_L1 <- bsh_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_CODE,
    period
  ) |>
  summarise(
    donor_n = n(),
    donor_first_year = min(YEAR),
    donor_last_year = max(YEAR),

    donor_target_mt =
      sum(rc_target, na.rm = TRUE),

    donor_bsh_mt =
      sum(rc_bsh, na.rm = TRUE),

    donor_ratio =
      donor_bsh_mt / donor_target_mt,

    .groups = "drop"
  ) |>
  filter(
    donor_n >= 2,
    donor_target_mt > 0
  ) |>
  rename(
    ratio_L1 = donor_ratio,
    n_L1 = donor_n,
    donor_first_year_L1 = donor_first_year,
    donor_last_year_L1 = donor_last_year
  )
# -----------------------------------------------------------------------------
# 5. Level 2: Fleet x Longline Group x Period
# -----------------------------------------------------------------------------

ratio_L2 <- bsh_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    period
  ) |>
  # summarise(
  #   donor_n = n(),
  #   donor_target_mt = sum(rc_target),
  #   donor_bsh_mt = sum(rc_bsh),

  #   donor_ratio =
  #     donor_bsh_mt / donor_target_mt,

  #   .groups = "drop"
  # ) |>
  # filter(
  #   donor_n >= 2,
  #   donor_target_mt > 0
  # ) |>
  # rename(
  #   ratio_L2 = donor_ratio,
  #   n_L2 = donor_n
  # )
summarise(
    donor_n = n(),
    donor_first_year = min(YEAR),
    donor_last_year = max(YEAR),

    donor_target_mt =
      sum(rc_target, na.rm = TRUE),

    donor_bsh_mt =
      sum(rc_bsh, na.rm = TRUE),

    donor_ratio =
      donor_bsh_mt / donor_target_mt,

    .groups = "drop"
  ) |>
  filter(
    donor_n >= 2,
    donor_target_mt > 0
  ) |>
  rename(
    ratio_L2 = donor_ratio,
    n_L2 = donor_n,
    donor_first_year_L2 = donor_first_year,
    donor_last_year_L2 = donor_last_year
  )

# -----------------------------------------------------------------------------
# 6. Level 3: Fishery x Period
# -----------------------------------------------------------------------------

ratio_L3 <- bsh_donors |>
  group_by(
    FISHERY_CODE,
    period
  ) |>
  # summarise(
  #   donor_n = n(),
  #   n_donor_fleets = n_distinct(FLEET_CODE),

  #   donor_target_mt = sum(
  #     rc_target,
  #     na.rm = TRUE
  #   ),

  #   donor_bsh_mt = sum(
  #     rc_bsh,
  #     na.rm = TRUE
  #   ),

  #   donor_ratio =
  #     donor_bsh_mt / donor_target_mt,

  #   .groups = "drop"
  # ) |>
  # filter(
  #   donor_n >= 3,
  #   donor_target_mt > 0
  # ) |>
  # rename(
  #   ratio_L3 = donor_ratio,
  #   n_L3 = donor_n,
  #   n_fleets_L3 = n_donor_fleets
  # )
summarise(
    donor_n = n(),
    n_donor_fleets = n_distinct(FLEET_CODE),
    donor_first_year = min(YEAR),
    donor_last_year = max(YEAR),

    donor_target_mt =
      sum(rc_target, na.rm = TRUE),

    donor_bsh_mt =
      sum(rc_bsh, na.rm = TRUE),

    donor_ratio =
      donor_bsh_mt / donor_target_mt,

    .groups = "drop"
  ) |>
  filter(
    donor_n >= 3,
    donor_target_mt > 0
  ) |>
  rename(
    ratio_L3 = donor_ratio,
    n_L3 = donor_n,
    n_fleets_L3 = n_donor_fleets,
    donor_first_year_L3 = donor_first_year,
    donor_last_year_L3 = donor_last_year
  )



 
# -----------------------------------------------------------------------------
# 7. Level 4: All longline within period
# -----------------------------------------------------------------------------

ratio_L4 <- bsh_donors |>
  group_by(period) |>
  summarise(
    donor_n = n(),
    n_donor_fleets = n_distinct(FLEET_CODE),

    donor_first_year = min(YEAR),
    donor_last_year = max(YEAR),

    donor_target_mt =
      sum(rc_target, na.rm = TRUE),

    donor_bsh_mt =
      sum(rc_bsh, na.rm = TRUE),

    donor_ratio =
      donor_bsh_mt / donor_target_mt,

    .groups = "drop"
  ) |>
  filter(
    donor_n >= 2,
    donor_target_mt > 0
  ) |>
  rename(
    ratio_L4 = donor_ratio,
    n_L4 = donor_n,
    n_fleets_L4 = n_donor_fleets,
    donor_first_year_L4 = donor_first_year,
    donor_last_year_L4 = donor_last_year
  )

# -----------------------------------------------------------------------------
#rename
# -----------------------------------------------------------------------------
ratio_L1 <- ratio_L1 |>
  rename(
    donor_target_mt_L1 = donor_target_mt,
    donor_bsh_mt_L1 = donor_bsh_mt
  )

ratio_L2 <- ratio_L2 |>
  rename(
    donor_target_mt_L2 = donor_target_mt,
    donor_bsh_mt_L2 = donor_bsh_mt
  )

ratio_L3 <- ratio_L3 |>
  rename(
    donor_target_mt_L3 = donor_target_mt,
    donor_bsh_mt_L3 = donor_bsh_mt
  )

ratio_L4 <- ratio_L4 |>
  rename(
    donor_target_mt_L4 = donor_target_mt,
    donor_bsh_mt_L4 = donor_bsh_mt
  )

glimpse(ratio_L4)
# -----------------------------------------------------------------------------
# 8. Candidate longline strata
# -----------------------------------------------------------------------------

bsh_ratio_candidates <- bsh_ratio_dat |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    bsh_reporting_state ==
      "R1_target_fishing_no_shark_reporting",
    rc_target > 0
  )


# -----------------------------------------------------------------------------
# 9. Attach donor ratios
# -----------------------------------------------------------------------------

bsh_ratio_candidates <- bsh_ratio_candidates |>

  # Level 1: same fleet x fishery x period
  left_join(
    ratio_L1,
    by = c(
      "FLEET_CODE",
      "FISHERY_CODE",
      "period"
    )
  ) |>

  # Level 2: same fleet x broad longline group x period
  left_join(
    ratio_L2,
    by = c(
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>

  # Level 3: same fishery x period
  left_join(
    ratio_L3,
    by = c(
      "FISHERY_CODE",
      "period"
    )
  ) |>

  # Level 4: all longline x period
  left_join(
    ratio_L4,
    by = "period"
  )

# -----------------------------------------------------------------------------
# 10. Select donor level and characterize evidential support
# -----------------------------------------------------------------------------

bsh_ratio_candidates <- bsh_ratio_candidates |>
  mutate(

    # -------------------------------------------------------------------------
    # Select most specific available donor ratio
    # -------------------------------------------------------------------------

    ratio_source = case_when(

      !is.na(ratio_L1) ~
        "fleet_fishery_period",

      !is.na(ratio_L2) ~
        "fleet_longline_period",

      !is.na(ratio_L3) ~
        "fishery_period",

      !is.na(ratio_L4) ~
        "longline_period",

      TRUE ~
        "no_donor"
    ),


    selected_ratio = coalesce(
      ratio_L1,
      ratio_L2,
      ratio_L3,
      ratio_L4
    ),


    # -------------------------------------------------------------------------
    # Evidential support
    # -------------------------------------------------------------------------

    ratio_support = case_when(

      ratio_source == "no_donor" ~
        "none",

      ratio_source == "fleet_fishery_period" ~
        "same_fleet_same_fishery",

      ratio_source == "fleet_longline_period" ~
        "same_fleet_broader_longline",

      ratio_source == "fishery_period" &
        coalesce(n_fleets_L3, 0L) == 1 ~
        "single_fleet_proxy",

      ratio_source == "fishery_period" &
        coalesce(n_fleets_L3, 0L) >= 2 ~
        "multi_fleet_same_fishery",

      ratio_source == "longline_period" &
        coalesce(n_fleets_L4, 0L) == 1 ~
        "single_fleet_proxy",

      ratio_source == "longline_period" &
        coalesce(n_fleets_L4, 0L) >= 2 ~
        "multi_fleet_longline",

      TRUE ~
        "other"
    ),


    # -------------------------------------------------------------------------
    # Years represented by selected donor pool
    # -------------------------------------------------------------------------

    selected_donor_first_year = case_when(

      ratio_source == "fleet_fishery_period" ~
        donor_first_year_L1,

      ratio_source == "fleet_longline_period" ~
        donor_first_year_L2,

      ratio_source == "fishery_period" ~
        donor_first_year_L3,

      ratio_source == "longline_period" ~
        donor_first_year_L4,

      TRUE ~
        NA_real_
    ),


    selected_donor_last_year = case_when(

      ratio_source == "fleet_fishery_period" ~
        donor_last_year_L1,

      ratio_source == "fleet_longline_period" ~
        donor_last_year_L2,

      ratio_source == "fishery_period" ~
        donor_last_year_L3,

      ratio_source == "longline_period" ~
        donor_last_year_L4,

      TRUE ~
        NA_real_
    ),


    # -------------------------------------------------------------------------
    # Temporal support
    #
    # Distinguishes interpolation within the observed donor years from
    # backcasting or forward extrapolation.
    # -------------------------------------------------------------------------

    temporal_support = case_when(

      ratio_source == "no_donor" ~
        "none",

      YEAR < selected_donor_first_year ~
        "backcast",

      YEAR > selected_donor_last_year ~
        "forward_extrapolation",

      TRUE ~
        "within_donor_period"
    ),


    # -------------------------------------------------------------------------
    # Reconstructed BSH
    # -------------------------------------------------------------------------

    bsh_ratio_estimate = if_else(
      !is.na(selected_ratio),
      rc_target * selected_ratio,
      NA_real_
    )
  )
#---------------------------------------------------------------

# diagnostic output
bsh_donors |>
  filter(YEAR <= 1980) |>
  count(
    FLEET_CODE,
    FLEET,
    FISHERY_CODE,
    sort = TRUE
  )

bsh_donors |>
  filter(YEAR <= 1980) |>
  group_by(
    FLEET_CODE,
    FLEET
  ) |>
  summarise(
    first_year = min(YEAR),
    last_year = max(YEAR),
    n_years = n_distinct(YEAR),
    target_mt = sum(rc_target),
    bsh_mt = sum(rc_bsh),
    pooled_ratio = bsh_mt / target_mt,
    .groups = "drop"
  ) |>
  arrange(desc(target_mt))

# other diagnostics
 
jpn_early_ratio <- reconstruction_strata |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    FLEET_CODE == "JPN",
    YEAR <= 1980,
    rc_target > 0
  ) |>
  group_by(YEAR) |>
  summarise(
    target_mt = sum(rc_target),
    bsh_mt = sum(rc_bsh),
    shark_nei_mt = sum(rc_shark_nei),
    bsh_target_ratio =
      if_else(
        target_mt > 0,
        bsh_mt / target_mt,
        NA_real_
      ),
    .groups = "drop"
  )

library(ggplot2)

ggplot(
  jpn_early_ratio,
  aes(YEAR, bsh_target_ratio)
) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Reported BSH / target catch",
    title = "Japan longline blue shark reporting ratio",
    subtitle = "1952–1980"
  ) +
  theme_bw()



bsh_ratio_candidates |>
  count(
    ratio_support,
    temporal_support,
    sort = TRUE
  )

 

bsh_ratio_support_summary <- bsh_ratio_candidates |>
  group_by(
    ratio_support,
    temporal_support
  ) |>
  summarise(
    n_strata = n(),
    target_mt = sum(rc_target, na.rm = TRUE),
    reconstructed_bsh_mt = sum(bsh_ratio_estimate, na.rm = TRUE),
    first_year = min(YEAR),
    last_year = max(YEAR),
    .groups = "drop"
  ) |>
  arrange(desc(reconstructed_bsh_mt))

print(bsh_ratio_support_summary, n = Inf)

bsh_ratio_early_support <- bsh_ratio_candidates |>
  filter(YEAR <= 1980) |>
  group_by(
    ratio_support,
    temporal_support
  ) |>
  summarise(
    n_strata = n(),
    target_mt = sum(rc_target, na.rm = TRUE),
    reconstructed_bsh_mt = sum(bsh_ratio_estimate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(reconstructed_bsh_mt))

print(bsh_ratio_early_support, n = Inf)
