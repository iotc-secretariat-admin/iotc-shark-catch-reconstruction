#R/reconstruction/bsh_nei_scenarios.R
# =============================================================================
# bsh_nei_scenarios.R
#
# Alternative methods for allocating generic shark catch (SKH) to blue shark.
#
# The purpose is to address under-reporting to species:
#
#     generic shark catch
#            ↓
#     estimate P(BSH | shark catch)
#            ↓
#     allocate a proportion of SKH to BSH
#
# The script produces three scenarios differing in the amount of information
# borrowed across fleets and fisheries.
#The scenarios are 
#  N1 Strict empirical: use same-fleet information where available; 
#  otherwise use a multi-fleet donor from the same fishery.
#  Do not allocate SKH when the only donor is a different single fleet.
#
#  N2 Hierarchical: your current approach—same fleet/fishery → same fleet/group → fishery/period → fishery-group/period.
#  This allows single-fleet proxies and is therefore more assumption-dependent.

#  N3 Broad fishery composition: use pooled multi-fleet fishery-group × period composition. 
#  This deliberately sacrifices fleet specificity for a larger donor base.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# =============================================================================
# 1. Start from species-composition donor data
# =============================================================================

# Assumes these objects already exist:
#
# shark_composition_strata
# bsh_comp_donors
#
# with:
#
# reported_bsh_mt
# identified_shark_mt
# generic_shark_mt
# period


# =============================================================================
# 2. Diagnostic: composition by fishery group and period
# =============================================================================

bsh_comp_group_period <- bsh_comp_donors |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period
  ) |>
  summarise(

    donor_n = n(),
    donor_fleets = n_distinct(FLEET_CODE),
    donor_years = n_distinct(YEAR),

    donor_bsh_mt =
      sum(reported_bsh_mt, na.rm = TRUE),

    donor_identified_shark_mt =
      sum(identified_shark_mt, na.rm = TRUE),

    p_bsh =
      donor_bsh_mt /
      donor_identified_shark_mt,

    .groups = "drop"
  ) |>
  arrange(
    FISHERY_GROUP_CODE,
    period
  )

print(bsh_comp_group_period, n = Inf)

# =============================================================================
# 3. Diagnostic: composition by fishery and period
# =============================================================================

bsh_comp_fishery_period <- bsh_comp_donors |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    FISHERY,
    period
  ) |>
  summarise(

    donor_n = n(),
    donor_fleets = n_distinct(FLEET_CODE),
    donor_years = n_distinct(YEAR),

    donor_bsh_mt =
      sum(reported_bsh_mt, na.rm = TRUE),

    donor_identified_shark_mt =
      sum(identified_shark_mt, na.rm = TRUE),

    p_bsh =
      donor_bsh_mt /
      donor_identified_shark_mt,

    .groups = "drop"
  ) |>
  arrange(
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    period
  )

print(bsh_comp_fishery_period, n = Inf)

# =============================================================================
# 4. Fleet x Fishery x Period donor
# =============================================================================

nei_L1 <- bsh_comp_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L1 = n(),
    years_L1 = n_distinct(YEAR),

    donor_first_year_L1 = min(YEAR),
    donor_last_year_L1  = max(YEAR),

    identified_shark_mt_L1 =
      sum(identified_shark_mt, na.rm = TRUE),

    bsh_mt_L1 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L1 =
      bsh_mt_L1 / identified_shark_mt_L1,

    .groups = "drop"
  ) |>
  filter(
    identified_shark_mt_L1 > 0,
    n_L1 >= 2
  )

# =============================================================================
# 5. Fleet x Fishery Group x Period donor
# =============================================================================

nei_L2 <- bsh_comp_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L2 = n(),
    years_L2 = n_distinct(YEAR),

    donor_first_year_L2 = min(YEAR),
    donor_last_year_L2  = max(YEAR),

    identified_shark_mt_L2 =
      sum(identified_shark_mt, na.rm = TRUE),

    bsh_mt_L2 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L2 =
      bsh_mt_L2 / identified_shark_mt_L2,

    .groups = "drop"
  ) |>
  filter(
    identified_shark_mt_L2 > 0,
    n_L2 >= 2
  )

# =============================================================================
# 6. Fishery x Period donor
# =============================================================================

nei_L3 <- bsh_comp_donors |>
  group_by(
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L3 = n(),
    n_fleets_L3 = n_distinct(FLEET_CODE),
    n_years_L3 = n_distinct(YEAR),

    donor_first_year_L3 = min(YEAR),
    donor_last_year_L3  = max(YEAR),

    identified_shark_mt_L3 =
      sum(identified_shark_mt, na.rm = TRUE),

    bsh_mt_L3 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L3 =
      bsh_mt_L3 / identified_shark_mt_L3,

    .groups = "drop"
  ) |>
  filter(
    identified_shark_mt_L3 > 0,
    n_L3 >= 3
  )

# =============================================================================
# 7. Fishery Group x Period donor
# =============================================================================

nei_L4 <- bsh_comp_donors |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L4 = n(),
    n_fleets_L4 = n_distinct(FLEET_CODE),
    n_years_L4 = n_distinct(YEAR),

    donor_first_year_L4 = min(YEAR),
    donor_last_year_L4  = max(YEAR),

    identified_shark_mt_L4 =
      sum(identified_shark_mt, na.rm = TRUE),

    bsh_mt_L4 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L4 =
      bsh_mt_L4 / identified_shark_mt_L4,

    .groups = "drop"
  ) |>
  filter(
    identified_shark_mt_L4 > 0
  )

# =============================================================================
# 8. Generic shark strata
# =============================================================================

nei_candidates <- shark_composition_strata |>
  filter(
    generic_shark_mt > 0
  ) |>

  left_join(
    nei_L1,
    by = c(
      "FLEET_CODE",
      "FISHERY_CODE",
      "period"
    )
  ) |>

  left_join(
    nei_L2,
    by = c(
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>

  left_join(
    nei_L3,
    by = c(
      "FISHERY_CODE",
      "period"
    )
  ) |>

  left_join(
    nei_L4,
    by = c(
      "FISHERY_GROUP_CODE",
      "period"
    )
  )

# =============================================================================
#  Define the different scenarios for borrowing historical information on species comp.
# the hierarchy is
#   use same fleet/fishery if available;
#   use same fleet/group if available;
#   otherwise use same-fishery donors only when at least two fleets contribute;
#   otherwise leave the generic shark catch unallocated. 
# 
# =============================================================================
# 9. N1: Strict empirical composition
# =============================================================================

nei_candidates <- nei_candidates |>
  mutate(

    nei_N1_source = case_when(

      !is.na(p_bsh_L1) ~
        "same_fleet_same_fishery",

      !is.na(p_bsh_L2) ~
        "same_fleet_same_group",

      !is.na(p_bsh_L3) &
        n_fleets_L3 >= 2 ~
        "multi_fleet_same_fishery",

      TRUE ~
        "insufficient_support"
    ),

    p_bsh_N1 = case_when(

      !is.na(p_bsh_L1) ~
        p_bsh_L1,

      !is.na(p_bsh_L2) ~
        p_bsh_L2,

      !is.na(p_bsh_L3) &
        n_fleets_L3 >= 2 ~
        p_bsh_L3,

      TRUE ~
        NA_real_
    ),

    bsh_from_nei_N1 =
      generic_shark_mt * p_bsh_N1
  )


# =============================================================================
# 10. N2: Hierarchical composition
# =============================================================================

nei_candidates <- nei_candidates |>
  mutate(

    nei_N2_source = case_when(

      !is.na(p_bsh_L1) ~
        "same_fleet_same_fishery",

      !is.na(p_bsh_L2) ~
        "same_fleet_same_group",

      !is.na(p_bsh_L3) ~
        "same_fishery",

      !is.na(p_bsh_L4) ~
        "same_fishery_group",

      TRUE ~
        "no_donor"
    ),

    p_bsh_N2 = coalesce(
      p_bsh_L1,
      p_bsh_L2,
      p_bsh_L3,
      p_bsh_L4
    ),

    bsh_from_nei_N2 =
      generic_shark_mt * p_bsh_N2
  )


# =============================================================================
# 11. N3: Broad fishery-group composition
# =============================================================================

nei_candidates <- nei_candidates |>
  mutate(

    p_bsh_N3 = if_else(
      n_fleets_L4 >= 2,
      p_bsh_L4,
      NA_real_
    ),

    nei_N3_source = if_else(
      !is.na(p_bsh_N3),
      "multi_fleet_fishery_group",
      "insufficient_support"
    ),

    bsh_from_nei_N3 =
      generic_shark_mt * p_bsh_N3
  )


# =============================================================================
# 12. quantify and t rack unallocated generic shark catch
# =============================================================================

nei_candidates <- nei_candidates |>
  mutate(

    generic_unallocated_N1 =
      if_else(
        is.na(p_bsh_N1),
        generic_shark_mt,
        0
      ),

    generic_unallocated_N2 =
      if_else(
        is.na(p_bsh_N2),
        generic_shark_mt,
        0
      ),

    generic_unallocated_N3 =
      if_else(
        is.na(p_bsh_N3),
        generic_shark_mt,
        0
      )
  )


# =============================================================================
# 13. Overall NEI scenario comparison
# =============================================================================

nei_scenario_summary <- nei_candidates |>
  summarise(

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_N3_mt =
      sum(bsh_from_nei_N3, na.rm = TRUE),

    unallocated_N1_mt =
      sum(generic_unallocated_N1, na.rm = TRUE),

    unallocated_N2_mt =
      sum(generic_unallocated_N2, na.rm = TRUE),

    unallocated_N3_mt =
      sum(generic_unallocated_N3, na.rm = TRUE)
  )

nei_scenario_summary


# =============================================================================
# 14. NEI scenarios by broad fishery group
# =============================================================================

nei_scenario_by_group <- nei_candidates |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP
  ) |>
  summarise(

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_N3_mt =
      sum(bsh_from_nei_N3, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(desc(generic_shark_mt))

print(nei_scenario_by_group, n = Inf)


# =============================================================================
# 15. Annual NEI reconstruction
# =============================================================================

nei_scenario_annual <- nei_candidates |>
  group_by(YEAR) |>
  summarise(

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_N3_mt =
      sum(bsh_from_nei_N3, na.rm = TRUE),

    .groups = "drop"
  )


# ==========================================================================
# Diagnostics and checks. 
#

nei_candidates |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_generic_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    median_p_N1 =
      median(p_bsh_N1, na.rm = TRUE),

    median_p_N2 =
      median(p_bsh_N2, na.rm = TRUE),

    p_group_N3 =
      first(p_bsh_L4),

    .groups = "drop"
  ) |>
  print(n = Inf)


nei_scenario_summary
nei_scenario_by_group
bsh_comp_group_period