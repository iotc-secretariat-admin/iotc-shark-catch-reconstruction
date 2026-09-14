# =============================================================================
# bsh_nei_historical_sensitivity.R
# =============================================================================

library(dplyr)

# Start only from unresolved N1 strata
nei_historical_candidates <- nei_candidates |>
  filter(
    generic_unallocated_N1 > 0
  )

# Earliest available composition by fishery group
earliest_group_composition <- bsh_comp_group_period |>
  filter(!is.na(p_bsh)) |>
  group_by(FISHERY_GROUP_CODE) |>
  slice_min(
    order_by = period,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  select(
    FISHERY_GROUP_CODE,
    historical_proxy_period = period,
    historical_proxy_p_bsh = p_bsh,
    donor_fleets,
    donor_years
  )


period_levels <- c(
  "1950_1963",
  "1964_1980",
  "1981_1994",
  "1995_2004",
  "2005_2014",
  "2015_2025"
)

earliest_group_composition <- bsh_comp_group_period |>
  mutate(
    period = factor(
      period,
      levels = period_levels,
      ordered = TRUE
    )
  ) |>
  filter(!is.na(p_bsh)) |>
  group_by(FISHERY_GROUP_CODE) |>
  slice_min(
    order_by = period,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()

nei_historical_candidates <- nei_historical_candidates |>
  left_join(
    earliest_group_composition |>
      select(
        FISHERY_GROUP_CODE,
        historical_proxy_period = period,
        historical_proxy_p_bsh = p_bsh,
        historical_proxy_fleets = donor_fleets,
        historical_proxy_years = donor_years
      ),
    by = "FISHERY_GROUP_CODE"
  )


nei_historical_candidates <- nei_historical_candidates |>
  mutate(

    # H0: leave unresolved
    bsh_hist_H0 = 0,

    # H1: earliest observed group-specific composition
    bsh_hist_H1 =
      generic_unallocated_N1 *
      historical_proxy_p_bsh,

    # Keep support metadata
    historical_support = case_when(
      is.na(historical_proxy_p_bsh) ~ "no_proxy",
      historical_proxy_fleets == 1 ~ "single_fleet_proxy",
      historical_proxy_fleets >= 2 ~ "multi_fleet_proxy",
      TRUE ~ "other"
    )
  )


nei_historical_summary <- nei_historical_candidates |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(
    unresolved_generic_mt =
      sum(generic_unallocated_N1, na.rm = TRUE),

    reconstructed_H1_mt =
      sum(bsh_hist_H1, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(desc(unresolved_generic_mt))

print(nei_historical_summary, n = Inf)
