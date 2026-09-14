#R/reconstruction/reconstruct_bsh_scenarios.R


# =============================================================================
# reconstruct_bsh_scenarios.R
#
# Assemble alternative blue shark catch reconstruction scenarios.
#
# This script does NOT estimate model parameters.
#
# Instead, it combines:
#
#   1. Species-specific reported BSH retained catch
#   2. BSH allocated from generic shark categories
#   3. BSH estimated for strata with no shark reporting
#
# General formulation:
#
# BSH_total =
#     reported_BSH
#   + BSH_from_generic_sharks
#   + BSH_from_non_reporting
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# =============================================================================
# 1. Define reconstruction keys
# =============================================================================

recon_keys <- c(
  "YEAR",
  "FLEET_CODE",
  "FISHERY_GROUP_CODE",
  "FISHERY_CODE"
)


# =============================================================================
# 2. Base table: reported BSH
# =============================================================================

bsh_scenario_base <- reconstruction_strata |>
  select(
    YEAR,
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY,

    rc_target,
    rc_bsh,
    rc_sma,
    rc_shark_nei,

    has_target,
    has_rc_bsh,
    has_shark_nei,
    bsh_reporting_state,
    bsh_problem
  ) |>
  rename(
    bsh_reported = rc_bsh
  )

# =============================================================================
# 3. Generic-shark disaggregation
# =============================================================================

bsh_nei_component <- bsh_nei_candidates |>
  select(
    all_of(recon_keys),

    generic_shark_mt,

    composition_source,
    selected_p_bsh,

    estimated_bsh_from_nei
  ) |>
  rename(
    bsh_nei_preferred =
      estimated_bsh_from_nei
  )

# =============================================================================
# 3. Generic-shark disaggregation
# =============================================================================

bsh_nei_component <- bsh_nei_candidates |>
  select(
    all_of(recon_keys),

    generic_shark_mt,

    composition_source,
    selected_p_bsh,

    estimated_bsh_from_nei
  ) |>
  rename(
    bsh_nei_preferred =
      estimated_bsh_from_nei
  )

# join
bsh_scenario_base <- bsh_scenario_base |>
  left_join(
    bsh_nei_component,
    by = recon_keys
  ) |>
  mutate(

    generic_shark_mt =
      coalesce(generic_shark_mt, 0),

    bsh_nei_preferred =
      coalesce(bsh_nei_preferred, 0)
  )


# =============================================================================
# 4. Ratio-based non-reporting reconstruction
# =============================================================================
glimpse(bsh_ratio_candidates)

bsh_ratio_component <- bsh_ratio_candidates |>
  select(
    all_of(recon_keys),

    selected_ratio,
    ratio_source,
    ratio_support,
    temporal_support,
    selected_donor_first_year,
    selected_donor_last_year,

    bsh_ratio_estimate
  ) |>
  rename(
    bsh_nonreport_ratio =
      bsh_ratio_estimate
  )

bsh_scenario_base <- bsh_scenario_base |>
  left_join(
    bsh_ratio_component,
    by = recon_keys
  ) |>
  mutate(

    bsh_nonreport_ratio =
      coalesce(
        bsh_nonreport_ratio,
        0
      )
  )

# =============================================================================
# 5. Double-counting checks
# =============================================================================

bsh_scenario_base <- bsh_scenario_base |>
  mutate(

    has_nei_reconstruction =
      bsh_nei_preferred > 0,

    has_nonreport_reconstruction =
      bsh_nonreport_ratio > 0,

    overlap_flag =
      has_nei_reconstruction &
      has_nonreport_reconstruction
  )

bsh_scenario_base |>
  count(overlap_flag)


# =============================================================================
# 6. Construct diagnostic scenarios
# =============================================================================

bsh_scenario_strata <- bsh_scenario_base |>
  mutate(

    # S00: current IOTC retained catch only
    bsh_s00_reported =
      bsh_reported,

    # S01: correct under-reporting to species only
    bsh_s01_nei =
      bsh_reported +
      bsh_nei_preferred,

    # S02: correct complete non-reporting only
    bsh_s02_ratio =
      bsh_reported +
      bsh_nonreport_ratio,

    # S03: both components
    bsh_s03_combined =
      bsh_reported +
      bsh_nei_preferred +
      bsh_nonreport_ratio,

    # Keep additions explicit
    added_bsh_nei =
      bsh_nei_preferred,

    added_bsh_nonreport =
      bsh_nonreport_ratio,

    added_bsh_total =
      added_bsh_nei +
      added_bsh_nonreport
  )


# =============================================================================
# 7. Annual scenarios
# =============================================================================

bsh_scenarios_annual <- bsh_scenario_strata |>
  group_by(YEAR) |>
  summarise(

    reported_bsh_mt =
      sum(bsh_s00_reported, na.rm = TRUE),

    reported_plus_nei_mt =
      sum(bsh_s01_nei, na.rm = TRUE),

    reported_plus_ratio_mt =
      sum(bsh_s02_ratio, na.rm = TRUE),

    combined_simple_mt =
      sum(bsh_s03_combined, na.rm = TRUE),

    nei_added_mt =
      sum(added_bsh_nei, na.rm = TRUE),

    nonreport_added_mt =
      sum(added_bsh_nonreport, na.rm = TRUE),

    .groups = "drop"
  )

print(bsh_scenarios_annual, n = Inf)



bsh_reconstruction_summary <- bsh_scenario_strata |>
  summarise(

    reported_bsh_mt =
      sum(bsh_reported, na.rm = TRUE),

    added_from_nei_mt =
      sum(added_bsh_nei, na.rm = TRUE),

    added_from_nonreport_mt =
      sum(added_bsh_nonreport, na.rm = TRUE),

    combined_simple_mt =
      sum(bsh_s03_combined, na.rm = TRUE)
  ) |>
  mutate(

    pct_increase_over_reported =
      100 *
      (combined_simple_mt - reported_bsh_mt) /
      reported_bsh_mt,

    pct_added_from_nei =
      100 *
      added_from_nei_mt /
      combined_simple_mt,

    pct_added_from_nonreport =
      100 *
      added_from_nonreport_mt /
      combined_simple_mt
  )

bsh_reconstruction_summary


bsh_early_ll_summary <- bsh_scenario_strata |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    YEAR <= 1980
  ) |>
  group_by(
    FLEET_CODE,
    FLEET
  ) |>
  summarise(

    first_year = min(YEAR),
    last_year = max(YEAR),

    target_mt =
      sum(rc_target, na.rm = TRUE),

    reported_bsh_mt =
      sum(bsh_reported, na.rm = TRUE),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_from_nei_mt =
      sum(added_bsh_nei, na.rm = TRUE),

    bsh_from_nonreport_mt =
      sum(added_bsh_nonreport, na.rm = TRUE),

    combined_bsh_mt =
      sum(bsh_s03_combined, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(desc(target_mt))

print(bsh_early_ll_summary, n = Inf)


ratio_source_summary <- bsh_scenario_strata |>
  filter(
    added_bsh_nonreport > 0
  ) |>
  group_by(
    ratio_source,
    ratio_support
  ) |>
  summarise(

    n_strata = n(),

    target_mt =
      sum(rc_target, na.rm = TRUE),

    reconstructed_bsh_mt =
      sum(added_bsh_nonreport, na.rm = TRUE),

    first_year = min(YEAR),
    last_year = max(YEAR),

    .groups = "drop"
  ) |>
  arrange(desc(reconstructed_bsh_mt))

ratio_source_summary
#
nei_source_summary <- bsh_scenario_strata |>
  filter(
    added_bsh_nei > 0
  ) |>
  group_by(
    composition_source
  ) |>
  summarise(

    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    reconstructed_bsh_mt =
      sum(added_bsh_nei, na.rm = TRUE),

    first_year = min(YEAR),
    last_year = max(YEAR),

    .groups = "drop"
  ) |>
  arrange(desc(reconstructed_bsh_mt))

nei_source_summary

############ plot 
bsh_scenarios_long <- bsh_scenarios_annual |>
  select(
    YEAR,
    reported_bsh_mt,
    reported_plus_nei_mt,
    reported_plus_ratio_mt,
    combined_simple_mt
  ) |>
  pivot_longer(
    cols = -YEAR,
    names_to = "scenario",
    values_to = "catch_mt"
  ) |>
  mutate(
    scenario = recode(
      scenario,
      reported_bsh_mt =
        "S00 Reported RC",
      reported_plus_nei_mt =
        "S01 Reported + generic shark allocation",
      reported_plus_ratio_mt =
        "S02 Reported + non-reporting ratio",
      combined_simple_mt =
        "S03 Combined simple reconstruction"
    )
  )
p_bsh_scenarios <- ggplot(
  bsh_scenarios_long,
  aes(
    x = YEAR,
    y = catch_mt / 1000,
    linetype = scenario, 
    color= scenario
  )
) +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Year",
    y = "BSH catch (1000 MT)",
    linetype = "Scenario",
    title = "Alternative blue shark catch scenarios",
    subtitle = "Simple reconstruction components"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_bsh_scenarios # this is old
#

#---------------------
# Join GAM allocation to the existing BSH annual scenario data
bsh_scenarios_annual2 <-
  bsh_scenarios_annual |>
  left_join(
    nei_scenario_annual |>
      select(
        YEAR,
        bsh_N1_mt,
        bsh_N2_mt
      ),
    by = "YEAR"
  ) |>
  left_join(
    bsh_gam_annual |>
      select(
        YEAR,
        bsh_GAM_mt
      ),
    by = "YEAR"
  ) |>
  mutate(
    final_HSCR_N1_mt =
      reported_bsh_mt +
      bsh_N1_mt +
      nonreport_added_mt,

    final_HSCR_N2_mt =
      reported_bsh_mt +
      bsh_N2_mt +
      nonreport_added_mt,

    final_HSCR_midpoint_mt =
      (final_HSCR_N1_mt + final_HSCR_N2_mt) / 2,

    final_GAM_mt =
      reported_bsh_mt +
      bsh_GAM_mt +
      nonreport_added_mt
  )

p_bsh_final_reconstruction <-
  ggplot(
    bsh_scenarios_annual2,
    aes(x = YEAR)
  ) +

  # HSCR N1-N2 range: primary reconstruction
  geom_ribbon(
    aes(
      ymin = pmin(final_HSCR_N1_mt, final_HSCR_N2_mt) / 1000,
      ymax = pmax(final_HSCR_N1_mt, final_HSCR_N2_mt) / 1000,
      fill = "HSCR N1-N2 range"
    ),
    alpha = 0.30
  ) +

  # HSCR midpoint
  geom_line(
    aes(
      y = final_HSCR_midpoint_mt / 1000,
      colour = "HSCR reconstruction"
    ),
    linewidth = 0.9
  ) +

  # GAM sensitivity
  geom_line(
    aes(
      y = final_GAM_mt / 1000,
      colour = "GAM sensitivity"
    ),
    linewidth = 0.9,
    linetype = "dashed"
  ) +

  # Reported catch for reference
  geom_line(
    aes(
      y = reported_bsh_mt / 1000,
      colour = "Reported BSH"
    ),
    linewidth = 0.8
  ) +

  scale_colour_viridis_d() +
  scale_fill_viridis_d() +

  labs(
    x = "Year",
    y = "BSH catch (1000 t)",
    colour = NULL,
    fill = NULL,
    title = "Reconstructed blue shark retained catch",
    subtitle =
      "HSCR N1-N2 primary reconstruction and GAM sensitivity"
  ) +

  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_bsh_final_reconstruction

ggsave(
  filename =
    "outputs/bsh/gam/figures/bsh_final_catch_comparison.png",
  plot =
    p_bsh_final_reconstruction,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

#-----------------

bsh_added_components <- bsh_scenarios_annual |>
  select(
    YEAR,
    nei_added_mt,
    nonreport_added_mt
  ) |>
  pivot_longer(
    cols = -YEAR,
    names_to = "component",
    values_to = "catch_mt"
  ) |>
  mutate(
    component = recode(
      component,
      nei_added_mt =
        "Under-reporting to species",
      nonreport_added_mt =
        "Complete non-reporting"
    )
  )

p_bsh_added_components <- ggplot(
  bsh_added_components,
  aes(
    x = YEAR,
    y = catch_mt / 1000,
    fill = component
  )
) +
  geom_col() +
  labs(
    x = "Year",
    y = "Additional estimated BSH (1000 MT)",
    fill = "Reconstruction component",
    title = "Sources of reconstructed blue shark catch"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_bsh_added_components

# summary tables
bsh_reconstruction_summary
bsh_early_ll_summary
ratio_source_summary
nei_source_summary


# # -------------------------------------------------------------------------
# Amount of generic catch by category
# -------------------------------------------------------------------------

generic_shark_inventory <- rc_dat |>
  filter(
    SPECIES_CODE %in% generic_shark_codes
  ) |>
  group_by(
    SPECIES_CODE,
    SPECIES
  ) |>
  summarise(
    catch_mt = sum(CATCH, na.rm = TRUE),
    first_year = min(YEAR),
    last_year = max(YEAR),
    n_years = n_distinct(YEAR),
    .groups = "drop"
  ) |>
  mutate(
    pct_generic_total =
      100 * catch_mt / sum(catch_mt)
  ) |>
  arrange(desc(catch_mt))

print(generic_shark_inventory, n = Inf)


# -------------------------------------------------------------------------
# Generic shark catch by broad fishery group
# -------------------------------------------------------------------------

generic_shark_by_fishery <- rc_dat |>
  filter(
    SPECIES_CODE %in% generic_shark_codes
  ) |>
  group_by(
    SPECIES_CODE,
    SPECIES,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP
  ) |>
  summarise(
    catch_mt = sum(CATCH, na.rm = TRUE),
    first_year = min(YEAR),
    last_year = max(YEAR),
    .groups = "drop"
  ) |>
  arrange(desc(catch_mt))

print(generic_shark_by_fishery, n = Inf)


nei_source_summary <- bsh_scenario_strata |>
  filter(
    added_bsh_nei > 0
  ) |>
  group_by(
    composition_source
  ) |>
  summarise(
    n_strata = n(),
    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),
    reconstructed_bsh_mt =
      sum(added_bsh_nei, na.rm = TRUE),
    first_year = min(YEAR),
    last_year = max(YEAR),
    .groups = "drop"
  ) |>
  arrange(desc(reconstructed_bsh_mt))

nei_source_summary


generic_core <- c(
  "SKH",
  "SHRK",
  "AG22"
)

generic_extended <- c(
  generic_core,
  "SKHS"
)

generic_broad <- c(
  generic_extended,
  "SKX"
)


 

write_csv(
  bsh_scenarios_annual,
  file.path(
    bsh_gam_output_dir,
    "bsh_scenarios_annual.csv"
  )
)


write_csv(
  bsh_scenarios_annual2,
  file.path(
    bsh_gam_output_dir,
    "bsh_scenarios_annual2.csv"
  )
)


bsh_scenarios_annual2