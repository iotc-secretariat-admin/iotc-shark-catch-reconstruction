# =============================================================================
# reconstruct_sma_scenarios.R
#
# Assemble complete shortfin mako reconstruction scenarios.
#
# Components:
#
#   1. Reported SMA
#   2. SMA allocated from MAK
#   3. SMA allocated from SKH
#   4. SMA reconstructed from complete non-reporting
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(readr)


# =============================================================================
# 1. Reconstruction keys
# =============================================================================

recon_keys <- c(
  "YEAR",
  "FLEET_CODE",
  "FISHERY_GROUP_CODE",
  "FISHERY_CODE"
)


# =============================================================================
# 2. Base RC SMA
# =============================================================================

sma_base <- sma_reporting_states |>
  select(
    all_of(
      recon_keys
    ),

    FLEET,
    FISHERY_GROUP,
    FISHERY,

    rc_sma,
    rc_mak,
    rc_shark_nei,
    rc_target,

    sma_reporting_state,
    sma_problem
  )


# =============================================================================
# 3. MAK component
# =============================================================================

mak_component <- sma_mak_candidates |>
  select(
    all_of(
      recon_keys
    ),

    p_sma_from_mak,
    sma_from_mak,
    mak_source
  )


# =============================================================================
# 4. SKH components
# =============================================================================

skh_component <- sma_skh_candidates |>
  select(
    all_of(
      recon_keys
    ),

    p_sma_N1,
    p_sma_N2,

    sma_from_skh_N1,
    sma_from_skh_N2,

    sma_N1_source,
    sma_N2_source
  )


# =============================================================================
# 5. Non-reporting component
# =============================================================================

nonreport_component <- sma_ratio_candidates |>
  select(
    all_of(
      recon_keys
    ),

    selected_ratio,

    ratio_source,
    ratio_support,
    temporal_support,

    sma_nonreport_estimate
  )


# =============================================================================
# 6. Join
# =============================================================================

sma_scenario_base <- sma_base |>

  left_join(
    mak_component,
    by = recon_keys
  ) |>

  left_join(
    skh_component,
    by = recon_keys
  ) |>

  left_join(
    nonreport_component,
    by = recon_keys
  )


 
# =============================================================================
# 7. Component exclusivity diagnostics
#
# MAK and SKH allocations may legitimately occur in the same stratum because
# they reconstruct SMA from distinct reported taxonomic catch categories.
#
# Complete non-reporting reconstruction must, however, be mutually exclusive
# from reported SMA, MAK, SKH, and their associated taxonomic reconstructions.
# =============================================================================

sma_scenario_base <- sma_scenario_base |>
  mutate(
    has_taxonomic_reconstruction =
      coalesce(sma_from_mak, 0) > 0 |
      coalesce(sma_from_skh_N1, 0) > 0,

    has_any_reported_shark_category =
      rc_sma > 0 |
      rc_mak > 0 |
      rc_shark_nei > 0,

    invalid_nonreport_overlap =
      coalesce(sma_nonreport_estimate, 0) > 0 &
      (
        has_any_reported_shark_category |
        has_taxonomic_reconstruction
      )
  )

if (
  any(
    sma_scenario_base$invalid_nonreport_overlap,
    na.rm = TRUE
  )
) {
  warning(
    "Complete non-reporting SMA reconstruction overlaps with reported or taxonomically reconstructed shark catch."
  )
}
 

# =============================================================================
# 8. Scenarios
# =============================================================================

sma_scenario_base <- sma_scenario_base |>
  mutate(

    # -------------------------------------------------------------
    # S00
    # Reported SMA only
    # -------------------------------------------------------------

    S00_reported =
      rc_sma,


    # -------------------------------------------------------------
    # S01
    # + genus-level MAK allocation
    # -------------------------------------------------------------

    S01_mak =
      rc_sma +
      coalesce(
        sma_from_mak,
        0
      ),


    # -------------------------------------------------------------
    # S02
    # + empirical HSCR allocation of SKH
    # -------------------------------------------------------------

    S02_hscr_empirical =
      rc_sma +
      coalesce(
        sma_from_mak,
        0
      ) +
      coalesce(
        sma_from_skh_N1,
        0
      ),


    # -------------------------------------------------------------
    # S03
    # + hierarchical SKH borrowing
    # -------------------------------------------------------------

    S03_hscr_hierarchical =
      rc_sma +
      coalesce(
        sma_from_mak,
        0
      ) +
      coalesce(
        sma_from_skh_N2,
        0
      ),


    # -------------------------------------------------------------
    # S04
    # + complete non-reporting
    #
    # Use N1 for SKH as the empirical species-composition baseline.
    # -------------------------------------------------------------

    S04_combined =
      rc_sma +
      coalesce(
        sma_from_mak,
        0
      ) +
      coalesce(
        sma_from_skh_N1,
        0
      ) +
      coalesce(
        sma_nonreport_estimate,
        0
      )
  )


# =============================================================================
# 9. Annual scenario series
# =============================================================================

sma_scenarios_annual <- sma_scenario_base |>
  group_by(YEAR) |>
  summarise(

    S00_reported =
      sum(
        S00_reported,
        na.rm = TRUE
      ),

    S01_mak =
      sum(
        S01_mak,
        na.rm = TRUE
      ),

    S02_hscr_empirical =
      sum(
        S02_hscr_empirical,
        na.rm = TRUE
      ),

    S03_hscr_hierarchical =
      sum(
        S03_hscr_hierarchical,
        na.rm = TRUE
      ),

    S04_combined =
      sum(
        S04_combined,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


# =============================================================================
# 10. Component summary
# =============================================================================

sma_reconstruction_summary <- sma_scenario_base |>
  summarise(

    reported_sma_mt =
      sum(
        rc_sma,
        na.rm = TRUE
      ),

    added_from_mak_mt =
      sum(
        sma_from_mak,
        na.rm = TRUE
      ),

    added_from_skh_N1_mt =
      sum(
        sma_from_skh_N1,
        na.rm = TRUE
      ),

    added_from_skh_N2_mt =
      sum(
        sma_from_skh_N2,
        na.rm = TRUE
      ),

    added_from_nonreport_mt =
      sum(
        sma_nonreport_estimate,
        na.rm = TRUE
      ),

    combined_N1_mt =
      reported_sma_mt +
      added_from_mak_mt +
      added_from_skh_N1_mt +
      added_from_nonreport_mt,

    combined_N2_mt =
      reported_sma_mt +
      added_from_mak_mt +
      added_from_skh_N2_mt +
      added_from_nonreport_mt
  )

print(
  sma_reconstruction_summary
)

# =============================================================================
# 11. Plot scenarios
# =============================================================================

sma_scenarios_long <- sma_scenarios_annual |>
  pivot_longer(
    cols = -YEAR,
    names_to = "scenario",
    values_to = "sma_mt"
  )


p_sma_scenarios <- ggplot(
  sma_scenarios_long,
  aes(
    x = YEAR,
    y = sma_mt / 1000,
    color = scenario
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  scale_color_viridis_d(
    option = "D"
  ) +
  labs(
    x = "Year",
    y = "SMA catch (1000 MT)",
    color = "Scenario",
    title = "Alternative shortfin mako catch scenarios",
    subtitle = "Reported catch and HSCR components"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_sma_scenarios


# =============================================================================
# 12. Reconstruction components through time
# =============================================================================

sma_components_annual <- sma_scenario_base |>
  group_by(YEAR) |>
  summarise(

    MAK_allocation =
      sum(
        sma_from_mak,
        na.rm = TRUE
      ),

    SKH_allocation =
      sum(
        sma_from_skh_N1,
        na.rm = TRUE
      ),

    Complete_non_reporting =
      sum(
        sma_nonreport_estimate,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  pivot_longer(
    cols = -YEAR,
    names_to = "component",
    values_to = "sma_mt"
  )

# plot  
p_sma_components <- ggplot(
  sma_components_annual,
  aes(
    x = YEAR,
    y = sma_mt / 1000,
    fill = component
  )
) +
  geom_col() +
  scale_fill_viridis_d() +
  labs(
    x = "Year",
    y = "Additional SMA catch (1000 MT)",
    fill = "Reconstruction component",
    title = "Sources of reconstructed shortfin mako catch"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_sma_components

 # =============================================================================
# 13. Save
# =============================================================================

saveRDS(
  sma_scenario_base,
  "outputs/sma/sma_scenario_base.rds"
)

write_csv(
  sma_scenarios_annual,
  "outputs/sma/sma_scenarios_annual.csv"
)

write_csv(
  sma_reconstruction_summary,
  "outputs/sma/sma_reconstruction_summary.csv"
)


#-------------- save plots
sma_output_dir <- "outputs/sma/fig" 
dir.create(sma_output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(
    sma_output_dir,
    "sma_sources_reconstruction.png"
  ),
  p_sma_components,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    sma_output_dir,
    "sma_scenarios.png"
  ),
  p_sma_scenarios,
  width = 9,
  height = 5.5,
  dpi = 300
)



#----------------------------

# check

sum(
  sma_ratio_candidates$sma_nonreport_estimate,
  na.rm = TRUE
)

sma_scenario_base |>
  summarise(
    added_from_nonreport_mt =
      sum(
        sma_nonreport_estimate,
        na.rm = TRUE
      )
  )

sma_reconstruction_summary


sma_component_summary <-
  sma_scenario_base |>
  summarise(

    reported_mt =
      sum(
        rc_sma,
        na.rm = TRUE
      ),

    mak_allocation_mt =
      sum(
        sma_from_mak,
        na.rm = TRUE
      ),

    skh_N1_mt =
      sum(
        sma_from_skh_N1,
        na.rm = TRUE
      ),

    skh_N2_mt =
      sum(
        sma_from_skh_N2,
        na.rm = TRUE
      ),

    nonreport_mt =
      sum(
        sma_nonreport_estimate,
        na.rm = TRUE
      ),

    total_N1_mt =
      sum(
        coalesce(rc_sma, 0) +
          coalesce(sma_from_mak, 0) +
          coalesce(sma_from_skh_N1, 0) +
          coalesce(sma_nonreport_estimate, 0)
      ),

    total_N2_mt =
      sum(
        coalesce(rc_sma, 0) +
          coalesce(sma_from_mak, 0) +
          coalesce(sma_from_skh_N2, 0) +
          coalesce(sma_nonreport_estimate, 0)
      )
  )

write_csv(
  sma_component_summary,
  "outputs/sma/sma_component_summary.csv"
)