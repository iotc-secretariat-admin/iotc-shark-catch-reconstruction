# =============================================================================
# 03_spatial_reporting_diagnostics.R
#
# Quantify temporal and spatial changes in shark reporting.
#
# Requires objects from:
#   01_decode_fishing_grounds.R
#   02_spatial_reporting_maps.R
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(scales)

# =============================================================================
# 1. Reporting state by period
# =============================================================================

reporting_period_summary <- spatial_reporting_simple |>
  group_by(
    period,
    reporting_class
  ) |>
  summarise(
    n_cells =
      n_distinct(FISHING_GROUND_CODE),

    shark_mt =
      sum(
        bsh_mt +
          sma_mt +
          shark_nei_mt,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  group_by(period) |>
  mutate(
    pct_cells =
      100 * n_cells / sum(n_cells),

    pct_shark_mt =
      100 * shark_mt / sum(shark_mt)
  ) |>
  ungroup()

print(reporting_period_summary, n = Inf)


# =============================================================================
# 2. First appearance of generic and species-specific reporting
# =============================================================================

cell_reporting_timing <- spatial_reporting_simple |>
  arrange(
    FISHING_GROUND_CODE,
    period
  ) |>
  group_by(FISHING_GROUND_CODE) |>
  summarise(

    first_period =
      first(period),

    first_generic_period =
      if (
        any(reporting_class == "Generic only")
      ) {
        min(period[
          reporting_class == "Generic only"
        ])
      } else {
        NA
      },

    first_species_period =
      if (
        any(
          reporting_class %in%
            c(
              "Species-specific",
              "Mixed generic and species-specific"
            )
        )
      ) {
        min(period[
          reporting_class %in%
            c(
              "Species-specific",
              "Mixed generic and species-specific"
            )
        ])
      } else {
        NA
      },

    ever_generic =
      any(
        reporting_class ==
          "Generic only"
      ),

    ever_species_specific =
      any(
        reporting_class %in%
          c(
            "Species-specific",
            "Mixed generic and species-specific"
          )
      ),

    .groups = "drop"
  ) |>
  mutate(
    transitioned =
      ever_generic &
        ever_species_specific
  )

cell_reporting_timing |>
  count(
    transitioned
  )

# =============================================================================
# 3. Reporting composition by fishery group
# =============================================================================

shark_reporting_by_fishery <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    CATCH_UNIT_CODE == "MT",
    SPECIES_CODE %in%
      c("BSH", "SMA", "SKH")
  ) |>
  group_by(
    period,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP
  ) |>
  summarise(

    bsh_mt =
      sum(
        CATCH[SPECIES_CODE == "BSH"],
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    shark_nei_mt =
      sum(
        CATCH[SPECIES_CODE == "SKH"],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(
    shark_total_mt =
      bsh_mt +
        sma_mt +
        shark_nei_mt,

    pct_species_specific =
      100 *
        (bsh_mt + sma_mt) /
        shark_total_mt,

    pct_generic =
      100 *
        shark_nei_mt /
        shark_total_mt
  )

print(
  shark_reporting_by_fishery,
  n = Inf
)

# =============================================================================
# 4. Reporting composition by fleet
# =============================================================================

shark_reporting_by_fleet <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    CATCH_UNIT_CODE == "MT",
    SPECIES_CODE %in%
      c("BSH", "SMA", "SKH")
  ) |>
  group_by(
    period,
    FLEET_CODE,
    FLEET
  ) |>
  summarise(

    bsh_mt =
      sum(
        CATCH[SPECIES_CODE == "BSH"],
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    shark_nei_mt =
      sum(
        CATCH[SPECIES_CODE == "SKH"],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(
    shark_total_mt =
      bsh_mt +
        sma_mt +
        shark_nei_mt,

    pct_species_specific =
      100 *
        (bsh_mt + sma_mt) /
        shark_total_mt
  ) |>
  arrange(
    period,
    desc(shark_total_mt)
  )


ggplot(
  shark_reporting_by_fishery,
  aes(
    x = period,
    y = pct_species_specific,
    group = FISHERY_GROUP_CODE,
    color = FISHERY_GROUP_CODE
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  labs(
    x = NULL,
    y = "Species-specific reporting (%)",
    color = "Fishery group",
    title = "Change in species-specific shark reporting by fishery",
    subtitle = "BSH + SMA as a percentage of BSH + SMA + Sharks NEI"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

