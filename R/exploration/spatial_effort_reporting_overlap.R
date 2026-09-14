#04_spatial_effort_reporting_overlap.R

# =============================================================================
# 04_spatial_effort_reporting_overlap.R
#
# Relate georeferenced fishing effort to shark reporting states.
#
# First implementation:
#   Longline effort measured in hooks.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(scales)

# =============================================================================
# 1. Longline hooks by period x spatial cell
# =============================================================================

ll_hooks_by_cell <- ef_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    FISHERY_GROUP_CODE == "LL",
    EFFORT_UNIT_CODE == "HOOKS"
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(
    hooks =
      sum(EFFORT, na.rm = TRUE),

    n_fleets =
      n_distinct(FLEET_CODE),

    .groups = "drop"
  )

# =============================================================================
# 2. Attach shark reporting state
# =============================================================================

ll_effort_reporting <- ll_hooks_by_cell |>
  left_join(
    spatial_reporting_simple |>
      select(
        period,
        FISHING_GROUND_CODE,
        reporting_class
      ),
    by = c(
      "period",
      "FISHING_GROUND_CODE"
    )
  ) |>
  mutate(
    reporting_class =
      replace_na(
        reporting_class,
        "No BSH/SMA/Sharks NEI reported"
      )
  )

# =============================================================================
# 3. Hooks by reporting state
# =============================================================================

ll_effort_reporting_summary <- ll_effort_reporting |>
  group_by(
    period,
    reporting_class
  ) |>
  summarise(
    hooks =
      sum(hooks, na.rm = TRUE),

    n_cells =
      n_distinct(FISHING_GROUND_CODE),

    .groups = "drop"
  ) |>
  group_by(period) |>
  mutate(
    pct_hooks =
      100 *
        hooks /
        sum(hooks)
  ) |>
  ungroup()

print(
  ll_effort_reporting_summary,
  n = Inf
)


ggplot(
  ll_effort_reporting_summary,
  aes(
    x = period,
    y = pct_hooks,
    fill = reporting_class
  )
) +
  geom_col() +
  scale_fill_viridis_d(
    option = "D",
    end = 0.9
  ) +
  labs(
    x = NULL,
    y = "Percentage of reported longline hooks",
    fill = "Shark reporting state",
    title = "Longline effort by shark reporting state",
    subtitle = "IOTC georeferenced effort and catch data"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


# =============================================================================
# 4. Fleet-specific longline effort overlap
# =============================================================================

ll_hooks_fleet_cell <- ef_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    FISHERY_GROUP_CODE == "LL",
    EFFORT_UNIT_CODE == "HOOKS"
  ) |>
  group_by(
    period,
    FLEET_CODE,
    FLEET,
    FISHING_GROUND_CODE
  ) |>
  summarise(
    hooks =
      sum(EFFORT, na.rm = TRUE),

    .groups = "drop"
  )


ll_shark_reporting_fleet_cell <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    FISHERY_GROUP_CODE == "LL",
    CATCH_UNIT_CODE == "MT",
    SPECIES_CODE %in%
      c("BSH", "SMA", "SKH")
  ) |>
  group_by(
    period,
    FLEET_CODE,
    FISHING_GROUND_CODE
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
    reporting_class = case_when(

      shark_nei_mt > 0 &
        bsh_mt == 0 &
        sma_mt == 0 ~
        "Generic only",

      shark_nei_mt > 0 &
        (bsh_mt > 0 | sma_mt > 0) ~
        "Mixed generic and species-specific",

      shark_nei_mt == 0 &
        (bsh_mt > 0 | sma_mt > 0) ~
        "Species-specific",

      TRUE ~
        "No focal shark reporting"
    )
  )

ll_fleet_effort_reporting <- ll_hooks_fleet_cell |>
  left_join(
    ll_shark_reporting_fleet_cell |>
      select(
        period,
        FLEET_CODE,
        FISHING_GROUND_CODE,
        reporting_class
      ),
    by = c(
      "period",
      "FLEET_CODE",
      "FISHING_GROUND_CODE"
    )
  ) |>
  mutate(
    reporting_class =
      replace_na(
        reporting_class,
        "No focal shark reporting"
      )
  )

ll_fleet_effort_reporting <- ll_hooks_fleet_cell |>
  left_join(
    ll_shark_reporting_fleet_cell |>
      select(
        period,
        FLEET_CODE,
        FISHING_GROUND_CODE,
        reporting_class
      ),
    by = c(
      "period",
      "FLEET_CODE",
      "FISHING_GROUND_CODE"
    )
  ) |>
  mutate(
    reporting_class =
      replace_na(
        reporting_class,
        "No focal shark reporting"
      )
  )



ll_fleet_reporting_summary <- ll_fleet_effort_reporting |>
  group_by(
    period,
    FLEET_CODE,
    FLEET,
    reporting_class
  ) |>
  summarise(
    hooks =
      sum(hooks, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(
    period,
    FLEET_CODE,
    FLEET
  ) |>
  mutate(
    pct_hooks =
      100 *
        hooks /
        sum(hooks)
  ) |>
  ungroup()
