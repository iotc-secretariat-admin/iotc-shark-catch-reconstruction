# =============================================================================
# 05_sma_data_inventory.R
#
# Initial data inventory for shortfin mako (SMA) catch reconstruction.
#
# PURPOSE
# -------
# Characterize:
#
#   1. Species-specific SMA reporting
#   2. Generic mako reporting (MAK = Isurus spp.)
#   3. Longfin mako reporting (LMA)
#   4. Generic shark reporting (SKH)
#   5. Fleet/fishery/temporal coverage
#   6. CA and EF availability
#   7. Spatial coverage where available
#
# This script is diagnostic only. It does not reconstruct catches.
#
# Requires:
#
#   rc_dat
#   ca_dat
#   ef_dat
#
# Optional:
#
#   ca_spatial
#   ef_spatial
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(readr)


# =============================================================================
# 1. Definitions
# =============================================================================

sma_code <- "SMA"

mako_codes <- c(
  "SMA",  # Shortfin mako
  "LMA",  # Longfin mako
  "MAK"   # Mako sharks, Isurus spp.
)

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

generic_shark_codes <- c(
  "SKH"
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
# 2. Species inventory
# =============================================================================

sma_species_inventory <- rc_dat |>
  filter(
    SPECIES_CATEGORY_CODE == "SHARKS" |
      SPECIES_CODE %in%
        c(
          identified_shark_codes,
          generic_shark_codes,
          "MAK"
        )
  ) |>
  distinct(
    SPECIES_CATEGORY_CODE,
    SPECIES_CATEGORY,
    SPECIES_CODE,
    SPECIES,
    SPECIES_SCIENTIFIC
  ) |>
  arrange(SPECIES_CODE)

print(
  sma_species_inventory,
  n = Inf
)


# =============================================================================
# 3. Annual RC catches of SMA, LMA, MAK and SKH
# =============================================================================

sma_rc_annual_species <- rc_dat |>
  filter(
    SPECIES_CODE %in%
      c(
        "SMA",
        "LMA",
        "MAK",
        "SKH"
      )
  ) |>
  group_by(
    YEAR,
    SPECIES_CODE
  ) |>
  summarise(
    catch_mt =
      sum(CATCH, na.rm = TRUE),
    .groups = "drop"
  ) |>
  complete(
    YEAR =
      min(rc_dat$YEAR, na.rm = TRUE):
      max(rc_dat$YEAR, na.rm = TRUE),

    SPECIES_CODE =
      c(
        "SMA",
        "LMA",
        "MAK",
        "SKH"
      ),

    fill = list(
      catch_mt = 0
    )
  )

print(
  sma_rc_annual_species,
  n = 30
)


# =============================================================================
# 4. Annual plot
# =============================================================================

p_sma_species_history <- ggplot(
  sma_rc_annual_species,
  aes(
    x = YEAR,
    y = catch_mt / 1000,
    color = SPECIES_CODE
  )
) +
  geom_line(linewidth = 0.8) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Reported catch (1000 MT)",
    color = "Reported category",
    title = "Reported mako and generic shark catches",
    subtitle = "IOTC retained catch data"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_sma_species_history


# =============================================================================
# 5. SMA retained catch by fishery group
# =============================================================================

sma_rc_by_group <- rc_dat |>
  add_reporting_period() |>
  filter(
    SPECIES_CODE == "SMA"
  ) |>
  group_by(
    period,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP
  ) |>
  summarise(
    sma_mt =
      sum(CATCH, na.rm = TRUE),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    n_years =
      n_distinct(YEAR),

    .groups = "drop"
  ) |>
  arrange(
    period,
    desc(sma_mt)
  )

print(
  sma_rc_by_group,
  n = Inf
)


# =============================================================================
# 6. SMA retained catch by fleet
# =============================================================================

sma_rc_by_fleet <- rc_dat |>
  filter(
    SPECIES_CODE == "SMA"
  ) |>
  group_by(
    FLEET_CODE,
    FLEET
  ) |>
  summarise(
    sma_mt =
      sum(CATCH, na.rm = TRUE),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    n_years =
      n_distinct(YEAR),

    .groups = "drop"
  ) |>
  arrange(
    desc(sma_mt)
  )

print(
  sma_rc_by_fleet,
  n = Inf
)


# =============================================================================
# 7. SMA retained catch by fleet x fishery
# =============================================================================

sma_rc_fleet_fishery <- rc_dat |>
  filter(
    SPECIES_CODE == "SMA"
  ) |>
  group_by(
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(
    sma_mt =
      sum(CATCH, na.rm = TRUE),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    n_years =
      n_distinct(YEAR),

    .groups = "drop"
  ) |>
  arrange(
    desc(sma_mt)
  )

print(
  sma_rc_fleet_fishery,
  n = 100
)


# =============================================================================
# 8. Generic MAK inventory
# =============================================================================

mak_rc_summary <- rc_dat |>
  filter(
    SPECIES_CODE == "MAK"
  ) |>
  group_by(
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(
    mak_mt =
      sum(CATCH, na.rm = TRUE),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    n_years =
      n_distinct(YEAR),

    .groups = "drop"
  ) |>
  arrange(
    desc(mak_mt)
  )

print(
  mak_rc_summary,
  n = Inf
)


# =============================================================================
# 9. Longfin mako inventory
# =============================================================================

lma_rc_summary <- rc_dat |>
  filter(
    SPECIES_CODE == "LMA"
  ) |>
  group_by(
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(
    lma_mt =
      sum(CATCH, na.rm = TRUE),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    n_years =
      n_distinct(YEAR),

    .groups = "drop"
  ) |>
  arrange(
    desc(lma_mt)
  )

print(
  lma_rc_summary,
  n = Inf
)


# =============================================================================
# 10. Period-level mako composition
#
# SMA / (SMA + LMA)
#
# This proportion will later be relevant for allocating MAK.
# =============================================================================

mako_identified_composition <- rc_dat |>
  add_reporting_period() |>
  filter(
    SPECIES_CODE %in%
      c(
        "SMA",
        "LMA"
      )
  ) |>
  group_by(
    period,
    FISHERY_GROUP_CODE
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

    identified_mako_mt =
      sma_mt + lma_mt,

    p_sma_given_mako =
      if_else(
        identified_mako_mt > 0,
        sma_mt / identified_mako_mt,
        NA_real_
      ),

    .groups = "drop"
  )

print(
  mako_identified_composition,
  n = Inf
)


# =============================================================================
# 11. Identified-shark SMA composition
#
# SMA / all identified shark species
#
# This is the quantity used to allocate SKH.
# =============================================================================

sma_identified_shark_composition <- rc_dat |>
  add_reporting_period() |>
  filter(
    SPECIES_CODE %in%
      identified_shark_codes
  ) |>
  group_by(
    period,
    FISHERY_GROUP_CODE
  ) |>
  summarise(
    sma_mt =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    identified_shark_mt =
      sum(CATCH, na.rm = TRUE),

    p_sma =
      if_else(
        identified_shark_mt > 0,
        sma_mt / identified_shark_mt,
        NA_real_
      ),

    .groups = "drop"
  )

print(
  sma_identified_shark_composition,
  n = Inf
)


# =============================================================================
# 12. CA inventory
# =============================================================================

sma_ca_summary <- ca_dat |>
  filter(
    SPECIES_CODE %in%
      c(
        "SMA",
        "LMA",
        "MAK",
        "SKH"
      )
  ) |>
  count(
    SPECIES_CODE,
    FISHERY_GROUP_CODE,
    GEAR_CODE,
    CATCH_UNIT_CODE,
    sort = TRUE
  )

print(
  sma_ca_summary,
  n = Inf
)


# =============================================================================
# 13. CA temporal coverage
# =============================================================================

sma_ca_temporal_coverage <- ca_dat |>
  filter(
    SPECIES_CODE %in%
      c(
        "SMA",
        "LMA",
        "MAK",
        "SKH"
      )
  ) |>
  group_by(
    SPECIES_CODE,
    FISHERY_GROUP_CODE,
    GEAR_CODE,
    CATCH_UNIT_CODE
  ) |>
  summarise(
    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    n_years =
      n_distinct(YEAR),

    n_records =
      n(),

    .groups = "drop"
  ) |>
  arrange(
    SPECIES_CODE,
    FISHERY_GROUP_CODE,
    first_year
  )

print(
  sma_ca_temporal_coverage,
  n = Inf
)


# =============================================================================
# 14. EF support for SMA-associated fisheries
# =============================================================================

sma_effort_summary <- ef_dat |>
  count(
    FISHERY_GROUP_CODE,
    GEAR_CODE,
    EFFORT_UNIT_CODE,
    sort = TRUE
  )

print(
  sma_effort_summary,
  n = 100
)


# =============================================================================
# 15. Optional spatial inventory
# =============================================================================

if (exists("ca_spatial")) {

  sma_spatial_summary <- ca_spatial |>
    add_reporting_period() |>
    filter(
      valid_grid_code,
      SPECIES_CODE %in%
        c(
          "SMA",
          "LMA",
          "MAK",
          "SKH"
        )
    ) |>
    group_by(
      period,
      SPECIES_CODE
    ) |>
    summarise(
      n_records =
        n(),

      n_cells =
        n_distinct(
          FISHING_GROUND_CODE
        ),

      n_fleets =
        n_distinct(
          FLEET_CODE
        ),

      .groups = "drop"
    )

  print(
    sma_spatial_summary,
    n = Inf
  )
}


# =============================================================================
# 16. Save core diagnostics
# =============================================================================

dir.create(
  "outputs/sma",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  sma_rc_by_group,
  "outputs/sma/sma_rc_by_group.csv"
)

write_csv(
  sma_rc_by_fleet,
  "outputs/sma/sma_rc_by_fleet.csv"
)

write_csv(
  mak_rc_summary,
  "outputs/sma/mak_rc_summary.csv"
)

write_csv(
  lma_rc_summary,
  "outputs/sma/lma_rc_summary.csv"
)

write_csv(
  mako_identified_composition,
  "outputs/sma/mako_identified_composition.csv"
)

write_csv(
  sma_identified_shark_composition,
  "outputs/sma/sma_identified_shark_composition.csv"
)
