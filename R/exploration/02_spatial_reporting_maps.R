# =============================================================================
# 02_spatial_reporting_maps.R
#
# Exploratory spatial maps of:
#
#   1. Fishing effort
#   2. Blue shark catch
#   3. Shortfin mako catch
#   4. Generic shark / Sharks NEI catch
#
# Maps are summarized by broad reporting periods in order to identify
# changes in spatial reporting coverage through time.
#
# Requires:
#
#   ca_spatial
#   ef_spatial
#
# from:
#
#   01_decode_fishing_grounds.R
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)


# =============================================================================
# 1. Reporting periods
# =============================================================================

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
# 2. Species groups
# =============================================================================

generic_shark_codes <- c(
  "SKH"
)

species_map_codes <- c(
  "BSH",
  "SMA"
)

# =============================================================================
# 3. Basemap
# =============================================================================

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)


indian_ocean_xlim <- c(20, 150)
indian_ocean_ylim <- c(-60, 30)

# =============================================================================
# 4. Blue shark catch by period
# =============================================================================

bsh_spatial_period <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    SPECIES_CODE == "BSH",
    CATCH_UNIT_CODE == "MT"
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(

    catch_mt =
      sum(CATCH, na.rm = TRUE),

    n_records =
      n(),

    n_fleets =
      n_distinct(FLEET_CODE),

    .groups = "drop"
  )

# =============================================================================
# 5. Shortfin mako catch by period
# =============================================================================

sma_spatial_period <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    SPECIES_CODE == "SMA",
    CATCH_UNIT_CODE == "MT"
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(

    catch_mt =
      sum(CATCH, na.rm = TRUE),

    n_records =
      n(),

    n_fleets =
      n_distinct(FLEET_CODE),

    .groups = "drop"
  )


# =============================================================================
# 6. Sharks NEI catch by period
# =============================================================================

shark_nei_spatial_period <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    SPECIES_CODE %in% generic_shark_codes,
    CATCH_UNIT_CODE == "MT"
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(

    catch_mt =
      sum(CATCH, na.rm = TRUE),

    n_records =
      n(),

    n_fleets =
      n_distinct(FLEET_CODE),

    .groups = "drop"
  )

#===========================================================================
# Function, Map of catch - 
#========================================================================== 
ptbg     <- rgb(red=105,green=105, blue=75, alpha=75, maxColorValue=255)  
myorange <- rgb(red=255, green=165, blue=0, alpha=75, maxColorValue=255) 
myred    <- rgb(red=255, green=2, blue=0, alpha=75, maxColorValue=255) 
mygrey   <-rgb(red=220,green=220, blue=220, alpha=75, maxColorValue=255)  
myblue   <-rgb(red=0,green=205, blue=255, alpha=175, maxColorValue=255)
mygold   <-rgb(red=255,green=215, blue=0, alpha=105, maxColorValue=255)



plot_spatial_catch <- function(
    dat,
    title,
    subtitle = NULL
) {

  ggplot() +

    geom_sf(
      data = world,
      fill = "grey85",
      linewidth = 0.2
    ) +

    geom_point(
      data = dat,
      aes(
        x = lon_mid,
        y = lat_mid,
        size = catch_mt
      ),
      alpha = 0.65, 
      color = myblue
    ) +

    facet_wrap(
      ~ period,
      ncol = 2
    ) +

    scale_size_continuous(
      trans = "sqrt",
      range = c(0.5, 7),
      labels = label_number(
        big.mark = ","
      )
    ) +

    coord_sf(
      xlim = indian_ocean_xlim,
      ylim = indian_ocean_ylim,
      expand = FALSE
    ) +

    labs(
      x = NULL,
      y = NULL,
      size = "Catch (MT)",
      title = title,
      subtitle = subtitle
    ) +

    theme_bw() +

    theme(
      panel.grid.major =
        element_line(linewidth = 0.2),

      panel.grid.minor =
        element_blank(),

      legend.position =
        "bottom"
    )
}

#============plot BSH reports by space and time period 
p_bsh_spatial <- plot_spatial_catch(
  bsh_spatial_period,
  title =
    "Spatial distribution of reported blue shark catch",
  subtitle =
    "IOTC geo-referenced catch data, metric tonnes"
)

p_bsh_spatial


p_sma_spatial <- plot_spatial_catch(
  sma_spatial_period,
  title =
    "Spatial distribution of reported shortfin mako catch",
  subtitle =
    "IOTC geo-referenced catch data, metric tonnes"
)

p_sma_spatial

# NEI -------------
p_shark_nei_spatial <- plot_spatial_catch(
  shark_nei_spatial_period,
  title =
    "Spatial distribution of reported Sharks NEI catch",
  subtitle =
    "IOTC geo-referenced catch data, SKH, metric tonnes"
)

p_shark_nei_spatial



# =============================================================================
# 8. Shortfin mako spatial catch
# =============================================================================

p_sma_spatial <- plot_spatial_catch(
  sma_spatial_period,
  title =
    "Spatial distribution of reported shortfin mako catch",
  subtitle =
    "IOTC geo-referenced catch data, metric tonnes"
)

p_sma_spatial


# =============================================================================
# 9. Sharks NEI spatial catch
# =============================================================================

p_shark_nei_spatial <- plot_spatial_catch(
  shark_nei_spatial_period,
  title =
    "Spatial distribution of reported Sharks NEI catch",
  subtitle =
    "IOTC geo-referenced catch data, SKH, metric tonnes"
)

p_shark_nei_spatial
# =============================================================================
# 10. Spatial distribution of longline effort
# =============================================================================

ll_effort_spatial_period <- ef_spatial |>
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

    effort_hooks =
      sum(
        EFFORT,
        na.rm = TRUE
      ),

    n_records =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    .groups = "drop"
  )


p_ll_effort_spatial <- ggplot() +

  geom_sf(
    data = world,
    fill = "grey85",
    linewidth = 0.2
  ) +

  geom_point(
    data = ll_effort_spatial_period,
    aes(
      x = lon_mid,
      y = lat_mid,
      size = effort_hooks
    ),
    alpha = 0.65,
    colour = "grey25"
  ) +

  facet_wrap(
    ~ period,
    ncol = 2
  ) +

  scale_size_continuous(
    trans = "sqrt",
    range = c(0.5, 7),
    labels = label_number(
      big.mark = ","
    )
  ) +

  coord_sf(
    xlim = indian_ocean_xlim,
    ylim = indian_ocean_ylim,
    expand = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    size = "Hooks",
    title =
      "Spatial distribution of reported longline effort",
    subtitle =
      "IOTC geo-referenced effort data"
  ) +

  theme_bw() +

  theme(
    panel.grid.major =
      element_line(
        linewidth = 0.2
      ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )


p_ll_effort_spatial

# =============================================================================
# 11. Spatial changes in shark reporting category
# =============================================================================

shark_reporting_spatial <- ca_spatial |>

  add_reporting_period() |>

  filter(
    valid_grid_code,
    CATCH_UNIT_CODE == "MT",
    SPECIES_CODE %in%
      c(
        "BSH",
        "SMA",
        "SKH"
      )
  ) |>

  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>

  summarise(

    bsh_mt =
      sum(
        CATCH[
          SPECIES_CODE == "BSH"
        ],
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        CATCH[
          SPECIES_CODE == "SMA"
        ],
        na.rm = TRUE
      ),

    skh_mt =
      sum(
        CATCH[
          SPECIES_CODE == "SKH"
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>

  mutate(

    has_bsh =
      bsh_mt > 0,

    has_sma =
      sma_mt > 0,

    has_skh =
      skh_mt > 0,

    reporting_category =
      case_when(

        has_bsh & has_sma & has_skh ~
          "BSH + SMA + Sharks NEI",

        has_bsh & has_sma ~
          "BSH + SMA",

        has_bsh & has_skh ~
          "BSH + Sharks NEI",

        has_sma & has_skh ~
          "SMA + Sharks NEI",

        has_bsh ~
          "BSH only",

        has_sma ~
          "SMA only",

        has_skh ~
          "Sharks NEI only",

        TRUE ~
          NA_character_
      )
  ) |>

  filter(
    !is.na(
      reporting_category
    )
  )


p_shark_reporting_spatial <- ggplot() +

  geom_sf(
    data = world,
    fill = "grey85",
    linewidth = 0.2
  ) +

  geom_point(
    data = shark_reporting_spatial,
    aes(
      x = lon_mid,
      y = lat_mid,
      colour = reporting_category
    ),
    size = 1.7,
    alpha = 0.75
  ) +

  facet_wrap(
    ~ period,
    ncol = 2
  ) +

  scale_colour_viridis_d(
    option = "D",
    end = 0.9
  ) +

  coord_sf(
    xlim = indian_ocean_xlim,
    ylim = indian_ocean_ylim,
    expand = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    colour = "Reported category",
    title =
      "Spatial changes in shark reporting",
    subtitle =
      "Blue shark, shortfin mako and Sharks NEI"
  ) +

  theme_bw() +

  theme(
    panel.grid.major =
      element_line(
        linewidth = 0.2
      ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )


p_shark_reporting_spatial

# =============================================================================
# 8. Gillnet effort
# =============================================================================

gn_trips_spatial_period <- ef_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    FISHERY_GROUP_CODE == "GN",
    EFFORT_UNIT_CODE == "TRIPS"
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(
    effort = sum(EFFORT, na.rm = TRUE),
    n_fleets = n_distinct(FLEET_CODE),
    .groups = "drop"
  )

gn_fdays_spatial_period <- ef_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    FISHERY_GROUP_CODE == "GN",
    EFFORT_UNIT_CODE == "FDAYS"
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(
    effort = sum(EFFORT, na.rm = TRUE),
    n_fleets = n_distinct(FLEET_CODE),
    .groups = "drop"
  )


# =============================================================================
# 9. Spatial shark reporting state
# =============================================================================

shark_reporting_spatial <- ca_spatial |>
  add_reporting_period() |>
  filter(
    valid_grid_code,
    CATCH_UNIT_CODE == "MT",
    SPECIES_CODE %in%
      c(
        "BSH",
        "SMA",
        "SKH"
      )
  ) |>
  group_by(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
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

    reporting_state = case_when(

      bsh_mt > 0 &
        sma_mt > 0 &
        shark_nei_mt > 0 ~
        "BSH + SMA + Sharks NEI",

      bsh_mt > 0 &
        shark_nei_mt > 0 ~
        "BSH + Sharks NEI",

      sma_mt > 0 &
        shark_nei_mt > 0 ~
        "SMA + Sharks NEI",

      bsh_mt > 0 &
        sma_mt > 0 ~
        "BSH + SMA",

      bsh_mt > 0 ~
        "BSH only",

      sma_mt > 0 ~
        "SMA only",

      shark_nei_mt > 0 ~
        "Sharks NEI only",

      TRUE ~
        "No shark reporting"
    )
  )

p_reporting_state <- ggplot() +

  geom_sf(
    data = world,
    fill = "grey85",
    linewidth = 0.2
  ) +

  geom_point(
    data = shark_reporting_spatial,
    aes(
      x = lon_mid,
      y = lat_mid,
      shape = reporting_state, 
      color = as.factor(reporting_state)
    ),
    size = 1.8,
    alpha = 0.7
  ) +

  facet_wrap(
    ~ period,
    ncol = 2
  ) +

  coord_sf(
    xlim = indian_ocean_xlim,
    ylim = indian_ocean_ylim,
    expand = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    shape = "Reported category",
    title =
      "Spatial changes in shark reporting",
    subtitle =
      "Blue shark, shortfin mako and Sharks NEI"
  ) +

  theme_bw() +

  theme(
    legend.position =
      "bottom"
  )

p_reporting_state


# =============================================================================
# Shark reporting composition by period
# =============================================================================

spatial_reporting_period_summary <- shark_reporting_spatial |>
  group_by(
    period,
    reporting_state
  ) |>
  summarise(

    n_cells =
      n_distinct(FISHING_GROUND_CODE),

    total_shark_mt =
      sum(
        bsh_mt +
        sma_mt +
        shark_nei_mt,
        na.rm = TRUE
      ),

    bsh_mt =
      sum(
        bsh_mt,
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        sma_mt,
        na.rm = TRUE
      ),

    shark_nei_mt =
      sum(
        shark_nei_mt,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  group_by(period) |>
  mutate(

    pct_cells =
      100 *
      n_cells /
      sum(n_cells),

    pct_shark_mt =
      100 *
      total_shark_mt /
      sum(total_shark_mt)
  ) |>
  ungroup()

print(
  spatial_reporting_period_summary,
  n = Inf
)

spatial_reporting_simple <- shark_reporting_spatial |>
  mutate(

    reporting_class = case_when(

      shark_nei_mt > 0 &
        bsh_mt == 0 &
        sma_mt == 0 ~
        "Generic only",

      (bsh_mt > 0 | sma_mt > 0) &
        shark_nei_mt > 0 ~
        "Mixed generic and species-specific",

      (bsh_mt > 0 | sma_mt > 0) &
        shark_nei_mt == 0 ~
        "Species-specific",

      TRUE ~
        "Other"
    )
  )

reporting_transition_summary <- spatial_reporting_simple |>
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
      100 *
      n_cells /
      sum(n_cells),

    pct_shark_mt =
      100 *
      shark_mt /
      sum(shark_mt)
  ) |>
  ungroup()


p_reporting_transition <- ggplot(
  reporting_transition_summary,
  aes(
    x = period,
    y = pct_shark_mt,
    fill = reporting_class
  )
) +
  geom_col() +
  labs(
    x = NULL,
    y = "Percentage of georeferenced shark catch",
    fill = "Reporting class",
    title = "Transition from generic to species-specific shark reporting",
    subtitle = "IOTC geo-referenced retained catch data"
  ) +
  theme_bw()  +
  scale_fill_viridis_d() +
  theme(
    legend.position = "bottom"
  )
p_reporting_transition

library(viridis)

ggplot(
  reporting_transition_summary,
  aes(
    x = period,
    y = pct_shark_mt,
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
    y = "Percentage of georeferenced shark catch",
    fill = "Reporting class",
    title = "Transition from generic to species-specific shark reporting",
    subtitle = "IOTC geo-referenced retained catch data"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


spatial_reporting_simple <- spatial_reporting_simple |>
  left_join(
    fishing_ground_lookup |>
      select(
        FISHING_GROUND_CODE,
        grid_height_deg,
        grid_width_deg 
      ),
    by = "FISHING_GROUND_CODE"
  ) |>
  mutate(

    approximate_cell_area =
      grid_height_deg *
      grid_width_deg *
      cos(
        lat_mid * pi / 180
      )
  )

reporting_footprint <- spatial_reporting_simple |>
  group_by(
    period,
    reporting_class
  ) |>
  summarise(

    n_cells =
      n_distinct(FISHING_GROUND_CODE),

    relative_area =
      sum(
        approximate_cell_area,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

cell_reporting_history <- spatial_reporting_simple |>
  select(
    period,
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid,
    reporting_class
  ) |>
  distinct()


cell_transition_summary <- cell_reporting_history |>
  group_by(
    FISHING_GROUND_CODE,
    lat_mid,
    lon_mid
  ) |>
  summarise(

    first_generic_period =
      if_else(
        any(reporting_class == "Generic only"),
        as.character(
          min(
            period[
              reporting_class ==
                "Generic only"
            ]
          )
        ),
        NA_character_
      ),

    first_species_period =
      if_else(
        any(
          reporting_class %in%
            c(
              "Species-specific",
              "Mixed generic and species-specific"
            )
        ),
        as.character(
          min(
            period[
              reporting_class %in%
                c(
                  "Species-specific",
                  "Mixed generic and species-specific"
                )
            ]
          )
        ),
        NA_character_
      ),

    ever_generic =
      any(
        reporting_class ==
          "Generic only"
      ),

    ever_species =
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
      ever_species
  )

# Map this 

p_reporting.changes <- ggplot() +

  geom_sf(
    data = world,
    fill = "grey90",
    linewidth = 0.2
  ) +

  geom_point(
    data = cell_transition_summary |>
      filter(transitioned),

    aes(
      x = lon_mid,
      y = lat_mid
    ),

    size = 1.5,
    alpha = 0.7
  ) +

  coord_sf(
    xlim = indian_ocean_xlim,
    ylim = indian_ocean_ylim,
    expand = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    title =
      "Areas transitioning from generic to species-specific shark reporting",
    subtitle =
      "Cells with Sharks NEI reporting in an earlier period and BSH/SMA reporting subsequently"
  ) +

  theme_bw()

p_reporting.changes


shark_reporting_by_fishery |>
  filter(
    FISHERY_GROUP_CODE %in% c("LL", "LI", "GN")
  ) |>
  ggplot(
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



# =============================================================================
# Output directory
# =============================================================================

spatial_figure_dir <- "outputs/spatial/figures"

dir.create(
  spatial_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 12. Save spatial figures
# =============================================================================

ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "bsh_spatial_reported_catch.png"
    ),
  plot = p_bsh_spatial,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "sma_spatial_reported_catch.png"
    ),
  plot = p_sma_spatial,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "skh_spatial_reported_catch.png"
    ),
  plot = p_shark_nei_spatial,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "longline_spatial_effort.png"
    ),
  plot = p_ll_effort_spatial,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "shark_spatial_reporting_categories.png"
    ),
  plot = p_shark_reporting_spatial,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)