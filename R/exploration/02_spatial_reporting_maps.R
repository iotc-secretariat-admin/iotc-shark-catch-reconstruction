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
