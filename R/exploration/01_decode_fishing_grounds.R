# =============================================================================
# 01_decode_fishing_grounds.R
#
# Decode IOTC numeric FISHING_GROUND_CODE values into geographic rectangles.
#
# IOTC regular-area coding convention:
#
#   [SIZE][QUADRANT][LATITUDE][LONGITUDE]
#
# Example:
#
#   6210115
#
#   SIZE      = 6   -> 5 x 5 degree rectangle
#   QUADRANT  = 2   -> southeast quadrant
#   LATITUDE  = 10
#   LONGITUDE = 115
#
# Therefore:
#
#   latitude  = 10 to 15 S
#   longitude = 115 to 120 E
#
# Quadrants:
#
#       4 | 1
#      ---+---
#       3 | 2
#
#   1 = Northeast: +latitude, +longitude
#   2 = Southeast: -latitude, +longitude
#   3 = Southwest: -latitude, -longitude
#   4 = Northwest: +latitude, -longitude
#
# For the Indian Ocean, quadrants 1 and 2 should dominate.
#
# OUTPUTS
# -------
# fishing_ground_lookup
#
# One row per unique numeric fishing-ground code, containing:
#
#   grid_size_code
#   quadrant
#   latitude_code
#   longitude_code
#   grid_height_deg
#   grid_width_deg
#   lat_min
#   lat_max
#   lon_min
#   lon_max
#   lat_mid
#   lon_mid
#   hemisphere_ns
#   hemisphere_ew
#   valid_grid_code
#
# =============================================================================


library(dplyr)
library(stringr)
library(tidyr)
library(readr)


# =============================================================================
# 1. IOTC rectangle-size lookup
# =============================================================================

iotc_grid_sizes <- tibble(
  grid_size_code = c(
    1L,
    2L,
    3L,
    4L,
    5L,
    6L
  ),

  grid_height_deg = c(
    5,
    10,
    10,
    20,
    1,
    5
  ),

  grid_width_deg = c(
    10,
    20,
    10,
    20,
    1,
    5
  ),

  grid_resolution = c(
    "5 x 10 degree",
    "10 x 20 degree",
    "10 x 10 degree",
    "20 x 20 degree",
    "1 x 1 degree",
    "5 x 5 degree"
  )
)

iotc_grid_sizes

# =============================================================================
# 2. Function to decode IOTC numeric fishing-ground codes
# =============================================================================

decode_iotc_fishing_ground <- function(x) {

  # Preserve original value
  code_original <- x

  # Convert safely to character
  code <- as.character(x) |>
    str_trim()

  # Remove decimal suffix if data were accidentally imported numerically
  # and then converted to character, e.g. "6210115.0"
  code <- str_remove(
    code,
    "\\.0+$"
  )

  # ---------------------------------------------------------------------------
  # Identify candidate regular-grid codes
  #
  # Expected structure:
  #
  # digit 1   = size
  # digit 2   = quadrant
  # digits 3-4 = absolute latitude
  # digits 5-7 = absolute longitude
  # ---------------------------------------------------------------------------

  is_numeric_grid <- str_detect(
    code,
    "^[1-6][1-4][0-9]{5}$"
  )

  out <- tibble(
    FISHING_GROUND_CODE = code_original,
    fishing_ground_code_chr = code,
    is_numeric_grid = is_numeric_grid
  ) |>

    mutate(

      grid_size_code = if_else(
        is_numeric_grid,
        as.integer(
          str_sub(
            fishing_ground_code_chr,
            1,
            1
          )
        ),
        NA_integer_
      ),

      quadrant = if_else(
        is_numeric_grid,
        as.integer(
          str_sub(
            fishing_ground_code_chr,
            2,
            2
          )
        ),
        NA_integer_
      ),

      latitude_code = if_else(
        is_numeric_grid,
        as.numeric(
          str_sub(
            fishing_ground_code_chr,
            3,
            4
          )
        ),
        NA_real_
      ),

      longitude_code = if_else(
        is_numeric_grid,
        as.numeric(
          str_sub(
            fishing_ground_code_chr,
            5,
            7
          )
        ),
        NA_real_
      )
    ) |>

    left_join(
      iotc_grid_sizes,
      by = "grid_size_code"
    ) |>

    mutate(

      # -----------------------------------------------------------------------
      # Hemisphere labels
      # -----------------------------------------------------------------------

      hemisphere_ns = case_when(
        quadrant %in% c(1L, 4L) ~ "N",
        quadrant %in% c(2L, 3L) ~ "S",
        TRUE ~ NA_character_
      ),

      hemisphere_ew = case_when(
        quadrant %in% c(1L, 2L) ~ "E",
        quadrant %in% c(3L, 4L) ~ "W",
        TRUE ~ NA_character_
      ),


      # -----------------------------------------------------------------------
      # Geographic boundaries
      #
      # Northern quadrants:
      #
      #   latitude_code = southern edge
      #
      # Southern quadrants:
      #
      #   latitude_code = absolute value of northern edge.
      #
      # Example:
      #
      #   quadrant 2, latitude 10, height 5
      #
      #   -15 <= latitude <= -10
      #
      # Eastern quadrants:
      #
      #   longitude_code = western edge
      #
      # Western quadrants:
      #
      #   longitude_code = absolute value of eastern edge.
      # -----------------------------------------------------------------------

      lat_min = case_when(

        quadrant %in% c(1L, 4L) ~
          latitude_code,

        quadrant %in% c(2L, 3L) ~
          -(latitude_code + grid_height_deg),

        TRUE ~
          NA_real_
      ),

      lat_max = case_when(

        quadrant %in% c(1L, 4L) ~
          latitude_code + grid_height_deg,

        quadrant %in% c(2L, 3L) ~
          -latitude_code,

        TRUE ~
          NA_real_
      ),

      lon_min = case_when(

        quadrant %in% c(1L, 2L) ~
          longitude_code,

        quadrant %in% c(3L, 4L) ~
          -(longitude_code + grid_width_deg),

        TRUE ~
          NA_real_
      ),

      lon_max = case_when(

        quadrant %in% c(1L, 2L) ~
          longitude_code + grid_width_deg,

        quadrant %in% c(3L, 4L) ~
          -longitude_code,

        TRUE ~
          NA_real_
      ),


      # -----------------------------------------------------------------------
      # Rectangle centroid
      # -----------------------------------------------------------------------

      lat_mid =
        (lat_min + lat_max) / 2,

      lon_mid =
        (lon_min + lon_max) / 2,


      # -----------------------------------------------------------------------
      # Basic validity checks
      # -----------------------------------------------------------------------

      valid_latitude =
        !is.na(lat_min) &
        !is.na(lat_max) &
        lat_min >= -90 &
        lat_max <= 90 &
        lat_min < lat_max,

      valid_longitude =
        !is.na(lon_min) &
        !is.na(lon_max) &
        lon_min >= -180 &
        lon_max <= 180 &
        lon_min < lon_max,

      valid_grid_code =
        is_numeric_grid &
        !is.na(grid_height_deg) &
        !is.na(grid_width_deg) &
        valid_latitude &
        valid_longitude
    )

  out
}

# =============================================================================
# 3. Test known IOTC example
# =============================================================================

test_6210115 <- decode_iotc_fishing_ground(
  6210115
)

t(test_6210115)


stopifnot(
  test_6210115$grid_size_code == 6,
  test_6210115$quadrant == 2,
  test_6210115$latitude_code == 10,
  test_6210115$longitude_code == 115,
  test_6210115$lat_min == -15,
  test_6210115$lat_max == -10,
  test_6210115$lon_min == 115,
  test_6210115$lon_max == 120,
  test_6210115$lat_mid == -12.5,
  test_6210115$lon_mid == 117.5
)


# =============================================================================
# 4. Inventory unique fishing-ground codes in CA and EF
# =============================================================================

ca_fishing_ground_codes <- ca_dat |>
  distinct(
    FISHING_GROUND_CODE
  ) |>
  mutate(
    source_CA = TRUE
  )


ef_fishing_ground_codes <- ef_dat |>
  distinct(
    FISHING_GROUND_CODE
  ) |>
  mutate(
    source_EF = TRUE
  )


all_fishing_ground_codes <- full_join(
  ca_fishing_ground_codes,
  ef_fishing_ground_codes,
  by = "FISHING_GROUND_CODE"
) |>
  mutate(
    source_CA =
      coalesce(
        source_CA,
        FALSE
      ),

    source_EF =
      coalesce(
        source_EF,
        FALSE
      )
  )

# =============================================================================
# 5. Decode all unique CA / EF fishing-ground codes
# =============================================================================

fishing_ground_lookup <- decode_iotc_fishing_ground(
  all_fishing_ground_codes$FISHING_GROUND_CODE
) |>
  left_join(
    all_fishing_ground_codes,
    by = "FISHING_GROUND_CODE"
  ) |>
  arrange(
    quadrant,
    latitude_code,
    longitude_code,
    grid_size_code
  )

glimpse(
  fishing_ground_lookup
)
fishing_ground_lookup <- fishing_ground_lookup |>
  mutate(
    spatial_code_type = case_when(
      valid_grid_code ~
        "regular_grid",

      str_detect(
        fishing_ground_code_chr,
        "^9"
      ) ~
        "non_regular_area",

      TRUE ~
        "unresolved"
    )
  )
# =============================================================================
# 6. QA: grid-code validity
# =============================================================================

fishing_ground_lookup |>
  count(
    spatial_code_type,
    sort = TRUE
  )

fishing_ground_lookup |>
  count(
    is_numeric_grid,
    valid_grid_code
  )

# Check 
fishing_ground_lookup |>
  filter(
    !valid_grid_code
  ) |>
  print(
    n = Inf
  )


# =============================================================================
# 7. QA: grid resolution inventory
# =============================================================================

fishing_ground_lookup |>
  filter(
    valid_grid_code
  ) |>
  count(
    grid_size_code,
    grid_resolution,
    grid_height_deg,
    grid_width_deg,
    sort = TRUE
  )


ca_spatial_resolution <- ca_dat |>
  left_join(
    fishing_ground_lookup |>
      select(
        FISHING_GROUND_CODE,
        grid_size_code,
        grid_resolution,
        valid_grid_code
      ),
    by = "FISHING_GROUND_CODE"
  ) |>
  count(
    grid_size_code,
    grid_resolution,
    valid_grid_code,
    name = "n_CA_records",
    sort = TRUE
  )

ca_spatial_resolution

ef_spatial_resolution <- ef_dat |>
  left_join(
    fishing_ground_lookup |>
      select(
        FISHING_GROUND_CODE,
        grid_size_code,
        grid_resolution,
        valid_grid_code
      ),
    by = "FISHING_GROUND_CODE"
  ) |>
  count(
    grid_size_code,
    grid_resolution,
    valid_grid_code,
    name = "n_EF_records",
    sort = TRUE
  )

ef_spatial_resolution


# =============================================================================
# 8. QA: geographic quadrants
# =============================================================================

fishing_ground_lookup |>
  filter(
    valid_grid_code
  ) |>
  count(
    quadrant,
    hemisphere_ns,
    hemisphere_ew,
    sort = TRUE
  )


# =============================================================================
# 9. QA: geographic extent
# =============================================================================

fishing_ground_lookup |>
  filter(
    valid_grid_code
  ) |>
  summarise(

    min_lat =
      min(
        lat_min,
        na.rm = TRUE
      ),

    max_lat =
      max(
        lat_max,
        na.rm = TRUE
      ),

    min_lon =
      min(
        lon_min,
        na.rm = TRUE
      ),

    max_lon =
      max(
        lon_max,
        na.rm = TRUE
      ),

    n_cells =
      n()
  )


# =============================================================================
# 10. Human-readable grid labels
# =============================================================================

fishing_ground_lookup <- fishing_ground_lookup |>
  mutate(

    grid_label = if_else(

      valid_grid_code,

      paste0(
        abs(lat_min),
        "-",
        abs(lat_max),
        "°",
        if_else(
          lat_mid >= 0,
          "N",
          "S"
        ),
        ", ",
        abs(lon_min),
        "-",
        abs(lon_max),
        "°",
        if_else(
          lon_mid >= 0,
          "E",
          "W"
        )
      ),

      NA_character_
    )
  )

# Check 
fishing_ground_lookup |>
  filter(
    fishing_ground_code_chr == "6210115"
  ) |>
  select(
    fishing_ground_code_chr,
    grid_resolution,
    grid_label,
    lat_mid,
    lon_mid
  )

# now fix 
fishing_ground_lookup <- fishing_ground_lookup |>
  mutate(

    lat_label = case_when(

      lat_mid >= 0 ~
        paste0(
          lat_min,
          "-",
          lat_max,
          "°N"
        ),

      TRUE ~
        paste0(
          abs(lat_max),
          "-",
          abs(lat_min),
          "°S"
        )
    ),

    lon_label = case_when(

      lon_mid >= 0 ~
        paste0(
          lon_min,
          "-",
          lon_max,
          "°E"
        ),

      TRUE ~
        paste0(
          abs(lon_max),
          "-",
          abs(lon_min),
          "°W"
        )
    ),

    grid_label =
      if_else(
        valid_grid_code,
        paste(
          lat_label,
          lon_label,
          sep = ", "
        ),
        NA_character_
      )
  )


# =============================================================================
# 11. Attach decoded geography to CA
# =============================================================================

ca_spatial <- ca_dat |>
  left_join(
    fishing_ground_lookup |>
      select(
        FISHING_GROUND_CODE,

        grid_size_code,
        grid_resolution,

        quadrant,

        latitude_code,
        longitude_code,

        lat_min,
        lat_max,
        lon_min,
        lon_max,

        lat_mid,
        lon_mid,

        hemisphere_ns,
        hemisphere_ew,

        grid_label,

        valid_grid_code
      ),
    by = "FISHING_GROUND_CODE"
  )


# =============================================================================
# 11. Attach decoded geography to CA
# =============================================================================

ca_spatial <- ca_dat |>
  left_join(
    fishing_ground_lookup |>
      select(
        FISHING_GROUND_CODE,

        grid_size_code,
        grid_resolution,

        quadrant,

        latitude_code,
        longitude_code,

        lat_min,
        lat_max,
        lon_min,
        lon_max,

        lat_mid,
        lon_mid,

        hemisphere_ns,
        hemisphere_ew,

        grid_label,

        valid_grid_code
      ),
    by = "FISHING_GROUND_CODE"
  )

# =============================================================================
# 12. Attach decoded geography to EF
# =============================================================================

ef_spatial <- ef_dat |>
  left_join(
    fishing_ground_lookup |>
      select(
        FISHING_GROUND_CODE,

        grid_size_code,
        grid_resolution,

        quadrant,

        latitude_code,
        longitude_code,

        lat_min,
        lat_max,
        lon_min,
        lon_max,

        lat_mid,
        lon_mid,

        hemisphere_ns,
        hemisphere_ew,

        grid_label,

        valid_grid_code
      ),
    by = "FISHING_GROUND_CODE"
  )
# =============================================================================
# 13. Join integrity checks
# =============================================================================

stopifnot(
  nrow(ca_spatial) ==
    nrow(ca_dat)
)

stopifnot(
  nrow(ef_spatial) ==
    nrow(ef_dat)
)

# make sure the lookupis unique
fishing_ground_lookup |>
  count(
    FISHING_GROUND_CODE
  ) |>
  filter(
    n > 1
  )


#  
# =============================================================================
# 14. Basic spatial coverage through time
# =============================================================================

ca_spatial_coverage_year <- ca_spatial |>
  group_by(YEAR) |>
  summarise(

    n_records =
      n(),

    n_spatial_records =
      sum(
        valid_grid_code,
        na.rm = TRUE
      ),

    pct_spatial_records =
      100 *
      n_spatial_records /
      n_records,

    n_unique_cells =
      n_distinct(
        FISHING_GROUND_CODE[
          valid_grid_code
        ]
      ),

    .groups = "drop"
  )

# effort 
ef_spatial_coverage_year <- ef_spatial |>
  group_by(YEAR) |>
  summarise(

    n_records =
      n(),

    n_spatial_records =
      sum(
        valid_grid_code,
        na.rm = TRUE
      ),

    pct_spatial_records =
      100 *
      n_spatial_records /
      n_records,

    n_unique_cells =
      n_distinct(
        FISHING_GROUND_CODE[
          valid_grid_code
        ]
      ),

    .groups = "drop"
  )

# =============================================================================
# 15. Save reusable spatial lookup and QA summaries
# =============================================================================

dir.create(
  "outputs/spatial",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  fishing_ground_lookup,
  "outputs/spatial/iotc_fishing_ground_lookup.rds"
)

write_csv(
  fishing_ground_lookup,
  "outputs/spatial/iotc_fishing_ground_lookup.csv"
)

write_csv(
  ca_spatial_coverage_year,
  "outputs/spatial/ca_spatial_coverage_by_year.csv"
)

write_csv(
  ef_spatial_coverage_year,
  "outputs/spatial/ef_spatial_coverage_by_year.csv"
)

# examine outputs 
fishing_ground_lookup |>
  count(
    grid_resolution,
    sort = TRUE
  )

fishing_ground_lookup |>
  count(
    quadrant,
    sort = TRUE
  )

fishing_ground_lookup |>
  filter(
    !valid_grid_code
  ) |>
  print(n = Inf)

fishing_ground_lookup |>
  filter(
    fishing_ground_code_chr == "6210115"
  )


