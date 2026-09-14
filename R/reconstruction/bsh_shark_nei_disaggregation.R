# =============================================================================
# bsh_shark_nei_disaggregation.R
#
# Disaggregate generic shark retained catch to blue shark (BSH).
#
# PURPOSE
# -------
# This script addresses under-reporting to species:
#
#   Shark catch was reported
#       ↓
#   Catch was reported as a generic shark category
#       ↓
#   Estimate the proportion attributable to BSH
#
# This is distinct from complete non-reporting, which is handled by the
# ratio/GAM reconstruction.
#
# IMPORTANT
# ---------
# The script does NOT add estimates for strata with no shark reporting.
# It only partitions generic shark catch that already exists in RC.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# -----------------------------------------------------------------------------
# 1. Generic shark categories
#
# These are the categories currently being treated as potentially containing
# blue shark.
# -----------------------------------------------------------------------------

generic_shark_codes <- c(
  "SKH",    # Sharks various nei
  "SHRK",   # Sharks various nei
  "AG22",   # Sharks various nei
  "SKX",    # Sharks, rays, skates, etc. nei
  "SKHS"    # Sharks finned
)


# Not currently included.
# Retained here explicitly for later sensitivity analyses.
excluded_generic_shark_codes <- c(
  "MSK",    # Sharks mackerel and porbeagles nei
  "AG20",   # Sharks mackerel, porbeagles nei
  "SHXX"    # Sharks nei other than oceanic whitetip shark and blue shark
)

 



# =============================================================================
# 2. First make an inventory of ALL shark-related species codes
#
# Do this before estimating species composition.
#
# We do not want to assume that BSH and SMA are the only identified sharks.
# =============================================================================

shark_species_inventory <- rc_dat |>
  filter(
    str_detect(
      coalesce(SPECIES_CATEGORY, ""),
      regex("shark", ignore_case = TRUE)
    ) |
      str_detect(
        coalesce(SPECIES, ""),
        regex("shark", ignore_case = TRUE)
      ) |
      SPECIES_CODE %in% c(
        "BSH",
        "SMA",
        generic_shark_codes,
        excluded_generic_shark_codes
      )
  ) |>
  distinct(
    SPECIES_CATEGORY_CODE,
    SPECIES_CATEGORY,
    SPECIES_CODE,
    SPECIES,
    SPECIES_SCIENTIFIC
  ) |>
  arrange(
    SPECIES_CATEGORY_CODE,
    SPECIES_CODE
  )

print(shark_species_inventory, n = Inf)

# =============================================================================
# 3. Define reconstruction periods
#
# These periods are provisional and represent reporting regimes rather than
# biological periods.
# =============================================================================

rc_shark_dat <- rc_dat |>
  mutate(
    period = case_when(
      YEAR <= 1963 ~ "1950_1963",
      YEAR <= 1980 ~ "1964_1980",
      YEAR <= 1994 ~ "1981_1994",
      YEAR <= 2004 ~ "1995_2004",
      YEAR <= 2014 ~ "2005_2014",
      TRUE         ~ "2015_2025"
    )
  )


# =============================================================================
# 4. Identified shark species
#
# IMPORTANT:
# Replace this vector after reviewing shark_species_inventory.
#
# Include species-specific shark codes.
# Do NOT include generic shark categories.
# =============================================================================

identified_shark_codes <- c(
  "ALV",  # Thresher
  "BSH",  # Blue shark
  "BTH",  # Bigeye thresher
  "FAL",  # Silky shark
  "LMA",  # Longfin mako
  "MAK",  # Mako sharks
  "OCS",  # Oceanic whitetip shark
  "POR",  # Porbeagle
  "PTH",  # Pelagic thresher
  "RHN",  # Whale shark
  "SMA",  # Shortfin mako
  "SPK",  # Great hammerhead
  "SPL",  # Scalloped hammerhead
  "SPN",  # Hammerhead sharks nei
  "SPY",  # Hammerheads etc. nei
  "SPZ",  # Smooth hammerhead
  "THR"   # Thresher sharks nei
)


# =============================================================================
# 5. Species-specific shark composition by Fleet x Fishery x Year
# =============================================================================

shark_composition_strata <- rc_shark_dat |>
  group_by(
    YEAR,
    period,
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY
  ) |>
  summarise(

    reported_bsh_mt = sum(
      CATCH[SPECIES_CODE == "BSH"],
      na.rm = TRUE
    ),

    reported_sma_mt = sum(
      CATCH[SPECIES_CODE == "SMA"],
      na.rm = TRUE
    ),

    identified_shark_mt = sum(
      CATCH[SPECIES_CODE %in% identified_shark_codes],
      na.rm = TRUE
    ),

    generic_shark_mt = sum(
      CATCH[SPECIES_CODE %in% generic_shark_codes],
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>
  mutate(

    p_bsh_observed = if_else(
      identified_shark_mt > 0,
      reported_bsh_mt / identified_shark_mt,
      NA_real_
    ), 
 
 
    p_sma_observed = if_else(
      identified_shark_mt > 0,
      reported_sma_mt / identified_shark_mt,
      NA_real_
    )
  )

# =============================================================================
# 6. Species-composition donors
# =============================================================================

bsh_comp_donors <- shark_composition_strata |>
  filter(
    identified_shark_mt > 0,
    !is.na(p_bsh_observed)
  )



# Build the hierarchical BSH composition estimates
# =============================================================================
# 7. Level 1:
#    Fleet x Fishery x Period
# =============================================================================

comp_L1 <- bsh_comp_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L1 = n(),
    n_years_L1 = n_distinct(YEAR),

    donor_identified_shark_mt_L1 =
      sum(identified_shark_mt, na.rm = TRUE),

    donor_bsh_mt_L1 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L1 =
      donor_bsh_mt_L1 /
      donor_identified_shark_mt_L1,

    .groups = "drop"
  ) |>
  filter(
    n_L1 >= 2,
    donor_identified_shark_mt_L1 > 0
  )

# =============================================================================
# 8. Level 2:
#    Fleet x Fishery Group x Period
# =============================================================================
comp_L2 <- bsh_comp_donors |>
  group_by(
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L2 = n(),

    donor_identified_shark_mt_L2 =
      sum(identified_shark_mt, na.rm = TRUE),

    donor_bsh_mt_L2 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L2 =
      donor_bsh_mt_L2 /
      donor_identified_shark_mt_L2,

    .groups = "drop"
  ) |>
  filter(
    n_L2 >= 2,
    donor_identified_shark_mt_L2 > 0
  )

# =============================================================================
# 9. Level 3:
#    Fishery x Period
# =============================================================================

comp_L3 <- bsh_comp_donors |>
  group_by(
    FISHERY_CODE,
    period
  ) |>
  summarise(

    n_L3 = n(),
    n_fleets_L3 = n_distinct(FLEET_CODE),

    donor_identified_shark_mt_L3 =
      sum(identified_shark_mt, na.rm = TRUE),

    donor_bsh_mt_L3 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L3 =
      donor_bsh_mt_L3 /
      donor_identified_shark_mt_L3,

    .groups = "drop"
  ) |>
  filter(
    n_L3 >= 3,
    donor_identified_shark_mt_L3 > 0
  )


# =============================================================================
# 10. Level 4:
#     Fishery Group x Period
# =============================================================================

comp_L4 <- bsh_comp_donors |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_L4 = n(),
    n_fleets_L4 = n_distinct(FLEET_CODE),

    donor_identified_shark_mt_L4 =
      sum(identified_shark_mt, na.rm = TRUE),

    donor_bsh_mt_L4 =
      sum(reported_bsh_mt, na.rm = TRUE),

    p_bsh_L4 =
      donor_bsh_mt_L4 /
      donor_identified_shark_mt_L4,

    .groups = "drop"
  ) |>
  filter(
    donor_identified_shark_mt_L4 > 0
  )

# =============================================================================
# 11. Candidate generic-shark strata
# =============================================================================

bsh_nei_candidates <- shark_composition_strata |>
  filter(generic_shark_mt > 0) |>
  left_join(
    comp_L1,
    by = c(
      "FLEET_CODE",
      "FISHERY_CODE",
      "period"
    )
  ) |>
  left_join(
    comp_L2,
    by = c(
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>
  left_join(
    comp_L3,
    by = c(
      "FISHERY_CODE",
      "period"
    )
  ) |>
  left_join(
    comp_L4,
    by = c(
      "FISHERY_GROUP_CODE",
      "period"
    )
  ) |>
  mutate(

    composition_source = case_when(
      !is.na(p_bsh_L1) ~ "fleet_fishery_period",
      !is.na(p_bsh_L2) ~ "fleet_group_period",
      !is.na(p_bsh_L3) ~ "fishery_period",
      !is.na(p_bsh_L4) ~ "fishery_group_period",
      TRUE ~ "no_donor"
    ),

    selected_p_bsh = coalesce(
      p_bsh_L1,
      p_bsh_L2,
      p_bsh_L3,
      p_bsh_L4
    ),

    estimated_bsh_from_nei = if_else(
      !is.na(selected_p_bsh),
      generic_shark_mt * selected_p_bsh,
      NA_real_
    ),

    bsh_species_adjusted =
      reported_bsh_mt +
      coalesce(estimated_bsh_from_nei, 0)
  )


# # ============================================================================= delete this 
# # 12. Attach BSH composition estimates
# # =============================================================================

# # bsh_nei_candidates <- bsh_nei_candidates |>

# #   left_join(
# #     comp_L1,
# #     by = c(
# #       "FLEET_CODE",
# #       "FISHERY_CODE",
# #       "period"
# #     )
# #   ) |>

# #   left_join(
# #     comp_L2,
# #     by = c(
# #       "FLEET_CODE",
# #       "FISHERY_GROUP_CODE",
# #       "period"
# #     )
# #   ) |>

# #   left_join(
# #     comp_L3,
# #     by = c(
# #       "FISHERY_CODE",
# #       "period"
# #     )
# #   ) |>

# #   left_join(
# #     comp_L4,
# #     by = c(
# #       "FISHERY_GROUP_CODE",
# #       "period"
# #     )
# #   )
# # =============================================================================
# # 13. Select composition donor
# # =============================================================================

# bsh_nei_candidates <- bsh_nei_candidates |>
#   mutate(

#     composition_source = case_when(

#       !is.na(p_bsh_L1) ~
#         "fleet_fishery_period",

#       !is.na(p_bsh_L2) ~
#         "fleet_group_period",

#       !is.na(p_bsh_L3) ~
#         "fishery_period",

#       !is.na(p_bsh_L4) ~
#         "fishery_group_period",

#       TRUE ~
#         "no_donor"
#     ),

#     selected_p_bsh = coalesce(
#       p_bsh_L1,
#       p_bsh_L2,
#       p_bsh_L3,
#       p_bsh_L4
#     ),

#     # Amount of generic shark catch assigned to BSH
#     bsh_from_nei = if_else(
#       !is.na(selected_p_bsh),
#       shark_nei_mt * selected_p_bsh,
#       NA_real_
#     )
#   )


# # =============================================================================
# # 14. Species-adjusted BSH
# #
# # This is NOT yet the complete reconstruction.
# #
# # It contains:
# #
# #     directly reported BSH
# #           +
# #     BSH allocated from generic shark catch
# #
# # No catch has yet been added for complete non-reporting.
# # =============================================================================

# bsh_nei_candidates <- bsh_nei_candidates |>
#   mutate(
#     bsh_species_adjusted =
#       bsh_mt +
#       coalesce(bsh_from_nei, 0)
#   )


# =============================================================================
# 12. Early longline generic-shark reconstruction
# =============================================================================

early_ll_nei_summary <- bsh_nei_candidates |>
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

    reported_bsh_mt =
      sum(reported_bsh_mt, na.rm = TRUE),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    estimated_bsh_from_nei =
      sum(estimated_bsh_from_nei, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(desc(generic_shark_mt))

print(early_ll_nei_summary, n = Inf)

# Which donor source was used for the early generic shark catch?

bsh_nei_candidates |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    YEAR <= 1980
  ) |>
  count(
    FLEET_CODE,
    FLEET,
    composition_source,
    sort = TRUE
  ) |>
  print(n = Inf)

#
# =============================================================================
# Japanese early longline shark species composition
# =============================================================================

jpn_early_shark_comp <- rc_dat |>
  filter(
    FLEET_CODE == "JPN",
    FISHERY_GROUP_CODE == "LL",
    YEAR >= 1964,
    YEAR <= 1980,
    SPECIES_CODE %in% identified_shark_codes
  ) |>
  group_by(
    SPECIES_CODE,
    SPECIES
  ) |>
  summarise(
    catch_mt = sum(CATCH, na.rm = TRUE),
    n_years = n_distinct(YEAR),
    first_year = min(YEAR),
    last_year = max(YEAR),
    .groups = "drop"
  ) |>
  mutate(
    prop_identified_sharks =
      catch_mt / sum(catch_mt)
  ) |>
  arrange(desc(catch_mt))

print(jpn_early_shark_comp, n = Inf)


jpn_early_shark_comp_year <- rc_dat |>
  filter(
    FLEET_CODE == "JPN",
    FISHERY_GROUP_CODE == "LL",
    YEAR >= 1964,
    YEAR <= 1980,
    SPECIES_CODE %in% identified_shark_codes
  ) |>
  group_by(
    YEAR,
    SPECIES_CODE,
    SPECIES
  ) |>
  summarise(
    catch_mt = sum(CATCH, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(YEAR) |>
  mutate(
    annual_identified_shark_mt = sum(catch_mt),
    proportion = catch_mt / annual_identified_shark_mt
  ) |>
  ungroup()


jpn_early_bsh_prop <- jpn_early_shark_comp_year |>
  filter(SPECIES_CODE == "BSH") |>
  select(
    YEAR,
    bsh_mt = catch_mt,
    identified_shark_mt = annual_identified_shark_mt,
    p_bsh = proportion
  )

print(jpn_early_bsh_prop, n = Inf)


library(ggplot2)

ggplot(
  jpn_early_bsh_prop,
  aes(YEAR, p_bsh)
) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = "Year",
    y = "BSH proportion of identified shark catch",
    title = "Blue shark proportion of identified shark catch",
    subtitle = "Japan longline, 1964–1980"
  ) +
  theme_bw()


