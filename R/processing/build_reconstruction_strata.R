# =============================================================================
# build_reconstruction_strata.R
#
# Build Fleet x Fishery x Year table used to identify strata requiring
# reconstruction of blue shark (BSH) and shortfin mako (SMA) catch.
#
# Inputs:
#   rc_dat  - IOTC retained catch data
#   ca_dat  - IOTC geo-referenced catch data
#   ef_dat  - IOTC geo-referenced effort data
#
# Output:
#   reconstruction_strata
#
# =============================================================================

library(dplyr)
library(tidyr)


# -----------------------------------------------------------------------------
# 1. Definitions
# -----------------------------------------------------------------------------

# Target species used in the 2025 ratio-based reconstruction
target_species <- c(
  "YFT",   # Yellowfin tuna
  "BET",   # Bigeye tuna
  "SKJ",   # Skipjack tuna
  "ALB",   # Albacore
  "SWO"    # Swordfish
)


# Generic shark categories currently used for reconstruction
generic_shark_codes <- c(
  "SKH",   # Sharks various nei
  "SHRK",  # Sharks various nei
  "AG22",  # Sharks various nei
  "SKX",   # Sharks, rays, skates, etc. nei
  "SKHS"   # Sharks finned
)


# Generic shark categories deliberately NOT included at this stage
excluded_generic_shark_codes <- c(
  "MSK",   # Sharks mackerel and porbeagles nei
  "AG20",  # Sharks mackerel, porbeagles nei
  "SHXX"   # Sharks nei other than OCS and BSH
)


# Reconstruction unit
recon_keys <- c(
  "YEAR",
  "FLEET_CODE",
  "FISHERY_GROUP_CODE",
  "FISHERY_CODE"
)


# -----------------------------------------------------------------------------
# 2. Summarise retained catch
# -----------------------------------------------------------------------------

rc_summary <- rc_dat |>
  group_by(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    # -------------------------------------------------------------
    # Total retained catch
    # -------------------------------------------------------------

    rc_total = sum(
      CATCH,
      na.rm = TRUE
    ),


    # -------------------------------------------------------------
    # Principal target species
    # -------------------------------------------------------------

    rc_target = sum(
      CATCH[SPECIES_CODE %in% target_species],
      na.rm = TRUE
    ),


    # Individual target species retained for diagnostics
    rc_yft = sum(
      CATCH[SPECIES_CODE == "YFT"],
      na.rm = TRUE
    ),

    rc_bet = sum(
      CATCH[SPECIES_CODE == "BET"],
      na.rm = TRUE
    ),

    rc_skj = sum(
      CATCH[SPECIES_CODE == "SKJ"],
      na.rm = TRUE
    ),

    rc_alb = sum(
      CATCH[SPECIES_CODE == "ALB"],
      na.rm = TRUE
    ),

    rc_swo = sum(
      CATCH[SPECIES_CODE == "SWO"],
      na.rm = TRUE
    ),


    # -------------------------------------------------------------
    # Species-specific shark catches
    # -------------------------------------------------------------

    rc_bsh = sum(
      CATCH[SPECIES_CODE == "BSH"],
      na.rm = TRUE
    ),

    rc_sma = sum(
      CATCH[SPECIES_CODE == "SMA"],
      na.rm = TRUE
    ),


    # -------------------------------------------------------------
    # Generic / unidentified shark catch
    # -------------------------------------------------------------

    rc_shark_nei = sum(
      CATCH[SPECIES_CODE %in% generic_shark_codes],
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 3. Add descriptive labels from RC
#
# Codes are used as keys; labels are descriptive only.
# -----------------------------------------------------------------------------

rc_labels <- rc_dat |>
  distinct(
    FLEET_CODE,
    FLEET,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY
  )


rc_summary <- rc_summary |>
  left_join(
    rc_labels,
    by = c(
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "FISHERY_CODE"
    )
  )


# -----------------------------------------------------------------------------
# 4. Summarise CA data
#
# Keep MT and NO separately.
# -----------------------------------------------------------------------------

ca_summary <- ca_dat |>
  filter(
    SPECIES_CODE %in% c("BSH", "SMA")
  ) |>
  group_by(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    ca_bsh_mt = sum(
      CATCH[
        SPECIES_CODE == "BSH" &
        CATCH_UNIT_CODE == "MT"
      ],
      na.rm = TRUE
    ),

    ca_bsh_no = sum(
      CATCH[
        SPECIES_CODE == "BSH" &
        CATCH_UNIT_CODE == "NO"
      ],
      na.rm = TRUE
    ),

    ca_sma_mt = sum(
      CATCH[
        SPECIES_CODE == "SMA" &
        CATCH_UNIT_CODE == "MT"
      ],
      na.rm = TRUE
    ),

    ca_sma_no = sum(
      CATCH[
        SPECIES_CODE == "SMA" &
        CATCH_UNIT_CODE == "NO"
      ],
      na.rm = TRUE
    ),

    ca_bsh_records = sum(SPECIES_CODE == "BSH"),

    ca_sma_records = sum(SPECIES_CODE == "SMA"),

    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 5. Summarise effort
#
# Do NOT sum across effort units.
# At this point we record which units and how much information are available.
# -----------------------------------------------------------------------------

ef_summary <- ef_dat |>
  group_by(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    effort_units = paste(
      sort(unique(EFFORT_UNIT_CODE)),
      collapse = ", "
    ),

    n_effort_units = n_distinct(
      EFFORT_UNIT_CODE
    ),

    effort_records = n(),

    has_hooks = any(
      EFFORT_UNIT_CODE == "HOOKS"
    ),

    has_trips = any(
      EFFORT_UNIT_CODE == "TRIPS"
    ),

    has_fdays = any(
      EFFORT_UNIT_CODE == "FDAYS"
    ),

    has_sets = any(
      EFFORT_UNIT_CODE == "SETS"
    ),

    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 6. Merge onto the RC universe
#
# RC defines the Fleet x Fishery x Year strata for which retained catch
# exists. CA and EF provide supporting information.
# -----------------------------------------------------------------------------

reconstruction_strata <- rc_summary |>
  left_join(
    ca_summary,
    by = recon_keys
  ) |>
  left_join(
    ef_summary,
    by = recon_keys
  )


# -----------------------------------------------------------------------------
# 7. Information indicators
#
# Important:
#   NA in CA/EF means that no corresponding information was located.
#   This is NOT converted into a biological zero.
# -----------------------------------------------------------------------------

reconstruction_strata <- reconstruction_strata |>
  mutate(

    # Fishing activity
    has_fishing = rc_total > 0,

    has_target = rc_target > 0,

    # Species-specific RC reporting
    has_rc_bsh = rc_bsh > 0,

    has_rc_sma = rc_sma > 0,

    # Generic shark reporting
    has_shark_nei = rc_shark_nei > 0,

    # Any currently relevant shark reporting
    has_any_shark_reporting =
      rc_bsh > 0 |
      rc_sma > 0 |
      rc_shark_nei > 0,

    # CA information
    has_ca_bsh =
      coalesce(ca_bsh_mt, 0) > 0 |
      coalesce(ca_bsh_no, 0) > 0,

    has_ca_sma =
      coalesce(ca_sma_mt, 0) > 0 |
      coalesce(ca_sma_no, 0) > 0,

    # Effort information
    has_effort =
      coalesce(effort_records, 0L) > 0
  )


# -----------------------------------------------------------------------------
# 8. Blue shark information state
#
# These states describe what is known.
# They do NOT yet decide how much catch should be added.
# -----------------------------------------------------------------------------

reconstruction_strata <- reconstruction_strata |>
  mutate(

    bsh_reporting_state = case_when(

      !has_fishing ~
        "R0_no_fishing",

      has_rc_bsh &
        has_ca_bsh &
        has_effort ~
        "R5_BSH_RC_CA_EF",

      has_rc_bsh &
        has_ca_bsh ~
        "R4_BSH_RC_CA",

      has_rc_bsh ~
        "R3_BSH_reported_RC",

      !has_rc_bsh &
        has_shark_nei ~
        "R2_generic_shark_reported",

      !has_rc_bsh &
        !has_shark_nei &
        has_target ~
        "R1_target_fishing_no_shark_reporting",

      TRUE ~
        "R6_other_fishing_no_BSH_information"
    )
  )


# -----------------------------------------------------------------------------
# 9. Candidate reconstruction problem
#
# This is intentionally descriptive rather than a final TRUE/FALSE decision.
# -----------------------------------------------------------------------------

reconstruction_strata <- reconstruction_strata |>
  mutate(

    bsh_problem = case_when(

      bsh_reporting_state == "R2_generic_shark_reported" ~
        "species_disaggregation",

      bsh_reporting_state == "R1_target_fishing_no_shark_reporting" ~
        "potential_non_reporting",

      bsh_reporting_state %in%
        c(
          "R3_BSH_reported_RC",
          "R4_BSH_RC_CA",
          "R5_BSH_RC_CA_EF"
        ) ~
        "reported_BSH",

      TRUE ~
        "insufficient_or_low_priority"
    )
  )


# -----------------------------------------------------------------------------
# 10. Simple BSH-to-target ratio
#
# Only calculated where both BSH and target catch are positive.
# This is NOT yet the reconstructed ratio estimator.
# -----------------------------------------------------------------------------

reconstruction_strata <- reconstruction_strata |>
  mutate(

    bsh_target_ratio = if_else(
      rc_bsh > 0 & rc_target > 0,
      rc_bsh / rc_target,
      NA_real_
    ),

    sma_target_ratio = if_else(
      rc_sma > 0 & rc_target > 0,
      rc_sma / rc_target,
      NA_real_
    )
  )


# -----------------------------------------------------------------------------
# 11. Arrange output
# -----------------------------------------------------------------------------

reconstruction_strata <- reconstruction_strata |>
  select(

    YEAR,

    FLEET_CODE,
    FLEET,

    FISHERY_GROUP_CODE,
    FISHERY_GROUP,

    FISHERY_CODE,
    FISHERY,

    # Fishing activity
    rc_total,
    rc_target,
    rc_yft,
    rc_bet,
    rc_skj,
    rc_alb,
    rc_swo,

    # Sharks
    rc_bsh,
    rc_sma,
    rc_shark_nei,

    # CA
    ca_bsh_mt,
    ca_bsh_no,
    ca_sma_mt,
    ca_sma_no,
    ca_bsh_records,
    ca_sma_records,

    # EF
    effort_units,
    n_effort_units,
    effort_records,
    has_hooks,
    has_trips,
    has_fdays,
    has_sets,

    # Information
    has_target,
    has_rc_bsh,
    has_rc_sma,
    has_shark_nei,
    has_any_shark_reporting,
    has_ca_bsh,
    has_ca_sma,
    has_effort,

    # Classification
    bsh_reporting_state,
    bsh_problem,

    # Initial diagnostic ratios
    bsh_target_ratio,
    sma_target_ratio
  ) |>
  arrange(
    YEAR,
    FISHERY_GROUP_CODE,
    FLEET_CODE,
    FISHERY_CODE
  )

reconstruction_strata |>
  count(
    bsh_reporting_state,
    sort = TRUE
  )

bsh_state_summary <- reconstruction_strata |>
  group_by(bsh_reporting_state) |>
  summarise(
    n_strata = n(),
    total_target_mt = sum(rc_target, na.rm = TRUE),
    total_bsh_mt = sum(rc_bsh, na.rm = TRUE),
    total_shark_nei_mt = sum(rc_shark_nei, na.rm = TRUE),
    first_year = min(YEAR),
    last_year = max(YEAR),
    .groups = "drop"
  ) |>
  arrange(desc(total_target_mt))

bsh_state_summary


early_ll <- reconstruction_strata |>
  filter(
    FISHERY_GROUP_CODE == "LL",
    YEAR <= 1980
  ) |>
  arrange(
    YEAR,
    FLEET_CODE,
    FISHERY_CODE
  )

early_ll_summary <- early_ll |>
  group_by(bsh_reporting_state) |>
  summarise(
    n_strata = n(),
    target_mt = sum(rc_target, na.rm = TRUE),
    bsh_mt = sum(rc_bsh, na.rm = TRUE),
    shark_nei_mt = sum(rc_shark_nei, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(target_mt))

early_ll_summary


early_ll_fleets <- early_ll |>
  group_by(
    FLEET_CODE,
    FLEET
  ) |>
  summarise(
    first_year = min(YEAR),
    last_year = max(YEAR),
    target_mt = sum(rc_target, na.rm = TRUE),
    bsh_mt = sum(rc_bsh, na.rm = TRUE),
    shark_nei_mt = sum(rc_shark_nei, na.rm = TRUE),
    n_years = n_distinct(YEAR),
    n_years_bsh = n_distinct(YEAR[rc_bsh > 0]),
    .groups = "drop"
  ) |>
  arrange(desc(target_mt))

print(early_ll_fleets, n = Inf)