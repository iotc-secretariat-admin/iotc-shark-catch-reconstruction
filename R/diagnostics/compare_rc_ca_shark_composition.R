# =============================================================================
# compare_rc_ca_shark_composition.R
#
# Compare shark species composition inferred from:
#
#   RC = annual retained catch
#   CA = geo-referenced retained catch
#
# PURPOSE
# -------
# Evaluate whether the species composition observed in the geo-referenced
# catch data is consistent with the composition in the retained catch data.
#
# This script is diagnostic only. It does NOT reconstruct catch.
#
# Initial comparison:
#   - identified shark species only
#   - CA records reported in metric tonnes (MT)
#   - BSH proportion of identified shark catch
#   - comparisons at several levels of aggregation
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


# =============================================================================
# 1. Species definitions
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
  "SPY",  # Hammerhead sharks, etc. nei
  "SPZ",  # Smooth hammerhead
  "THR"   # Thresher sharks nei
)

generic_shark_codes <- c(
  "SKH",
  "SHRK",
  "AG22",
  "SKX",
  "SKHS"
)
# =============================================================================
# 2. Reporting periods
# =============================================================================

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
      )
    )
}


# =============================================================================
# 3. RC species-specific shark catch
# =============================================================================

rc_shark_species <- rc_dat |>
  filter(
    SPECIES_CODE %in% identified_shark_codes
  ) |>
  add_reporting_period() |>
  select(
    YEAR,
    period,

    FLEET_CODE,
    FLEET,

    FISHERY_GROUP_CODE,
    FISHERY_GROUP,

    FISHERY_CODE,
    FISHERY,

    GEAR_CODE,
    GEAR,

    SPECIES_CODE,
    SPECIES,

    catch_mt = CATCH
  )

# =============================================================================
# 4. CA species-specific shark catch
#
# Use metric-tonne records only.
# =============================================================================

ca_shark_species <- ca_dat |>
  filter(
    SPECIES_CODE %in% identified_shark_codes,
    CATCH_UNIT_CODE == "MT"
  ) |>
  add_reporting_period() |>
  select(
    YEAR,
    period,

    FLEET_CODE,
    FLEET,

    FISHERY_GROUP_CODE,
    FISHERY_GROUP,

    FISHERY_CODE,
    FISHERY,

    GEAR_CODE,
    GEAR,

    SPECIES_CODE,
    SPECIES,

    catch_mt = CATCH
  )


# =============================================================================
# 5. Composition by Fishery Group x Period
# =============================================================================

rc_group_period <- rc_shark_species |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period,
    SPECIES_CODE,
    SPECIES
  ) |>
  summarise(
    catch_mt = sum(catch_mt, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period
  ) |>
  mutate(
    identified_shark_mt =
      sum(catch_mt, na.rm = TRUE),

    proportion =
      catch_mt / identified_shark_mt
  ) |>
  ungroup()


ca_group_period <- ca_shark_species |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period,
    SPECIES_CODE,
    SPECIES
  ) |>
  summarise(
    catch_mt = sum(catch_mt, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period
  ) |>
  mutate(
    identified_shark_mt =
      sum(catch_mt, na.rm = TRUE),

    proportion =
      catch_mt / identified_shark_mt
  ) |>
  ungroup()


# =============================================================================
# 6. Compare BSH fraction directly
# =============================================================================

rc_bsh_group_period <- rc_group_period |>
  filter(SPECIES_CODE == "BSH") |>
  select(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period,

    rc_identified_shark_mt =
      identified_shark_mt,

    rc_bsh_mt =
      catch_mt,

    rc_p_bsh =
      proportion
  )


ca_bsh_group_period <- ca_group_period |>
  filter(SPECIES_CODE == "BSH") |>
  select(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period,

    ca_identified_shark_mt =
      identified_shark_mt,

    ca_bsh_mt =
      catch_mt,

    ca_p_bsh =
      proportion
  )
 

# Join


bsh_group_period_comparison <-
  full_join(
    rc_bsh_group_period,
    ca_bsh_group_period,
    by = c(
      "FISHERY_GROUP_CODE",
      "FISHERY_GROUP",
      "period"
    )
  ) |>
  mutate(

    difference_p_bsh =
      ca_p_bsh - rc_p_bsh,

    abs_difference_p_bsh =
      abs(difference_p_bsh),

    ratio_ca_rc =
      if_else(
        rc_p_bsh > 0,
        ca_p_bsh / rc_p_bsh,
        NA_real_
      )
  ) |>
  arrange(
    FISHERY_GROUP_CODE,
    period
  )

print(
  bsh_group_period_comparison,
  n = Inf
)


# =============================================================================
# 7. Fleet x Fishery x Year compositions
# =============================================================================

rc_year_comp <- rc_shark_species |>
  group_by(
    YEAR,
    period,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    identified_shark_mt =
      sum(catch_mt, na.rm = TRUE),

    bsh_mt =
      sum(
        catch_mt[SPECIES_CODE == "BSH"],
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        catch_mt[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(
    p_bsh =
      if_else(
        identified_shark_mt > 0,
        bsh_mt / identified_shark_mt,
        NA_real_
      ),

    p_sma =
      if_else(
        identified_shark_mt > 0,
        sma_mt / identified_shark_mt,
        NA_real_
      )
  )


ca_year_comp <- ca_shark_species |>
  group_by(
    YEAR,
    period,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE
  ) |>
  summarise(

    identified_shark_mt =
      sum(catch_mt, na.rm = TRUE),

    bsh_mt =
      sum(
        catch_mt[SPECIES_CODE == "BSH"],
        na.rm = TRUE
      ),

    sma_mt =
      sum(
        catch_mt[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(
    p_bsh =
      if_else(
        identified_shark_mt > 0,
        bsh_mt / identified_shark_mt,
        NA_real_
      ),

    p_sma =
      if_else(
        identified_shark_mt > 0,
        sma_mt / identified_shark_mt,
        NA_real_
      )
  )
# inner join to check 
# =============================================================================
# 8. Matched RC–CA strata
# note that identified_shark_coverage variable is useful, but it   is simply:
# CA_identified sharks / RC_ Identified sharks 
# =============================================================================

rc_ca_matched <- rc_year_comp |>
  inner_join(
    ca_year_comp,
    by = c(
      "YEAR",
      "period",
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "FISHERY_CODE"
    ),
    suffix = c("_rc", "_ca")
  ) |>
  mutate(

    p_bsh_difference =
      p_bsh_ca - p_bsh_rc,

    abs_p_bsh_difference =
      abs(p_bsh_difference),

    identified_shark_coverage =
      if_else(
        identified_shark_mt_rc > 0,
        identified_shark_mt_ca /
          identified_shark_mt_rc,
        NA_real_
      )
  )


# =============================================================================
# 9. Agreement by broad fishery group
# =============================================================================

rc_ca_agreement_by_group <- rc_ca_matched |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    n_matched_strata = n(),

    n_fleets =
      n_distinct(FLEET_CODE),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    median_rc_p_bsh =
      median(
        p_bsh_rc,
        na.rm = TRUE
      ),

    median_ca_p_bsh =
      median(
        p_bsh_ca,
        na.rm = TRUE
      ),

    median_difference =
      median(
        p_bsh_difference,
        na.rm = TRUE
      ),

    median_absolute_difference =
      median(
        abs_p_bsh_difference,
        na.rm = TRUE
      ),

    correlation =
      if_else(
        sum(
          complete.cases(
            p_bsh_rc,
            p_bsh_ca
          )
        ) >= 3,

        cor(
          p_bsh_rc,
          p_bsh_ca,
          use = "complete.obs",
          method = "spearman"
        ),

        NA_real_
      ),

    .groups = "drop"
  ) |>
  arrange(
    desc(n_matched_strata)
  )

print(
  rc_ca_agreement_by_group,
  n = Inf
)



# =============================================================================
# 10. Agreement by fishery group x period
# =============================================================================

rc_ca_agreement_by_group_period <- rc_ca_matched |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_matched_strata = n(),

    n_fleets =
      n_distinct(FLEET_CODE),

    median_rc_p_bsh =
      median(
        p_bsh_rc,
        na.rm = TRUE
      ),

    median_ca_p_bsh =
      median(
        p_bsh_ca,
        na.rm = TRUE
      ),

    median_difference =
      median(
        p_bsh_ca - p_bsh_rc,
        na.rm = TRUE
      ),

    median_absolute_difference =
      median(
        abs(p_bsh_ca - p_bsh_rc),
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  arrange(
    FISHERY_GROUP_CODE,
    period
  )

print(
  rc_ca_agreement_by_group_period,
  n = Inf
)


# =============================================================================
# 11. RC vs CA BSH composition
# =============================================================================

p_rc_ca_bsh <- ggplot(
  rc_ca_matched,
  aes(
    x = p_bsh_rc,
    y = p_bsh_ca
  )
) +
  geom_point(
    alpha = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  facet_wrap(
    ~ FISHERY_GROUP_CODE
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    x = "BSH proportion of identified shark catch: RC",
    y = "BSH proportion of identified shark catch: CA",
    title = "Comparison of blue shark species composition",
    subtitle = "Matched fleet × fishery × year strata"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )

p_rc_ca_bsh


# =============================================================================
# 12. RC–CA composition differences through time
# =============================================================================

p_rc_ca_difference <- ggplot(
  rc_ca_matched,
  aes(
    x = YEAR,
    y = p_bsh_difference
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point(
    alpha = 0.5
  ) +
  facet_wrap(
    ~ FISHERY_GROUP_CODE,
    scales = "free_x"
  ) +
  labs(
    x = "Year",
    y = "CA BSH proportion - RC BSH proportion",
    title = "Difference between CA and RC blue shark composition",
    subtitle = "Positive values indicate a higher BSH proportion in CA"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )

p_rc_ca_difference


# =============================================================================
# 13. Full shark species-composition comparison
# =============================================================================

rc_full_comp <- rc_group_period |>
  select(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period,
    SPECIES_CODE,
    SPECIES,
    rc_proportion = proportion,
    rc_catch_mt = catch_mt
  )


ca_full_comp <- ca_group_period |>
  select(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    period,
    SPECIES_CODE,
    SPECIES,
    ca_proportion = proportion,
    ca_catch_mt = catch_mt
  )


rc_ca_full_composition <- full_join(
  rc_full_comp,
  ca_full_comp,
  by = c(
    "FISHERY_GROUP_CODE",
    "FISHERY_GROUP",
    "period",
    "SPECIES_CODE",
    "SPECIES"
  )
) |>
  mutate(

    rc_proportion =
      coalesce(rc_proportion, 0),

    ca_proportion =
      coalesce(ca_proportion, 0),

    difference =
      ca_proportion - rc_proportion
  )

# check ll 
rc_ca_full_composition |>
  filter(
    FISHERY_GROUP_CODE == "LL"
  ) |>
  arrange(
    period,
    desc(rc_proportion)
  ) |>
  print(n = Inf)
#
# Check gillnet
rc_ca_full_composition |>
  filter(
    FISHERY_GROUP_CODE == "GN"
  ) |>
  arrange(
    period,
    desc(rc_proportion)
  ) |>
  print(n = Inf)

rc_ca_agreement_by_group
rc_ca_agreement_by_group_period
bsh_group_period_comparison
