# =============================================================================
# compare_rc_ca_sma_composition.R
#
# Compare shortfin mako species composition between:
#
#   RC = retained catch
#   CA = georeferenced catch
#
# Two composition metrics are evaluated:
#
#   1. SMA proportion of identified shark catch
#
#        SMA / identified sharks
#
#      Used for SKH -> SMA allocation.
#
#   2. SMA proportion of identified mako catch
#
#        SMA / (SMA + LMA)
#
#      Used for MAK -> SMA allocation.
#
# Comparison is performed on matched:
#
#   Fleet x Fishery x Year
#
# strata.
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


safe_cor <- function(x, y) {

  ok <- complete.cases(x, y)

  x <- x[ok]
  y <- y[ok]

  if (
    length(x) < 3 ||
    sd(x) == 0 ||
    sd(y) == 0
  ) {
    return(NA_real_)
  }

  cor(
    x,
    y,
    method = "spearman"
  )
}

# =============================================================================
# 2. RC Fleet x Fishery x Year composition
# =============================================================================

sma_rc_composition <- rc_dat |>
  add_reporting_period() |>
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

    rc_sma_mt =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    rc_lma_mt =
      sum(
        CATCH[SPECIES_CODE == "LMA"],
        na.rm = TRUE
      ),

    rc_mak_mt =
      sum(
        CATCH[SPECIES_CODE == "MAK"],
        na.rm = TRUE
      ),

    rc_skh_mt =
      sum(
        CATCH[SPECIES_CODE == "SKH"],
        na.rm = TRUE
      ),

    rc_identified_shark_mt =
      sum(
        CATCH[
          SPECIES_CODE %in%
            identified_shark_codes
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(

    rc_identified_mako_mt =
      rc_sma_mt +
      rc_lma_mt,

    rc_p_sma_shark =
      if_else(
        rc_identified_shark_mt > 0,
        rc_sma_mt /
          rc_identified_shark_mt,
        NA_real_
      ),

    rc_p_sma_mako =
      if_else(
        rc_identified_mako_mt > 0,
        rc_sma_mt /
          rc_identified_mako_mt,
        NA_real_
      )
  )


# =============================================================================
# 3. CA Fleet x Fishery x Year composition
#
# Use MT only so the denominator is directly comparable with RC.
# =============================================================================

sma_ca_composition <- ca_dat |>
  add_reporting_period() |>
  filter(
    CATCH_UNIT_CODE == "MT"
  ) |>
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

    ca_sma_mt =
      sum(
        CATCH[SPECIES_CODE == "SMA"],
        na.rm = TRUE
      ),

    ca_lma_mt =
      sum(
        CATCH[SPECIES_CODE == "LMA"],
        na.rm = TRUE
      ),

    ca_mak_mt =
      sum(
        CATCH[SPECIES_CODE == "MAK"],
        na.rm = TRUE
      ),

    ca_skh_mt =
      sum(
        CATCH[SPECIES_CODE == "SKH"],
        na.rm = TRUE
      ),

    ca_identified_shark_mt =
      sum(
        CATCH[
          SPECIES_CODE %in%
            identified_shark_codes
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  mutate(

    ca_identified_mako_mt =
      ca_sma_mt +
      ca_lma_mt,

    ca_p_sma_shark =
      if_else(
        ca_identified_shark_mt > 0,
        ca_sma_mt /
          ca_identified_shark_mt,
        NA_real_
      ),

    ca_p_sma_mako =
      if_else(
        ca_identified_mako_mt > 0,
        ca_sma_mt /
          ca_identified_mako_mt,
        NA_real_
      )
  )


# =============================================================================
# 4. Match RC and CA strata
# =============================================================================

comparison_keys <- c(
  "YEAR",
  "FLEET_CODE",
  "FISHERY_GROUP_CODE",
  "FISHERY_CODE"
)


sma_rc_ca_comparison <- sma_rc_composition |>

  select(
    all_of(comparison_keys),

    period,
    FLEET,
    FISHERY_GROUP,
    FISHERY,

    rc_sma_mt,
    rc_lma_mt,
    rc_mak_mt,
    rc_skh_mt,

    rc_identified_shark_mt,
    rc_identified_mako_mt,

    rc_p_sma_shark,
    rc_p_sma_mako
  ) |>

  inner_join(

    sma_ca_composition |>
      select(
        all_of(comparison_keys),

        ca_sma_mt,
        ca_lma_mt,
        ca_mak_mt,
        ca_skh_mt,

        ca_identified_shark_mt,
        ca_identified_mako_mt,

        ca_p_sma_shark,
        ca_p_sma_mako
      ),

    by = comparison_keys
  )


# =============================================================================
# 5. Shark-composition comparison
# =============================================================================

sma_shark_comp_comparison <- sma_rc_ca_comparison |>
  filter(
    !is.na(rc_p_sma_shark),
    !is.na(ca_p_sma_shark)
  ) |>
  mutate(

    difference_p_sma_shark =
      ca_p_sma_shark -
      rc_p_sma_shark,

    abs_difference_p_sma_shark =
      abs(
        difference_p_sma_shark
      )
  )


# =============================================================================
# 6. Mako-composition comparison
# =============================================================================

sma_mako_comp_comparison <- sma_rc_ca_comparison |>
  filter(
    !is.na(rc_p_sma_mako),
    !is.na(ca_p_sma_mako)
  ) |>
  mutate(

    difference_p_sma_mako =
      ca_p_sma_mako -
      rc_p_sma_mako,

    abs_difference_p_sma_mako =
      abs(
        difference_p_sma_mako
      )
  )


# =============================================================================
# 7. Agreement by fishery group: shark composition
# =============================================================================

sma_shark_agreement_by_group <-
  sma_shark_comp_comparison |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    n_matched_strata =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    median_rc_p_sma =
      median(
        rc_p_sma_shark,
        na.rm = TRUE
      ),

    median_ca_p_sma =
      median(
        ca_p_sma_shark,
        na.rm = TRUE
      ),

    median_difference =
      median(
        difference_p_sma_shark,
        na.rm = TRUE
      ),

    median_absolute_difference =
      median(
        abs_difference_p_sma_shark,
        na.rm = TRUE
      ),

    correlation =
  safe_cor(
    rc_p_sma_mako,
    ca_p_sma_mako
  ),

    .groups = "drop"
  )


print(
  sma_shark_agreement_by_group,
  n = Inf
)


# =============================================================================
# 8. Agreement by fishery group: mako composition
# =============================================================================

sma_mako_agreement_by_group <-
  sma_mako_comp_comparison |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    n_matched_strata =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    first_year =
      min(YEAR),

    last_year =
      max(YEAR),

    median_rc_p_sma =
      median(
        rc_p_sma_mako,
        na.rm = TRUE
      ),

    median_ca_p_sma =
      median(
        ca_p_sma_mako,
        na.rm = TRUE
      ),

    median_difference =
      median(
        difference_p_sma_mako,
        na.rm = TRUE
      ),

    median_absolute_difference =
      median(
        abs_difference_p_sma_mako,
        na.rm = TRUE
      ),

  correlation =
  safe_cor(
    rc_p_sma_mako,
    ca_p_sma_mako
  ),

    .groups = "drop"
  )


print(
  sma_mako_agreement_by_group,
  n = Inf
)


# =============================================================================
# 9. Agreement by group x period
# =============================================================================

sma_shark_agreement_by_group_period <-
  sma_shark_comp_comparison |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_matched_strata =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    median_rc_p_sma =
      median(
        rc_p_sma_shark,
        na.rm = TRUE
      ),

    median_ca_p_sma =
      median(
        ca_p_sma_shark,
        na.rm = TRUE
      ),

    median_difference =
      median(
        difference_p_sma_shark,
        na.rm = TRUE
      ),

    median_absolute_difference =
      median(
        abs_difference_p_sma_shark,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


sma_mako_agreement_by_group_period <-
  sma_mako_comp_comparison |>
  group_by(
    FISHERY_GROUP_CODE,
    period
  ) |>
  summarise(

    n_matched_strata =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    median_rc_p_sma =
      median(
        rc_p_sma_mako,
        na.rm = TRUE
      ),

    median_ca_p_sma =
      median(
        ca_p_sma_mako,
        na.rm = TRUE
      ),

    median_difference =
      median(
        difference_p_sma_mako,
        na.rm = TRUE
      ),

    median_absolute_difference =
      median(
        abs_difference_p_sma_mako,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


# =============================================================================
# 10. RC vs CA: SMA proportion of identified sharks
# =============================================================================

p_sma_rc_ca_sharks <- ggplot(
  sma_shark_comp_comparison,
  aes(
    x = rc_p_sma_shark,
    y = ca_p_sma_shark
  )
) +
  geom_point(
    alpha = 0.45
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
    x =
      "SMA proportion of identified shark catch: RC",

    y =
      "SMA proportion of identified shark catch: CA",

    title =
      "Comparison of shortfin mako shark composition",

    subtitle =
      "Matched fleet × fishery × year strata"
  ) +
  theme_bw()

p_sma_rc_ca_sharks


# =============================================================================
# 11. RC vs CA: SMA proportion of identified makos
# =============================================================================

p_sma_rc_ca_makos <- ggplot(
  sma_mako_comp_comparison,
  aes(
    x = rc_p_sma_mako,
    y = ca_p_sma_mako
  )
) +
  geom_point(
    alpha = 0.45
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
    x =
      "SMA proportion of identified mako catch: RC",

    y =
      "SMA proportion of identified mako catch: CA",

    title =
      "Comparison of shortfin mako composition within identified makos",

    subtitle =
      "Matched fleet × fishery × year strata"
  ) +
  theme_bw()

p_sma_rc_ca_makos


# =============================================================================
# 12. Difference through time: shark composition
# =============================================================================

p_sma_shark_difference <- ggplot(
  sma_shark_comp_comparison,
  aes(
    x = YEAR,
    y = difference_p_sma_shark
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point(
    alpha = 0.45
  ) +
  facet_wrap(
    ~ FISHERY_GROUP_CODE
  ) +
  labs(
    x = "Year",
    y = "CA - RC SMA proportion",
    title =
      "Difference between CA and RC shortfin mako shark composition",
    subtitle =
      "Positive values indicate a higher SMA proportion in CA"
  ) +
  theme_bw()

p_sma_shark_difference


# =============================================================================
# 13. Difference through time: mako composition
# =============================================================================

p_sma_mako_difference <- ggplot(
  sma_mako_comp_comparison,
  aes(
    x = YEAR,
    y = difference_p_sma_mako
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point(
    alpha = 0.45
  ) +
  facet_wrap(
    ~ FISHERY_GROUP_CODE
  ) +
  labs(
    x = "Year",
    y = "CA - RC SMA proportion within makos",
    title =
      "Difference between CA and RC mako composition",
    subtitle =
      "Positive values indicate a higher SMA fraction in CA"
  ) +
  theme_bw()

p_sma_mako_difference


# =============================================================================
# 14. Longline-class diagnostic
#
# Particularly important for SMA because LLO differs substantially from
# deep-freezing and fresh longline fisheries.
# =============================================================================

add_longline_class <- function(dat) {

  dat |>
    mutate(

      longline_class = case_when(

        FISHERY_GROUP_CODE != "LL" ~
          NA_character_,

        FISHERY_CODE == "LLD" ~
          "LL_deep_freezing",

        FISHERY_CODE == "LLF" ~
          "LL_fresh",

        FISHERY_CODE == "LLO" ~
          "LL_other",

        TRUE ~
          "LL_other_unspecified"
      )
    )
}


sma_ll_rc_ca <- sma_rc_ca_comparison |>
  add_longline_class() |>
  filter(
    FISHERY_GROUP_CODE == "LL"
  )


sma_ll_composition_summary <- sma_ll_rc_ca |>
  group_by(
    longline_class,
    period
  ) |>
  summarise(

    n_strata =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    rc_sma_mt =
      sum(
        rc_sma_mt,
        na.rm = TRUE
      ),

    ca_sma_mt =
      sum(
        ca_sma_mt,
        na.rm = TRUE
      ),

    median_rc_p_sma_shark =
      median(
        rc_p_sma_shark,
        na.rm = TRUE
      ),

    median_ca_p_sma_shark =
      median(
        ca_p_sma_shark,
        na.rm = TRUE
      ),

    median_rc_p_sma_mako =
      median(
        rc_p_sma_mako,
        na.rm = TRUE
      ),

    median_ca_p_sma_mako =
      median(
        ca_p_sma_mako,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

print(
  sma_ll_composition_summary,
  n = Inf
)


# =============================================================================
# 15. Save
# =============================================================================

dir.create(
  "outputs/sma",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  sma_shark_agreement_by_group,
  "outputs/sma/sma_rc_ca_shark_agreement_by_group.csv"
)

write_csv(
  sma_mako_agreement_by_group,
  "outputs/sma/sma_rc_ca_mako_agreement_by_group.csv"
)

write_csv(
  sma_shark_agreement_by_group_period,
  "outputs/sma/sma_rc_ca_shark_agreement_by_group_period.csv"
)

write_csv(
  sma_mako_agreement_by_group_period,
  "outputs/sma/sma_rc_ca_mako_agreement_by_group_period.csv"
)

write_csv(
  sma_ll_composition_summary,
  "outputs/sma/sma_rc_ca_longline_class_summary.csv"
)