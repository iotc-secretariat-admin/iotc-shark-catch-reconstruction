# =============================================================================
# sma_gam_species_composition.R
#
# Statistical species-composition reconstruction for shortfin mako (SMA)
#
# Two distinct taxonomic allocation problems are modelled:
#
#   A. SKH -> SMA
#
#      Estimate:
#
#          P(SMA | identified shark catch)
#
#      This provides a statistical alternative to the HSCR allocation of
#      Sharks NEI (SKH).
#
#
#   B. MAK -> SMA
#
#      Estimate:
#
#          P(SMA | SMA + LMA)
#
#      This provides a statistical alternative to the HSCR allocation of
#      Mako sharks NEI (MAK).
#
#
# IMPORTANT
# ---------
# Complete non-reporting reconstruction is NOT handled here.
#
# That component remains in:
#
#   sma_hscr_non_reporting.R
#
# and is based on SMA / target-catch donor ratios.
#
# =============================================================================


# =============================================================================
# 1. Packages
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(readr)
library(purrr)
library(tibble)
library(viridis)


# =============================================================================
# 2. Output directories
# =============================================================================

dir.create(
  "outputs/sma/gam",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/sma/gam/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. Definitions
# =============================================================================

period_levels <- c(
  "1950_1963",
  "1964_1980",
  "1981_1994",
  "1995_2004",
  "2005_2014",
  "2015_2025"
)


# ---------------------------------------------------------------------------
# Identified shark species
#
# IMPORTANT:
# Keep this list synchronized with the list used in the BSH reconstruction.
# Generic categories such as SKH and MAK are intentionally excluded.
# ---------------------------------------------------------------------------

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


target_species <- c(
  "YFT",
  "BET",
  "SKJ",
  "ALB",
  "SWO",
  "SBF"
)


# =============================================================================
# 4. Helper functions
# =============================================================================


# ---------------------------------------------------------------------------
# Reporting period
# ---------------------------------------------------------------------------

add_reporting_period <- function(dat) {

  dat |>
    mutate(

      period = case_when(

        YEAR <= 1963 ~
          "1950_1963",

        YEAR <= 1980 ~
          "1964_1980",

        YEAR <= 1994 ~
          "1981_1994",

        YEAR <= 2004 ~
          "1995_2004",

        YEAR <= 2014 ~
          "2005_2014",

        TRUE ~
          "2015_2025"
      ),

      period = factor(
        period,
        levels = period_levels,
        ordered = TRUE
      )
    )
}


# ---------------------------------------------------------------------------
# Longline class
#
# This mirrors the classes now used in the deterministic SMA reconstruction.
# ---------------------------------------------------------------------------

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

        FISHERY_CODE == "SLL" ~
          "LL_shark_target",

        FISHERY_CODE == "LLO" ~
          "LL_other",

        TRUE ~
          "LL_other_unspecified"
      )
    )
}


# ---------------------------------------------------------------------------
# Combined fishery modelling class
#
# LL is split into the longline classes established during the deterministic
# reconstruction. Other fisheries retain their broad fishery-group identity.
# ---------------------------------------------------------------------------

add_model_fishery_class <- function(dat) {

  dat |>
    add_longline_class() |>
    mutate(

      model_fishery_class = case_when(

        FISHERY_GROUP_CODE == "LL" ~
          longline_class,

        TRUE ~
          FISHERY_GROUP_CODE
      ),

      model_fishery_class =
        factor(
          model_fishery_class
        )
    )
}


# ---------------------------------------------------------------------------
# Move proportions away from exact 0 and 1 for beta regression.
#
# The transformation is deliberately small. Exact boundaries remain important
# biologically and are retained in the raw response variables.
# ---------------------------------------------------------------------------

squeeze_beta <- function(x, eps = 0.001) {

  pmin(
    pmax(
      x,
      eps
    ),
    1 - eps
  )
}


# ---------------------------------------------------------------------------
# Safe correlation for CV diagnostics
# ---------------------------------------------------------------------------

safe_cor <- function(x, y) {

  ok <- complete.cases(
    x,
    y
  )

  x <- x[ok]
  y <- y[ok]

  if (
    length(x) < 3 ||
    sd(x) == 0 ||
    sd(y) == 0
  ) {
    return(
      NA_real_
    )
  }

  cor(
    x,
    y,
    method = "spearman"
  )
}


# =============================================================================
# 5. Construct RC composition data
# =============================================================================

sma_gam_base <- rc_dat |>

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

    # -------------------------------------------------------------------------
    # Shortfin and longfin mako
    # -------------------------------------------------------------------------

    sma_mt =
      sum(
        CATCH[
          SPECIES_CODE == "SMA"
        ],
        na.rm = TRUE
      ),

    lma_mt =
      sum(
        CATCH[
          SPECIES_CODE == "LMA"
        ],
        na.rm = TRUE
      ),

    mak_mt =
      sum(
        CATCH[
          SPECIES_CODE == "MAK"
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


    # -------------------------------------------------------------------------
    # Identified shark catch
    # -------------------------------------------------------------------------

    identified_shark_mt =
      sum(
        CATCH[
          SPECIES_CODE %in%
            identified_shark_codes
        ],
        na.rm = TRUE
      ),


    # -------------------------------------------------------------------------
    # Target catch
    # -------------------------------------------------------------------------

    tropical_mt =
      sum(
        CATCH[
          SPECIES_CODE %in%
            c(
              "YFT",
              "BET",
              "SKJ"
            )
        ],
        na.rm = TRUE
      ),

    alb_mt =
      sum(
        CATCH[
          SPECIES_CODE == "ALB"
        ],
        na.rm = TRUE
      ),

    swo_mt =
      sum(
        CATCH[
          SPECIES_CODE == "SWO"
        ],
        na.rm = TRUE
      ),

    sbf_mt =
      sum(
        CATCH[
          SPECIES_CODE == "SBF"
        ],
        na.rm = TRUE
      ),

    target_mt =
      sum(
        CATCH[
          SPECIES_CODE %in%
            target_species
        ],
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>

  mutate(

    identified_mako_mt =
      sma_mt +
      lma_mt,


    # -------------------------------------------------------------------------
    # Response A:
    # SMA proportion of identified shark catch
    # -------------------------------------------------------------------------

    p_sma_shark = case_when(

      identified_shark_mt > 0 ~
        sma_mt /
          identified_shark_mt,

      TRUE ~
        NA_real_
    ),


    # -------------------------------------------------------------------------
    # Response B:
    # SMA proportion of identified mako catch
    # -------------------------------------------------------------------------

    p_sma_mako = case_when(

      identified_mako_mt > 0 ~
        sma_mt /
          identified_mako_mt,

      TRUE ~
        NA_real_
    ),


    # -------------------------------------------------------------------------
    # Composition of target catch
    #
    # These predictors describe broad targeting characteristics rather than
    # using SMA catch itself to identify shark-oriented fisheries.
    # -------------------------------------------------------------------------

    prop_tropical = case_when(

      target_mt > 0 ~
        tropical_mt /
          target_mt,

      TRUE ~
        0
    ),

    prop_swo = case_when(

      target_mt > 0 ~
        swo_mt /
          target_mt,

      TRUE ~
        0
    ),

    prop_alb = case_when(

      target_mt > 0 ~
        alb_mt /
          target_mt,

      TRUE ~
        0
    ),

    prop_sbf = case_when(

      target_mt > 0 ~
        sbf_mt /
          target_mt,

      TRUE ~
        0
    )
  ) |>

  add_model_fishery_class()


# =============================================================================
# 6. Inspect available training information
# =============================================================================

sma_gam_inventory <- sma_gam_base |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    n_strata =
      n(),

    n_shark_composition =
      sum(
        !is.na(
          p_sma_shark
        )
      ),

    n_mako_composition =
      sum(
        !is.na(
          p_sma_mako
        )
      ),

    n_fleets_shark =
      n_distinct(
        FLEET_CODE[
          !is.na(
            p_sma_shark
          )
        ]
      ),

    n_fleets_mako =
      n_distinct(
        FLEET_CODE[
          !is.na(
            p_sma_mako
          )
        ]
      ),

    first_year =
      min(
        YEAR
      ),

    last_year =
      max(
        YEAR
      ),

    .groups = "drop"
  )


print(
  sma_gam_inventory,
  n = Inf
)


write_csv(
  sma_gam_inventory,
  "outputs/sma/gam/sma_gam_training_inventory.csv"
)


# =============================================================================
# PART A
#
# SKH -> SMA
#
# Model:
#
#     P(SMA | identified sharks)
#
# =============================================================================


# =============================================================================
# 7. SKH-model training data
# =============================================================================

sma_shark_training <- sma_gam_base |>

  filter(

    !is.na(
      p_sma_shark
    ),

    identified_shark_mt > 0,

    # Keep groups with reasonable support.
    #
    # Initially mirror the BSH approach and exclude OT and BB.
    FISHERY_GROUP_CODE %in%
      c(
        "GN",
        "LI",
        "LL",
        "PS"
      )
  ) |>

  mutate(

    p_sma_shark_beta =
      squeeze_beta(
        p_sma_shark
      ),

    FLEET_CODE =
      factor(
        FLEET_CODE
      ),

    FISHERY_GROUP_CODE =
      factor(
        FISHERY_GROUP_CODE
      )
  ) |>

  droplevels()


# =============================================================================
# 8. Diagnostics for response support
# =============================================================================

sma_shark_response_summary <-
  sma_shark_training |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    n =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    n_zero =
      sum(
        p_sma_shark == 0
      ),

    n_one =
      sum(
        p_sma_shark == 1
      ),

    proportion_zero =
      mean(
        p_sma_shark == 0
      ),

    proportion_one =
      mean(
        p_sma_shark == 1
      ),

    median_p =
      median(
        p_sma_shark
      ),

    mean_p =
      mean(
        p_sma_shark
      ),

    .groups = "drop"
  )


print(
  sma_shark_response_summary,
  n = Inf
)


write_csv(
  sma_shark_response_summary,
  "outputs/sma/gam/sma_shark_response_summary.csv"
)


# =============================================================================
# 9. Candidate SKH -> SMA GAMs
#
# Begin with a sequence analogous to the BSH model development.
#
# The final model should be chosen based on:
#
#   - convergence;
#   - explanatory performance;
#   - leave-one-fleet-out prediction;
#   - plausible extrapolation;
#   - stability of reconstructed SKH totals.
# =============================================================================


# ---------------------------------------------------------------------------
# GAM S1
#
# Fishery class + common year smooth
# ---------------------------------------------------------------------------

sma_shark_gam1 <- gam(

  p_sma_shark_beta ~

    model_fishery_class +

    s(
      YEAR,
      k = 8,
      bs = "tp"
    ),

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_shark_training
)


# ---------------------------------------------------------------------------
# GAM S2
#
# Fishery-specific temporal patterns
# ---------------------------------------------------------------------------

sma_shark_gam2 <- gam(

  p_sma_shark_beta ~

    model_fishery_class +

    s(
      YEAR,
      by = FISHERY_GROUP_CODE,
      k = 8,
      bs = "tp"
    ),

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_shark_training
)


# ---------------------------------------------------------------------------
# GAM S3
#
# Add fleet random effect
# ---------------------------------------------------------------------------

sma_shark_gam3 <- gam(

  p_sma_shark_beta ~

    model_fishery_class +

    s(
      YEAR,
      by = FISHERY_GROUP_CODE,
      k = 8,
      bs = "tp"
    ) +

    s(
      FLEET_CODE,
      bs = "re"
    ),

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_shark_training
)


# ---------------------------------------------------------------------------
# GAM S4
#
# Add fishery-targeting composition predictors.
#
# This is the SMA analogue of the final BSH candidate model.
# ---------------------------------------------------------------------------

sma_shark_gam4 <- gam(

  p_sma_shark_beta ~

    model_fishery_class +

    s(
      YEAR,
      by = FISHERY_GROUP_CODE,
      k = 8,
      bs = "tp"
    ) +

    s(
      FLEET_CODE,
      bs = "re"
    ) +

    s(
      prop_tropical,
      k = 5
    ) +

    s(
      prop_swo,
      k = 5
    ) +

    s(
      prop_alb,
      k = 5
    ) +

    s(
      prop_sbf,
      k = 5
    ),

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_shark_training
)


# =============================================================================
# 10. SKH-model comparison
# =============================================================================

sma_shark_model_comparison <- tibble(

  model =
    c(
      "S1",
      "S2",
      "S3",
      "S4"
    ),

  n =
    c(
      nobs(
        sma_shark_gam1
      ),
      nobs(
        sma_shark_gam2
      ),
      nobs(
        sma_shark_gam3
      ),
      nobs(
        sma_shark_gam4
      )
    ),

  deviance_explained =
    c(
      summary(
        sma_shark_gam1
      )$dev.expl,

      summary(
        sma_shark_gam2
      )$dev.expl,

      summary(
        sma_shark_gam3
      )$dev.expl,

      summary(
        sma_shark_gam4
      )$dev.expl
    ),

  adjusted_r2 =
    c(
      summary(
        sma_shark_gam1
      )$r.sq,

      summary(
        sma_shark_gam2
      )$r.sq,

      summary(
        sma_shark_gam3
      )$r.sq,

      summary(
        sma_shark_gam4
      )$r.sq
    ),

  reml_score =
    c(
      sma_shark_gam1$gcv.ubre,
      sma_shark_gam2$gcv.ubre,
      sma_shark_gam3$gcv.ubre,
      sma_shark_gam4$gcv.ubre
    )
)


print(
  sma_shark_model_comparison
)


write_csv(
  sma_shark_model_comparison,
  "outputs/sma/gam/sma_shark_model_comparison.csv"
)


# =============================================================================
# 11. Diagnostic checks for S4
# =============================================================================

summary(
  sma_shark_gam4
)

gam.check(
  sma_shark_gam4
)


# =============================================================================
# 12. Leave-one-fleet-out CV function
#
# IMPORTANT:
# The fleet random effect cannot be meaningfully estimated for the omitted
# fleet. Therefore the CV model excludes the fleet RE and evaluates whether
# the fixed/smooth structure transfers among fleets.
# =============================================================================

fit_predict_sma_shark_lofo <- function(
    fleet_out,
    dat
) {

  train_dat <- dat |>
    filter(
      as.character(
        FLEET_CODE
      ) != fleet_out
    ) |>
    droplevels()


  test_dat <- dat |>
    filter(
      as.character(
        FLEET_CODE
      ) == fleet_out
    )


  if (
    nrow(
      train_dat
    ) < 50 ||
    nrow(
      test_dat
    ) == 0
  ) {

    return(
      NULL
    )
  }


  # Remove levels that are absent from training.

  valid_classes <-
    levels(
      train_dat$model_fishery_class
    )


  test_dat <- test_dat |>
    filter(
      as.character(
        model_fishery_class
      ) %in%
        valid_classes
    )


  if (
    nrow(
      test_dat
    ) == 0
  ) {

    return(
      NULL
    )
  }


  cv_model <- try(

    gam(

      p_sma_shark_beta ~

        model_fishery_class +

        s(
          YEAR,
          by = FISHERY_GROUP_CODE,
          k = 8,
          bs = "tp"
        ) +

        s(
          prop_tropical,
          k = 5
        ) +

        s(
          prop_swo,
          k = 5
        ) +

        s(
          prop_alb,
          k = 5
        ) +

        s(
          prop_sbf,
          k = 5
        ),

      family =
        betar(
          link = "logit"
        ),

      method =
        "REML",

      data =
        train_dat
    ),

    silent = TRUE
  )


  if (
    inherits(
      cv_model,
      "try-error"
    )
  ) {

    return(
      NULL
    )
  }


  test_dat |>
    mutate(

      predicted =
        predict(
          cv_model,
          newdata = test_dat,
          type = "response"
        ),

      observed =
        p_sma_shark,

      residual =
        predicted -
        observed
    )
}


# =============================================================================
# 13. SKH-model leave-one-fleet-out CV
# =============================================================================

sma_shark_lofo <- map_dfr(

  unique(
    as.character(
      sma_shark_training$FLEET_CODE
    )
  ),

  fit_predict_sma_shark_lofo,

  dat =
    sma_shark_training
)


sma_shark_lofo_summary <- sma_shark_lofo |>

  summarise(

    n =
      n(),

    mae =
      mean(
        abs(
          predicted -
            observed
        ),
        na.rm = TRUE
      ),

    rmse =
      sqrt(
        mean(
          (
            predicted -
              observed
          )^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted -
          observed,
        na.rm = TRUE
      ),

    rho =
      safe_cor(
        observed,
        predicted
      )
  )


sma_shark_lofo_by_group <- sma_shark_lofo |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    n =
      n(),

    mae =
      mean(
        abs(
          predicted -
            observed
        ),
        na.rm = TRUE
      ),

    rmse =
      sqrt(
        mean(
          (
            predicted -
              observed
          )^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted -
          observed,
        na.rm = TRUE
      ),

    rho =
      safe_cor(
        observed,
        predicted
      ),

    .groups = "drop"
  )


print(
  sma_shark_lofo_summary
)

print(
  sma_shark_lofo_by_group,
  n = Inf
)


write_csv(
  sma_shark_lofo_summary,
  "outputs/sma/gam/sma_shark_lofo_summary.csv"
)

write_csv(
  sma_shark_lofo_by_group,
  "outputs/sma/gam/sma_shark_lofo_by_group.csv"
)


# =============================================================================
# 14. Plot SKH-model LOFO performance
# =============================================================================

p_sma_shark_lofo <- ggplot(

  sma_shark_lofo,

  aes(
    x = observed,
    y = predicted,
    colour = FISHERY_GROUP_CODE
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

  coord_equal(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +

  scale_colour_viridis_d(
    end = 0.9
  ) +

  labs(

    x =
      "Observed SMA proportion of identified sharks",

    y =
      "Predicted SMA proportion",

    colour =
      "Fishery group",

    title =
      "SMA species-composition GAM",

    subtitle =
      "Leave-one-fleet-out cross-validation"
  ) +

  theme_bw()


p_sma_shark_lofo


ggsave(

  "outputs/sma/gam/figures/sma_shark_gam_lofo.png",

  p_sma_shark_lofo,

  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# PART B
#
# MAK -> SMA
#
# Model:
#
#     P(SMA | SMA + LMA)
#
# =============================================================================


# =============================================================================
# 15. Mako-model training data
# =============================================================================

sma_mako_training <- sma_gam_base |>

  filter(

    !is.na(
      p_sma_mako
    ),

    identified_mako_mt > 0,

    FISHERY_GROUP_CODE %in%
      c(
        "GN",
        "LI",
        "LL",
        "PS"
      )
  ) |>

  mutate(

    p_sma_mako_beta =
      squeeze_beta(
        p_sma_mako
      ),

    FLEET_CODE =
      factor(
        FLEET_CODE
      ),

    FISHERY_GROUP_CODE =
      factor(
        FISHERY_GROUP_CODE
      )
  ) |>

  droplevels()


# =============================================================================
# 16. Mako response diagnostics
# =============================================================================

sma_mako_response_summary <-
  sma_mako_training |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    n =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    n_zero =
      sum(
        p_sma_mako == 0
      ),

    n_one =
      sum(
        p_sma_mako == 1
      ),

    proportion_zero =
      mean(
        p_sma_mako == 0
      ),

    proportion_one =
      mean(
        p_sma_mako == 1
      ),

    median_p =
      median(
        p_sma_mako
      ),

    mean_p =
      mean(
        p_sma_mako
      ),

    .groups = "drop"
  )


print(
  sma_mako_response_summary,
  n = Inf
)


write_csv(
  sma_mako_response_summary,
  "outputs/sma/gam/sma_mako_response_summary.csv"
)


# =============================================================================
# 17. Simple empirical benchmark for MAK allocation
#
# This benchmark is important because SMA appears to dominate identified
# makos in many fisheries. A complex model should improve upon, or at least
# meaningfully modify, this simpler estimator before it is preferred.
# =============================================================================

sma_mako_benchmark <- sma_mako_training |>

  group_by(
    model_fishery_class,
    period
  ) |>

  summarise(

    sma_mt =
      sum(
        sma_mt,
        na.rm = TRUE
      ),

    lma_mt =
      sum(
        lma_mt,
        na.rm = TRUE
      ),

    identified_mako_mt =
      sma_mt +
      lma_mt,

    pooled_p_sma =
      sma_mt /
        identified_mako_mt,

    n =
      n(),

    n_fleets =
      n_distinct(
        FLEET_CODE
      ),

    .groups = "drop"
  )


print(
  sma_mako_benchmark,
  n = Inf
)


write_csv(
  sma_mako_benchmark,
  "outputs/sma/gam/sma_mako_empirical_benchmark.csv"
)


# =============================================================================
# 18. Candidate MAK -> SMA GAMs
#
# Keep these deliberately simpler than the SKH model.
# =============================================================================


# ---------------------------------------------------------------------------
# M1: fishery class only
# ---------------------------------------------------------------------------

sma_mako_gam1 <- gam(

  p_sma_mako_beta ~

    model_fishery_class,

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_mako_training
)


# ---------------------------------------------------------------------------
# M2: fishery class + year
# ---------------------------------------------------------------------------

sma_mako_gam2 <- gam(

  p_sma_mako_beta ~

    model_fishery_class +

    s(
      YEAR,
      k = 6,
      bs = "tp"
    ),

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_mako_training
)


# ---------------------------------------------------------------------------
# M3: fishery class + year + fleet random effect
# ---------------------------------------------------------------------------

sma_mako_gam3 <- gam(

  p_sma_mako_beta ~

    model_fishery_class +

    s(
      YEAR,
      k = 6,
      bs = "tp"
    ) +

    s(
      FLEET_CODE,
      bs = "re"
    ),

  family =
    betar(
      link = "logit"
    ),

  method =
    "REML",

  data =
    sma_mako_training
)


# =============================================================================
# 19. MAK-model comparison
# =============================================================================

sma_mako_model_comparison <- tibble(

  model =
    c(
      "M1",
      "M2",
      "M3"
    ),

  n =
    c(
      nobs(
        sma_mako_gam1
      ),
      nobs(
        sma_mako_gam2
      ),
      nobs(
        sma_mako_gam3
      )
    ),

  deviance_explained =
    c(
      summary(
        sma_mako_gam1
      )$dev.expl,

      summary(
        sma_mako_gam2
      )$dev.expl,

      summary(
        sma_mako_gam3
      )$dev.expl
    ),

  adjusted_r2 =
    c(
      summary(
        sma_mako_gam1
      )$r.sq,

      summary(
        sma_mako_gam2
      )$r.sq,

      summary(
        sma_mako_gam3
      )$r.sq
    ),

  reml_score =
    c(
      sma_mako_gam1$gcv.ubre,
      sma_mako_gam2$gcv.ubre,
      sma_mako_gam3$gcv.ubre
    )
)


print(
  sma_mako_model_comparison
)


write_csv(
  sma_mako_model_comparison,
  "outputs/sma/gam/sma_mako_model_comparison.csv"
)


summary(
  sma_mako_gam3
)

gam.check(
  sma_mako_gam3
)


# =============================================================================
# 20. Leave-one-fleet-out CV for MAK model
#
# Again omit the fleet random effect from the CV model.
# =============================================================================

fit_predict_sma_mako_lofo <- function(
    fleet_out,
    dat
) {

  train_dat <- dat |>
    filter(
      as.character(
        FLEET_CODE
      ) != fleet_out
    ) |>
    droplevels()


  test_dat <- dat |>
    filter(
      as.character(
        FLEET_CODE
      ) == fleet_out
    )


  if (
    nrow(
      train_dat
    ) < 30 ||
    nrow(
      test_dat
    ) == 0
  ) {

    return(
      NULL
    )
  }


  valid_classes <-
    levels(
      train_dat$model_fishery_class
    )


  test_dat <- test_dat |>
    filter(
      as.character(
        model_fishery_class
      ) %in%
        valid_classes
    )


  if (
    nrow(
      test_dat
    ) == 0
  ) {

    return(
      NULL
    )
  }


  cv_model <- try(

    gam(

      p_sma_mako_beta ~

        model_fishery_class +

        s(
          YEAR,
          k = 6,
          bs = "tp"
        ),

      family =
        betar(
          link = "logit"
        ),

      method =
        "REML",

      data =
        train_dat
    ),

    silent = TRUE
  )


  if (
    inherits(
      cv_model,
      "try-error"
    )
  ) {

    return(
      NULL
    )
  }


  test_dat |>
    mutate(

      predicted =
        predict(
          cv_model,
          newdata = test_dat,
          type = "response"
        ),

      observed =
        p_sma_mako,

      residual =
        predicted -
        observed
    )
}


# =============================================================================
# 21. MAK-model LOFO
# =============================================================================

sma_mako_lofo <- map_dfr(

  unique(
    as.character(
      sma_mako_training$FLEET_CODE
    )
  ),

  fit_predict_sma_mako_lofo,

  dat =
    sma_mako_training
)


sma_mako_lofo_summary <- sma_mako_lofo |>

  summarise(

    n =
      n(),

    mae =
      mean(
        abs(
          predicted -
            observed
        ),
        na.rm = TRUE
      ),

    rmse =
      sqrt(
        mean(
          (
            predicted -
              observed
          )^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted -
          observed,
        na.rm = TRUE
      ),

    rho =
      safe_cor(
        observed,
        predicted
      )
  )


sma_mako_lofo_by_group <- sma_mako_lofo |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    n =
      n(),

    mae =
      mean(
        abs(
          predicted -
            observed
        ),
        na.rm = TRUE
      ),

    rmse =
      sqrt(
        mean(
          (
            predicted -
              observed
          )^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted -
          observed,
        na.rm = TRUE
      ),

    rho =
      safe_cor(
        observed,
        predicted
      ),

    .groups = "drop"
  )


print(
  sma_mako_lofo_summary
)

print(
  sma_mako_lofo_by_group,
  n = Inf
)


write_csv(
  sma_mako_lofo_summary,
  "outputs/sma/gam/sma_mako_lofo_summary.csv"
)

write_csv(
  sma_mako_lofo_by_group,
  "outputs/sma/gam/sma_mako_lofo_by_group.csv"
)


# =============================================================================
# 22. Plot MAK-model LOFO
# =============================================================================

p_sma_mako_lofo <- ggplot(

  sma_mako_lofo,

  aes(
    x = observed,
    y = predicted,
    colour = FISHERY_GROUP_CODE
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

  coord_equal(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +

  scale_colour_viridis_d(
    end = 0.9
  ) +

  labs(

    x =
      "Observed SMA proportion of identified makos",

    y =
      "Predicted SMA proportion",

    colour =
      "Fishery group",

    title =
      "SMA composition within identified makos",

    subtitle =
      "Leave-one-fleet-out cross-validation"
  ) +

  theme_bw()


p_sma_mako_lofo


ggsave(

  "outputs/sma/gam/figures/sma_mako_gam_lofo.png",

  p_sma_mako_lofo,

  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# PART C
#
# Apply fitted models to unresolved taxonomic catch
# =============================================================================


# =============================================================================
# 23. Build prediction data
#
# Use the same Fleet x Fishery x Year covariates as the training data.
# =============================================================================

sma_prediction_base <- sma_gam_base |>

  mutate(

    FLEET_CODE =
      factor(
        FLEET_CODE,
        levels =
          levels(
            sma_shark_training$FLEET_CODE
          )
      ),

    FISHERY_GROUP_CODE =
      factor(
        FISHERY_GROUP_CODE,
        levels =
          levels(
            sma_shark_training$FISHERY_GROUP_CODE
          )
      )
  )


# =============================================================================
# 24. Predict SKH -> SMA
#
# Only groups represented by the fitted model are eligible.
# =============================================================================

sma_skh_gam_prediction <- sma_prediction_base |>

  filter(
    skh_mt > 0
  ) |>

  mutate(

    gam_supported =
      FISHERY_GROUP_CODE %in%
        levels(
          sma_shark_training$FISHERY_GROUP_CODE
        ) &

      model_fishery_class %in%
        levels(
          sma_shark_training$model_fishery_class
        )
  )


# ---------------------------------------------------------------------------
# Prediction from S4
#
# exclude = "s(FLEET_CODE)" means prediction is based on the population-level
# fishery/year/targeting relationship rather than requiring a known fleet RE.
#
# This also allows application to fleets not represented in training.
# ---------------------------------------------------------------------------

sma_skh_gam_prediction <-
  sma_skh_gam_prediction |>

  mutate(

    p_sma_skh_gam =
      NA_real_
  )


supported_rows <-
  which(
    sma_skh_gam_prediction$gam_supported
  )


if (
  length(
    supported_rows
  ) > 0
) {

  sma_skh_gam_prediction$p_sma_skh_gam[
    supported_rows
  ] <-

    predict(

      sma_shark_gam4,

      newdata =
        sma_skh_gam_prediction[
          supported_rows,
        ],

      type =
        "response",

      exclude =
        "s(FLEET_CODE)"
    )
}


sma_skh_gam_prediction <-
  sma_skh_gam_prediction |>

  mutate(

    sma_from_skh_gam_sensitivity = case_when(

      gam_supported &
        !is.na(
          p_sma_skh_gam
        ) ~
        skh_mt *
          p_sma_skh_gam,

      TRUE ~
        NA_real_
    )
  )


# =============================================================================
# 25. Interpretation and use of SMA species-composition GAMs
#
# SKH -> SMA
# ----------
# The SKH GAM is retained as a statistical sensitivity analysis only.
#
# Although the final model explained substantial in-sample variation,
# leave-one-fleet-out predictive performance was modest. In addition:
#
#   * predictions for historical GN SKH depend strongly on temporal
#     backcasting beyond observed species-composition support;
#
#   * the response contains many exact zeros, especially in LI and PS;
#
#   * the beta-regression formulation requires squeezing zeros away from
#     the boundary and produces substantial positive SMA predictions in
#     some periods where empirical SMA composition is zero or near zero;
#
#   * predicted compositions therefore should not be interpreted as the
#     preferred reconstruction of unresolved SKH.
#
# The HSCR N1/N2 estimates remain the primary SKH -> SMA reconstruction.
# GAM predictions are retained for sensitivity analyses and diagnostics.
#
# MAK -> SMA
# ----------
# The MAK GAM is retained as a diagnostic only.
#
# Identified mako catch is strongly dominated by SMA and the empirical
# hierarchical SMA / (SMA + LMA) estimator has a direct taxonomic
# interpretation. The GAM showed poor leave-one-fleet-out transferability,
# including negative rank correlations across the principal fishery groups.
#
# The hierarchical empirical MAK allocation therefore remains the primary
# MAK -> SMA reconstruction.
# =============================================================================

use_skh_gam_as_primary <- FALSE
use_skh_gam_as_sensitivity <- TRUE

use_mako_gam_for_reconstruction <- FALSE
use_mako_gam_as_diagnostic <- TRUE



# =============================================================================
# 26. MAK GAM diagnostic predictions
# =============================================================================

sma_mak_prediction_base <- sma_gam_base |>

  filter(
    mak_mt > 0
  ) |>

  mutate(

    FLEET_CODE =
      factor(
        FLEET_CODE,
        levels =
          levels(
            sma_mako_training$FLEET_CODE
          )
      ),

    FISHERY_GROUP_CODE =
      factor(
        FISHERY_GROUP_CODE,
        levels =
          levels(
            sma_mako_training$FISHERY_GROUP_CODE
          )
      ),

    gam_supported =
      FISHERY_GROUP_CODE %in%
        levels(
          sma_mako_training$FISHERY_GROUP_CODE
        ) &

      model_fishery_class %in%
        levels(
          sma_mako_training$model_fishery_class
        ),

    p_sma_mak_gam =
      NA_real_
  )


supported_rows <-
  which(
    sma_mak_prediction_base$gam_supported
  )


if (
  length(
    supported_rows
  ) > 0
) {

  sma_mak_prediction_base$p_sma_mak_gam[
    supported_rows
  ] <-

    predict(

      sma_mako_gam3,

      newdata =
        sma_mak_prediction_base[
          supported_rows,
        ],

      type =
        "response",

      exclude =
        "s(FLEET_CODE)"
    )
}


sma_mak_gam_diagnostic_prediction <-
  sma_mak_prediction_base |>

  mutate(

    sma_from_mak_gam_diagnostic =
      case_when(

        gam_supported &
          !is.na(
            p_sma_mak_gam
          ) ~
          mak_mt *
            p_sma_mak_gam,

        TRUE ~
          NA_real_
      )
  )

# =============================================================================
# PART D
#
# GAM sensitivity and diagnostic summaries
# =============================================================================


# =============================================================================
# 27. SKH GAM sensitivity summary
# =============================================================================

sma_skh_gam_sensitivity_summary <-
  sma_skh_gam_prediction |>
  summarise(

    total_skh_mt =
      sum(
        skh_mt,
        na.rm = TRUE
      ),

    gam_supported_skh_mt =
      sum(
        skh_mt[
          gam_supported
        ],
        na.rm = TRUE
      ),

    gam_supported_fraction =
      gam_supported_skh_mt /
        total_skh_mt,

    sma_from_skh_gam_sensitivity_mt =
      sum(
        sma_from_skh_gam,
        na.rm = TRUE
      )
  )


print(
  sma_skh_gam_sensitivity_summary
)

# =============================================================================
# 28. MAK GAM diagnostic summary
# =============================================================================

sma_mak_gam_diagnostic_summary <-
  sma_mak_gam_diagnostic_prediction |>
  summarise(

    total_mak_mt =
      sum(
        mak_mt,
        na.rm = TRUE
      ),

    gam_supported_mak_mt =
      sum(
        mak_mt[
          gam_supported
        ],
        na.rm = TRUE
      ),

    gam_supported_fraction =
      gam_supported_mak_mt /
        total_mak_mt,

    sma_from_mak_gam_diagnostic_mt =
      sum(
        sma_from_mak_gam_diagnostic,
        na.rm = TRUE
      )
  )


print(
  sma_mak_gam_diagnostic_summary
)

#========================================
write_csv(
  sma_skh_gam_summary,
  "outputs/sma/gam/sma_skh_gam_summary.csv"
)

write_csv(
  sma_mak_gam_summary,
  "outputs/sma/gam/sma_mak_gam_summary.csv"
)


# =============================================================================
# 28. Reconstruction by fishery group
# =============================================================================

sma_skh_gam_by_group <-
  sma_skh_gam_prediction |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    skh_mt =
      sum(
        skh_mt,
        na.rm = TRUE
      ),

    supported_skh_mt =
      sum(
        skh_mt[
          gam_supported
        ],
        na.rm = TRUE
      ),

    sma_from_skh_gam_mt =
      sum(
        sma_from_skh_gam,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


sma_mak_gam_by_group <-
  sma_mak_gam_prediction |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    mak_mt =
      sum(
        mak_mt,
        na.rm = TRUE
      ),

    supported_mak_mt =
      sum(
        mak_mt[
          gam_supported
        ],
        na.rm = TRUE
      ),

    sma_from_mak_gam_mt =
      sum(
        sma_from_mak_gam,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


print(
  sma_skh_gam_by_group,
  n = Inf
)

print(
  sma_mak_gam_by_group,
  n = Inf
)


write_csv(
  sma_skh_gam_by_group,
  "outputs/sma/gam/sma_skh_gam_by_group.csv"
)

write_csv(
  sma_mak_gam_by_group,
  "outputs/sma/gam/sma_mak_gam_by_group.csv"
)


# =============================================================================
# 29. Annual GAM reconstructions
# =============================================================================

sma_skh_gam_annual <-
  sma_skh_gam_prediction |>

  group_by(
    YEAR
  ) |>

  summarise(

    skh_mt =
      sum(
        skh_mt,
        na.rm = TRUE
      ),

    sma_from_skh_gam_mt =
      sum(
        sma_from_skh_gam,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


sma_mak_gam_annual <-
  sma_mak_gam_prediction |>

  group_by(
    YEAR
  ) |>

  summarise(

    mak_mt =
      sum(
        mak_mt,
        na.rm = TRUE
      ),

    sma_from_mak_gam_mt =
      sum(
        sma_from_mak_gam,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


write_csv(
  sma_skh_gam_annual,
  "outputs/sma/gam/sma_skh_gam_annual.csv"
)

write_csv(
  sma_mak_gam_annual,
  "outputs/sma/gam/sma_mak_gam_annual.csv"
)


# =============================================================================
# 30. Save model objects
# =============================================================================

saveRDS(
  sma_shark_gam4,
  "outputs/sma/gam/sma_shark_gam4.rds"
)

saveRDS(
  sma_mako_gam3,
  "outputs/sma/gam/sma_mako_gam3.rds"
)


# =============================================================================
# 31. Objects retained in workspace
#
# Key objects:
#
#   sma_shark_training
#   sma_shark_gam1
#   sma_shark_gam2
#   sma_shark_gam3
#   sma_shark_gam4
#   sma_shark_model_comparison
#   sma_shark_lofo
#   sma_shark_lofo_summary
#   sma_shark_lofo_by_group
#
#   sma_mako_training
#   sma_mako_gam1
#   sma_mako_gam2
#   sma_mako_gam3
#   sma_mako_model_comparison
#   sma_mako_lofo
#   sma_mako_lofo_summary
#   sma_mako_lofo_by_group
#
#   sma_skh_gam_prediction
#   sma_mak_gam_prediction
#   sma_skh_gam_summary
#   sma_mak_gam_summary
#
# =============================================================================