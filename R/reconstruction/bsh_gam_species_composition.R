# =============================================================================
# bsh_gam_species_composition.R
#
# GAM-based reconstruction of blue shark species composition
#
# PURPOSE
# -------
#   Evaluate a beta-regression GAM as a sensitivity to the empirical
#   HSCR allocation of generic Sharks NEI (SKH) to blue shark. 
# Inputs:
#   - BSH reconstruction strata
#   - identified shark species composition
#   - target catch composition
#
# Main outputs:
#   - bsh_gam_hscr_comparison
#   - bsh_gam_hscr_summary
#   - bsh_gam_hscr_by_group
#   - gam_cv_final
#   - gam_cv_by_group
#
# Estimation Method:
#
#   P(BSH | identified shark catch)
#
# for Fleet x Fishery x Year strata with incomplete species reporting.
#
# This provides a statistical alternative to the
# Hierarchical species-composition reconstruction (HSCR).
#
# The model is fitted only to strata containing species-specific shark
# information. It is then used to predict the BSH proportion in strata
# containing generic shark catch.
#
# IMPORTANT
# ---------
# This script addresses UNDER-REPORTING TO SPECIES.#
# It does NOT estimate completely unreported shark catch.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(purrr)
library(tibble)


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
  "SKH"
)

# Other generic categories retained for later sensitivity testing
excluded_generic_shark_codes <- c(
  "MSK",
  "AG20",
  "SHXX",
  "SHRK",
  "AG22",
  "SKX",
  "SKHS"
)


# =============================================================================
# 2. Reporting periods
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
# 3. Build Fleet x Fishery x Year shark-composition dataset
# =============================================================================

gam_shark_data <- rc_dat |>
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

    # -------------------------------------------------------------
    # Species-specific sharks
    # -------------------------------------------------------------

    bsh_mt = sum(
      CATCH[SPECIES_CODE == "BSH"],
      na.rm = TRUE
    ),

    sma_mt = sum(
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

    # -------------------------------------------------------------
    # Target species
    # -------------------------------------------------------------

    yft_mt = sum(
      CATCH[SPECIES_CODE == "YFT"],
      na.rm = TRUE
    ),

    bet_mt = sum(
      CATCH[SPECIES_CODE == "BET"],
      na.rm = TRUE
    ),

    skj_mt = sum(
      CATCH[SPECIES_CODE == "SKJ"],
      na.rm = TRUE
    ),

    alb_mt = sum(
      CATCH[SPECIES_CODE == "ALB"],
      na.rm = TRUE
    ),

    swo_mt = sum(
      CATCH[SPECIES_CODE == "SWO"],
      na.rm = TRUE
    ),

    sbf_mt = sum(
      CATCH[SPECIES_CODE == "SBF"],
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# =============================================================================
# 4. Derive explanatory variables
# =============================================================================

gam_shark_data <- gam_shark_data |>
  mutate(

    target_primary_mt =
      yft_mt +
      bet_mt +
      skj_mt +
      alb_mt +
      swo_mt,

    target_plus_sbf_mt =
      target_primary_mt +
      sbf_mt,

    tropical_tuna_mt =
      yft_mt +
      bet_mt +
      skj_mt,

    temperate_target_mt =
      alb_mt +
      swo_mt,

    # -------------------------------------------------------------
    # Target composition
    # -------------------------------------------------------------

    prop_tropical =
      if_else(
        target_plus_sbf_mt > 0,
        tropical_tuna_mt / target_plus_sbf_mt,
        0
      ),

    prop_alb =
      if_else(
        target_plus_sbf_mt > 0,
        alb_mt / target_plus_sbf_mt,
        0
      ),

    prop_swo =
      if_else(
        target_plus_sbf_mt > 0,
        swo_mt / target_plus_sbf_mt,
        0
      ),

    prop_sbf =
      if_else(
        target_plus_sbf_mt > 0,
        sbf_mt / target_plus_sbf_mt,
        0
      ),

    # -------------------------------------------------------------
    # Response
    # -------------------------------------------------------------

    p_bsh =
      if_else(
        identified_shark_mt > 0,
        bsh_mt / identified_shark_mt,
        NA_real_
      )
  )

# =============================================================================
# 5. Prepare beta-regression response mgcv::betar() requires values strictly between 0 and 1.
# =============================================================================

eps <- 0.001

gam_shark_data <- gam_shark_data |>
  mutate(

    p_bsh_beta = case_when(

      is.na(p_bsh) ~
        NA_real_,

      p_bsh <= 0 ~
        eps,

      p_bsh >= 1 ~
        1 - eps,

      TRUE ~
        p_bsh
    )
  )


# =============================================================================
# 6. Training data
# =============================================================================

min_identified_shark_mt <- 1

gam_supported_groups <- gam_shark_data |>
  filter(
    identified_shark_mt >= min_identified_shark_mt,
    !is.na(p_bsh_beta),
    target_plus_sbf_mt > 0
  ) |>
  group_by(FISHERY_GROUP_CODE) |>
  summarise(
    n_strata = n(),
    n_fleets = n_distinct(FLEET_CODE),
    .groups = "drop"
  ) |>
  filter(
    n_strata >= 10,
    n_fleets >= 2
  )

gam_supported_groups

# Check
gam_train <- gam_shark_data |>
  semi_join(
    gam_supported_groups,
    by = "FISHERY_GROUP_CODE"
  ) |>
  filter(
    identified_shark_mt >= min_identified_shark_mt,
    !is.na(p_bsh_beta),
    target_plus_sbf_mt > 0
  ) |>
  mutate(
    FLEET_CODE = factor(FLEET_CODE),
    FISHERY_CODE = factor(FISHERY_CODE),
    FISHERY_GROUP_CODE = factor(FISHERY_GROUP_CODE)
  )

gam_train |>
  count(
    FISHERY_GROUP_CODE,
    sort = TRUE
  )


# =============================================================================
# 7. Response diagnostics
# =============================================================================

ggplot(
  gam_train,
  aes(p_bsh)
) +
  geom_histogram(
    bins = 40
  ) +
  facet_wrap(
    ~ FISHERY_GROUP_CODE,
    scales = "free_y"
  ) +
  labs(
    x = "BSH proportion of identified shark catch",
    y = "Number of strata",
    title = "Observed blue shark species composition",
    subtitle = "Fleet x fishery x year strata"
  ) +
  theme_bw()

# =============================================================================
# 8a. GAM 0: broad fishery structure
# =============================================================================

gam0 <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE,

  family = betar(link = "logit"),

  method = "REML",

  data = gam_train
)

# =============================================================================
# 8b. GAM 1: fishery structure + temporal trend
# =============================================================================

gam1 <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE +
    s(
      YEAR,
      k = 10,
      bs = "tp"
    ),

  family = betar(link = "logit"),

  method = "REML",

  data = gam_train
)
# =============================================================================
# 8c. GAM 2: fishery-specific temporal patterns
# =============================================================================

gam2 <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE +

    s(
      YEAR,
      by = FISHERY_GROUP_CODE,
      k = 8,
      bs = "tp"
    ),

  family = betar(link = "logit"),

  method = "REML",

  data = gam_train
)

# =============================================================================
# 8d. GAM 3: fishery-specific temporal trends + fleet effects
# =============================================================================

gam3 <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE +

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

  family = betar(link = "logit"),

  method = "REML",

  data = gam_train
)

# =============================================================================
# 8e. GAM 4: target assemblage
# =============================================================================

gam4 <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE +

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

  family = betar(link = "logit"),

  method = "REML",

  select = TRUE,

  data = gam_train
)


# =============================================================================
# 8f. GAM 5: fleet + fishery random effects
# =============================================================================

gam5 <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE +

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
      FISHERY_CODE,
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

  family = betar(link = "logit"),

  method = "REML",

  select = TRUE,

  data = gam_train
)


# =============================================================================
# 9. Model summaries
# =============================================================================

summary(gam0)
summary(gam1)
summary(gam2)
summary(gam3)
summary(gam4)
summary(gam5)

gam.check(gam3)
gam.check(gam4)
gam.check(gam5)

# =============================================================================
# 10. Candidate model comparison
# =============================================================================

gam_models <- list(
  GAM0 = gam0,
  GAM1 = gam1,
  GAM2 = gam2,
  GAM3 = gam3,
  GAM4 = gam4,
  GAM5 = gam5
)

gam_model_comparison <- imap_dfr(
  gam_models,
  ~ tibble(
    model = .y,

    AIC =
      AIC(.x),

    deviance_explained =
      summary(.x)$dev.expl,

    adjusted_r2 =
      summary(.x)$r.sq,

    n =
      nobs(.x)
  )
) |>
  arrange(AIC)

gam_model_comparison

# =============================================================================
# 11. In-sample diagnostics
# =============================================================================

gam_train_pred <- gam_train |>
  mutate(

    pred_gam0 =
      predict(
        gam0,
        newdata = gam_train,
        type = "response"
      ),

    pred_gam3 =
      predict(
        gam3,
        newdata = gam_train,
        type = "response"
      ),

    pred_gam4 =
      predict(
        gam4,
        newdata = gam_train,
        type = "response"
      ),

    pred_gam5 =
      predict(
        gam5,
        newdata = gam_train,
        type = "response"
      )
  )


ggplot(
  gam_train_pred,
  aes(
    x = p_bsh,
    y = pred_gam4
  )
) +
  geom_point(
    alpha = 0.35
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
    x = "Observed BSH proportion",
    y = "Predicted BSH proportion",
    title = "Observed versus GAM-predicted shark composition",
    subtitle = "Candidate GAM 4"
  ) +
  theme_bw()


# Check what the Cross Validatoin excluedes 
cv_support_by_group <- gam_train |>
  group_by(FISHERY_GROUP_CODE) |>
  summarise(
    n_strata = n(),
    n_fleets = n_distinct(FLEET_CODE),
    .groups = "drop"
  )

cv_support_by_group

# =============================================================================
# 12. Leave-one-fleet-out cross-validation
# start with Model 4. w/o random effect. the unseen fleet anc on have an estiamted random effect. 
#
#
# Predictions are made only for fishery groups represented in the
# corresponding training fold.
# =============================================================================

cv_fleets <- sort(
  unique(as.character(gam_train$FLEET_CODE))
)

cv_results <- map_dfr(

  cv_fleets,

  function(test_fleet) {

    train_i <- gam_train |>
      filter(
        as.character(FLEET_CODE) != test_fleet
      )

    test_i <- gam_train |>
      filter(
        as.character(FLEET_CODE) == test_fleet
      )

    if (
      nrow(train_i) < 30 ||
      nrow(test_i) == 0
    ) {
      return(NULL)
    }

    # -----------------------------------------------------------------------
    # Restrict test observations to factor levels represented in training.
    #
    # This prevents a held-out fleet from introducing an entirely unseen
    # fishery group.
    # -----------------------------------------------------------------------

    supported_groups <- unique(
      as.character(train_i$FISHERY_GROUP_CODE)
    )

    test_i <- test_i |>
      filter(
        as.character(FISHERY_GROUP_CODE) %in%
          supported_groups
      )

    if (nrow(test_i) == 0) {
      return(NULL)
    }

    # -----------------------------------------------------------------------
    # Re-factor training and test data using identical levels
    # -----------------------------------------------------------------------

    group_levels <- levels(gam_train$FISHERY_GROUP_CODE)

    train_i <- train_i |>
      mutate(
        FISHERY_GROUP_CODE = factor(
          FISHERY_GROUP_CODE,
          levels = group_levels
        )
      )

    test_i <- test_i |>
      mutate(
        FISHERY_GROUP_CODE = factor(
          FISHERY_GROUP_CODE,
          levels = group_levels
        )
      )

    # -----------------------------------------------------------------------
    # Fit fold
    # -----------------------------------------------------------------------

    model_i <- gam(
      gam_cv_formula,

      family =
        betar(link = "logit"),

      method =
        "REML",

      select =
        TRUE,

      data =
        train_i,

      drop.unused.levels =
        FALSE
    )

    # -----------------------------------------------------------------------
    # Predict
    # -----------------------------------------------------------------------

    pred_i <- predict(
      model_i,
      newdata = test_i,
      type = "response"
    )

    tibble(

      YEAR =
        test_i$YEAR,

      FLEET_CODE =
        as.character(test_i$FLEET_CODE),

      FISHERY_GROUP_CODE =
        as.character(test_i$FISHERY_GROUP_CODE),

      observed =
        test_i$p_bsh,

      predicted =
        as.numeric(pred_i)
    )
  }
)





cv_coverage <- gam_train |>
  count(
    FISHERY_GROUP_CODE,
    name = "n_training"
  ) |>
  left_join(
    cv_results |>
      count(
        FISHERY_GROUP_CODE,
        name = "n_cv_predictions"
      ),
    by = "FISHERY_GROUP_CODE"
  ) |>
  mutate(
    n_cv_predictions =
      coalesce(n_cv_predictions, 0L),

    pct_cv_coverage =
      100 *
      n_cv_predictions /
      n_training
  )

cv_coverage

# =============================================================================
# 13. Cross-validation metrics ::
# =============================================================================

gam_cv_metrics <- cv_results |>
  summarise(

    n = n(),

    MAE =
      mean(
        abs(predicted - observed),
        na.rm = TRUE
      ),

    RMSE =
      sqrt(
        mean(
          (predicted - observed)^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted - observed,
        na.rm = TRUE
      ),

    correlation =
      cor(
        predicted,
        observed,
        method = "spearman",
        use = "complete.obs"
      )
  )

gam_cv_metrics

# now do by fishery group
gam_cv_by_group <- cv_results |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    n = n(),

    MAE =
      mean(
        abs(predicted - observed),
        na.rm = TRUE
      ),

    RMSE =
      sqrt(
        mean(
          (predicted - observed)^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted - observed,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

gam_cv_by_group

# 
gam_cv_metrics
gam_cv_by_group

# Overall observed vs predicted CV plot

p_lofo_cv_bsh<- ggplot(
  cv_results,
  aes(
    x = observed,
    y = predicted
  )
) +
  geom_point(alpha = 0.4) +
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
    x = "Observed BSH proportion",
    y = "Cross-validated predicted BSH proportion",
    title = "Leave-one-fleet-out cross-validation",
    subtitle = "GAM species-composition model"
  ) +
  theme_bw()

p_lofo_cv_bsh



#-  ==================================================================================
# compare gam4 and gam5

cv_gam_model <- function(
    formula,
    data,
    fleet_var = "FLEET_CODE"
) {

  fleets <- sort(
    unique(
      as.character(
        data[[fleet_var]]
      )
    )
  )

  map_dfr(
    fleets,
    function(test_fleet) {

      train_i <- data |>
        filter(
          as.character(.data[[fleet_var]]) != test_fleet
        )

      test_i <- data |>
        filter(
          as.character(.data[[fleet_var]]) == test_fleet
        )

      supported_groups <- unique(
        as.character(
          train_i$FISHERY_GROUP_CODE
        )
      )

      test_i <- test_i |>
        filter(
          as.character(FISHERY_GROUP_CODE) %in%
            supported_groups
        )

      if (
        nrow(train_i) < 30 ||
        nrow(test_i) == 0
      ) {
        return(NULL)
      }

      model_i <- gam(
        formula,
        family = betar(link = "logit"),
        method = "REML",
        select = TRUE,
        data = train_i,
        drop.unused.levels = FALSE
      )

      pred_i <- predict(
        model_i,
        newdata = test_i,
        type = "response"
      )

      tibble(
        YEAR = test_i$YEAR,
        FLEET_CODE =
          as.character(test_i$FLEET_CODE),
        FISHERY_GROUP_CODE =
          as.character(test_i$FISHERY_GROUP_CODE),
        observed = test_i$p_bsh,
        predicted = as.numeric(pred_i)
      )
    }
  ) 
}

# lets compare modesl 3, 4,  & 5
# check 
cv3 <- cv_gam_model(
  formula(gam3),
  gam_train
)

cv4 <- cv_gam_model(
  formula(gam4),
  gam_train
)

cv5 <- cv_gam_model(
  formula(gam5),
  gam_train
)
# summarize 
summarise_cv <- function(x, name) {

  x |>
    summarise(
      model = name,
      n = n(),
      MAE =
        mean(abs(predicted - observed)),
      RMSE =
        sqrt(
          mean(
            (predicted - observed)^2
          )
        ),
      bias =
        mean(predicted - observed),
      rho =
        cor(
          predicted,
          observed,
          method = "spearman"
        )
    )
}

bind_rows(
  summarise_cv(cv3, "GAM3"),
  summarise_cv(cv4, "GAM4"),
  summarise_cv(cv5, "GAM5")
)


# =============================================================================
# 14. Generic-shark prediction strata
#
# IMPORTANT:
# Preserve original fleet and fishery codes.
#
# GAM4 does not use FLEET_CODE or FISHERY_CODE directly in prediction.
# They are retained as identifiers for comparison with HSCR.
# =============================================================================

gam_nei_prediction_data <- gam_shark_data |>
  filter(
    generic_shark_mt > 0
  ) |>
  mutate(

    # Preserve identifiers as character strings
    FLEET_CODE =
      as.character(FLEET_CODE),

    FISHERY_CODE =
      as.character(FISHERY_CODE),

    # FISHERY_GROUP_CODE IS used by the model, so its factor levels
    # must correspond to those used during fitting.
    FISHERY_GROUP_CODE =
      factor(
        as.character(FISHERY_GROUP_CODE),
        levels = levels(gam_train$FISHERY_GROUP_CODE)
      )
  )

#identify prediction rows whose broad fishery group is unsupported:
gam_nei_prediction_data <- gam_nei_prediction_data |>
  mutate(

    gam_group_supported =
      !is.na(FISHERY_GROUP_CODE)
  )


gam_nei_prediction_data <- gam_shark_data |>
  filter(
    generic_shark_mt > 0
  ) |>
  mutate(

    FLEET_CODE =
      as.character(FLEET_CODE),

    FISHERY_CODE =
      as.character(FISHERY_CODE),

    FISHERY_GROUP_CODE_ORIG =
      as.character(FISHERY_GROUP_CODE),

    FISHERY_GROUP_CODE =
      factor(
        FISHERY_GROUP_CODE_ORIG,
        levels = levels(gam_train$FISHERY_GROUP_CODE)
      ),

    gam_group_supported =
      !is.na(FISHERY_GROUP_CODE)
  )


gam_nei_supported <- gam_nei_prediction_data |>
  filter(
    gam_group_supported
  )

# =============================================================================
# 15. Production GAM for prediction
# =============================================================================

gam_production <- gam(
  gam_cv_formula,

  family =
    betar(link = "logit"),

  method =
    "REML",

  select =
    TRUE,

  data =
    gam_train
)


# =============================================================================
# 16. Predict BSH composition in generic-shark strata
# =============================================================================
gam_pred_link <- predict(
  gam_production,
  newdata = gam_nei_supported,
  type = "link",
  se.fit = TRUE
)


inv_logit <- function(x) {
  plogis(x)
}

gam_nei_predictions <- gam_nei_supported |>
  mutate(

    gam_link =
      as.numeric(gam_pred_link$fit),

    gam_se =
      as.numeric(gam_pred_link$se.fit),

    p_bsh_gam =
      plogis(gam_link),

    p_bsh_gam_lower =
      plogis(
        gam_link - 1.96 * gam_se
      ),

    p_bsh_gam_upper =
      plogis(
        gam_link + 1.96 * gam_se
      ),

    bsh_from_nei_gam =
      generic_shark_mt *
      p_bsh_gam,

    bsh_from_nei_gam_lower =
      generic_shark_mt *
      p_bsh_gam_lower,

    bsh_from_nei_gam_upper =
      generic_shark_mt *
      p_bsh_gam_upper
  ) |>
  mutate(

    # Convert back for joins and reporting
    FISHERY_GROUP_CODE =
      FISHERY_GROUP_CODE_ORIG
  )


# =============================================================================
# 17. Determine temporal support for GAM predictions
# =============================================================================

gam_training_support <- gam_train |>
  mutate(
    FISHERY_GROUP_CODE =
      as.character(FISHERY_GROUP_CODE)
  ) |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    donor_first_year =
      min(YEAR, na.rm = TRUE),

    donor_last_year =
      max(YEAR, na.rm = TRUE),

    donor_n =
      n(),

    donor_fleets =
      n_distinct(FLEET_CODE),

    .groups = "drop"
  )


gam_nei_predictions <- gam_nei_predictions |>
  left_join(
    gam_training_support,
    by = "FISHERY_GROUP_CODE"
  ) |>
  mutate(

    gam_temporal_support = case_when(

      is.na(donor_first_year) ~
        "unsupported_fishery_group",

      YEAR < donor_first_year ~
        "backcast",

      YEAR > donor_last_year ~
        "forward_extrapolation",

      TRUE ~
        "within_training_period"
    )
  )
 
# dbug
gam_nei_predictions |>
  count(
    FISHERY_GROUP_CODE,
    gam_temporal_support,
    sort = TRUE
  ) |>
  print(n = Inf)
#
names(gam_nei_predictions)[
  grepl(
    "support|donor",
    names(gam_nei_predictions)
  )
]
# =============================================================================
# 18. Compare GAM against HSCR
# =============================================================================

comparison_keys <- c(
  "YEAR",
  "FLEET_CODE",
  "FISHERY_GROUP_CODE",
  "FISHERY_CODE"
)

gam_comparison_dat <- gam_nei_predictions |>
  mutate(
    FLEET_CODE = as.character(FLEET_CODE),
    FISHERY_GROUP_CODE = as.character(FISHERY_GROUP_CODE),
    FISHERY_CODE = as.character(FISHERY_CODE)
  )

hscr_comparison_dat <- nei_candidates |>
  mutate(
    FLEET_CODE = as.character(FLEET_CODE),
    FISHERY_GROUP_CODE = as.character(FISHERY_GROUP_CODE),
    FISHERY_CODE = as.character(FISHERY_CODE)
  )

bsh_gam_hscr_comparison <- gam_comparison_dat |>
  select(
    all_of(comparison_keys),

    generic_shark_mt,

    p_bsh_gam,
    p_bsh_gam_lower,
    p_bsh_gam_upper,

    bsh_from_nei_gam,
    bsh_from_nei_gam_lower,
    bsh_from_nei_gam_upper,

    donor_first_year,
    donor_last_year,
    donor_n,
    donor_fleets,

    gam_temporal_support,
    fleet_seen_in_training
  ) |>
  left_join(
    hscr_comparison_dat |>
      select(
        all_of(comparison_keys),

        p_bsh_N1,
        p_bsh_N2,

        bsh_from_nei_N1,
        bsh_from_nei_N2,

        nei_N1_source,
        nei_N2_source
      ),
    by = comparison_keys
  )


# did I join it correctly  ? THese should be 0 rows for each object
gam_comparison_dat |>
  count(
    across(all_of(comparison_keys))
  ) |>
  filter(n > 1)

hscr_comparison_dat |>
  count(
    across(all_of(comparison_keys))
  ) |>
  filter(n > 1)

# how many generic shark strata are left out there 
training_fleets <- gam_train |>
  distinct(
    FLEET_CODE
  ) |>
  pull(FLEET_CODE) |>
  as.character()


gam_nei_predictions <- gam_nei_predictions |>
  mutate(

    fleet_seen_in_training =
      FLEET_CODE %in% training_fleets
  )

# glimpse(gam_nei_predictions)


gam_nei_predictions |>
  group_by(
    FISHERY_GROUP_CODE,
    fleet_seen_in_training,
    gam_temporal_support
  ) |>
  summarise(

    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    gam_bsh_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(
    desc(generic_shark_mt)
  ) |>
  print(n = Inf)

# =============================================================================
# Annual BSH HSCR and GAM sensitivity comparison
#
# HSCR:
#   canonical annual scenarios from nei_scenario_annual
#
# GAM:
#   annual predictions aggregated from gam_hscr_comparison
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Canonical annual HSCR series
# -----------------------------------------------------------------------------

bsh_hscr_annual <-
  nei_scenario_annual |>

  transmute(

    YEAR,

    generic_shark_mt,

    bsh_HSCR_N1_mt =
      bsh_N1_mt,

    bsh_HSCR_N2_mt =
      bsh_N2_mt
  )


# -----------------------------------------------------------------------------
# 2. Annual GAM sensitivity series
# -----------------------------------------------------------------------------

bsh_gam_annual <-
  bsh_gam_hscr_comparison |>

  group_by(
    YEAR
  ) |>

  summarise(

    bsh_GAM_mt =
      sum(
        bsh_from_nei_gam,
        na.rm = TRUE
      ),

    bsh_GAM_lower_mt =
      sum(
        bsh_from_nei_gam_lower,
        na.rm = TRUE
      ),

    bsh_GAM_upper_mt =
      sum(
        bsh_from_nei_gam_upper,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 3. Combine HSCR and GAM annual series
# -----------------------------------------------------------------------------

bsh_skh_method_annual <-
  bsh_hscr_annual |>

  left_join(
    bsh_gam_annual,
    by = "YEAR"
  ) |>

  mutate(

    hscr_lower_mt =
      pmin(
        bsh_HSCR_N1_mt,
        bsh_HSCR_N2_mt
      ),

    hscr_upper_mt =
      pmax(
        bsh_HSCR_N1_mt,
        bsh_HSCR_N2_mt
      ),

    hscr_midpoint_mt =
      (
        bsh_HSCR_N1_mt +
        bsh_HSCR_N2_mt
      ) / 2,

    gam_minus_N1_mt =
      bsh_GAM_mt -
      bsh_HSCR_N1_mt,

    gam_minus_N2_mt =
      bsh_GAM_mt -
      bsh_HSCR_N2_mt
  )


# =============================================================================
# 19. Overall GAM vs HSCR comparison
# =============================================================================

bsh_gam_hscr_summary <-
  bsh_gam_hscr_comparison |>
  summarise(
    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE)
  ) |>
  mutate(
    difference_GAM_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt,

    difference_GAM_N2_mt =
      bsh_GAM_mt - bsh_HSCR_N2_mt
  )

glimpse(bsh_gam_hscr_summary)
# =============================================================================
# 20. GAM vs HSCR by fishery group
# =============================================================================

bsh_gam_hscr_by_group <- bsh_gam_hscr_comparison |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    difference_GAM_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt,

    difference_GAM_N2_mt =
      bsh_GAM_mt - bsh_HSCR_N2_mt,

    .groups = "drop"
  ) |>
  arrange(
    desc(generic_shark_mt)
  )

# =============================================================================
# 21. GAM extrapolation support
# =============================================================================

gam_extrapolation_summary <- bsh_gam_hscr_comparison |>
  group_by(
    FISHERY_GROUP_CODE,
    gam_temporal_support
  ) |>
  summarise(

    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(
    FISHERY_GROUP_CODE,
    gam_temporal_support
  )

#  fleet vs temporal extrapolation 
glimpse(bsh_gam_hscr_comparison)


gam_support_summary <- bsh_gam_hscr_comparison |>
  mutate(
    support_class = case_when(
      gam_temporal_support == "backcast" &
        !fleet_seen_in_training ~
        "backcast_unseen_fleet",

      gam_temporal_support == "backcast" &
        fleet_seen_in_training ~
        "backcast_seen_fleet",

      gam_temporal_support == "within_training_period" &
        !fleet_seen_in_training ~
        "within_period_unseen_fleet",

      gam_temporal_support == "within_training_period" &
        fleet_seen_in_training ~
        "within_period_seen_fleet",

      TRUE ~
        "other"
    )
  ) |>
  group_by(
    FISHERY_GROUP_CODE,
    support_class
  ) |>
  summarise(

    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    .groups = "drop"
  ) |>
  arrange(
    desc(generic_shark_mt)
  )

print(gam_support_summary, n = Inf)



#===========================================


# output
print(gam_extrapolation_summary, n = Inf)
print(bsh_gam_hscr_by_group, n = Inf)
print(bsh_gam_hscr_summary)


# 
plot(
  gam4,
  pages = 1,
  shade = TRUE
)

# I have to drop the OT group. way too much of the generic shark is ending up as bsh
# =============================================================================
# Final fishery groups retained for GAM species-composition modelling
# =============================================================================

gam_supported_groups_final <- c(
  "GN",
  "LI",
  "LL",
  "PS"
)

gam_train <- gam_shark_data |>
  filter(
    FISHERY_GROUP_CODE %in% gam_supported_groups_final,
    identified_shark_mt >= min_identified_shark_mt,
    !is.na(p_bsh_beta),
    target_plus_sbf_mt > 0
  ) |>
  mutate(
    FLEET_CODE =
      factor(FLEET_CODE),

    FISHERY_CODE =
      factor(FISHERY_CODE),

    FISHERY_GROUP_CODE =
      factor(FISHERY_GROUP_CODE)
  )

gam_train |>
  count(
    FISHERY_GROUP_CODE,
    sort = TRUE
  )
# =============================================================================
# Refit GAM4 without OT
# =============================================================================

gam4_final <- gam(
  p_bsh_beta ~
    FISHERY_GROUP_CODE +

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
    betar(link = "logit"),

  method =
    "REML",

  select =
    TRUE,

  data =
    gam_train
)

summary(gam4_final)
gam.check(gam4_final)


# =============================================================================
# Prediction formula for unseen fleets
# =============================================================================

gam_predictive_formula <-
  p_bsh_beta ~
    FISHERY_GROUP_CODE +

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
    )

# run the Cross validation. 
cv_final <- cv_gam_model(
  gam_predictive_formula,
  gam_train
)

gam_cv_final <- cv_final |>
  summarise(

    n = n(),

    MAE =
      mean(
        abs(predicted - observed),
        na.rm = TRUE
      ),

    RMSE =
      sqrt(
        mean(
          (predicted - observed)^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted - observed,
        na.rm = TRUE
      ),

    rho =
      cor(
        predicted,
        observed,
        method = "spearman",
        use = "complete.obs"
      )
  )

gam_cv_final
 


gam_cv_final_by_group <- cv_final |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    n = n(),

    MAE =
      mean(
        abs(predicted - observed),
        na.rm = TRUE
      ),

    RMSE =
      sqrt(
        mean(
          (predicted - observed)^2,
          na.rm = TRUE
        )
      ),

    bias =
      mean(
        predicted - observed,
        na.rm = TRUE
      ),

    rho =
      cor(
        predicted,
        observed,
        method = "spearman",
        use = "complete.obs"
      ),

    .groups = "drop"
  )

gam_cv_final_by_group


# =============================================================================
# Final predictive GAM
#
# Used to predict composition for both observed and previously unseen fleets.
# this is the GAM species-composition reconstruction model, while gam4_final is the explanatory/fitted model.

# =============================================================================

gam_production <- gam(
  gam_predictive_formula,

  family =
    betar(link = "logit"),

  method =
    "REML",

  select =
    TRUE,

  data =
    gam_train
)

summary(gam_production)
gam.check(gam_production)


# =============================================================================
# Generic-shark prediction universe
# =============================================================================

gam_nei_prediction_data <- gam_shark_data |>
  filter(
    generic_shark_mt > 0
  ) |>
  mutate(

    FLEET_CODE =
      as.character(FLEET_CODE),

    FISHERY_CODE =
      as.character(FISHERY_CODE),

    FISHERY_GROUP_CODE_ORIG =
      as.character(FISHERY_GROUP_CODE),

    gam_group_supported =
      FISHERY_GROUP_CODE_ORIG %in%
        gam_supported_groups_final
  )

#========================================================================
# This is what we can predict. 
gam_nei_supported <- gam_nei_prediction_data |>
  filter(
    gam_group_supported
  ) |>
  mutate(

    FISHERY_GROUP_CODE =
      factor(
        FISHERY_GROUP_CODE_ORIG,
        levels =
          levels(gam_train$FISHERY_GROUP_CODE)
      )
  )


gam_pred_link <- predict(
  gam_production,
  newdata = gam_nei_supported,
  type = "link",
  se.fit = TRUE
)

gam_nei_supported <- gam_nei_supported |>
  mutate(

    gam_link =
      as.numeric(
        gam_pred_link$fit
      ),

    gam_se =
      as.numeric(
        gam_pred_link$se.fit
      ),

    p_bsh_gam =
      plogis(gam_link),

    p_bsh_gam_lower =
      plogis(
        gam_link -
          1.96 * gam_se
      ),

    p_bsh_gam_upper =
      plogis(
        gam_link +
          1.96 * gam_se
      ),

    bsh_from_nei_gam =
      generic_shark_mt *
      p_bsh_gam,

    bsh_from_nei_gam_lower =
      generic_shark_mt *
      p_bsh_gam_lower,

    bsh_from_nei_gam_upper =
      generic_shark_mt *
      p_bsh_gam_upper,

    FISHERY_GROUP_CODE =
      FISHERY_GROUP_CODE_ORIG
  )

gam_nei_unsupported <- gam_nei_prediction_data |>
  filter(
    !gam_group_supported
  ) |>
  mutate(

    FISHERY_GROUP_CODE =
      FISHERY_GROUP_CODE_ORIG,

    gam_link =
      NA_real_,

    gam_se =
      NA_real_,

    p_bsh_gam =
      NA_real_,

    p_bsh_gam_lower =
      NA_real_,

    p_bsh_gam_upper =
      NA_real_,

    bsh_from_nei_gam =
      NA_real_,

    bsh_from_nei_gam_lower =
      NA_real_,

    bsh_from_nei_gam_upper =
      NA_real_
  )

gam_nei_predictions <- bind_rows(
  gam_nei_supported,
  gam_nei_unsupported
)



gam_training_support <- gam_train |>
  mutate(
    FISHERY_GROUP_CODE =
      as.character(FISHERY_GROUP_CODE)
  ) |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(

    donor_first_year =
      min(YEAR),

    donor_last_year =
      max(YEAR),

    donor_n =
      n(),

    donor_fleets =
      n_distinct(FLEET_CODE),

    .groups = "drop"
  )


gam_nei_predictions <- gam_nei_predictions |>
  left_join(
    gam_training_support,
    by = "FISHERY_GROUP_CODE"
  ) |>
  mutate(

    gam_temporal_support = case_when(

      !gam_group_supported ~
        "unsupported_fishery_group",

      YEAR < donor_first_year ~
        "backcast",

      YEAR > donor_last_year ~
        "forward_extrapolation",

      TRUE ~
        "within_training_period"
    )
  )


training_fleets <- gam_train |>
  distinct(
    FLEET_CODE
  ) |>
  pull(
    FLEET_CODE
  ) |>
  as.character()


gam_nei_predictions <- gam_nei_predictions |>
  mutate(

    fleet_seen_in_training =
      FLEET_CODE %in%
        training_fleets
  )


gam_coverage_summary <- gam_nei_predictions |>
  summarise(

    total_generic_shark_mt =
      sum(
        generic_shark_mt,
        na.rm = TRUE
      ),

    gam_supported_generic_mt =
      sum(
        generic_shark_mt[
          gam_group_supported
        ],
        na.rm = TRUE
      ),

    gam_unsupported_generic_mt =
      sum(
        generic_shark_mt[
          !gam_group_supported
        ],
        na.rm = TRUE
      ),

    pct_generic_supported =
      100 *
      gam_supported_generic_mt /
      total_generic_shark_mt,

    gam_bsh_estimated_mt =
      sum(
        bsh_from_nei_gam,
        na.rm = TRUE
      )
  )

gam_coverage_summary

gam_cv_final
gam_cv_final_by_group
gam_coverage_summary
bsh_gam_hscr_by_group

# some comparison plotting 
# =============================================================================
# Rebuild BSH GAM-HSCR comparison by fishery group
# =============================================================================

bsh_gam_hscr_by_group <-
  bsh_gam_hscr_comparison |>
  group_by(
    FISHERY_GROUP_CODE
  ) |>
  summarise(
    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    difference_GAM_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt,

    difference_GAM_N2_mt =
      bsh_GAM_mt - bsh_HSCR_N2_mt
  )

bsh_gam_hscr_by_group

# =============================================================================
# BSH GAM allocation by temporal support
# =============================================================================

bsh_gam_temporal_support <-
  bsh_gam_hscr_comparison |>
  group_by(
    FISHERY_GROUP_CODE,
    gam_temporal_support
  ) |>
  summarise(
    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    GAM_minus_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt,

    GAM_minus_N2_mt =
      bsh_GAM_mt - bsh_HSCR_N2_mt
  )

bsh_gam_temporal_support
# =============================================================================
# BSH GAM allocation by temporal support
# =============================================================================

bsh_gam_temporal_support <-
  bsh_gam_hscr_comparison |>
  group_by(
    FISHERY_GROUP_CODE,
    gam_temporal_support
  ) |>
  summarise(
    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    GAM_minus_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt,

    GAM_minus_N2_mt =
      bsh_GAM_mt - bsh_HSCR_N2_mt
  )

bsh_gam_temporal_support

bsh_gam_temporal_support_total <-
  bsh_gam_hscr_comparison |>
  group_by(
    gam_temporal_support
  ) |>
  summarise(
    n_strata = n(),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    proportion_GAM =
      bsh_GAM_mt /
      sum(bsh_GAM_mt),

    GAM_minus_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt
  )

bsh_gam_temporal_support_total

# =============================================================================
# Detailed OT diagnostic
# =============================================================================

bsh_ot_gam_diagnostic <-
  bsh_gam_hscr_comparison |>
  filter(
    FISHERY_GROUP_CODE == "OT"
  ) |>
  group_by(
    gam_temporal_support,
    fleet_seen_in_training
  ) |>
  summarise(
    n_strata = n(),

    first_year =
      min(YEAR, na.rm = TRUE),

    last_year =
      max(YEAR, na.rm = TRUE),

    mean_p_bsh_gam =
      weighted.mean(
        p_bsh_gam,
        w = generic_shark_mt,
        na.rm = TRUE
      ),

    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    .groups = "drop"
  ) |>
  mutate(
    GAM_minus_N1_mt =
      bsh_GAM_mt - bsh_HSCR_N1_mt,

    GAM_minus_N2_mt =
      bsh_GAM_mt - bsh_HSCR_N2_mt
  )

bsh_ot_gam_diagnostic


bsh_ot_top_strata <-
  bsh_gam_hscr_comparison |>
  filter(
    FISHERY_GROUP_CODE == "OT"
  ) |>
  mutate(
    difference_N1 =
      bsh_from_nei_gam -
      bsh_from_nei_N1
  ) |>
  arrange(
    desc(abs(difference_N1))
  ) |>
  select(
    YEAR,
    FLEET_CODE,
    FISHERY_CODE,
    generic_shark_mt,

    p_bsh_N1,
    p_bsh_N2,
    p_bsh_gam,

    bsh_from_nei_N1,
    bsh_from_nei_N2,
    bsh_from_nei_gam,

    difference_N1,

    gam_temporal_support,
    fleet_seen_in_training,

    donor_first_year,
    donor_last_year,
    donor_n,
    donor_fleets
  )

print(
  bsh_ot_top_strata,
  n = 30
)

bsh_ot_annual <-
  bsh_gam_hscr_comparison |>
  filter(
    FISHERY_GROUP_CODE == "OT"
  ) |>
  group_by(
    YEAR
  ) |>
  summarise(
    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    HSCR_N1 =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    HSCR_N2 =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    GAM =
      sum(bsh_from_nei_gam, na.rm = TRUE),

    .groups = "drop"
  )
# =============================================================================
# Annual BSH HSCR and GAM comparison
# =============================================================================

bsh_skh_method_annual <-
  bsh_gam_hscr_comparison |>

  group_by(
    YEAR
  ) |>

  summarise(

    generic_shark_mt =
      sum(
        generic_shark_mt,
        na.rm = TRUE
      ),

    bsh_HSCR_N1_mt =
      sum(
        bsh_from_nei_N1,
        na.rm = TRUE
      ),

    bsh_HSCR_N2_mt =
      sum(
        bsh_from_nei_N2,
        na.rm = TRUE
      ),

    bsh_GAM_mt =
      sum(
        bsh_from_nei_gam,
        na.rm = TRUE
      ),

    bsh_GAM_lower_mt =
      sum(
        bsh_from_nei_gam_lower,
        na.rm = TRUE
      ),

    bsh_GAM_upper_mt =
      sum(
        bsh_from_nei_gam_upper,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>

  mutate(

    hscr_lower_mt =
      pmin(
        bsh_HSCR_N1_mt,
        bsh_HSCR_N2_mt
      ),

    hscr_upper_mt =
      pmax(
        bsh_HSCR_N1_mt,
        bsh_HSCR_N2_mt
      ),

    hscr_midpoint_mt =
      (
        bsh_HSCR_N1_mt +
        bsh_HSCR_N2_mt
      ) / 2
  )

#plot
p_bsh_skh_annual_comparison <-
  ggplot(
    bsh_skh_method_annual,
    aes(
      x = YEAR
    )
  ) +

  geom_ribbon(
    aes(
      ymin = hscr_lower_mt,
      ymax = hscr_upper_mt,
      fill = "HSCR N1–N2"
    ),
    alpha = 0.25
  ) +

  geom_line(
    aes(
      y = hscr_midpoint_mt,
      colour = "HSCR midpoint"
    ),
    linewidth = 0.9
  ) +

geom_line(
  aes(
    y = bsh_GAM_mt,
    colour = "GAM sensitivity"
  ),
  linewidth = 0.9
) +

  scale_colour_viridis_d(
    option = "D",
    end = 0.85
  ) +

  scale_fill_viridis_d(
    option = "D",
    end = 0.85
  ) +

  scale_y_continuous(
    labels =
      scales::label_number(
        big.mark = ","
      )
  ) +

  labs(
    x = "Year",
    y = "Blue shark allocated from Sharks NEI (t)",
    colour = NULL,
    fill = NULL,
    title =
      "Blue shark allocation from generic shark catch",
    subtitle =
      "HSCR and GAM species-composition alternatives"
  ) +

  theme_bw() +

  theme(
    legend.position =
      "bottom"
  )


p_bsh_skh_annual_comparison
dir.create(
  "outputs/bsh/gam/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename =
    "outputs/bsh/gam/figures/bsh_skh_hscr_gam_annual_comparison.png",
  plot =
    p_bsh_skh_annual_comparison,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

write_csv(
  bsh_skh_annual_comparison,
  "outputs/bsh/gam/bsh_skh_hscr_gam_annual_comparison.csv"
)


# =============================================================================
# BSH HSCR and GAM comparison by fishery group
# =============================================================================

bsh_skh_by_group <-
  bsh_skh_gam_prediction |>

  group_by(
    FISHERY_GROUP_CODE
  ) |>

  summarise(

    skh_mt =
      sum(
        skh_mt,
        na.rm = TRUE
      ),

    hscr_N1_mt =
      sum(
        bsh_from_skh_N1,
        na.rm = TRUE
      ),

    hscr_N2_mt =
      sum(
        bsh_from_skh_N2,
        na.rm = TRUE
      ),

    gam_sensitivity_mt =
      sum(
        bsh_from_skh_gam_sensitivity,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>

  mutate(

    gam_minus_N1_mt =
      gam_sensitivity_mt -
      hscr_N1_mt,

    gam_minus_N2_mt =
      gam_sensitivity_mt -
      hscr_N2_mt
  )


print(
  bsh_skh_by_group,
  n = Inf
)

write_csv(
  bsh_skh_by_group,
  "outputs/bsh/gam/bsh_skh_hscr_gam_group_comparison.csv"
)


# =============================================================================
# BSH manuscript output directories
# =============================================================================

bsh_gam_output_dir <-
  "outputs/bsh/gam"

bsh_gam_figure_dir <-
  file.path(
    bsh_gam_output_dir,
    "figures"
  )

dir.create(
  bsh_gam_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  bsh_gam_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# Save BSH HSCR-GAM comparison tables
# =============================================================================

write_csv(
  bsh_gam_hscr_summary,
  file.path(
    bsh_gam_output_dir,
    "bsh_skh_hscr_gam_summary.csv"
  )
)

write_csv(
  bsh_gam_hscr_by_group,
  file.path(
    bsh_gam_output_dir,
    "bsh_skh_hscr_gam_by_group.csv"
  )
)

write_csv(
  bsh_group_period_comparison,
  file.path(
    bsh_gam_output_dir,
    "bsh_skh_hscr_gam_group_period_comparison.csv"
  )
)

write_csv(
  gam_cv_by_group,
  file.path(
    bsh_gam_output_dir,
    "bsh_skh_gam_cv_by_group.csv"
  )
)

write_csv(
  bsh_gam_hscr_comparison,
  file.path(
    bsh_gam_output_dir,
    "bsh_skh_hscr_gam_comparison.csv"
  )
)



bsh_gam_cv_by_group_long <-
  gam_cv_by_group |>
  dplyr::select(
    FISHERY_GROUP_CODE,
    MAE,
    RMSE
  ) |>
  tidyr::pivot_longer(
    cols = c(MAE, RMSE),
    names_to = "metric",
    values_to = "value"
  )

p_bsh_skh_gam_cv_by_group <- ggplot(
  bsh_gam_cv_by_group_long,
  aes(
    x = FISHERY_GROUP_CODE,
    y = value,
    fill = metric
  )
) +
  geom_col(
    position = "dodge"
  ) +
  scale_fill_viridis_d(
    option = "D"
  ) +
  labs(
    x = "Fishery group",
    y = "Prediction error",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


ggsave(
  filename =
    "outputs/bsh/gam/figures/bsh_skh_gam_cv_by_group.png",
  plot =
    p_bsh_skh_gam_cv_by_group,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)


 
ggsave(
  filename =
    "outputs/bsh/gam/figures/bsh_lofo_cv_bsh.png",
  plot =
    p_lofo_cv_bsh,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

# =============================================================================
# Identify canonical BSH HSCR annual source
# =============================================================================

candidate_objects <-
  c(
    "bsh_nei_candidates",
    "bsh_nei_component",
    "nei_scenario_annual",
    "nei_scenario_by_group",
    "nei_scenario_summary",
    "hscr_comparison_dat",
    "hscr_gam_comparison",
    "gam_hscr_comparison"
  )

for (obj in candidate_objects) {

  cat(
    "\n\n========================================\n",
    obj,
    "\n========================================\n"
  )

  x <- get(obj)

  print(
    names(x)
  )
}


# =============================================================================
# Consistency checks against canonical totals
# =============================================================================

bsh_skh_method_total_check <-
  bsh_skh_method_annual |>

  summarise(

    generic_shark_mt =
      sum(
        generic_shark_mt,
        na.rm = TRUE
      ),

    bsh_HSCR_N1_mt =
      sum(
        bsh_HSCR_N1_mt,
        na.rm = TRUE
      ),

    bsh_HSCR_N2_mt =
      sum(
        bsh_HSCR_N2_mt,
        na.rm = TRUE
      ),

    bsh_GAM_mt =
      sum(
        bsh_GAM_mt,
        na.rm = TRUE
      )
  )


print(
  bsh_skh_method_total_check
)





# check 
# helper
check_hscr_totals <- function(dat) {

  nm <- names(dat)

  cat("\nColumns containing N1/N2:\n")
  print(
    nm[grepl("N1|N2", nm)]
  )

  invisible(dat)
}

check_hscr_totals(hscr_gam_matched)
check_hscr_totals(gam_nei_supported)
check_hscr_totals(gam_nei_unsupported)
check_hscr_totals(supported_rows)



# ============================================================
# BSH-only object audit
# ============================================================

bsh_objects <- c(
  "bsh_nei_candidates",
  "bsh_nei_component",
  "nei_scenario_annual",
  "nei_scenario_by_group",
  "nei_scenario_summary",
  "hscr_comparison_dat",
  "gam_hscr_comparison",
  "bsh_gam_hscr_summary",
  "bsh_gam_hscr_by_group",
  "gam_nei_supported",
  "gam_nei_unsupported"
)

for (obj in bsh_objects) {

  cat(
    "\n\n========================================\n",
    obj,
    "\n========================================\n"
  )

  x <- get(obj)

  cat("Class:\n")
  print(class(x))

  cat("\nDimensions:\n")
  if (is.data.frame(x)) {
    print(dim(x))
  }

  cat("\nRelevant columns:\n")
  if (is.data.frame(x)) {
    print(
      names(x)[
        grepl(
          "bsh|N1|N2|generic|gam|hscr",
          names(x),
          ignore.case = TRUE
        )
      ]
    )
  }
}


#  these are the bsh objects
# # Full HSCR reconstruction
# hscr_comparison_dat
# nei_scenario_annual
# nei_scenario_by_group
# nei_scenario_summary

# # GAM reconstruction / comparison
# gam_nei_supported
# gam_nei_unsupported
# gam_hscr_comparison
# gam_hscr_summary
# gam_hscr_by_group

# ============================================================
# Identify HSCR strata absent from GAM comparison
# ============================================================

bsh_missing_from_gam <-
  hscr_comparison_dat |>
  anti_join(
    gam_hscr_comparison,
    by = c(
      "YEAR",
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "FISHERY_CODE"
    )
  )

bsh_missing_from_gam |>
  select(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    generic_shark_mt,
    p_bsh_N1,
    p_bsh_N2,
    bsh_from_nei_N1,
    bsh_from_nei_N2,
    nei_N1_source,
    nei_N2_source
  )


bsh_hscr_full <-
  hscr_comparison_dat |>
  summarise(
    N1 =
      sum(bsh_from_nei_N1, na.rm = TRUE),
    N2 =
      sum(bsh_from_nei_N2, na.rm = TRUE)
  )

bsh_hscr_gam_subset <-
  gam_hscr_comparison |>
  summarise(
    N1 =
      sum(bsh_from_nei_N1, na.rm = TRUE),
    N2 =
      sum(bsh_from_nei_N2, na.rm = TRUE)
  )

bsh_hscr_full
bsh_hscr_gam_subset



# =============================================================================
# Compare HSCR values before and after construction of gam_hscr_comparison
# =============================================================================

bsh_hscr_value_check <-
  hscr_comparison_dat |>
  select(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    generic_shark_mt,
    p_bsh_N1_original = p_bsh_N1,
    p_bsh_N2_original = p_bsh_N2,
    bsh_N1_original = bsh_from_nei_N1,
    bsh_N2_original = bsh_from_nei_N2
  ) |>

  inner_join(
    gam_hscr_comparison |>
      select(
        YEAR,
        FLEET_CODE,
        FISHERY_GROUP_CODE,
        FISHERY_CODE,
        p_bsh_N1_gam_object = p_bsh_N1,
        p_bsh_N2_gam_object = p_bsh_N2,
        bsh_N1_gam_object = bsh_from_nei_N1,
        bsh_N2_gam_object = bsh_from_nei_N2
      ),
    by = c(
      "YEAR",
      "FLEET_CODE",
      "FISHERY_GROUP_CODE",
      "FISHERY_CODE"
    )
  ) |>

  mutate(
    diff_p_N1 =
      p_bsh_N1_gam_object - p_bsh_N1_original,

    diff_p_N2 =
      p_bsh_N2_gam_object - p_bsh_N2_original,

    diff_mt_N1 =
      bsh_N1_gam_object - bsh_N1_original,

    diff_mt_N2 =
      bsh_N2_gam_object - bsh_N2_original
  )


bsh_hscr_value_check |>
  summarise(
    n = n(),

    original_N1 =
      sum(bsh_N1_original, na.rm = TRUE),

    gam_object_N1 =
      sum(bsh_N1_gam_object, na.rm = TRUE),

    difference_N1 =
      sum(diff_mt_N1, na.rm = TRUE),

    original_N2 =
      sum(bsh_N2_original, na.rm = TRUE),

    gam_object_N2 =
      sum(bsh_N2_gam_object, na.rm = TRUE),

    difference_N2 =
      sum(diff_mt_N2, na.rm = TRUE),

    n_different_N1 =
      sum(
        abs(diff_mt_N1) > 1e-8,
        na.rm = TRUE
      ),

    n_different_N2 =
      sum(
        abs(diff_mt_N2) > 1e-8,
        na.rm = TRUE
      )
  )


bsh_hscr_value_check |>
  filter(
    abs(diff_mt_N1) > 1e-8 |
      abs(diff_mt_N2) > 1e-8
  ) |>
  arrange(
    desc(
      pmax(
        abs(diff_mt_N1),
        abs(diff_mt_N2),
        na.rm = TRUE
      )
    )
  ) |>
  select(
    YEAR,
    FLEET_CODE,
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    generic_shark_mt,
    p_bsh_N1_original,
    p_bsh_N1_gam_object,
    p_bsh_N2_original,
    p_bsh_N2_gam_object,
    bsh_N1_original,
    bsh_N1_gam_object,
    diff_mt_N1,
    bsh_N2_original,
    bsh_N2_gam_object,
    diff_mt_N2
  ) |>
  print(n = 30)


# so where is the differnce 
gam_hscr_comparison |>
  summarise(
    generic_shark_mt =
      sum(generic_shark_mt, na.rm = TRUE),

    bsh_HSCR_N1_mt =
      sum(bsh_from_nei_N1, na.rm = TRUE),

    bsh_HSCR_N2_mt =
      sum(bsh_from_nei_N2, na.rm = TRUE),

    bsh_GAM_mt =
      sum(bsh_from_nei_gam, na.rm = TRUE)
  )


read_csv(
  file.path(
    "outputs/bsh/gam",
    "bsh_skh_hscr_gam_summary.csv"
  ),
  show_col_types = FALSE
)

