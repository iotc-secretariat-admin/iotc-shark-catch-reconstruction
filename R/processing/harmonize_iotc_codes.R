# R/processing/harmonize_iotc_codes.R

library(dplyr)
library(tidyr)

# -------------------------------------------------------------------
# Extract unique fishery / gear definitions from each dataset
# -------------------------------------------------------------------

rc_codes <- rc_dat |>
  distinct(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY,
    GEAR_CODE,
    GEAR
  ) |>
  mutate(DATASET = "RC")

ca_codes <- ca_dat |>
  distinct(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY,
    GEAR_CODE,
    GEAR
  ) |>
  mutate(DATASET = "CA")

ef_codes <- ef_dat |>
  distinct(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY,
    GEAR_CODE,
    GEAR
  ) |>
  mutate(DATASET = "EF")


# -------------------------------------------------------------------
# Combined code inventory
# -------------------------------------------------------------------

fishery_code_inventory <- bind_rows(
  rc_codes,
  ca_codes,
  ef_codes
) |>
  arrange(
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    GEAR_CODE,
    DATASET
  )


# -------------------------------------------------------------------
# Identify combinations that do not occur in all three datasets
# -------------------------------------------------------------------

fishery_code_comparison <- fishery_code_inventory |>
  distinct(
    DATASET,
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    GEAR_CODE
  ) |>
  mutate(present = TRUE) |>
  pivot_wider(
    names_from = DATASET,
    values_from = present,
    values_fill = FALSE
  ) |>
  arrange(
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    GEAR_CODE
  )


# Show mismatches
fishery_code_mismatches <- fishery_code_comparison |>
  filter(!(RC & CA & EF))

print(fishery_code_mismatches, n = Inf)