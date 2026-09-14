# R/processing/harmonize_iotc_codes.R

library(dplyr)
library(tidyr)

# -------------------------------------------------------------------
# Extract unique fishery / gear definitions from each dataset
# -------------------------------------------------------------------

# read in the RC, CA and EF data, 

rc_file <- "data/raw/rc/IOTC-DATASETS-2026-08-13-RC-RAW-1950-2025.csv"
rc_dat <- readr::read_csv(rc_file )

# Catch
ca_file <- "data/raw/ca/IOTC-DATASETS-2026-08-13-CA-1952-2025.csv"
ca_dat <- readr::read_csv(ca_file )

# Effort 
ef_file <- "data/raw/ef/IOTC-DATASETS-2026-08-13-EF-1952-2025.csv"
ef_dat <- readr::read_csv(ef_file )



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


label_check <- bind_rows(
  rc_codes,
  ca_codes,
  ef_codes
) |>
  group_by(
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    GEAR_CODE
  ) |>
  summarise(
    n_group_labels   = n_distinct(FISHERY_GROUP),
    n_fishery_labels = n_distinct(FISHERY),
    n_gear_labels    = n_distinct(GEAR),
    datasets         = paste(sort(unique(DATASET)), collapse = ", "),
    .groups = "drop"
  ) |>
  filter(
    n_group_labels > 1 |
    n_fishery_labels > 1 |
    n_gear_labels > 1
  )

print(label_check, n = Inf)
