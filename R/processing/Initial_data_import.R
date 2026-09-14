# Joel Rice
# Setup code
# 


getwd()
 # "c:/Users/JSR/Dropbox/JoelWork/IOTC/BSH_SMA_CATCH_WPEB_22_2026/IOTC_Shark_Catch_Reconstruction"

# -------------------------------------------------------------------------
# Create project folder structure
# -------------------------------------------------------------------------

# Folder names
r_dirs <- c(
  "R/functions",
  "R/import",
  "R/processing",
  "R/reconstruction",
  "R/diagnostics",
  "R/plotting",
  "R/tests"
)

# Create folders if they do not already exist
invisible(lapply(r_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Add .gitkeep files so Git tracks empty folders
invisible(
  lapply(
    r_dirs,
    function(x) {
      gitkeep <- file.path(x, ".gitkeep")
      if (!file.exists(gitkeep))
        file.create(gitkeep)
    }
  )
)

cat("R project folders successfully created.\n")



# read in the RC, CA and EF data, 

rc_file <- "data/raw/rc/IOTC-DATASETS-2026-08-13-RC-RAW-1950-2025.csv"
rc_dat <- readr::read_csv(rc_file )

# Catch
ca_file <- "data/raw/ca/IOTC-DATASETS-2026-08-13-CA-1952-2025.csv"
ca_dat <- readr::read_csv(ca_file )

# Effort 
ef_file <- "data/raw/ef/IOTC-DATASETS-2026-08-13-EF-1952-2025.csv"
ef_dat <- readr::read_csv(ef_file )


write.table(head(ca_dat, n=15), file="clipboard", sep = "\t", row.names = FALSE)

write.table(head(ef_dat, n=15), file="clipboard", sep = "\t", row.names = FALSE)


unique(ca_dat$CATCH_UNIT_CODE)
unique(ef_dat$EFFORT_UNIT_CODE)

unique(ca_dat$FISHERY_GROUP_CODE)
unique(ef_dat$FISHERY_GROUP_CODE)

range(ca_dat$YEAR)
range(ef_dat$YEAR)

table(ca_dat$SPECIES_CODE %in% c("BSH", "SMA"))
table(ef_dat$GEAR_CODE)

ca_dat |>
  dplyr::filter(SPECIES_CODE %in% c("BSH", "SMA")) |>
  dplyr::count(SPECIES_CODE, CATCH_UNIT_CODE, sort = TRUE)

print(n=40,ef_dat |>
  dplyr::count(FISHERY_GROUP_CODE, GEAR_CODE, EFFORT_UNIT_CODE, sort = TRUE) )


shark_ca_summary <-
  ca_dat |>
  dplyr::filter(SPECIES_CODE %in% c("BSH", "SMA")) |>
  dplyr::count(
    SPECIES_CODE,
    FISHERY_GROUP_CODE,
    GEAR_CODE,
    CATCH_UNIT_CODE,
    sort = TRUE
  )

print(shark_ca_summary, n = 100)

# effor summary 
shark_fisheries <-
  ca_dat |>
  dplyr::filter(SPECIES_CODE %in% c("BSH", "SMA")) |>
  dplyr::distinct(
    FISHERY_GROUP_CODE,
    GEAR_CODE
  )

shark_effort_summary <-
  ef_dat |>
  dplyr::semi_join(
    shark_fisheries,
    by = c("FISHERY_GROUP_CODE", "GEAR_CODE")
  ) |>
  dplyr::count(
    FISHERY_GROUP_CODE,
    GEAR_CODE,
    EFFORT_UNIT_CODE,
    sort = TRUE
  )

print(shark_effort_summary, n = 100)


shark_temporal_coverage <-
  ca_dat |>
  dplyr::filter(SPECIES_CODE %in% c("BSH", "SMA")) |>
  dplyr::group_by(
    SPECIES_CODE,
    FISHERY_GROUP_CODE,
    GEAR_CODE,
    CATCH_UNIT_CODE
  ) |>
  dplyr::summarise(
    first_year = min(YEAR, na.rm = TRUE),
    last_year  = max(YEAR, na.rm = TRUE),
    n_years    = dplyr::n_distinct(YEAR),
    n_records  = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    SPECIES_CODE,
    dplyr::desc(n_records)
  )

print(shark_temporal_coverage, n = 100)

 rc_geatr_group_tbl <- rc_dat |>
  dplyr::distinct(
    FISHERY_GROUP_CODE,
    FISHERY_GROUP,
    FISHERY_CODE,
    FISHERY,
    GEAR_CODE,
    GEAR
  ) |>
  dplyr::arrange(
    FISHERY_GROUP_CODE,
    FISHERY_CODE,
    GEAR_CODE
  )

print(rc_geatr_group_tbl, n = 100)
