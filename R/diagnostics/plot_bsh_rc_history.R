# Diagnostic script to reproduce previous catch plots. 


library(dplyr)
library(ggplot2)
library(scales)

# -------------------------------------------------------------------
# Blue shark retained catch by broad IOTC fishery group
# 2026 retained catch dataset
# -------------------------------------------------------------------

bsh_rc_by_group <- rc_dat |>
  filter(SPECIES_CODE == "BSH") |>
  group_by(
    YEAR,
    FISHERY_GROUP_CODE,
    FISHERY_GROUP
  ) |>
  summarise(
    catch_mt = sum(CATCH, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------
# Check annual totals
# -------------------------------------------------------------------

bsh_rc_annual <- bsh_rc_by_group |>
  group_by(YEAR) |>
  summarise(
    catch_mt = sum(catch_mt, na.rm = TRUE),
    .groups = "drop"
  )

print(bsh_rc_annual, n = Inf)


# -------------------------------------------------------------------
# Stacked annual retained catch
# -------------------------------------------------------------------

p_bsh_rc_history <- ggplot(
  bsh_rc_by_group,
  aes(
    x = YEAR,
    y = catch_mt / 1000,
    fill = FISHERY_GROUP
  )
) +
  geom_col(width = 0.9) +
  scale_x_continuous(
    breaks = seq(
      floor(min(bsh_rc_by_group$YEAR) / 10) * 10,
      max(bsh_rc_by_group$YEAR),
      by = 10
    ),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Year",
    y = "Retained BSH catch (1000 MT)",
    fill = "Fishery group",
    title = "Indian Ocean blue shark retained catch",
    subtitle = "IOTC retained catch dataset, 2026 release"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p_bsh_rc_history


bsh_rc_by_group_code <- rc_dat |>
  filter(SPECIES_CODE == "BSH") |>
  group_by(
    YEAR,
    FISHERY_GROUP_CODE
  ) |>
  summarise(
    catch_mt = sum(CATCH, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    FISHERY_GROUP_CODE = factor(
      FISHERY_GROUP_CODE,
      levels = c("GN", "LI", "LL", "PS", "BB", "OT")
    )
  )

p_bsh_rc_history_code <- ggplot(
  bsh_rc_by_group_code,
  aes(
    x = YEAR,
    y = catch_mt / 1000,
    fill = FISHERY_GROUP_CODE
  )
) +
  geom_col(width = 0.9) +
  scale_x_continuous(
    breaks = seq(1950, 2025, by = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Year",
    y = "Retained BSH catch (1000 MT)",
    fill = "Fishery group",
    title = "Indian Ocean blue shark retained catch",
    subtitle = "IOTC-DATASETS-2026-08-13-RC-RAW-1950-2025"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p_bsh_rc_history_code

# Total BSH retained catch by broad fishery group
bsh_group_totals <- bsh_rc_by_group_code |>
  group_by(FISHERY_GROUP_CODE) |>
  summarise(
    total_catch_mt = sum(catch_mt),
    first_year = min(YEAR[catch_mt > 0]),
    last_year = max(YEAR[catch_mt > 0]),
    .groups = "drop"
  ) |>
  arrange(desc(total_catch_mt))

print(bsh_group_totals)

# Annual contribution by group as proportions
bsh_group_prop <- bsh_rc_by_group_code |>
  group_by(YEAR) |>
  mutate(
    annual_total_mt = sum(catch_mt),
    proportion = catch_mt / annual_total_mt
  ) |>
  ungroup()

# Find years with the largest year-to-year changes in total retained catch
bsh_annual_change <- bsh_rc_annual |>
  arrange(YEAR) |>
  mutate(
    change_mt = catch_mt - lag(catch_mt),
    pct_change = 100 * (catch_mt / lag(catch_mt) - 1)
  ) |>
  arrange(desc(abs(change_mt)))

print(bsh_annual_change, n = 20)
