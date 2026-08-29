


# =============================================================================
# RC reporting history: SKH, BSH, SMA and MAK
# =============================================================================

rc_shark_reporting_annual <-
  rc_dat |>

  filter(
    SPECIES_CODE %in%
      c(
        "SKH",
        "BSH",
        "SMA",
        "MAK"
      )
  ) |>

  group_by(
    YEAR,
    SPECIES_CODE
  ) |>

  summarise(
    catch_mt =
      sum(
        CATCH,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>

  mutate(
    reported_category =
      recode(
        SPECIES_CODE,
        SKH = "Sharks NEI",
        BSH = "Blue shark",
        SMA = "Shortfin mako",
        MAK = "Mako sharks NEI"
      ),

    reported_category =
      factor(
        reported_category,
        levels = c(
          "Sharks NEI",
          "Mako sharks NEI",
          "Blue shark",
          "Shortfin mako"
        )
      )
  )



p_rc_shark_reporting_annual <-
  ggplot(
    rc_shark_reporting_annual,
    aes(
      x = YEAR,
      y = catch_mt,
      fill = reported_category
    )
  ) +

  geom_col(
    width = 0.9
  ) +

  facet_wrap(
    ~ reported_category,
    ncol = 1,
    scales = "free_y"
  ) +

  scale_fill_viridis_d(
    option = "D",
    end = 0.9,
    guide = "none"
  ) +

  labs(
    x = "Year",
    y = "Reported retained catch (t)",
    title =
      "Historical reporting of shark catch categories",
    subtitle =
      "IOTC retained catch data"
  ) +

  theme_bw()


p_rc_shark_reporting_annual


rc_shark_reporting_period <-
  rc_dat |>

  add_reporting_period() |>

  filter(
    SPECIES_CODE %in%
      c(
        "SKH",
        "BSH",
        "SMA",
        "MAK"
      )
  ) |>

  group_by(
    period,
    SPECIES_CODE
  ) |>

  summarise(
    catch_mt =
      sum(
        CATCH,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>

  mutate(
    reported_category =
      recode(
        SPECIES_CODE,
        SKH = "Sharks NEI",
        BSH = "Blue shark",
        SMA = "Shortfin mako",
        MAK = "Mako sharks NEI"
      )
  )


p_rc_shark_reporting_period <-
  ggplot(
    rc_shark_reporting_period,
    aes(
      x = period,
      y = catch_mt,
      fill = reported_category
    )
  ) +

  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +

  scale_fill_viridis_d(
    option = "D",
    end = 0.9
  ) +

  scale_y_continuous(
    labels =
      label_number(
        big.mark = ","
      )
  ) +

  labs(
    x = "Reporting period",
    y = "Reported retained catch (t)",
    fill = "Reported category",
    title =
      "Reported shark catch by reporting period"
  ) +

  theme_bw() +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    legend.position =
      "bottom"
  )


p_rc_shark_reporting_period


ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "rc_shark_reporting_annual.png"
    ),
  plot =
    p_rc_shark_reporting_annual,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename =
    file.path(
      spatial_figure_dir,
      "rc_shark_reporting_by_period.png"
    ),
  plot =
    p_rc_shark_reporting_period,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)