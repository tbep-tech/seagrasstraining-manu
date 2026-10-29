library(tidyverse)
library(here)

funcsurl <- 'https://raw.githubusercontent.com/tbep-tech/seagrasstransect-training-reports/main/R/funcs.R'
funcsfl <- paste(tempdir(), 'funcs.R', sep = '/')
download.file(funcsurl, destfile = funcsfl, quiet = TRUE)
source(funcsfl)

load(file = here('data/trndat.rda'))

yrs <- sort(unique(trndat$yr))
metvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')
scr_metric <- c(Abundance = 'abs', `Blade Length` = 'abs', `Short Shoot Density` = 'pct')

cal <- calibrate_scr_fun(trndat, metric = scr_metric)

# one row per group x year x metric, for both the calibrated score and the
# raw (pre-calibration) absolute deviation -- every participant that year
# counts (no group filtering), consistent with R/tabs.R's training summary
# table: this is a descriptive look at the full data, not a trend model
alldat <- purrr::map(yrs, function(yr) {

  tv <- truvar_fun(trndat, yr)

  scr <- allgrpscr_fun(trndat, yr, tv, cal = cal, raw = TRUE, metric = scr_metric) |>
    dplyr::select(-Total) |>
    tidyr::pivot_longer(dplyr::all_of(metvars), names_to = 'var', values_to = 'val') |>
    dplyr::mutate(type = 'Score')

  dev <- allgrpscr_fun(trndat, yr, tv, raw_diff = TRUE, metric = scr_metric) |>
    tidyr::pivot_longer(dplyr::all_of(metvars), names_to = 'var', values_to = 'val') |>
    dplyr::mutate(type = 'Raw deviation')

  dplyr::bind_rows(scr, dev) |>
    dplyr::mutate(yr = yr)

}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    var  = factor(var, levels = metvars),
    type = factor(type, levels = c('Score', 'Raw deviation'))
  )

# mean +/- 95% CI (t-based, since n is as small as 4-5 groups in a year)
# per type x metric x year
summdat <- alldat |>
  dplyr::summarise(
    n    = dplyr::n(),
    mean = mean(val, na.rm = TRUE),
    se   = sd(val, na.rm = TRUE) / sqrt(n),
    .by  = c(type, var, yr)
  ) |>
  dplyr::mutate(
    lwr = mean - se * qt(0.975, n - 1),
    upr = mean + se * qt(0.975, n - 1)
  )

# 2 x 3 facet: rows are Score / Raw deviation, columns are metric. jittered
# raw values sit underneath (smaller, translucent), the mean +/- 95% CI
# summary is drawn on top. each panel gets its own y-axis -- scores are a
# common 50-100 scale, but the three metrics' raw deviations are on very
# different units (category positions, cm, percent difference)
p <- ggplot(summdat, aes(x = yr, y = mean)) +
  geom_jitter(
    data = alldat, aes(x = yr, y = val),
    width = 0.15, height = 0, size = 0.9, alpha = 0.3, color = 'grey40',
    inherit.aes = FALSE
  ) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), color = '#004F7E', linewidth = 0.6, size = 0.5) +
  facet_grid(type ~ var, scales = 'free_y') +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p
