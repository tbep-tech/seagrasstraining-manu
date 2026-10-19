library(tidyverse)
library(here)
library(patchwork)
library(lme4)
library(lmerTest)
library(tbeptools)

funcsurl <- 'https://raw.githubusercontent.com/tbep-tech/seagrasstransect-training-reports/main/R/funcs.R'
funcsfl <- paste(tempdir(), 'funcs.R', sep = '/')
download.file(funcsurl, destfile = funcsfl, quiet = TRUE)
source(funcsfl)

load(file = here('data/trndat.rda'))
data(file = 'trnlns', package = 'tbeptools')

# leave one site out sensitivity -----------------------------------------

scr_metric <- c(Abundance = 'abs', `Blade Length` = 'abs', `Short Shoot Density` = 'pct')
yrs <- sort(unique(trndat$yr))

devvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')
scrvars <- devvars

devfun <- function(sub, yr) {
  tv <- truvar_fun(sub, yr)
  allgrpscr_fun(sub, yr, tv, raw_diff = TRUE, metric = scr_metric)
}

# fixed historical baseline, built once from the full dataset (same as
# production use) -- not rebuilt per Site dropped, since cal is meant to
# stay fixed regardless of which Sites a given focal year's own data has
cal <- calibrate_scr_fun(trndat, metric = scr_metric)

scrfun <- function(sub, yr) {
  tv <- truvar_fun(sub, yr)
  allgrpscr_fun(sub, yr, tv, cal = cal, raw = TRUE, metric = scr_metric)
}

# a dropped Site can occasionally remove all data for one of the metrics
# entirely, which errors inside allgrpscr_fun -- skip that iteration rather
# than fail the whole run
devfun_safe <- purrr::possibly(devfun, otherwise = NULL)
scrfun_safe <- purrr::possibly(scrfun, otherwise = NULL)

# runs the leave-one-Site-out loop for a given scoring function, returning
# one row per group x year x metric x Site left out, with the full-sample
# and leave-one-out values side by side
run_loo <- function(scorefun_safe, varnames) {

  purrr::map(yrs, function(yr) {

    message('year ', yr)

    yrdat <- trndat |> dplyr::filter(yr == !!yr)
    sites <- sort(unique(yrdat$Site))

    full <- scorefun_safe(yrdat, yr)
    if (is.null(full)) return(NULL)

    full_long <- full |>
      tidyr::pivot_longer(dplyr::all_of(varnames), names_to = 'metric', values_to = 'full')

    purrr::map(sites, function(s) {

      sub <- yrdat |> dplyr::filter(Site != s)
      loo <- scorefun_safe(sub, yr)
      if (is.null(loo)) return(NULL)

      loo |>
        tidyr::pivot_longer(dplyr::all_of(varnames), names_to = 'metric', values_to = 'loo') |>
        dplyr::left_join(full_long, by = c('grpact', 'metric')) |>
        dplyr::mutate(site_out = s, yr = yr)

    }) |>
      purrr::compact() |>
      dplyr::bind_rows()

  }) |>
    purrr::compact() |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      grp    = gsub('^.*:\\s(.*?)\\s\\(.*$', '\\1', grpact),
      metric = factor(metric, levels = varnames)
    ) |>
    dplyr::filter(!is.na(full), !is.na(loo)) |>
    dplyr::filter(!grp %in% 'NA') |>
    dplyr::filter(grp %in% tbeptools::trnlns$MonAgency)

}

# how much does dropping any one Site move a group's value? correlation
# between full-sample and leave-one-Site-out values per metric (groups are
# on different scales by metric, same reasoning as the other tables), plus
# the median/max absolute difference as a plain-units influence measure.
# also builds the full-sample vs leave-one-out scatter, one panel per
# metric, dashed 1:1 line
summarize_loo <- function(loodat) {
  loocor <- loodat |>
    dplyr::summarise(r = cor(full, loo, method = 'pearson'), .by = metric)

  loodat |>
    dplyr::mutate(absdiff = abs(full - loo)) |>
    dplyr::summarise(
      absdiff_median = median(absdiff, na.rm = TRUE),
      absdiff_max    = max(absdiff, na.rm = TRUE),
      .by = metric
    ) |>
    dplyr::left_join(loocor, by = 'metric') |>
    print()

  loocor
}

# formal test for a systematic bias between full-sample and leave-one-Site-
# out values, separately per metric (different scales, same reasoning as
# summarize_loo()). testing whether an ordinary regression's slope differs
# from 1 isn't used here -- OLS assumes the x-axis is measured without
# error, but full and loo are both noisy measurements of the same group x
# year quantity, so OLS is biased toward a shallower slope purely from
# measurement noise (regression dilution), which would flag "significant"
# deviation even for a perfectly unbiased method. instead this tests
# whether the mean paired difference (full - loo) differs from zero, with
# random intercepts for yr and grp accounting for the same group x year
# appearing once per Site dropped that year (not independent observations)
diff_test_loo <- function(loodat) {
  purrr::map(levels(loodat$metric), function(m) {
    dat <- loodat |>
      dplyr::filter(metric == m) |>
      dplyr::mutate(diff = full - loo)
    mod <- lmer(diff ~ 1 + (1 | yr) + (1 | grp), data = dat)
    cfs <- summary(mod)$coefficients
    tibble::tibble(
      metric   = m,
      estimate = cfs['(Intercept)', 'Estimate'],
      se       = cfs['(Intercept)', 'Std. Error'],
      pval     = cfs['(Intercept)', 'Pr(>|t|)']
    )
  }) |>
    dplyr::bind_rows()
}

# one full-sample vs leave-one-Site-out panel for a single metric, used to
# build the combined 2 x 3 figure below (rows: raw deviation, score;
# columns: metric) -- one ggplot object per metric rather than facet_wrap,
# so the two rows (different datasets) can be stacked with patchwork
make_loo_panel <- function(loodat, loocor, m, xlab, ylab, title) {

  datm <- loodat |> dplyr::filter(metric == m)
  r    <- loocor$r[loocor$metric == m]
  lab  <- sprintf('r = %.2f', r)
  rng  <- range(c(datm$full, datm$loo))

  ggplot(datm, aes(x = full, y = loo)) +
    geom_point(aes(color = factor(yr)), alpha = 0.8, size = 1.3) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey30') +
    annotate('text', x = rng[1], y = rng[2], label = lab, hjust = 0, vjust = 1, size = 3.2) +
    scale_color_viridis_d(name = 'Year') +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal(base_size = 11) + 
    theme(
      panel.grid.minor = element_blank()
    )

}

leaveoutdev <- run_loo(devfun_safe, devvars)
leaveoutdevcor <- summarize_loo(leaveoutdev)
diff_test_loo(leaveoutdev) |> print()

leaveoutscr <- run_loo(scrfun_safe, scrvars)
leaveoutscrcor <- summarize_loo(leaveoutscr)
diff_test_loo(leaveoutscr) |> print()

# combined 2 x 3 figure: top row raw deviation, bottom row calibrated
# score, columns = metric. y-axis label dropped on the 2nd/3rd columns
# (redundant once shown on the left-most panel), same convention as
# figs.r's combined score/deviation figure
devpanels <- purrr::map(devvars, function(m) {
  make_loo_panel(
    leaveoutdev, leaveoutdevcor, m,
    xlab = NULL, ylab = 'Leave-one-Site-out\nraw deviation', title = m
  )
})
scrpanels <- purrr::map(scrvars, function(m) {
  make_loo_panel(
    leaveoutscr, leaveoutscrcor, m,
    xlab = 'Full-sample value', ylab = 'Leave-one-Site-out\nscore', title = NULL
  )
})

devpanels[[2]] <- devpanels[[2]] + labs(y = NULL)
devpanels[[3]] <- devpanels[[3]] + labs(y = NULL)
scrpanels[[2]] <- scrpanels[[2]] + labs(y = NULL)
scrpanels[[3]] <- scrpanels[[3]] + labs(y = NULL)

# both rows share the same Year color scale (unlike figs.r's combined
# figure, where each row had its own distinct legend), so all 6 panels are
# combined into one wrap_plots() call with a single collected legend,
# rather than collecting legends separately per row
p <- wrap_plots(c(devpanels, scrpanels), ncol = 3, guides = 'collect') &
  theme(legend.position = 'bottom') &
  guides(color = guide_legend(nrow = 1))

png(here('figs/leaveoneout.png'), height = 6, width = 8, units = 'in', res = 300)
print(p)
dev.off()

# floor sensitivity to calibration years ---------------------------------

scr_metric <- c(Abundance = 'abs', `Blade Length` = 'abs', `Short Shoot Density` = 'pct')
yrs <- sort(unique(trndat$yr))

# one floor per metric x year, for each candidate calibration window length
calsens <- purrr::map(c(2, 3, 4, 5, 6), function(n_cal_yrs) {

  cal <- calibrate_scr_fun(trndat, metric = scr_metric, n_cal_yrs = n_cal_yrs)

  cal_lkp <- tibble::tibble(
    metric  = c('Abundance', 'Blade Length', 'Short Shoot Density'),
    mean_sd = c(cal$mean_sd$Abundance, cal$mean_sd$`Blade Length`,
                cal$mean_sd$`Short Shoot Density`)
  )

  purrr::map(yrs, function(yr) {
    tv <- truvar_fun(trndat, yr)
    allgrpscr_fun(trndat, yr, tv, raw_diff = TRUE, metric = scr_metric) |>
      dplyr::summarize(dplyr::across(Abundance:`Short Shoot Density`,
                                      ~ sd(.x, na.rm = TRUE))) |>
      dplyr::mutate(yr = yr)
  }) |>
    dplyr::bind_rows() |>
    tidyr::pivot_longer(-yr, names_to = 'metric', values_to = 'sd_avediff') |>
    dplyr::left_join(cal_lkp, by = 'metric') |>
    dplyr::mutate(
      ratio     = sd_avediff / mean_sd,
      floor_adj = as.integer(pmax(50, 50 + (1 - ratio) * 50)),
      metric    = factor(metric, levels = c('Abundance', 'Blade Length', 'Short Shoot Density')),
      n_cal_yrs = n_cal_yrs
    )
  }) |>
  dplyr::bind_rows()

# does calibration window length change the floors? -------------------------
# repeated-measures design: each yr x metric combo has a floor_adj under
# every candidate n_cal_yrs, so window length is a within-unit factor,
# with (1 | yr) and (1 | metric) absorbing year- and metric-level baseline
# differences. floor_adj is on the same 50-100 scale for all three metrics
# (see calibrate_scr_fun()), so no need for separate per-metric models here.
calsensmod <- lmer(floor_adj ~ factor(n_cal_yrs) + (1 | yr) + (1 | metric), data = calsens)

# omnibus F-test for the n_cal_yrs factor as a whole (not just pairwise
# contrasts against the reference level) -- non-significant means floors
# are not sensitive to calibration window length
anova(calsensmod)

calsensrng <- calsens |>
  summarise(
    floor_mean = mean(floor_adj),
    .by = c(yr, metric)
  )

senscols <- c(
  setNames(viridisLite::viridis(5), sort(unique(calsens$n_cal_yrs))),
  Mean = 'grey20'
)

p <- ggplot(calsens, aes(x = yr, y = floor_adj)) +
  geom_line(
    data = calsensrng,
    aes(x = yr, y = floor_mean, color = 'Mean'),
    inherit.aes = FALSE, linewidth = 0.9
  ) +
  geom_point(
    aes(color = factor(n_cal_yrs)), size = 2,
    position = position_jitter(width = 0.15, height = 0, seed = 123)
  ) +
  scale_color_manual(
    values = senscols, name = 'Calibration window (yrs)',
    breaks = c(sort(unique(calsens$n_cal_yrs)), 'Mean')
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_y_continuous(limits = c(50, 100)) +
  facet_wrap(~ metric, ncol = 1) +
  labs(
    x = NULL, y = 'Score floor',
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = 'bottom', 
    panel.grid.minor = element_blank()
  )

png(here('figs/calsens.png'), width = 6, height = 6, units = 'in', res = 300)
print(p)
dev.off()

# mixeff model for scores over tim ---------------------------------------

load(file = here('data/allyrscrs.RData'))
load(file = here('data/scrmods.RData'))

# two-entry legend: individual group trends vs. the overall group mean
linecols <- c('Groups' = 'grey', 'Group mean' = '#0b0b0b')

# plot the population-level trend (thick) against each group's fitted
# trend (thin) for one score category, using the pre-fit random-intercept model
fit_and_plot <- function(varnm, alldat, scrmods) {

  vardat <- alldat |> filter(var == varnm)
  mod <- scrmods[[varnm]]

  cfs <- summary(mod)$coefficients
  est <- cfs['yrctr', 'Estimate']
  pval <- cfs['yrctr', 'Pr(>|t|)']
  civ <- confint(mod, parm = 'yrctr', quiet = TRUE)

  fmt1 <- function(x) formatC(x, format = 'f', digits = 1)
  estlab <- paste0(fmt1(est), ' (', fmt1(civ[1, 1]), ', ', fmt1(civ[1, 2]), ')')
  plab <- if (pval >= 0.05) 'ns' else if (pval < 0.001) 'p < 0.001' else paste('p =', formatC(pval, format = 'f', digits = 3))

  prdgrd <- expand.grid(
    yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50),
    grp = sort(unique(vardat$grp))
  )
  prdgrd$yr <- prdgrd$yrctr + min(vardat$yr)
  prdgrd$fit <- predict(mod, newdata = prdgrd, re.form = ~(1 | grp))

  ovrfit <- data.frame(yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50))
  ovrfit$yr <- ovrfit$yrctr + min(vardat$yr)
  ovrfit$fit <- predict(mod, newdata = ovrfit, re.form = NA)

  plot <- ggplot() +
    geom_line(data = prdgrd, aes(x = yr, y = fit, color = 'Groups', group = grp), linewidth = 0.6, linetype = 'dashed') +
    geom_line(data = ovrfit, aes(x = yr, y = fit, color = 'Group mean'), linewidth = 2) +
    scale_color_manual(values = linecols, name = NULL) +
    labs(
      x = NULL, y = 'Score', title = varnm,
      subtitle = paste0('Chg yr⁻¹: ', estlab, ', ', plab)
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = 'bottom')

  yrng <- range(c(prdgrd$fit, ovrfit$fit))

  list(plot = plot, yrng = yrng)
}

# FDEP and HC-ES dropped: they have the fewest years of data of the eight groups
alldat <- allyrscrs |>
  filter(!grp %in% c('FDEP', 'HC-ES')) |>
  mutate(yrctr = yr - min(yr))

# combined score + deviation trend figure, 2 x 3 ---------------------------
# top row: calibrated score ~ yr by metric (same models/style as above,
# Total dropped to keep a clean 3-column row to match the bottom row)
# bottom row: absolute weighted-mean deviation ~ yr by metric, from the
# species-level models (abs(devval) ~ yrctr + Species + (1 | grp)). Species
# is a fixed effect there (only 3 levels), not part of the random-intercept
# structure the top row uses, so each group's fitted trend is really one
# line per group x species combination. Color is used for Species (the new
# aesthetic the top row doesn't need); groups are shown the same way as the
# top row's thin grey lines, just one per species color instead of grey,
# so the "many thin lines vs. one thick population mean" visual language
# carries over between rows

load(file = here('data/rawdiffscrs.RData'))
load(file = here('data/spmods.RData'))

spvars <- names(spmods)

spdat <- rawdiffscrs |>
  mutate(
    absdev  = abs(devval),
    yrctr   = yr - min(yr),
    Species = factor(Species)
  ) |>
  filter(!is.na(absdev))

devcols <- c(
  Halodule = '#1b9e77', Syringodium = '#d95f02', Thalassia = '#7570b3',
  'Group/species or overall group mean' = '#0b0b0b'
)

# same per-group (thin) logic as fit_and_plot above, but the per-group
# fitted trend now varies by Species too (color). Three reference lines
# per panel: the thin dashed lines are per-group x species conditional
# (BLUP) predictions (re.form = NULL, both random intercepts included);
# the solid colored lines are per-species means, averaged over grp only
# (re.form = ~(1 | Species)), guaranteed to sit centered among that
# species' own thin lines since they only average within that species;
# the thick black "Group mean" line is the single fixed-effect-only
# prediction (re.form = NA, both random intercepts excluded), a "typical
# group, typical species" value, matching fit_and_plot's "Group mean"
# line above. That black line need not fall anywhere near any one
# species' own cluster, e.g. for Blade Length it happens to sit right at
# Halodule's single highest group (a coincidence of that group's and that
# species' estimated effects nearly canceling, not a bug). Short Shoot
# Density's pooled model estimates exactly zero variance for both grp and
# Species (a singular fit, see spmods summary), so all of its lines
# legitimately collapse onto one flat line, that's a real finding (no
# detectable group or species difference for that metric), not a
# plotting artifact
dev_fit_and_plot <- function(varnm, spdat, spmods) {

  vardat <- spdat |> filter(var == varnm)
  mod <- spmods[[varnm]]

  cfs <- summary(mod)$coefficients
  est <- cfs['yrctr', 'Estimate']
  pval <- cfs['yrctr', 'Pr(>|t|)']
  civ <- confint(mod, parm = 'yrctr', quiet = TRUE)

  fmt1 <- function(x) formatC(x, format = 'f', digits = 3)
  estlab <- paste0(fmt1(est), ' (', fmt1(civ[1, 1]), ', ', fmt1(civ[1, 2]), ')')
  plab <- if (pval >= 0.05) 'ns' else if (pval < 0.001) 'p < 0.001' else paste('p =', formatC(pval, format = 'f', digits = 3))

  prdgrd <- expand.grid(
    yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50),
    grp = sort(unique(vardat$grp)),
    Species = levels(vardat$Species)
  )
  prdgrd$yr <- prdgrd$yrctr + min(vardat$yr)
  prdgrd$fit <- predict(mod, newdata = prdgrd, re.form = NULL)

  sppmean <- expand.grid(
    yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50),
    Species = levels(vardat$Species)
  )
  sppmean$yr <- sppmean$yrctr + min(vardat$yr)
  sppmean$fit <- predict(mod, newdata = sppmean, re.form = ~(1 | Species))

  ovrfit <- data.frame(yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50))
  ovrfit$yr <- ovrfit$yrctr + min(vardat$yr)
  ovrfit$fit <- predict(mod, newdata = ovrfit, re.form = NA)

  ggplot() +
    geom_line(
      data = prdgrd,
      aes(x = yr, y = fit, color = Species, group = interaction(grp, Species)),
      linewidth = 0.6, alpha = 0.7, linetype = 'dashed'
    ) +
    geom_line(data = sppmean, aes(x = yr, y = fit, color = Species, group = Species), linewidth = 1.3) +
    geom_line(data = ovrfit, aes(x = yr, y = fit, color = 'Group/species or overall group mean'), linewidth = 2) +
    scale_color_manual(values = devcols, name = NULL) +
    labs(
      x = NULL, y = 'Absolute deviation', title = NULL,
      subtitle = paste0('Chg yr⁻¹: ', estlab, ', ', plab)
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
}

scrres <- lapply(setdiff(names(scrmods), 'Total'), fit_and_plot, alldat = alldat, scrmods = scrmods)

# shared y-axis within the score row only, same logic as the original
# scrmods.png panel, the deviation row keeps each metric's own scale
# since the three metrics aren't on comparable units (category positions,
# cm, percent difference)
scryrng <- range(unlist(lapply(scrres, `[[`, 'yrng')))
scrplts <- lapply(scrres, function(x) x$plot + coord_cartesian(ylim = scryrng))

devplts <- lapply(spvars, dev_fit_and_plot, spdat = spdat, spmods = spmods)
# each metric keeps its own y-axis scale (different units), but the axis
# title is redundant once shown on the left-most panel
devplts[[2]] <- devplts[[2]] + labs(y = NULL)
devplts[[3]] <- devplts[[3]] + labs(y = NULL)

scrrow <- wrap_plots(scrplts, ncol = 3, guides = 'collect', axes = 'collect', axis_titles = 'collect') &
  theme(legend.position = 'bottom')
devrow <- wrap_plots(devplts, ncol = 3, guides = 'collect') &
  theme(legend.position = 'bottom')

p2 <- scrrow / devrow

png(here('figs/scrdevmods.png'), width = 11, height = 6.5, units = 'in', res = 300)
print(p2)
dev.off()