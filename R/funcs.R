# leave-one-unit-out sensitivity --------------------------------------------

# raw (pre-calibration) deviation scoring function for one Site/group subset
devfun <- function(sub, yr, scr_metric) {
  tv <- truvar_fun(sub, yr)
  allgrpscr_fun(sub, yr, tv, raw_diff = TRUE, metric = scr_metric)
}

# calibrated (0-100) scoring function for one Site/group subset, using a
# fixed calibration baseline (cal) built once from the full dataset -- not
# rebuilt per Site/group dropped, since cal is meant to stay fixed
# regardless of which Sites/groups a given focal year's own data has
scrfun <- function(sub, yr, cal, scr_metric) {
  tv <- truvar_fun(sub, yr)
  allgrpscr_fun(sub, yr, tv, cal = cal, raw = TRUE, metric = scr_metric)
}

# runs the leave-one-Site-out loop for a given scoring function, returning
# one row per group x year x metric x Site left out, with the full-sample
# and leave-one-out values side by side
run_loo <- function(trndat, yrs, scorefun_safe, varnames) {

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
      grp    = extract_grp(grpact),
      metric = factor(metric, levels = varnames)
    ) |>
    dplyr::filter(!is.na(full), !is.na(loo)) |>
    dplyr::filter(!grp %in% 'NA')

}

# same as run_loo() above, but drops one Group's data entirely each
# iteration instead of one Site. a dropped group isn't scored itself (it has
# no data left in sub), but truvar_fun()/truespp_fun() inside devfun/scrfun
# rebuild the cross-group consensus from the remaining groups each time, so
# every other group's consensus-dependent score can shift. the left_join +
# !is.na(loo) filter already used downstream naturally drops the removed
# group's own (now-missing) row, so the comparison is exactly "how much does
# removing one group from the consensus move the scores of the groups that
# are left," complementing the Site-level diagnostic above with a
# group-level one (Site variation vs. group/consensus variation)
run_loo_grp <- function(trndat, yrs, scorefun_safe, varnames) {

  purrr::map(yrs, function(yr) {

    message('year ', yr)

    yrdat <- trndat |> dplyr::filter(yr == !!yr)
    grps  <- sort(unique(yrdat$grpact))

    full <- scorefun_safe(yrdat, yr)
    if (is.null(full)) return(NULL)

    full_long <- full |>
      tidyr::pivot_longer(dplyr::all_of(varnames), names_to = 'metric', values_to = 'full')

    purrr::map(grps, function(g) {

      sub <- yrdat |> dplyr::filter(grpact != g)
      loo <- scorefun_safe(sub, yr)
      if (is.null(loo)) return(NULL)

      loo |>
        tidyr::pivot_longer(dplyr::all_of(varnames), names_to = 'metric', values_to = 'loo') |>
        dplyr::left_join(full_long, by = c('grpact', 'metric')) |>
        dplyr::mutate(grp_out = g, yr = yr)

    }) |>
      purrr::compact() |>
      dplyr::bind_rows()

  }) |>
    purrr::compact() |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      grp    = extract_grp(grpact),
      metric = factor(metric, levels = varnames)
    ) |>
    dplyr::filter(!is.na(full), !is.na(loo)) |>
    dplyr::filter(!grp %in% 'NA')

}

# how much does dropping any one Site/group move a group's value? correlation
# between full-sample and leave-one-out values per metric (groups are on
# different scales by metric), plus the median/max absolute difference as a
# plain-units influence measure
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

# formal test for a systematic bias between full-sample and leave-one-out
# values, separately per metric (different scales, same reasoning as
# summarize_loo()). testing whether an ordinary regression's slope differs
# from 1 isn't used here -- OLS assumes the x-axis is measured without
# error, but full and loo are both noisy measurements of the same group x
# year quantity, so OLS is biased toward a shallower slope purely from
# measurement noise (regression dilution), which would flag "significant"
# deviation even for a perfectly unbiased method. instead this tests
# whether the mean paired difference (full - loo) differs from zero, with
# random intercepts for yr and grp accounting for the same group x year
# appearing once per Site/group dropped that year (not independent
# observations)
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

# plotting -------------------------------------------------------------------

# one full-sample vs leave-one-out panel for a single metric (Site- or
# Group-level, raw deviation or score, depending on which loodat/loocor are
# passed in) -- one ggplot object per metric rather than facet_wrap, so
# the dev/score blocks (different datasets) can be stacked with patchwork
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

  # two-entry legend: individual group trends vs. the overall group mean
  linecols <- c('Groups' = 'grey', 'Group mean' = '#0b0b0b')

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

  devcols <- c(
    Halodule = '#1b9e77', Syringodium = '#d95f02', Thalassia = '#7570b3',
    'Group/species or overall group mean' = '#0b0b0b'
  )

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

# one Score or Raw deviation panel for a single metric, year on x, group
# values jittered underneath a mean +/- 95% CI summary -- built as individual
# ggplot objects (not facet_wrap) so the Score/Raw deviation blocks can be
# stacked with patchwork and labeled (a)/(b), same pattern as
# make_loo_panel() above
make_scoredev_panel <- function(summdat, alldat, m, typ, varlab, ylim = NULL) {

  summ <- summdat |> dplyr::filter(type == typ, var == m)
  jit  <- alldat  |> dplyr::filter(type == typ, var == m)

  plot <- ggplot(summ, aes(x = yr, y = mean)) +
    geom_jitter(
      data = jit, aes(x = yr, y = val),
      width = 0.15, height = 0, size = 0.9, alpha = 0.3, color = 'grey10',
      inherit.aes = FALSE
    ) +
    geom_pointrange(aes(ymin = lwr, ymax = upr), color = '#004F7E', linewidth = 0.6, size = 0.5) +
    labs(x = NULL, y = varlab) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  if (!is.null(ylim)) plot <- plot + coord_cartesian(ylim = ylim)

  plot
}

# table-building helpers ------------------------------------------------------

# extracts the short group code from a grpact string (e.g. "2021: TBW (RR
# MG TK PL)" -> "TBW"), then normalizes known synonyms to their canonical
# code -- "TBW" is this training's own abbreviation for Tampa Bay Watch,
# confirmed to be the same group as "TBWI" (the code used elsewhere in
# tbeptools for its separate, ongoing transect monitoring program)
extract_grp <- function(grpact) {
  grp <- gsub('^.*:\\s(.*?)\\s\\(.*$', '\\1', grpact)
  dplyr::recode(grp, TBW = 'TBWI')
}

# applies extract_grp() without dropping anything -- every row that
# participated that year counts, including the literal grp == 'NA'
# placeholder (an entry with no affiliation given, e.g. "2021: NA (UT)"):
# that participant's data still fed into that year's baseline variation,
# species weights, and consensus calculations like everyone else's, and
# this helper is only ever used for per-year descriptive summaries (not
# for tracking a group's identity across years, where an anonymous entry
# can't be meaningfully attributed to the same group from one year to the
# next -- see trendgrps in R/dat_proc.R for that separate exclusion).
# this also does NOT filter against tbeptools::trnlns$MonAgency: that list
# is for the separate, ongoing transect monitoring program, not this
# training event, and has no entry at all for several real training
# groups (e.g. one-off participants like "USF" or "SWK")
validgrp <- function(dat) {
  dat |>
    dplyr::mutate(grp = extract_grp(grpact))
}

# one cell as "mean (min, max)" (or "mean (min)" if show_max = FALSE)
# across a vector, e.g. a metric's values across groups in a given year
fmtavgrng <- function(x, digits, show_max = TRUE) {
  fmt <- paste0('%.', digits, 'f')
  if (show_max) {
    sprintf(paste0(fmt, ' (', fmt, ', ', fmt, ')'), mean(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  } else {
    sprintf(paste0(fmt, ' (', fmt, ')'), mean(x, na.rm = TRUE), min(x, na.rm = TRUE))
  }
}

# a blank full-width row with `label` in `label_col` and NA elsewhere, used
# to insert a section-header row directly into a wide summary table
# (merged across columns afterward with flextable's merge_at()) rather
# than relying on flextable's as_grouped_data()
blank_row <- function(label, cols, label_col = 'Metric') {
  row <- stats::setNames(as.list(rep(NA_character_, length(cols))), cols)
  tibble::as_tibble(c(stats::setNames(list(label), label_col), row))
}
