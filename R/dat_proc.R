library(tbeptools)
library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(glmmTMB)

source(here('R/funcs.R'))

# all training data from report card repo --------------------------------

dataurl <- 'https://github.com/tbep-tech/seagrasstransect-training-reports/raw/refs/heads/main/data/trndat.rda'
fl <- paste(tempdir(), basename(dataurl), sep = "/")
utils::download.file(dataurl, destfile = fl, quiet = TRUE)
load(file = fl)
save(trndat, file = here('data/trndat.rda'), compress = 'xz')

# the report repo's own scoring functions (truvar_fun(), allgrpscr_fun(),
# calibrate_scr_fun(), allyrscr_fun(), etc.) -- sourced early since
# everything below needs them
funcsurl <- 'https://raw.githubusercontent.com/tbep-tech/seagrasstransect-training-reports/main/R/funcs.R'
funcsfl <- paste(tempdir(), 'funcs.R', sep = '/')
download.file(funcsurl, destfile = funcsfl, quiet = TRUE)
source(funcsfl)

yrs <- sort(unique(trndat$yr))
scr_metric <- c(Abundance = 'abs', `Blade Length` = 'abs', `Short Shoot Density` = 'pct')

# fixed historical baseline, built once from the full dataset (same as
# production use)
cal <- calibrate_scr_fun(trndat, metric = scr_metric)

# all score data, computed locally rather than downloaded -----------------
# the report repo's own allyrscr_fun() is reused here directly, but with
# usemon = FALSE. Its default (usemon = TRUE, used when the report repo
# builds its own app/data/allyrscrs.RData) filters groups against
# tbeptools::trnlns$MonAgency, a code list for the *separate*, ongoing
# transect monitoring program (not this training event) -- it's missing
# "TBW" (this training's actual code for Tampa Bay Watch, every year
# 2020-2026) and has no entry for one-off participants like "USF" (2020)
# or "SWK" (2026), so usemon = TRUE silently drops those groups' data from
# every trend model, slope table, and figure built on top of it. The only
# other post-hoc filter allyrscr_fun() applies (na.rm = TRUE, its default)
# drops the literal grp == 'NA' placeholder for rows with no recorded
# affiliation, which is kept here since that's not a real participating
# group

allyrscrs <- allyrscr_fun(trndat, usemon = FALSE, cal = cal, metric = scr_metric) |>
  mutate(grp = dplyr::recode(grp, TBW = 'TBWI'))
save(allyrscrs, file = here('data/allyrscrs.Rdata'), compress = 'xz')

# groups qualifying for trend analyses -------------------------------------
# every figure/table above and below uses every group that participated,
# but an analysis that models CHANGE OVER TIME is only meaningful for a
# group attending enough trainings to show one -- restricted to groups
# with more than three years of participation, aligning with the
# hypothesis that groups attending more than a few years of training
# might show measurable improvement. this restriction applies only to the
# trend models and whatever is built from them: scrmods, spmods, slopetab
# (R/tabs.R), and the scrdevmods figure (R/figs.R) -- not to allyrscrs/
# rawdiffscrs themselves, which stay complete, or to anything else
# (leave-one-out, floor calibration, the training summary table)
trendgrps <- allyrscrs |>
  distinct(grp, yr) |>
  count(grp, name = 'n_yrs') |>
  filter(n_yrs > 3) |>
  pull(grp)

save(trendgrps, file = here('data/trendgrps.RData'), compress = 'xz')

# random-intercept models of score trend by metric ------------------------

load(file = here('data/allyrscrs.RData'))

alldat <- allyrscrs |>
  filter(grp %in% trendgrps) |>
  mutate(yrctr = yr - min(yr))

scrvars <- c('Abundance', 'Blade Length', 'Short Shoot Density', 'Total')
scrmods <- lapply(scrvars, function(varnm){
  vardat <- alldat |> filter(var == varnm)
  lmer(scr ~ yrctr + (1 | grp), data = vardat)
})
names(scrmods) <- scrvars

save(scrmods, file = here('data/scrmods.RData'), compress = 'xz')

# deviation data by species, metric, group, year -------------------------
# kept complete (every group, not just trendgrps) -- used as-is by
# R/tabs.R's training summary table, and restricted to trendgrps only
# where a consumer specifically needs it (slopetab, the scrdevmods figure)

rawdiffscrs <- tibble(yr = yrs) |>
  group_nest(yr) |>
  mutate(
    data = purrr::map(yr, function(x){
      truvar <- truvar_fun(trndat, x)
      allgrpscr_fun(trndat, x, truvar, spp_diff = T)
    })
  ) |>
  unnest('data') |>
  select(-aveval, -truval, -avediff, -aveperc) |>
  mutate(grp = extract_grp(grpact)) |>
  filter(!grp %in% 'NA') |>
  filter(Species %in% c('Halodule', 'Syringodium', 'Thalassia'))

save(rawdiffscrs, file = here('data/rawdiffscrs.RData'), compress = 'xz')

# species, metric mods ---------------------------------------------------
# rawdiffscrs holds the Step 6 per-species weighted-mean deviations
# (devval), so species differences are collapsed first by averaging their
# absolute values within each group x year x metric -- the same Step 7
# collapse allgrpscr_fun(raw_diff = TRUE) performs, and the same one used
# for the scoredevtime figure's "Raw deviation" panels above -- leaving
# one row per group x year x metric.
#
# unlike scrmods (a bounded 0-100 score, reasonably close to normal), the
# absolute deviations are non-negative, right-skewed, and their spread
# scales with their mean -- a Gaussian LMM on the raw scale is a poor fit
# (heavy heteroscedasticity, a singular group-variance fit for two of the
# three metrics) and can predict physically impossible negative
# deviations. A Gamma GLMM (log link), fit via glmmTMB since that's the
# engine with a well-supported Gamma family and a straightforward
# dispersion-model extension if needed later, fixes both: variance scales
# with the mean^2 as expected for this kind of data, and all three metrics
# converge to a non-singular group-variance estimate instead of two of the
# three collapsing to zero under the Gaussian LMM. Gamma requires a
# strictly positive response, so the one exact-zero group/year/metric mean
# (Abundance) gets a small continuity correction (half the smallest
# nonzero value for that metric, a common convention)

spvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')

spdat <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  mutate(absdev = abs(devval)) |>
  summarise(absdev = mean(absdev, na.rm = TRUE), .by = c(var, grp, yr)) |>
  filter(!is.nan(absdev)) |>
  mutate(yrctr = yr - min(yr))

spmods <- lapply(spvars, function(vr){
  vardat <- spdat |> filter(var == vr)
  eps <- min(vardat$absdev[vardat$absdev > 0]) / 2
  vardat <- vardat |> mutate(absdev = absdev + eps * (absdev == 0))
  glmmTMB(absdev ~ yrctr + (1 | grp), family = Gamma(link = 'log'), data = vardat)
})
names(spmods) <- spvars

save(spmods, file = here('data/spmods.RData'))

# per-group x species slopes of raw deviation over time ---------------------
# descriptive (not a model p-value) slope of abs(devval) ~ yr for each
# group x species x metric combination, restricted to trendgrps like
# spmods above -- used by R/tabs.R's slopetab and for inline reporting in
# the manuscript. group identity is also anonymized here (grp_anon: A, B,
# C...), the same convention slopetab displayed on its own before -- the
# original grp is kept alongside it for anything that needs the real
# identity (e.g. inline reporting in the manuscript text)

slopedat <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  mutate(absdev = abs(devval)) |>
  filter(!is.na(absdev)) |>
  summarise(
    slope = coef(lm(absdev ~ yr))[['yr']],
    .by = c(var, grp, Species)
  )

slopegrps <- sort(unique(slopedat$grp))
grplabs <- setNames(LETTERS[seq_along(slopegrps)], slopegrps)
message('slopedat group key: ', paste(names(grplabs), grplabs, sep = ' = ', collapse = ', '))

slopedat <- slopedat |>
  mutate(grp_anon = unname(grplabs[grp]))

save(slopedat, file = here('data/slopedat.RData'), compress = 'xz')

# variance (Kendall trend) tests -------------------------------------------
# a different question than the mean-trend models above: is the SPREAD
# (across-group SD within a year) of the metric scores and raw absolute
# deviations shrinking over time? computed directly from alldat/spdat
# (the same trendgrps-restricted data scrmods/spmods are fit to, not the
# per-species or per-group values those models operate on), one SD per
# metric x year, then a simple non-parametric Kendall trend test per
# metric x type across the 7 yearly SD values. one-tailed (alternative =
# 'less', i.e. tau < 0) since the hypothesis is specifically that
# variance is DECREASING over time, not just changing in either direction

yearly_sd <- dplyr::bind_rows(
  alldat |>
    dplyr::filter(var %in% spvars) |>
    dplyr::summarise(sd = sd(scr, na.rm = TRUE), .by = c(var, yr)) |>
    dplyr::mutate(type = 'Metric scores'),
  spdat |>
    dplyr::summarise(sd = sd(absdev, na.rm = TRUE), .by = c(var, yr)) |>
    dplyr::mutate(type = 'Raw metric deviations')
)

kendalltests <- yearly_sd |>
  dplyr::mutate(
    var  = factor(var, levels = spvars),
    type = factor(type, levels = c('Metric scores', 'Raw metric deviations'))
  ) |>
  dplyr::group_nest(type, var) |>
  dplyr::mutate(
    kt   = purrr::map(data, ~ cor.test(.x$yr, .x$sd, method = 'kendall', alternative = 'less')),
    tau  = purrr::map_dbl(kt, ~ unname(.x$estimate)),
    pval = purrr::map_dbl(kt, ~ .x$p.value)
  ) |>
  dplyr::arrange(type, var) |>
  dplyr::select(type, metric = var, tau, pval)

save(kendalltests, file = here('data/kendalltests.RData'), compress = 'xz')

# leave-one-Site-out / leave-one-group-out sensitivity ---------------------
# does any single Site or group disproportionately drive a group's score?
# for each year, drops one Site (or one group) at a time, rescores from
# what's left, and compares to the full-sample value. see run_loo()/
# run_loo_grp() in R/funcs.R for the mechanics and rationale

# a dropped Site/group can occasionally remove all data for one of the
# metrics entirely, which errors inside allgrpscr_fun -- skip that
# iteration rather than fail the whole run
devfun_safe <- purrr::possibly(function(sub, yr) devfun(sub, yr, scr_metric), otherwise = NULL)
scrfun_safe <- purrr::possibly(function(sub, yr) scrfun(sub, yr, cal, scr_metric), otherwise = NULL)

leaveoutdev <- run_loo(trndat, yrs, devfun_safe, spvars)
leaveoutdevcor <- summarize_loo(leaveoutdev)
leaveoutdevtest <- diff_test_loo(leaveoutdev)
print(leaveoutdevtest)

leaveoutscr <- run_loo(trndat, yrs, scrfun_safe, spvars)
leaveoutscrcor <- summarize_loo(leaveoutscr)
leaveoutscrtest <- diff_test_loo(leaveoutscr)
print(leaveoutscrtest)

leaveoutdev_grp <- run_loo_grp(trndat, yrs, devfun_safe, spvars)
leaveoutdevcor_grp <- summarize_loo(leaveoutdev_grp)
leaveoutdevtest_grp <- diff_test_loo(leaveoutdev_grp)
print(leaveoutdevtest_grp)

leaveoutscr_grp <- run_loo_grp(trndat, yrs, scrfun_safe, spvars)
leaveoutscrcor_grp <- summarize_loo(leaveoutscr_grp)
leaveoutscrtest_grp <- diff_test_loo(leaveoutscr_grp)
print(leaveoutscrtest_grp)

save(
  leaveoutdev, leaveoutdevcor, leaveoutdevtest,
  leaveoutscr, leaveoutscrcor, leaveoutscrtest,
  leaveoutdev_grp, leaveoutdevcor_grp, leaveoutdevtest_grp,
  leaveoutscr_grp, leaveoutscrcor_grp, leaveoutscrtest_grp,
  file = here('data/leaveoneout.RData'), compress = 'xz'
)

# full-sample vs. leave-one-out Pearson correlations, for inline reporting
# in the manuscript text (one row per test x metric, the same r values
# annotated on each panel of the leaveoneout figure in R/figs.R)
leaveoneoutcors <- dplyr::bind_rows(
  leaveoutscrcor |> dplyr::mutate(test = 'site_scr'),
  leaveoutdevcor |> dplyr::mutate(test = 'site_dev'),
  leaveoutscrcor_grp |> dplyr::mutate(test = 'grp_scr'),
  leaveoutdevcor_grp |> dplyr::mutate(test = 'grp_dev')
) |>
  dplyr::rename(corr = r) |>
  tidyr::separate(test, c('test', 'calc')) |> 
  dplyr::select(test, calc, metric, corr)

save(leaveoneoutcors, file = here('data/leaveoneoutcors.RData'), compress = 'xz')

# floor sensitivity to calibration years ------------------------------------

# one floor per metric x year, for each candidate calibration window length
calsens <- purrr::map(c(2, 3, 4, 5, 6), function(n_cal_yrs) {

  calsub <- calibrate_scr_fun(trndat, metric = scr_metric, n_cal_yrs = n_cal_yrs)

  cal_lkp <- tibble::tibble(
    metric  = spvars,
    mean_sd = c(calsub$mean_sd$Abundance, calsub$mean_sd$`Blade Length`,
                calsub$mean_sd$`Short Shoot Density`)
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
      metric    = factor(metric, levels = spvars),
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
print(anova(calsensmod))

save(calsens, calsensmod, file = here('data/calsens.RData'), compress = 'xz')
