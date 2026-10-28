library(tbeptools)
library(tidyverse)
library(here)
library(lme4)
library(lmerTest)

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

spvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')

spdat <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  mutate(
    absdev  = abs(devval),
    yrctr   = yr - min(yr),
    Species = factor(Species)
  ) |>
  filter(!is.na(absdev))

spmods <- lapply(spvars, function(vr){
  vardat <- spdat |> filter(var == vr)
  lmer(absdev ~ yrctr + (1 | Species) + (1 | grp), data = vardat)
})
names(spmods) <- spvars

save(spmods, file = here('data/spmods.RData'))

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
