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

# all score data from report card repo -----------------------------------

allyrscrs <- tbepreport::util_rdataload('https://github.com/tbep-tech/seagrasstransect-training-reports/raw/refs/heads/main/app/data/allyrscrs.RData')
save(allyrscrs, file = here('data/allyrscrs.Rdata'), compress = 'xz')

# random-intercept models of score trend by metric ------------------------

load(file = here('data/allyrscrs.RData'))

# FDEP and HC-ES dropped: they have the fewest years of data of the eight groups
alldat <- allyrscrs |>
  filter(!grp %in% c('FDEP', 'HC-ES')) |>
  mutate(yrctr = yr - min(yr))

scrvars <- c('Abundance', 'Blade Length', 'Short Shoot Density', 'Total')
scrmods <- lapply(scrvars, function(varnm){
  vardat <- alldat |> filter(var == varnm)
  lmer(scr ~ yrctr + (1 | grp), data = vardat)
})
names(scrmods) <- scrvars

save(scrmods, file = here('data/scrmods.RData'), compress = 'xz')

# deviation data by species, metric, group, year -------------------------

funcsurl <- 'https://raw.githubusercontent.com/tbep-tech/seagrasstransect-training-reports/main/R/funcs.R'
funcsfl <- paste(tempdir(), 'funcs.R', sep = '/')
download.file(funcsurl, destfile = funcsfl, quiet = TRUE)
source(funcsfl)

data(file = 'trnlns', package = 'tbeptools')

yrs <- sort(unique(trndat$yr))

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
  mutate(
    grp = gsub('^.*:\\s(.*?)\\s\\(.*$', '\\1', grpact)
  ) |>
  filter(!grp %in% 'NA') |>
  filter(grp %in% tbeptools::trnlns$MonAgency) |> 
  filter(Species %in% c('Halodule', 'Syringodium', 'Thalassia'))

save(rawdiffscrs, file = here('data/rawdiffscrs.RData'), compress = 'xz')

# species, metric mods ---------------------------------------------------


spvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')

spdat <- rawdiffscrs |>
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

scr_metric <- c(Abundance = 'abs', `Blade Length` = 'abs', `Short Shoot Density` = 'pct')

# fixed historical baseline, built once from the full dataset (same as
# production use) -- not rebuilt per Site/group dropped, since cal is meant
# to stay fixed regardless of which Sites a given focal year's own data has
cal <- calibrate_scr_fun(trndat, metric = scr_metric)

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