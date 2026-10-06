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