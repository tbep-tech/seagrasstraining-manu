library(tidyverse)
library(here)

load(here("data", "trndat.rda"))
load(here("data", "allyrscrs.RData"))

tmp <- allyrscrs |> 
  filter(!var %in% c('Total')) |> 
  summarise(
    std = sd(scr, na.rm = T), 
    .by = c(yr, var)
  )

ggplot(tmp, aes(x = yr, y = std)) +
  geom_col() +
  scale_x_continuous(breaks = unique(tmp$yr)) +
  facet_wrap(~var, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Standard Deviation"
  )

tmp2 <- allyrscrs |> 
  filter(!var %in% c('Total'))

ggplot(tmp2, aes(x = yr, y = scr, group = yr)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_x_continuous(breaks = unique(tmp2$yr)) +
  facet_wrap(~var, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Score"
  )

# raw (pre-rescale) differences instead of the calibrated 0-100 scr --------
# scr is a within-year min-max rescale (worst group that year always lands
# at the floor, best always at 100), so it reflects a group's rank among
# its peers that year, not its absolute accuracy, an absolute "have groups
# gotten more accurate" trend can be invisible to scr even when it's real.
# The sections below use avediff/aveperc from allgrpscr_fun(raw_diff = TRUE)
# instead: the weighted-mean deviation from truth used right before that
# rescale happens, still on each metric's own scale (cm for Blade Length, a
# proportion for Short Shoot Density, an ordinal-ish measure for Abundance)

library(lme4)
library(lmerTest)
library(patchwork)

# how data/rawdiffscrs.RData was built --------------------------------------
# kept here for reproducibility, the models below just load the already-
# built file. Pulls R/funcs.R from the report repo the same way trndat.rda/
# allyrscrs.RData are pulled in dat_proc.R, rather than assuming a local
# checkout of that repo

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
      allgrpscr_fun(trndat, x, truvar, raw_diff = TRUE)
    })
  ) |>
  unnest('data') |>
  pivot_longer(cols = c(Abundance, `Blade Length`, `Short Shoot Density`),
               names_to = 'var', values_to = 'avediff') |>
  mutate(
    grp = gsub('^.*:\\s(.*?)\\s\\(.*$', '\\1', grpact)
  ) |>
  summarise(
    avediff = mean(avediff, na.rm = TRUE),
    .by = c(yr, grp, var)
  ) |>
  filter(!grp %in% 'NA') |>
  filter(grp %in% tbeptools::trnlns$MonAgency)

save(rawdiffscrs, file = here('data/rawdiffscrs.RData'), compress = 'xz')

# mixed model of raw (pre-rescale) differences over time -------------------
# avediff is strictly positive and right-skewed (a distance from truth), so
# a Gamma GLMM with a log link is used instead of a Gaussian lmer on the raw
# scale, this also means the yrctr coefficient is a proportional (not
# additive) change per year, and a *negative* coefficient is improvement
# (smaller deviation from truth), opposite of the scr-based models

load(file = here('data/rawdiffscrs.RData'))

# FDEP and HC-ES dropped: they have the fewest years of data of the eight groups
alldat <- rawdiffscrs |>
  filter(!grp %in% c('FDEP', 'HC-ES')) |>
  mutate(yrctr = yr - min(yr))

scrvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')
scrmods <- lapply(scrvars, function(varnm){
  vardat <- alldat |> filter(var == varnm)
  glmer(avediff ~ yrctr + (1 | grp), family = Gamma(link = 'log'), data = vardat)
})
names(scrmods) <- scrvars

# population-level trend (thick) against each group's fitted trend (thin),
# on the response (raw difference) scale
fit_and_plot <- function(varnm, alldat, scrmods) {

  vardat <- alldat |> filter(var == varnm)
  mod <- scrmods[[varnm]]

  cfs <- summary(mod)$coefficients
  est <- cfs['yrctr', 'Estimate']
  pval <- cfs['yrctr', 'Pr(>|z|)']
  civ <- suppressMessages(confint(mod, parm = 'yrctr', method = 'Wald'))

  # exponentiate: proportional change in raw difference per year
  fmt1 <- function(x) formatC(100 * (exp(x) - 1), format = 'f', digits = 1)
  estlab <- paste0(fmt1(est), '% (', fmt1(civ[1, 1]), ', ', fmt1(civ[1, 2]), ')')
  plab <- if (pval >= 0.05) 'ns' else if (pval < 0.001) 'p < 0.001' else paste('p =', formatC(pval, format = 'f', digits = 3))

  prdgrd <- expand.grid(
    yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50),
    grp = sort(unique(vardat$grp))
  )
  prdgrd$yr <- prdgrd$yrctr + min(vardat$yr)
  prdgrd$fit <- predict(mod, newdata = prdgrd, re.form = ~(1 | grp), type = 'response')

  ovrfit <- data.frame(yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50))
  ovrfit$yr <- ovrfit$yrctr + min(vardat$yr)
  ovrfit$fit <- predict(mod, newdata = ovrfit, re.form = NA, type = 'response')

  linecols <- c('Groups' = 'grey', 'Group mean' = '#0b0b0b')

  plot <- ggplot() +
    geom_point(data = vardat, aes(x = yr, y = avediff), color = 'grey60', alpha = 0.6) +
    geom_line(data = prdgrd, aes(x = yr, y = fit, color = 'Groups', group = grp), linewidth = 0.6) +
    geom_line(data = ovrfit, aes(x = yr, y = fit, color = 'Group mean'), linewidth = 2) +
    scale_color_manual(values = linecols, name = NULL) +
    labs(
      x = NULL, y = 'Raw difference from truth', title = varnm,
      subtitle = paste0('Chg yr⁻¹: ', estlab, ', ', plab)
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = 'bottom')

  yrng <- range(c(vardat$avediff, prdgrd$fit, ovrfit$fit))

  list(plot = plot, yrng = yrng)
}

res <- lapply(names(scrmods), fit_and_plot, alldat = alldat, scrmods = scrmods)

# each metric keeps its own y-axis scale here (different units/basis), unlike
# the scr-based version where a common 0-100 scale makes a shared axis sensible
plts <- lapply(res, `[[`, 'plot')

p <- wrap_plots(plts, ncol = 2, guides = 'collect') &
  theme(legend.position = 'bottom')

print(p)

# within-year SD of raw differences over time -------------------------------
# a second, distinct kind of improvement: the raw-difference models above ask
# whether groups are getting closer to the truth, this asks whether groups
# are getting closer to EACH OTHER (this is the same quantity
# calibrate_scr_fun uses for the floor, but there only over the fixed 5-year
# baseline window, here tracked across all years as its own outcome). all
# eight groups (not just the six used above) are kept here since this is a
# per-year summary stat, not a per-group trend, so the FDEP/HC-ES short
# history isn't a problem, one metric-year with too few groups is
sddat <- rawdiffscrs |>
  summarise(
    sdval = sd(avediff, na.rm = TRUE),
    n = dplyr::n(),
    .by = c(yr, var)
  ) |>
  filter(n > 1) |>
  mutate(yrctr = yr - min(yr))

# one GLM per metric instead of a single pooled model: there's no group
# structure left at this point (one sdval per year per metric), so a mixed
# model doesn't add anything, and fitting each metric on its own scale
# sidesteps the cm vs proportion vs abundance unit mismatch entirely instead
# of needing to standardize first
sdmods <- lapply(scrvars, function(varnm){
  vardat <- sddat |> filter(var == varnm)
  glm(sdval ~ yrctr, family = Gamma(link = 'log'), data = vardat)
})
names(sdmods) <- scrvars

sd_fit_and_plot <- function(varnm, sddat, sdmods) {

  vardat <- sddat |> filter(var == varnm)
  mod <- sdmods[[varnm]]

  cfs <- summary(mod)$coefficients
  est <- cfs['yrctr', 'Estimate']
  pval <- cfs['yrctr', 'Pr(>|t|)']
  civ <- confint.default(mod, parm = 'yrctr')

  # exponentiate: proportional change in SD per year
  fmt1 <- function(x) formatC(100 * (exp(x) - 1), format = 'f', digits = 1)
  estlab <- paste0(fmt1(est), '% (', fmt1(civ[1, 1]), ', ', fmt1(civ[1, 2]), ')')
  plab <- if (pval >= 0.05) 'ns' else if (pval < 0.001) 'p < 0.001' else paste('p =', formatC(pval, format = 'f', digits = 3))

  prdgrd <- data.frame(yrctr = seq(min(vardat$yrctr), max(vardat$yrctr), length.out = 50))
  prdgrd$yr <- prdgrd$yrctr + min(vardat$yr)
  prlink <- predict(mod, newdata = prdgrd, type = 'link', se.fit = TRUE)
  prdgrd$fit <- exp(prlink$fit)
  prdgrd$lwr <- exp(prlink$fit - 1.96 * prlink$se.fit)
  prdgrd$upr <- exp(prlink$fit + 1.96 * prlink$se.fit)

  ggplot() +
    geom_ribbon(data = prdgrd, aes(x = yr, ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_point(data = vardat, aes(x = yr, y = sdval), color = 'grey30', size = 2) +
    geom_line(data = prdgrd, aes(x = yr, y = fit), linewidth = 1) +
    labs(
      x = NULL, y = 'Within-year SD of raw differences\n(agreement across groups)', title = varnm,
      subtitle = paste0('Chg yr⁻¹: ', estlab, ', ', plab)
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

sdplts <- lapply(scrvars, sd_fit_and_plot, sddat = sddat, sdmods = sdmods)

sdp <- wrap_plots(sdplts, ncol = 2)

print(sdp)
