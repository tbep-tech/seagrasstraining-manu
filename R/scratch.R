library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(performance)

source(here('R/funcs.R'))

# rebuild the same collapsed-across-species absolute deviation data used
# for spmods in R/dat_proc.R (not saved on its own, so recomputed here
# from the saved intermediates)
load(file = here('data/rawdiffscrs.RData'))
load(file = here('data/trendgrps.RData'))

spvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')

spdat <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  mutate(absdev = abs(devval)) |>
  summarise(absdev = mean(absdev, na.rm = TRUE), .by = c(var, grp, yr)) |>
  filter(!is.nan(absdev)) |>
  mutate(yrctr = yr - min(yr))

# any exact zeros? Gamma/log-normal need strictly positive values
spdat |>
  summarise(n = n(), n_zero = sum(absdev == 0), min = min(absdev), .by = var) |>
  print()

# three candidate model families per metric --------------------------------
# 1) current approach: Gaussian LMM on the raw absolute deviation
# 2) Gaussian LMM on log(absdev) -- trend becomes a %/yr change, variance
#    no longer forced constant on the raw scale
# 3) Gamma GLMM (log link) -- models the mean on the log scale directly,
#    variance scales with the mean^2 (Gamma's natural mean-variance
#    relationship), no need to log-transform the response itself

fits <- map(spvars, function(vr) {

  vardat <- spdat |> filter(var == vr)

  # Abundance has one exact-zero group/year mean (printed above) -- log()
  # and Gamma() both require strictly positive y, so a small continuity
  # correction (half the smallest nonzero value for that metric, a common
  # convention) is added only for those two model types, leaving the
  # Gaussian LMM on the untransformed scale as-is
  eps <- min(vardat$absdev[vardat$absdev > 0]) / 2
  vardat <- vardat |> mutate(absdev_adj = absdev + eps * (absdev == 0))

  gauss <- lmer(absdev ~ yrctr + (1 | grp), data = vardat)
  gauss_log <- lmer(log(absdev_adj) ~ yrctr + (1 | grp), data = vardat)
  gamma_glmm <- glmer(absdev_adj ~ yrctr + (1 | grp), data = vardat, family = Gamma(link = 'log'))

  list(gauss = gauss, gauss_log = gauss_log, gamma_glmm = gamma_glmm)

})
names(fits) <- spvars

# comparison table: trend estimate/p-value (on each model's own scale),
# AIC, and singular-fit flag, one row per metric x model
comptab <- map(spvars, function(vr) {

  map(names(fits[[vr]]), function(modnm) {

    mod <- fits[[vr]][[modnm]]
    cfs <- summary(mod)$coefficients

    tibble(
      metric = vr,
      model = modnm,
      estimate = cfs['yrctr', 'Estimate'],
      pval = cfs['yrctr', grepl('Pr', colnames(cfs))],
      AIC = AIC(mod),
      singular = isSingular(mod)
    )

  }) |>
    bind_rows()

}) |>
  bind_rows()

print(comptab, n = Inf)

# residual diagnostics -- do the Gaussian LMM residuals actually look
# heteroscedastic / skewed relative to the alternatives? saved to a
# scratch png rather than relying on an interactive plot window
png(here('figs/scratch_devmod_diagnostics.png'), width = 10, height = 10, units = 'in', res = 150)
par(mfrow = c(3, 3))
for (vr in spvars) {
  for (modnm in c('gauss', 'gauss_log', 'gamma_glmm')) {
    mod <- fits[[vr]][[modnm]]
    plot(fitted(mod), resid(mod, type = 'pearson'),
         main = paste(vr, modnm), xlab = 'Fitted', ylab = 'Pearson resid')
    abline(h = 0, lty = 2, col = 'red')
  }
}
dev.off()

message('wrote figs/scratch_devmod_diagnostics.png -- fitted vs. Pearson residuals, 3 metrics x 3 model types')

# does the VARIANCE of the raw (signed) deviation decrease over time? ------
# variance is a property of the signed difference, not its absolute value,
# so these two approaches both use the signed, species-collapsed devval
# (same collapse as spdat above, just without abs())

spdat_signed <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  summarise(devval = mean(devval, na.rm = TRUE), .by = c(var, grp, yr)) |>
  filter(!is.nan(devval)) |>
  mutate(yrctr = yr - min(yr))

# 1) simple approach: SD across groups within each year, then a basic
# monotonic trend test (Kendall's tau) on the 7 yearly SD values per metric
yearly_sd <- spdat_signed |>
  summarise(sd = sd(devval, na.rm = TRUE), n = n(), .by = c(var, yr))

print(yearly_sd, n = Inf)

sd_trend <- map(spvars, function(vr) {
  dat <- yearly_sd |> filter(var == vr) |> arrange(yr)
  kt <- cor.test(dat$yr, dat$sd, method = 'kendall')
  tibble(metric = vr, tau = unname(kt$estimate), pval = kt$p.value)
}) |>
  bind_rows()

print(sd_trend)

png(here('figs/scratch_yearly_sd.png'), width = 9, height = 3.5, units = 'in', res = 150)
par(mfrow = c(1, 3))
for (vr in spvars) {
  dat <- yearly_sd |> filter(var == vr) |> arrange(yr)
  plot(dat$yr, dat$sd, type = 'b', main = vr, xlab = 'Year', ylab = 'SD across groups')
}
dev.off()

# 2) glmmTMB dispersion-model approach: same mean structure as the current
# spmods (devval ~ yrctr + (1 | grp)), but with an explicit model for how
# the residual dispersion itself changes with year (dispformula = ~
# yrctr). A significant dispformula yrctr term is a direct test of
# "variance changes with year" -- gaussian() dispersion uses a log link,
# so a negative estimate means variance shrinking over time
library(glmmTMB)

disp_fits <- map(spvars, function(vr) {
  vardat <- spdat_signed |> filter(var == vr)
  glmmTMB(devval ~ yrctr + (1 | grp), dispformula = ~ yrctr, family = gaussian(), data = vardat)
})
names(disp_fits) <- spvars

for (vr in spvars) {
  cat('\n==', vr, '==\n')
  print(summary(disp_fits[[vr]]))
}
