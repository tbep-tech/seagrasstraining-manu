library(tidyverse)
library(here)
library(patchwork)
library(lme4)
library(lmerTest)

source(here('R/funcs.R'))

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