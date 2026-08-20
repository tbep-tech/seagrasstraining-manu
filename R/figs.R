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
    geom_line(data = prdgrd, aes(x = yr, y = fit, color = 'Groups', group = grp), linewidth = 0.6) +
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

res <- lapply(names(scrmods), fit_and_plot, alldat = alldat, scrmods = scrmods)

# common y-axis scale across all four panels
yrng <- range(unlist(lapply(res, `[[`, 'yrng')))
plts <- lapply(res, function(x) x$plot + coord_cartesian(ylim = yrng))

p <- wrap_plots(plts, ncol = 2, guides = 'collect', axes = 'collect', axis_titles = 'collect') &
  theme(legend.position = 'bottom')

png(here('figs/scrmods.png'), width = 6.5, height = 6.5, units = 'in', res = 300)
print(p)
dev.off()