library(tidyverse)
library(here)
library(patchwork)
library(lme4)
library(lmerTest)

source(here('R/funcs.R'))

load(file = here('data/leaveoneout.RData'))
load(file = here('data/calsens.RData'))
load(file = here('data/allyrscrs.RData'))
load(file = here('data/scrmods.RData'))
load(file = here('data/rawdiffscrs.RData'))
load(file = here('data/spmods.RData'))

# combined leave-one-Site-out (a) / leave-one-group-out (b) figure --------
# stacks two 2 x 3 blocks into one 4 x 3 figure (Site block on top, Group
# block on bottom): no x-axis titles in the Site block (metric identity and
# axis meaning are already established there for both blocks), no
# per-metric panel titles in the Group block (redundant with the Site
# block's titles directly above), and short shared y-axis titles ("Raw
# deviation" / "Score"). (a)/(b) labels are added via labs(subtitle = ...)
# on each block's top-left panel. combining all 12 panels into a single
# wrap_plots() call (rather than stacking two separately-collected blocks
# with /) keeps one shared, collected legend

devvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')
scrvars <- devvars

devpanels_comb <- purrr::map(devvars, function(m) {
  make_loo_panel(
    leaveoutdev, leaveoutdevcor, m,
    xlab = NULL, ylab = 'Raw deviation', title = m
  )
})
scrpanels_comb <- purrr::map(scrvars, function(m) {
  make_loo_panel(
    leaveoutscr, leaveoutscrcor, m,
    xlab = NULL, ylab = 'Score', title = NULL
  )
})
devpanels_grp_comb <- purrr::map(devvars, function(m) {
  make_loo_panel(
    leaveoutdev_grp, leaveoutdevcor_grp, m,
    xlab = NULL, ylab = 'Raw deviation', title = NULL
  )
})
scrpanels_grp_comb <- purrr::map(scrvars, function(m) {
  make_loo_panel(
    leaveoutscr_grp, leaveoutscrcor_grp, m,
    xlab = 'Full-sample value', ylab = 'Score', title = NULL
  )
})

devpanels_comb[[2]] <- devpanels_comb[[2]] + labs(y = NULL)
devpanels_comb[[3]] <- devpanels_comb[[3]] + labs(y = NULL)
scrpanels_comb[[2]] <- scrpanels_comb[[2]] + labs(y = NULL)
scrpanels_comb[[3]] <- scrpanels_comb[[3]] + labs(y = NULL)
devpanels_grp_comb[[2]] <- devpanels_grp_comb[[2]] + labs(y = NULL)
devpanels_grp_comb[[3]] <- devpanels_grp_comb[[3]] + labs(y = NULL)
scrpanels_grp_comb[[2]] <- scrpanels_grp_comb[[2]] + labs(y = NULL)
scrpanels_grp_comb[[3]] <- scrpanels_grp_comb[[3]] + labs(y = NULL)

devpanels_comb[[1]] <- devpanels_comb[[1]] + labs(subtitle = '(a) Leave-one-site-out')
devpanels_grp_comb[[1]] <- devpanels_grp_comb[[1]] + labs(subtitle = '(b) Leave-one-group-out')

p_all <- wrap_plots(
  c(devpanels_comb, scrpanels_comb, devpanels_grp_comb, scrpanels_grp_comb),
  ncol = 3, guides = 'collect'
) &
  theme(legend.position = 'bottom') &
  guides(color = guide_legend(nrow = 1))

png(here('figs/leaveoneout.png'), height = 11, width = 8, units = 'in', res = 300)
print(p_all)
dev.off()

# floor sensitivity to calibration years figure -----------------------------

senscols <- c(
  setNames(viridisLite::viridis(5), sort(unique(calsens$n_cal_yrs))),
  Mean = 'grey20'
)

calsensrng <- calsens |>
  summarise(
    floor_mean = mean(floor_adj),
    .by = c(yr, metric)
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

# group score / deviation trend over time figure ----------------------------
# two-entry legend: individual group trends vs. the overall group mean (see
# fit_and_plot()/dev_fit_and_plot() in R/funcs.R for plot mechanics)

# FDEP and HC-ES dropped: they have the fewest years of data of the eight groups
alldat <- allyrscrs |>
  filter(!grp %in% c('FDEP', 'HC-ES')) |>
  mutate(yrctr = yr - min(yr))

spvars <- names(spmods)

spdat <- rawdiffscrs |>
  mutate(
    absdev  = abs(devval),
    yrctr   = yr - min(yr),
    Species = factor(Species)
  ) |>
  filter(!is.na(absdev))

# combined score + deviation trend figure, 2 x 3: top row calibrated score ~
# yr by metric (Total dropped to keep a clean 3-column row to match the
# bottom row), bottom row absolute weighted-mean deviation ~ yr by metric
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
