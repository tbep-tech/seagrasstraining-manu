library(tidyverse)
library(here)
library(flextable)
library(lme4)
library(lmerTest)
library(performance)

source(here('R/funcs.R'))

load(file = here('data/scrmods.RData'))

fmt1 <- function(x) formatC(x, format = 'f', digits = 1)
fmt2 <- function(x) formatC(x, format = 'f', digits = 2)
fmt3 <- function(x) formatC(x, format = 'f', digits = 3)

# table of score trend model summaries by metric --------------------------

# N and number of groups are the same across all models
allN <- unique(sapply(scrmods, nobs))
allGroups <- unique(sapply(scrmods, function(mod) unname(summary(mod)$ngrps)))
stopifnot(length(allN) == 1, length(allGroups) == 1)
totab <- c('Abundance', 'Blade Length', 'Short Shoot Density')

sumdat <- lapply(totab, function(varnm) {

  mod <- scrmods[[varnm]]
  cfs <- summary(mod)$coefficients
  civ <- confint(mod, parm = 'yrctr', quiet = TRUE)

  # random-intercept (among-group) variance, from the fitted model
  vc <- as.data.frame(VarCorr(mod))
  grpvar <- vc$vcov[vc$grp == 'grp']

  r2 <- r2_nakagawa(mod)

  data.frame(
    Metric = varnm,
    Estimate = cfs['yrctr', 'Estimate'],
    lwr = civ[1, 1],
    upr = civ[1, 2],
    pval = cfs['yrctr', 'Pr(>|t|)'],
    grpsd = sqrt(grpvar),
    r2marg = r2$R2_marginal,
    r2cond = r2$R2_conditional,
    AIC = AIC(mod)
  )

}) |>
  bind_rows() |>
  mutate(
    sig = case_when(pval < 0.001 ~ '**', pval < 0.05 ~ '*', TRUE ~ ''),
    `Change yr⁻¹ (95% CI)` = paste0(fmt1(Estimate), sig, ' (', fmt1(lwr), ', ', fmt1(upr), ')'),
    `Group SD` = fmt1(grpsd),
    `R² (marg/cond)` = paste0(fmt2(r2marg), ' / ', fmt2(r2cond)),
    AIC = fmt1(AIC)
  ) |>
  select(Metric, `Change yr⁻¹ (95% CI)`, `Group SD`, `R² (marg/cond)`, AIC)

scrmodstab <- flextable(sumdat) |>
  autofit() |>
  align(align = 'center', part = 'all') |>
  align(j = 'Metric', align = 'left', part = 'all') |>
  add_footer_lines(values = paste0(
    'N = ', allN, ', Groups = ', allGroups, ' for all models.'
  ))

save(scrmodstab, file = here('tabs/scrmodstab.RData'), compress = 'xz')

# table of absolute deviation trend model summaries by metric -------------
# unlike scrmods above, these models have two crossed random intercepts
# (grp and Species, see dat_proc.R), so both get their own SD column.
# N varies slightly by metric (a few missing deviations per metric, see
# dat_proc.R), so it's reported per row instead of a single footer value;
# Groups and Species counts are constant across models, same convention
# as scrmods' footer

load(file = here('data/spmods.RData'))

allGroups2 <- unique(sapply(spmods, function(mod) unname(summary(mod)$ngrps['grp'])))
allSpecies <- unique(sapply(spmods, function(mod) unname(summary(mod)$ngrps['Species'])))
stopifnot(length(allGroups2) == 1, length(allSpecies) == 1)

sumdat2 <- lapply(names(spmods), function(varnm) {

  mod <- spmods[[varnm]]
  cfs <- summary(mod)$coefficients
  civ <- confint(mod, parm = 'yrctr', quiet = TRUE)

  # random-intercept (among-group, among-species) variances, from the fitted model
  vc <- as.data.frame(VarCorr(mod))
  grpvar <- vc$vcov[vc$grp == 'grp']
  sppvar <- vc$vcov[vc$grp == 'Species']

  r2 <- r2_nakagawa(mod)

  data.frame(
    Metric = varnm,
    N = nobs(mod),
    Estimate = cfs['yrctr', 'Estimate'],
    lwr = civ[1, 1],
    upr = civ[1, 2],
    pval = cfs['yrctr', 'Pr(>|t|)'],
    grpsd = sqrt(grpvar),
    sppsd = sqrt(sppvar),
    r2marg = r2$R2_marginal,
    r2cond = r2$R2_conditional,
    AIC = AIC(mod)
  )

}) |>
  bind_rows() |>
  mutate(
    sig = case_when(pval < 0.001 ~ '**', pval < 0.05 ~ '*', TRUE ~ ''),
    `Change yr⁻¹ (95% CI)` = paste0(fmt3(Estimate), sig, ' (', fmt3(lwr), ', ', fmt3(upr), ')'),
    `Group SD` = fmt2(grpsd),
    `Species SD` = fmt2(sppsd),
    `R² (marg/cond)` = paste0(fmt2(r2marg), ' / ', fmt2(r2cond)),
    AIC = fmt1(AIC)
  ) |>
  select(Metric, N, `Change yr⁻¹ (95% CI)`, `Group SD`, `Species SD`, `R² (marg/cond)`, AIC)

spmodstab <- flextable(sumdat2) |>
  fontsize(size = 9, part = 'all') |>
  autofit() |>
  align(align = 'center', part = 'all') |>
  align(j = 'Metric', align = 'left', part = 'all') |>
  add_footer_lines(values = paste0(
    'Groups = ', allGroups2, ', Species = ', allSpecies,
    ' for all models. '
  ))

save(spmodstab, file = here('tabs/spmodstab.RData'), compress = 'xz')

# matrix table: slope of abs(devval) ~ yr per group x species, one combined --
# flextable with a row block per metric (rows = species within each block,
# columns = group), rather than three separate tables. p-values are ignored
# (each regression has only 4-9 points, one per year that group/species pair
# has data), same six groups as the per-group models in R/scratch.R
# (FDEP/HC-ES dropped for too few years). Cells are colored by slope,
# diverging red (increasing deviation, getting worse) to green (decreasing,
# getting better), white at a slope of exactly 0, matching the direction of
# spmodstab above (negative = improvement). A single shared color domain
# across all three metrics wouldn't make sense (Blade Length slopes run
# -1.9 to 1.7 cm/yr, Abundance only -0.16 to 0.25 category-positions/yr),
# but flextable's bg() accepts a color function scoped to a row/column
# subset, so each metric's block of rows gets its own domain within the
# one table instead of needing three separate table objects

load(file = here('data/rawdiffscrs.RData'))

spvars <- names(spmods)

spdat <- rawdiffscrs |>
  mutate(
    absdev  = abs(devval),
    yrctr   = yr - min(yr),
    Species = factor(Species)
  ) |>
  filter(!is.na(absdev))

# same six groups as the per-group models in R/scratch.R
grpdat <- spdat |> filter(!grp %in% c('FDEP', 'HC-ES'))

slopedat <- grpdat |>
  summarise(
    slope = coef(lm(absdev ~ yr))[['yr']],
    .by = c(var, grp, Species)
  )

slopegrpcols <- sort(unique(slopedat$grp))

# anonymize group identity as A, B, C... in this table only; the mapping is
# printed below for reference but not saved anywhere in slopetab itself
grplabs <- setNames(LETTERS[seq_along(slopegrpcols)], slopegrpcols)
message('slopetab group key: ', paste(names(grplabs), grplabs, sep = ' = ', collapse = ', '))

slopewide <- slopedat |>
  mutate(
    var = factor(var, levels = spvars),
    grp = grplabs[grp]
  ) |>
  tidyr::pivot_wider(names_from = grp, values_from = slope) |>
  arrange(var, Species) |>
  select(Metric = var, Species, dplyr::all_of(sort(unname(grplabs))))

slopegrpcols <- sort(unname(grplabs))

slopetab <- flextable(slopewide) |>
  fontsize(size = 10, part = 'all') |>
  colformat_double(j = slopegrpcols, digits = 2) |>
  merge_v(j = 'Metric') |>
  valign(j = 'Metric', valign = 'top') |>
  fix_border_issues() |>
  add_header_row(values = c('', 'Groups'), colwidths = c(2, length(slopegrpcols)), top = TRUE) |>
  fix_border_issues()

for (vr in spvars) {
  rows <- which(slopewide$Metric == vr)
  rng <- max(abs(unlist(slopewide[rows, slopegrpcols])), na.rm = TRUE)
  slopetab <- bg(
    slopetab, i = rows, j = slopegrpcols,
    bg = scales::col_numeric(palette = c('#1a9850', '#f7f7f7', '#d73027'), domain = c(-rng, rng))
  )
}

slopetab <- slopetab |>
  align(align = 'center', part = 'all') |>
  align(j = c('Metric', 'Species'), align = 'left', part = 'all') |>
  bold(part = 'header') |>
  autofit()

save(slopetab, file = here('tabs/slopetab.RData'), compress = 'xz')

