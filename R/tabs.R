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

# table of score and deviation trend model summaries, one combined --------

load(file = here('data/spmods.RData'))

allN <- unique(c(sapply(scrmods, nobs), sapply(spmods, nobs)))
allGroups <- unique(c(
  sapply(scrmods, function(mod) unname(summary(mod)$ngrps)),
  sapply(spmods, function(mod) unname(summary(mod)$ngrps['grp']))
))
stopifnot(length(allN) == 1, length(allGroups) == 1)
totab <- c('Abundance', 'Blade Length', 'Short Shoot Density')

modsumm <- function(mod, fmtfun) {

  cfs <- summary(mod)$coefficients
  civ <- confint(mod, parm = 'yrctr', quiet = TRUE)

  # random-intercept (among-group) variance, from the fitted model
  vc <- as.data.frame(VarCorr(mod))
  grpvar <- vc$vcov[vc$grp == 'grp']

  r2 <- r2_nakagawa(mod)

  est <- cfs['yrctr', 'Estimate']
  lwr <- civ[1, 1]
  upr <- civ[1, 2]
  pval <- cfs['yrctr', 'Pr(>|t|)']
  sig <- case_when(pval < 0.001 ~ '**', pval < 0.05 ~ '*', TRUE ~ '')

  data.frame(
    `Change yr⁻¹ (95% CI)` = paste0(fmtfun(est), sig, ' (', fmtfun(lwr), ', ', fmtfun(upr), ')'),
    `Group SD` = fmt2(sqrt(grpvar)),
    `R² (marg/cond)` = paste0(fmt2(r2$R2_marginal), ' / ', fmt2(as.numeric(r2$R2_conditional))),
    AIC = fmt1(AIC(mod)),
    check.names = FALSE
  )

}

scrsumm <- lapply(totab, function(varnm) cbind(Metric = varnm, modsumm(scrmods[[varnm]], fmt1))) |>
  bind_rows() |>
  mutate(block = 'Metric scores')
spsumm <- lapply(totab, function(varnm) cbind(Metric = varnm, modsumm(spmods[[varnm]], fmt3))) |>
  bind_rows() |>
  mutate(block = 'Raw metric deviations')

sumdat <- bind_rows(scrsumm, spsumm)

# as_grouped_data() inserts a full-width label row ahead of each block
grouped <- as_grouped_data(sumdat, groups = 'block')
blockrows <- which(is.na(grouped$Metric))

modstab <- as_flextable(grouped, hide_grouplabel = TRUE) |>
  autofit() |>
  align(align = 'center', part = 'all') |>
  align(j = 'Metric', align = 'left', part = 'all') |>
  bold(i = blockrows, part = 'body') |>
  add_footer_lines(values = paste0(
    'N = ', allN, ', Groups = ', allGroups, ' for each model.'
  ))

save(modstab, file = here('tabs/modstab.RData'), compress = 'xz')

# matrix table: slope of abs(devval) ~ yr per group x species, one combined --

load(file = here('data/rawdiffscrs.RData'))
load(file = here('data/trendgrps.RData'))

spvars <- names(spmods)

spdat <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  mutate(
    absdev  = abs(devval),
    yrctr   = yr - min(yr),
    Species = factor(Species)
  ) |>
  filter(!is.na(absdev))

slopedat <- spdat |>
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

# explicit column widths (not autofit(), which only sizes to natural
# content width and can run past the page) keep the font size fixed and
# the table within a standard 6.5in printable page width -- a group
# column's slope values are short enough to not need it, but a long
# Species name can still wrap onto a second line
page_width <- 6.5
metric_width <- 0.9
species_width <- 0.8
grp_width <- (page_width - metric_width - species_width) / length(slopegrpcols)

slopetab <- slopetab |>
  align(align = 'center', part = 'all') |>
  align(j = c('Metric', 'Species'), align = 'left', part = 'all') |>
  bold(part = 'header') |>
  width(j = 'Metric', width = metric_width) |>
  width(j = 'Species', width = species_width) |>
  width(j = slopegrpcols, width = grp_width) |>
  set_table_properties(layout = 'fixed')

save(slopetab, file = here('tabs/slopetab.RData'), compress = 'xz')

# training summary table: groups, metric scores, raw deviations, and the --
# consensus species list by year. Metric score and raw deviation rows reuse
# the same metric names as their Metric label, so a 'type' column (not the
# label text) is what keeps those two blocks apart during assembly, same
# reasoning as the per-year summary rows below

funcsurl <- 'https://raw.githubusercontent.com/tbep-tech/seagrasstransect-training-reports/main/R/funcs.R'
funcsfl <- paste(tempdir(), 'funcs.R', sep = '/')
download.file(funcsurl, destfile = funcsfl, quiet = TRUE)
source(funcsfl)

load(file = here('data/trndat.rda'))

yrs <- sort(unique(trndat$yr))
metvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')
scr_metric <- c(Abundance = 'abs', `Blade Length` = 'abs', `Short Shoot Density` = 'pct')

# fixed historical baseline, built once from the full dataset (same as
# production use)
cal <- calibrate_scr_fun(trndat, metric = scr_metric)

# per-year summary, kept long with a 'type' column distinguishing the four
# row categories: number of groups, metric scores and raw (pre-
# calibration) deviations (both averaged across groups, with (min, max) in
# one cell -- deviations are summarized in absolute terms since the raw
# differences are signed, scores are not since they're already a bounded
# 0-100 scale), and how many Sites each consensus species (truespp_fun(),
# >= 2 groups reporting) was found at
yrsumm <- purrr::map(yrs, function(yr) {

  yrdat <- trndat |> dplyr::filter(yr == !!yr)
  tv <- truvar_fun(trndat, yr)

  n_groups <- yrdat |>
    validgrp() |>
    dplyr::pull(grp) |>
    unique() |>
    length()

  scrsumm <- allgrpscr_fun(trndat, yr, tv, cal = cal, raw = TRUE, metric = scr_metric) |>
    validgrp() |>
    tidyr::pivot_longer(dplyr::all_of(metvars), names_to = 'Species', values_to = 'val') |>
    dplyr::summarise(val = fmtavgrng(val, 0, show_max = FALSE), .by = 'Species') |>
    dplyr::mutate(type = 'score')

  devsumm <- allgrpscr_fun(trndat, yr, tv, raw_diff = TRUE, metric = scr_metric) |>
    validgrp() |>
    tidyr::pivot_longer(dplyr::all_of(metvars), names_to = 'Species', values_to = 'val') |>
    dplyr::summarise(val = fmtavgrng(abs(val), 1), .by = 'Species') |>
    dplyr::mutate(type = 'deviation')

  truespp <- truespp_fun(trndat, yr)
  sppcnt <- truespp |>
    dplyr::count(Species, name = 'val') |>
    dplyr::mutate(val = as.character(val), type = 'species')

  dplyr::bind_rows(
    tibble::tibble(Species = 'Number of groups', val = as.character(n_groups), type = 'groups'),
    scrsumm, devsumm, sppcnt
  ) |>
    dplyr::mutate(yr = yr)

}) |>
  dplyr::bind_rows()

# macroalgae species codes share these prefixes (same convention used in
# the scoring repo's docs, see 'non' in docs/scoring.qmd), anything else is
# SAV -- SAV listed first, each block alphabetical within itself
non <- '^AA|^DA|^DB|^DG|^DR|^CYAN'

allspp <- yrsumm |>
  dplyr::filter(type == 'species') |>
  dplyr::pull(Species) |>
  unique()

savspp <- sort(allspp[!grepl(non, allspp)])
algspp <- sort(allspp[grepl(non, allspp)])

yrcols <- as.character(yrs)

yrwide_all <- yrsumm |>
  tidyr::pivot_wider(id_cols = c(Species, type), names_from = yr, values_from = val, values_fill = '-') |>
  dplyr::rename(Metric = Species)

yrwide <- dplyr::bind_rows(
  yrwide_all |> dplyr::filter(type == 'groups'),
  blank_row('Metric scores', yrcols),
  yrwide_all |> dplyr::filter(type == 'score') |> dplyr::arrange(factor(Metric, levels = metvars)),
  blank_row('Raw metric deviations', yrcols),
  yrwide_all |> dplyr::filter(type == 'deviation') |> dplyr::arrange(factor(Metric, levels = metvars)),
  blank_row('Seagrass', yrcols),
  yrwide_all |> dplyr::filter(type == 'species', Metric %in% savspp) |> dplyr::arrange(Metric),
  blank_row('Macroalgae', yrcols),
  yrwide_all |> dplyr::filter(type == 'species', Metric %in% algspp) |> dplyr::arrange(Metric)
) |>
  dplyr::select(-type)

hdr_i <- which(yrwide$Metric %in% c('Metric scores', 'Raw metric deviations', 'Seagrass', 'Macroalgae'))
spp_i <- which(yrwide$Metric %in% c(savspp, algspp))

# compact layout for the Word render: fixed font size and tighter row
# padding (default padding leaves a lot of vertical whitespace per row
# given how many rows this table has). Column widths are set explicitly
# (not autofit()/fit_to_width(), which shrinks font size -- see
# flextable:::fit_to_width -- to force content to fit instead of wrapping
# it) and layout = 'fixed' tells Word to honor those widths rather than
# resize the table to fit its content; a cell whose text is wider than its
# column just wraps onto more lines, which is fine here
page_width <- 6.5
metric_width <- 1.2
yr_width <- (page_width - metric_width) / length(yrcols)

trainsumtab <- flextable(yrwide) |>
  fontsize(size = 10, part = 'all') |>
  padding(padding.top = 1, padding.bottom = 1, part = 'all') |>
  italic(i = spp_i) |>
  bold(i = c(1, hdr_i)) |>
  bold(part = 'header') |>
  align(align = 'center', part = 'all') |>
  align(j = 'Metric', align = 'left', part = 'all') |>
  hline(i = c(1, hdr_i[2] - 1, hdr_i[3] - 1, hdr_i[4] - 1)) |>
  set_header_labels(Metric = '') |>
  width(j = 'Metric', width = metric_width) |>
  width(j = yrcols, width = yr_width) |>
  set_table_properties(layout = 'fixed')

for (i in hdr_i) {
  trainsumtab <- merge_at(trainsumtab, i = i, j = seq_len(length(yrcols) + 1), part = 'body')
}

save(trainsumtab, file = here('tabs/trainsumtab.RData'), compress = 'xz')

