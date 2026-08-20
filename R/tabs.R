library(tidyverse)
library(here)
library(flextable)
library(lme4)
library(lmerTest)
library(performance)

source(here('R/funcs.R'))

# table of score trend model summaries by metric --------------------------

load(file = here('data/scrmods.RData'))

fmt1 <- function(x) formatC(x, format = 'f', digits = 1)
fmt2 <- function(x) formatC(x, format = 'f', digits = 2)

# N and number of groups are the same across all models
allN <- unique(sapply(scrmods, nobs))
allGroups <- unique(sapply(scrmods, function(mod) unname(summary(mod)$ngrps)))
stopifnot(length(allN) == 1, length(allGroups) == 1)

sumdat <- lapply(names(scrmods), function(varnm) {

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
  set_caption('Random-intercept model summaries of score trend by metric.') |>
  add_footer_lines(values = paste0(
    '* p < 0.05, ** p < 0.001; N = ', allN, ', Groups = ', allGroups, ' for all models.'
  ))

save(scrmodstab, file = here('tabs/scrmodstab.RData'), compress = 'xz')
