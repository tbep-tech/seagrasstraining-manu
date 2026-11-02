library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(performance)

source(here('R/funcs.R'))

# rebuild the same trendgrps-restricted data scrmods/spmods are fit to
# (not saved on their own, so recomputed here from the saved intermediates)
load(file = here('data/allyrscrs.RData'))
load(file = here('data/rawdiffscrs.RData'))
load(file = here('data/trendgrps.RData'))

alldat <- allyrscrs |>
  filter(grp %in% trendgrps) |>
  mutate(yrctr = yr - min(yr))

spvars <- c('Abundance', 'Blade Length', 'Short Shoot Density')

spdat <- rawdiffscrs |>
  filter(grp %in% trendgrps) |>
  mutate(absdev = abs(devval)) |>
  summarise(absdev = mean(absdev, na.rm = TRUE), .by = c(var, grp, yr)) |>
  filter(!is.nan(absdev)) |>
  mutate(yrctr = yr - min(yr))

# does a random SLOPE for group (in addition to the random intercept) ------
# change anything? random-intercept-only is the production choice (see
# R/dat_proc.R); this compares it against (1 + yrctr | grp) for both the
# score model (Gaussian LMM) and the deviation model (Gamma GLMM), one
# random-intercept vs. random-intercept+slope pair per metric

# scrmods (metric scores, Gaussian LMM) -------------------------------------
scr_fits <- map(spvars, function(vr) {
  vardat <- alldat |> filter(var == vr)
  ri <- lmer(scr ~ yrctr + (1 | grp), data = vardat)
  rs <- lmer(scr ~ yrctr + (1 + yrctr | grp), data = vardat)
  list(ri = ri, rs = rs)
})
names(scr_fits) <- spvars

# spmods (raw deviations, Gamma GLMM via glmmTMB) ----------------------------
dev_fits <- map(spvars, function(vr) {
  vardat <- spdat |> filter(var == vr)
  eps <- min(vardat$absdev[vardat$absdev > 0]) / 2
  vardat <- vardat |> mutate(absdev = absdev + eps * (absdev == 0))
  ri <- glmmTMB(absdev ~ yrctr + (1 | grp), family = Gamma(link = 'log'), data = vardat)
  rs <- glmmTMB(absdev ~ yrctr + (1 + yrctr | grp), family = Gamma(link = 'log'), data = vardat)
  list(ri = ri, rs = rs)
})
names(dev_fits) <- spvars

# comparison table: trend estimate/p-value under each random-effects
# structure, AIC, convergence/singularity, and an LRT comparing the two
# (same fixed effects, nested random effects -- a valid likelihood ratio
# test either way)
comptab <- map(c('scr', 'dev'), function(mtype) {

  fits <- if (mtype == 'scr') scr_fits else dev_fits

  map(spvars, function(vr) {

    ri <- fits[[vr]]$ri
    rs <- fits[[vr]]$rs

    lrt <- tryCatch(anova(ri, rs), error = function(e) NULL)
    lrt_p <- if (!is.null(lrt)) lrt[['Pr(>Chisq)']][2] else NA

    get_row <- function(mod, label) {
      cfs <- mod_coefs(mod)
      tibble(
        type     = mtype,
        metric   = vr,
        model    = label,
        estimate = cfs['yrctr', 'Estimate'],
        pval     = cfs['yrctr', grepl('^Pr', colnames(cfs))],
        AIC      = AIC(mod),
        singular = if (inherits(mod, 'glmmTMB')) check_singularity(mod) else isSingular(mod)
      )
    }

    bind_rows(
      get_row(ri, 'rand intercept'),
      get_row(rs, 'rand intercept + slope')
    ) |>
      mutate(lrt_pval = lrt_p)

  }) |>
    bind_rows()

}) |>
  bind_rows()

print(comptab, n = Inf)

message(
  '\nLRT p-values (same for both rows of a metric x type, testing whether\n',
  'adding a random slope improves fit): ', 'see lrt_pval column above.'
)
