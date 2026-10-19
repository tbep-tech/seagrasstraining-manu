library(tidyverse)
library(here)
library(tbeptools)
library(lme4)
library(lmerTest)
library(patchwork)

# leave-one-site-out sensitivity --------------------------------------------
# is any single Site disproportionately driving a group's score? for each
# year, drop one Site at a time, recompute every group's score from the
# remaining Sites, and compare it to the group's full-sample score (all
# Sites). unlike the 50/50 split-half test, removing just one Site barely
# changes the total sample -- almost every Site x Species consensus pair
# used by truespp_fun()/truvar_fun() survives each iteration, since only
# one Site is missing rather than half of them. so this tests whether a
# group's score depends too heavily on any one Site, rather than testing
# reliability under a materially smaller sample -- the same idea as a
# leverage/influence diagnostic (e.g. Cook's distance) rather than
# classical split-half reliability. exhaustive (one iteration per Site per
# year), not random, so no repeated draws or seed are needed.
#
# run twice: once on the raw (pre-calibration) weighted-mean deviations
# (what calibrate_scr_fun() itself works with), and once on the final
# calibrated 0-100 scores. the two shouldn't necessarily move in lockstep:
# the calibrated score also depends on that focal year's own cross-group
# spread (used to set the score floor, see allgrpscr_fun()), which itself
# shifts slightly when a Site is dropped, on top of whatever the raw
# deviation does -- so the calibrated version is a noisier, compound
# version of the same question, not a strict repeat of it.

funcsurl <- 'https://raw.githubusercontent.com/tbep-tech/seagrasstransect-training-reports/main/R/funcs.R'
funcsfl <- paste(tempdir(), 'funcs.R', sep = '/')
download.file(funcsurl, destfile = funcsfl, quiet = TRUE)
source(funcsfl)

load(file = here('data/trndat.rda'))
data(file = 'trnlns', package = 'tbeptools')


ggsave(here('figs/leaveoneout.png'), leaveoutplot, width = 11, height = 6.5, units = 'in', dpi = 300)
