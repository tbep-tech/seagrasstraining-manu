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
