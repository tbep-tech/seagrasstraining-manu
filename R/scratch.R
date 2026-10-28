ggplot(allyrscrs, aes(x = factor(yr), y = scr)) + 
  geom_boxplot() + 
  facet_wrap(~var, ncol = 1)
