
#--------------------------------------------------------
# Packages, global variables
#--------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)
theme_set(theme_bw()) # ggplot theme
library(lmerTest)
library(xtable)
library(cowplot)

#------------------------------------------------------------
# Helper functions
#------------------------------------------------------------

# conf_summary()
# partially based on: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
# group data, and calculate: mean, standard error, 95% confidence interval
conf_summary <- function(data, measurevar, ...) {
  d %>% 
    group_by(...) %>%
    summarise(mean = mean({{measurevar}}, na.rm = T),
              sd = sd({{measurevar}}, na.rm = T),
              N = n(),
              se = sd / sqrt(N),
              ci = se * qt(0.95/2 + .5, N-1)) %>% 
    ungroup()
}

# save plots with cairo
ggsave.alt <- function(plot.lst = list(last_plot()), custom.name = "lastplot",
                       height = fh, width = fw, figscale = 1,
                       figdir = '../plots/',
                       device = cairo_pdf, device_suffix = 'pdf', ...) {
  if (is.null(names(plot.lst))) names(plot.lst) <- custom.name
  plot.names <- paste0(figdir, names(plot.lst), '.', device_suffix)
  pmap(list(plot.names, plot.lst), ggsave,
       device = device, scale = figscale,
       height = height, width = width, ...)
}

#--------------------------------------------------------
# Load data
#--------------------------------------------------------

# semantic labels for conditions
condlabels <- tibble(condition = c(1, 2, 3),
                     condition_label = c('Cond. 1: constant IOI & ISI',
                                         'Cond. 2: variable ISI',
                                         'Cond. 3: variable IOI'))

# read response data
d <- read_csv('../data/seh2_response_data.csv') %>% 
  left_join(condlabels) %>% 
  mutate(condition = as_factor(condition))

#--------------------------------------------------------
# Distributions
#--------------------------------------------------------

# Conditions were conducted on different machines:
# Conditions 1 and 2: Dell XPS M1330, Ubuntu 12.04.
# Condition 3: Lenovo T440s, Ubuntu 15.10.
# RTs are consistently shorter for condition 3.

ggplot(d) + aes(x = relRT, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5)

ggplot(d) + aes(x = zRT, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5)

conf_summary(d, relRT, condition_label)
conf_summary(d, zRT, condition_label)

#--------------------------------------------------------
# RTs by probe
#--------------------------------------------------------

## Summary by probe
dsum2 <- conf_summary(d, zRT, probe)
ggplot(dsum2) + aes(x = factor(probe), y = mean, group = 1) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = "Position of probe within sequence",
       y = "Normalised reaction time (z)")
# ggsave("dsum_all.pdf", scale = 2)

## Summary by probe and condition
dsum <- conf_summary(d, zRT, probe, condition_label) %>% 
  rename(condition = condition_label)
ggplot(dsum) + aes(x = factor(probe), y = mean,
                   group = condition, color = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = "Position of probe within sequence",
       y = "Normalised reaction time (z)")
# ggsave("dsum_conditions.pdf", scale = 2)

ggplot(dsum) + aes(x = factor(probe), y = mean, group = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = "Position of probe within sequence",
       y = "Normalised reaction time (z)") +
  theme(legend.position="none") +
  theme_cowplot()
# ggsave(paste0(figdir, "dsum_conditions2.pdf"), width = 10, height = 6)

# summary statistics (for plot overlay)
dsum.stat <- dsum %>% 
  group_by(condition) %>% 
  summarise(r = cor(probe, mean) %>% round(2)) %>% 
  mutate(highlight = condition)

# with hihglight
dsum %>% 
  crossing(highlight = unique(.$condition)) %>%
  ggplot() +
  aes(x = factor(probe), y = mean, group = condition,
      alpha = highlight == condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = "Position of probe within sequence",
       y = "Normalised reaction time (z)") +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.1, 1)) +
  geom_text(aes(x = 7, y = 1.1, label = paste0("italic(r) == ", r)), parse = T, 
            data = dsum.stat) +
  theme_cowplot() +
  theme(legend.position="none")

ggsave.alt(custom.name = "dsum_conditions3", width = 7.5, height = 5)
ggsave.alt(custom.name = "dsum_conditions3", width = 7.5, height = 5,
           device = 'jpeg', device_suffix = 'jpg')

#--------------------------------------------------------
# mixed models
#--------------------------------------------------------

# probe + probe.dist
mm1 <- lmer(zRT ~ probe + probe.dist + (1 + probe | subjectID),
            d, REML = F)
summary(mm1)
mm1.null <- lmer(zRT ~ probe.dist + (1 + probe | subjectID),
                 d, REML = T)
anova(mm1, mm1.null)

# probe + probe.dist + condition
mm2 <- lmer(zRT ~ probe + probe.dist + condition + (1 + probe | subjectID),
            d, REML = F)
summary(mm2)
mm2.null <- lmer(zRT ~ probe + probe.dist + (1 + probe | subjectID),
                 d, REML = T)
anova(mm2, mm2.null)

# probe + probe.dist + condition*probe
mm3 <- lmer(zRT ~ probe.dist + probe + condition + (probe*condition) +
              (1 + probe | subjectID),
            d, REML = F)
summary(mm3)
mm3.null <- lmer(zRT ~ probe.dist + probe + condition + (1 + probe | subjectID),
                 d, REML = T)
anova(mm3, mm3.null)

#--------------------------------------------------------

print(xtable(summary(mm3)$coefficients, digits = 4),
      scalebox = '.8', include.rownames = T, booktabs = T)

print(xtable(anova(mm3), digits = 4),
      scalebox = '.8', include.rownames = T, booktabs = T)
