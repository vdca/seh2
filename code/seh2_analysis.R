
#--------------------------------------------------------
# Packages, global variables
#--------------------------------------------------------

# remove previous objects from workspace
rm(list=ls())

library(tidyverse)
library(cowplot)
library(lmerTest)
library(broom.mixed)
library(gtools)
library(kableExtra)

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

# save plots with cairo pdf device
ggsave_pdf <- function(custom_name = "lastplot", plot_lst = list(last_plot()),
                       height = 5, width = 7.5, figscale = 1,
                       figdir = '../plots/',
                       device = cairo_pdf, device_suffix = 'pdf', ...) {
  if (is.null(names(plot_lst))) names(plot_lst) <- custom_name
  plot_names <- paste0(figdir, names(plot_lst), '.', device_suffix)
  pmap(list(plot_names, plot_lst), ggsave,
       device = device, scale = figscale,
       height = height, width = width, ...)
}

ggsave_jpg <- function(custom_name = "lastplot", ...) {
  ggsave_pdf(device = 'jpeg', device_suffix = 'jpg', custom_name, ...)
}

# print fixed effects of mixed model
print_model <- function(mx) {
  mx %>% 
    tidy(effects = 'fixed') %>% 
    mutate(signif = stars.pval(p.value),
           t.value = statistic,
           p.value = format.pval(p.value, digits = 2),
           term = str_replace_all(term, 'probe.dist', 'probe.distance')) %>% 
    mutate_if(is.numeric, ~ round(., digits = 2)) %>% 
    select(term, estimate, std.error, t.value, p.value, signif) %>% 
    kable(format = 'latex', booktabs = T, linesep = "")
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
  mutate(condition = as_factor(condition)) %>% 
  group_by(condition) %>% 
  mutate(zRT = scale(logRT),
         stdRT = scale(relRT)) %>% 
  ungroup()

#--------------------------------------------------------
# Distributions
#--------------------------------------------------------

# themes for distribution plots.
# legend.justification = inside reference point for legend.position 
theme_dist <- list(theme_cowplot(),
                   theme(legend.title = element_blank(),
                         legend.position = c(1, 1),
                         legend.justification = c('right', 'top')))
theme_dist_nolegend <- list(theme_cowplot(),
                            theme(legend.position = 'none'))

# Conditions were conducted on different machines:
# Condition 1: Dell XPS M1330, Ubuntu 12.04.
# Condition 2: Dell XPS M1330, Ubuntu 12.04.
# Condition 3: Lenovo T440s, Ubuntu 15.10.
# RTs are consistently shorter for condition 3.
# Hence, z-normalise RTs to make conditions comparable.

dist_rt <- ggplot(d) +
  aes(x = relRT, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5) +
  theme_dist + xlab('reaction time (ms)')
dist_rt

dist_logrt <- ggplot(d) +
  aes(x = logRT, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5) +
  theme_dist_nolegend + xlab('reaction time (log-transformed)')
dist_logrt

dist_stdrt <- ggplot(d) +
  aes(x = stdRT, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5) +
  theme_dist_nolegend + xlab('reaction time (z-normalised)')
dist_stdrt

dist_zrt <- ggplot(d) +
  aes(x = zRT, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5) +
  theme_dist_nolegend + xlab('reaction time (log-transformed, z-normalised)')
dist_zrt

plot_grid(dist_rt, dist_zrt, ncol = 1, labels = 'AUTO')
ggsave_jpg('rt_distributions_2', height = 9)

plot_grid(dist_rt, dist_logrt, dist_zrt, ncol = 1, labels = 'AUTO')
ggsave_jpg('rt_distributions_3', height = 9)

conf_summary(d, relRT, condition_label)
conf_summary(d, logRT, condition_label)
conf_summary(d, zRT, condition_label)

#--------------------------------------------------------
# RTs by probe
#--------------------------------------------------------

# axis labels
probe_lab <- 'position of probe within sequence'
zrt_lab <- 'reaction time (log-transformed, z-normalised)'
lrt_lab <- 'reaction time (log-transformed)'

## Summary by probe
dsum2 <- conf_summary(d, zRT, probe)
ggplot(dsum2) + aes(x = factor(probe), y = mean, group = 1) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = probe_lab, y = zrt_lab)
# ggsave("dsum_all.pdf", scale = 2)

## Summary by probe and condition
dsum <- conf_summary(d, zRT, probe, condition_label) %>% 
  rename(condition = condition_label)
ggplot(dsum) + aes(x = factor(probe), y = mean,
                   group = condition, color = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = probe_lab, y = zrt_lab) +
  theme_cowplot() +
  theme(legend.title = element_blank())
# ggsave("dsum_conditions.pdf", scale = 2)

#--------------------------------------------------------
# log-transformed, z-normalised RTs by probe (with overlay)
#--------------------------------------------------------

# summary statistics (for plot overlay)
dsum.stat <- dsum %>% 
  group_by(condition) %>% 
  summarise(r = cor(probe, mean) %>% round(3)) %>% 
  mutate(highlight = condition)

# no highlight
dsum %>% 
  ggplot() + aes(x = factor(probe), y = mean, group = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = probe_lab, y = zrt_lab) +
  theme(legend.position="none") +
  geom_text(aes(x = 7, y = 1.1, label = paste0("italic(r) == ", r)), parse = T, 
            data = dsum.stat, inherit.aes = F) +
  theme_cowplot() +
  background_grid(major = 'y')

ggsave_jpg(custom_name = "zrt_conditions", width = 9)

# with hihglight
dsum %>% 
  crossing(highlight = unique(.$condition)) %>%
  ggplot() +
  aes(x = factor(probe), y = mean, group = condition,
      alpha = highlight == condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = probe_lab, y = zrt_lab) +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.1, 1)) +
  geom_text(aes(x = 7, y = 1.1, label = paste0("italic(r) == ", r)), parse = T, 
            data = dsum.stat, inherit.aes = F) +
  theme_cowplot() +
  background_grid(major = 'y') +
  theme(legend.position = 'none')

ggsave_jpg(custom_name = "zrt_conditions_hi", width = 9)

#--------------------------------------------------------
# log-transformed RTs by probe (with overlay)
#--------------------------------------------------------

## summary by probe and condition
log_sum <- conf_summary(d, logRT, probe, condition_label) %>% 
  rename(condition = condition_label)

# summary statistics (for plot overlay)
log_sum_stat <- log_sum %>% 
  group_by(condition) %>% 
  summarise(r = cor(probe, mean) %>% round(3)) %>% 
  mutate(highlight = condition)

# no highlight
ggplot(log_sum) + aes(x = factor(probe), y = mean, group = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_text(aes(x = 7, y = 6.5, label = paste0("italic(r) == ", r)),
            parse = T, data = log_sum_stat) +
  facet_wrap(~ condition) +
  labs(x = probe_lab, y = lrt_lab) +
  theme_cowplot() +
  background_grid(major = 'y')

ggsave_jpg(custom_name = "logrt_conditions", width = 9)

# with hihglight
log_sum %>% 
  crossing(highlight = unique(.$condition)) %>%
  ggplot() +
  aes(x = factor(probe), y = mean, group = condition,
      alpha = highlight == condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = probe_lab, y = lrt_lab) +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.1, 1)) +
  geom_text(aes(x = 7, y = 6.5, label = paste0("italic(r) == ", r)),
            parse = T, data = log_sum_stat) +
  theme_cowplot() +
  background_grid(major = 'y') +
  theme(legend.position="none")

ggsave_jpg(custom_name = "logrt_conditions_hi", width = 9)

#--------------------------------------------------------
# regression models (zRT)
#--------------------------------------------------------

# saturated model
mm_sat <- lmer(zRT ~ probe.dist + probe + condition +
              probe:condition + probe.dist:condition +
              (1|subjectID), d, REML = F)
summary(mm_sat)

# model selection
step(mm_sat)

# step-selected final model
mm_final <- lmer(zRT ~ probe.dist + probe +
                   (1|subjectID), d, REML = F)
summary(mm_final)
print_model(mm_final)

#--------------------------------------------------------
# regression models (logRT)
#--------------------------------------------------------

# saturated model
mm_sat_log <- lmer(logRT ~ probe.dist + probe + condition +
                 probe:condition + probe.dist:condition +
                 (1|subjectID), d, REML = F)
summary(mm_sat_log)

# model selection
step(mm_sat_log)

# final model
mm_final_log <- lmer(logRT ~ probe.dist + probe + condition + probe:condition +
                   (1|subjectID), d, REML = F)
summary(mm_final_log)
print_model(mm_final_log)



