
#--------------------------------------------------------
# packages, global variables
#--------------------------------------------------------

# remove previous objects from workspace
rm(list=ls())

library(tidyverse)
library(skimr)
library(cowplot)
theme_set(theme_cowplot())
library(ggpubr)
library(ggeffects)
library(lmerTest)
library(LMERConvenienceFunctions)
library(broom.mixed)
library(gtools)
library(kableExtra)
library(rstatix) # only for pvalue formatting (without leading zeroes)

#------------------------------------------------------------
# helper functions
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

ggsave_png <- function(custom_name = "lastplot", ...) {
  ggsave_pdf(device = 'png', device_suffix = 'png', custom_name, ...)
}

# print table in latex format
ltx_tbl <- function(mx) {
  mx %>% 
    mutate(across(where(is.numeric), ~ round(., digits = 2))) %>% 
    kable(format = 'latex', booktabs = T, linesep = "", escape = T)
}

# print fixed effects of mixed model
print_model <- function(mx) {
  
  # 95% confidence intervals
  mx_ci <- confint(mx, method = 'Wald') %>% 
    as.data.frame() %>% 
    slice_tail(n = nrow(.)-2) %>% 
    rownames_to_column('term')
  colnames(mx_ci) <- c('term', 'lower.CI', 'upper.CI')
  
  mx %>% 
    tidy(effects = 'fixed') %>% 
    left_join(mx_ci) %>% 
    mutate(signif = stars.pval(p.value),
           t.value = statistic,
           p.value = p_format(p.value, digits = 1,
                              accuracy = .001, leading.zero = F), # APA guidelines
           term = str_replace_all(term, 'deviant_probability', 'dev.probability'),
           term = str_replace_all(term, 'preceding_stds', 'dev.distance'),
           term = str_replace_all(term, 'deviant', 'dev.position')) %>% 
    select(term, estimate, lower.CI, upper.CI,
           std.error, t.value, p.value, signif) %>% 
    ltx_tbl() %>% 
    str_replace_all('signif', '') %>% 
    str_replace_all('<', '$<$') %>% # less-than symbol is not well escaped
    cat()
}

# themes for distribution plots
# legend.justification = inside reference point for legend.position 
theme_dist <- list(theme(legend.title = element_blank(),
                         legend.position = c(1, 1),
                         legend.justification = c('right', 'top')))
theme_dist_nolegend <- list(theme(legend.position = 'none'))

#--------------------------------------------------------
# load data, add condition-level features
#--------------------------------------------------------

# semantic labels for conditions.
# machine (hard/soft-ware) for conditions.
condlabels <- tibble(condition = c(1, 2, 3),
                     condition_label = c('1: constant ISI & ITI',
                                         '2: variable ITI',
                                         '3: variable ISI'),
                     machine = c('dell_ubuntu12',
                                 'dell_ubuntu12',
                                 'lenovo_ubuntu15'))

# read all response data (including misses, false alarms, etc.)
alld <- read_tsv('../data/seh2_response_data.tsv') %>% 
  left_join(condlabels) %>% 
  mutate(condition = as_factor(condition))

#--------------------------------------------------------
# exclude block-initial trials
#--------------------------------------------------------

# remove all block-initial trials,
# since the experiment didn't control that each block started
# with a filler trial (i.e. with no deviants in it).
# trials are numbered from 1 to 76 (col=itemID):
# 1:4 belong to training block,
# 5:40 (n=36) belong to experimental block 1
# 41:76 (n=36) belong to experimental block 2.
# hence, remove trials 5 and 41 from dataset.
alld %>% 
  group_by(section) %>% 
  filter(itemID == min(itemID)) %>% 
  count(itemID)

alld <- alld %>% 
  group_by(section) %>% 
  filter(itemID != min(itemID)) %>% 
  ungroup()

#--------------------------------------------------------
# overall accuracy; exclude low-performance participants
#--------------------------------------------------------

# compute hits, overreactions, misses, false alarms.
# coding:
# n_reactions = 0 = missed deviants
# n_reactions = 1 = hits
# n_reactions > 1 = overreactions, multiple reactions to single deviant
# n_reactions = -1 (arbitrary code; when deviant absent) = false alarm
alld <- alld %>% 
  group_by(subjectID, itemID) %>% 
  nest() %>% 
  mutate(n_reactions = map_dbl(data, nrow)) %>%
  unnest(data) %>% 
  ungroup() %>%
  mutate(n_reactions = if_else(is.na(relRT), 0, n_reactions),
         n_reactions = if_else(is.na(deviant), -1, n_reactions), # arbitrary code for false alarms
         response_type = if_else(n_reactions==-1, 'false_alarm', 'undefined!'),
         response_type = if_else(n_reactions>1, 'overreaction', response_type),
         response_type = if_else(n_reactions==1, 'hit', response_type),
         response_type = if_else(n_reactions==0, 'miss', response_type))

# address overreactions.
# multiple key presses may be correct responses to a deviant.
# in any case, keep only the first response for each trial.
# after doing so, mark all negative RTs as false_alarms,
# and re-code kept overreactions as legit hits.
alld <- alld %>% 
  arrange(subjectID, itemID, absRT) %>% 
  group_by(subjectID, itemID) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(response_type = if_else(relRT <= 0, 'false_alarm', response_type, response_type),
         response_type = if_else(response_type == 'overreaction', 'hit', response_type))

# summary of hits, misses, false alarms.
alld %>% 
  count(response_type) %>% 
  mutate(total = sum(n),
         p = (n/total*100) %>% round(0))

# proportion of response type by subject
hits_by_subj <- alld %>% 
  group_by(condition, subjectID) %>% 
  count(response_type) %>% 
  mutate(total = sum(n),
         p = (n/total*100) %>% round(0)) %>%
  select(subjectID, response_type, p) %>% 
  pivot_wider(names_from = response_type, values_from = p)

# exclude participants whose proportion of hits is below 90%
low_performance <- hits_by_subj %>% filter(hit < 90)
low_performance

d_perform <- alld %>% 
  filter(!subjectID %in% low_performance$subjectID)

# n of participants per condition (after data filtering)
d_perform %>% 
  count(condition, subjectID) %>% 
  count(condition)

# n of participants per condition (before data filtering)
alld %>% 
  count(condition, subjectID) %>% 
  count(condition)

#--------------------------------------------------------
# accuracy per condition
#--------------------------------------------------------

# accuracy (=proportion of hits) by subject
acc_subj <- d_perform %>%
  group_by(condition_label, subjectID, itemID, response_type) %>% 
  nest() %>% 
  group_by(condition_label, subjectID) %>% 
  select(-data) %>% 
  count(response_type) %>% 
  mutate(total = sum(n),
         p = n/total*100) %>%
  ungroup() %>% 
  filter(response_type %in% c('hit'))

# accuracy by condition
acc_cond <- acc_subj %>%
  group_by(condition_label) %>% 
  summarise(accuracy_mean = mean(p),
            accuracy_sd = sd(p)) %>% 
  ungroup()

acc_cond

#--------------------------------------------------------
# exclude all responses except hits
#--------------------------------------------------------

d_perform %>% 
  count(response_type) %>% 
  mutate(total = sum(n),
         p = (n/total*100) %>% round(0))

# dataset for analyses:
# only keep hits, i.e. deviant is present and participant reacts once;
# exclude misses, false_alarms and trials with overreactions.
d_hits <- d_perform %>% 
  filter(response_type == 'hit')

#--------------------------------------------------------
# exlude impossible/extreme RTs
#--------------------------------------------------------

# negative RTs are equivalent to false alarms (these have been removed previously);
# extremely low RTs too;
# define extreme (=outlier) as less/more than 2 SDs from mean;
# define composite threshold:
#   - by-subject threshold: subject mean +/- SDs
#   - by-position threshold: position mean +/- SDs.
# only exclude if RT exceeds both thresholds.
# otherwise, too many RTs may be excluded from initial positions.
# credit: rev4

# remove extreme
d_hits <- d_hits %>% 
  group_by(subjectID) %>% 
  mutate(extreme_lo = relRT < mean(relRT) - 2*sd(relRT),
         extreme_hi = relRT > mean(relRT) + 2*sd(relRT)) %>% 
  group_by(deviant) %>%
  mutate(extreme_lo_pos = relRT < mean(relRT) - 2*sd(relRT),
         extreme_hi_pos = relRT > mean(relRT) + 2*sd(relRT)) %>% 
  ungroup() %>% 
  filter((extreme_lo == F) | (extreme_lo_pos == F),
         (extreme_hi == F) | (extreme_hi_pos == F)) %>%
  select(-starts_with('extreme'))

#--------------------------------------------------------
# covariate: deviant probability
#--------------------------------------------------------

# following anonymous reviewer #4:
# introduce deviant probability as a covariate.
# given that 2/3s of the sequences contain a deviant,
# the overall probability of a sequence containing a deviant is of 2/3.
# further, the deviant probability for each position within a sequence is:
# 1st position = 2/3 * 1/8; 2nd = 2/3 * 1/7; 8th = 2/3 * 1/1.
# that is, the deviant probability increases along the sequence.
dev_prob <- tibble(deviant = seq(8),
                   deviant_probability = 2/3 * 1/seq(8, 1))

#--------------------------------------------------------
# final dataset, add transformed RT features
#--------------------------------------------------------

# scale() = center values, & normalise them (dividing by SD)
d <- d_hits %>% 
  mutate(logRT = log(relRT)) %>% 
  group_by(machine) %>% 
  # group_by(subjectID) %>%
  mutate(logRT_z = scale(logRT),
         relRT_z = scale(relRT)) %>% 
  ungroup() %>% 
  left_join(dev_prob)

# d %>% write_tsv('../data/seh2_processed_data.tsv')
# d <- read_tsv('../data/seh2_processed_data.tsv') %>% 
#   mutate(condition = as_factor(condition))

#--------------------------------------------------------
# descriptive stats for RTs by condition
#--------------------------------------------------------

# frontiers: Table S1
d %>%
  group_by(condition_label) %>% 
  summarise(RT_mean = mean(relRT) %>% round(0),
            RT_sd = sd(relRT) %>% round(0),
            participant_n = n_distinct(subjectID),
            sample_size = n()) %>% 
  ungroup() %>% 
  # summarise(across(where(is.numeric), sum)) %>% 
  left_join(acc_cond) %>% 
  mutate(accuracy_mean = round(accuracy_mean, 1)) %>% 
  ltx_tbl()

#--------------------------------------------------------
# distributions
#--------------------------------------------------------

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
  aes(x = relRT_z, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5) +
  theme_dist_nolegend + xlab('reaction time (z-normalised)')
dist_stdrt

dist_zrt <- ggplot(d) +
  aes(x = logRT_z, fill = condition_label, colour = condition_label) +
  geom_density(alpha = .5) +
  theme_dist_nolegend + xlab('reaction time (log-transformed, z-normalised)')
dist_zrt

plot_grid(dist_rt, dist_zrt, ncol = 1, labels = 'AUTO')
# ggsave_png('rt_distributions_2', height = 9)
# frontiers: Figure S1

plot_grid(dist_rt, dist_logrt, dist_zrt, ncol = 1, labels = 'AUTO')
# ggsave_png('rt_distributions_3', height = 9)

conf_summary(d, relRT, condition_label)
conf_summary(d, logRT, condition_label)
conf_summary(d, logRT_z, condition_label)

#--------------------------------------------------------
# preceding standards
#--------------------------------------------------------

# distribution of preceding standards, i.e. distance since preceding deviant
d %>%
  mutate(condition = as_factor(condition)) %>% 
  ggplot() +
  aes(x = preceding_stds,
      group = condition,
      fill = condition,
      colour = condition) +
  geom_density(alpha = .5) +
  theme_dist

d %>%
  mutate(deviant = as_factor(deviant)) %>%
  ggplot() +
  aes(x = deviant, y = preceding_stds) +
  geom_boxplot()

# d %>% 
#   group_by(deviant) %>% 
#   summarise(across(preceding_stds, list(mean=mean, sd=sd)))

d %>% 
  select(deviant, preceding_stds) %>% 
  group_by(deviant) %>% 
  skim()

d %>% 
  select(preceding_stds) %>% 
  skim()

# correlation between deviant position and preceding_stds
cor_test(d, deviant, preceding_stds)

# correlation between preceding_stds and RT
cor_test(d, preceding_stds, logRT_z)

#--------------------------------------------------------
# RTs by position of deviant
#--------------------------------------------------------

# axis labels
deviant_lab <- 'position of deviant within sequence'
zrt_lab <- 'reaction time (log-transformed, z-normalised)'
lrt_lab <- 'reaction time (log-transformed)'

## Summary by deviant position
dsum2 <- conf_summary(d, logRT_z, deviant)
ggplot(dsum2) + aes(x = factor(deviant), y = mean, group = 1) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = deviant_lab, y = zrt_lab)
# ggsave("dsum_all.pdf", scale = 2)

## Summary by deviant and condition
dsum <- conf_summary(d, logRT_z, deviant, condition_label) %>% 
  rename(condition = condition_label)
ggplot(dsum) + aes(x = factor(deviant), y = mean,
                   group = condition, color = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = deviant_lab, y = zrt_lab) +
  theme_dist
# ggsave("dsum_conditions.pdf", scale = 2)

#--------------------------------------------------------
# log-transformed, z-normalised RTs by deviant (with overlay)
#--------------------------------------------------------

# summary statistics (for plot overlay)
dsum.stat <- dsum %>% 
  group_by(condition) %>% 
  summarise(r = cor(deviant, mean) %>% round(2)) %>% 
  mutate(highlight = condition)

# no highlight
dsum %>% 
  ggplot() + aes(x = factor(deviant), y = mean, group = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = deviant_lab, y = zrt_lab) +
  stat_cor(label.x = 6, label.y = 1.1, cor.coef.name = 'r',
           aes(label = ..r.label..)) +
  theme_pubr(legend = 'none') +
  background_grid(major = 'y')

# frontiers: figure2
# ggsave_png(custom_name = "zrt_conditions", width = 9)

# with hihglight
dsum %>% 
  crossing(highlight = unique(.$condition)) %>%
  ggplot() +
  aes(x = factor(deviant), y = mean, group = condition,
      alpha = highlight == condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = deviant_lab, y = zrt_lab) +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.1, 1)) +
  geom_text(aes(x = 7, y = 1.1, label = paste0("italic(r) == ", r)), parse = T, 
            data = dsum.stat, inherit.aes = F) +
  background_grid(major = 'y') +
  theme(legend.position = 'none')

# ggsave_png(custom_name = "zrt_conditions_hi", width = 9)

#--------------------------------------------------------
# log-transformed RTs (not normalised by machine) by deviant (with overlay)
#--------------------------------------------------------

## summary by deviant and condition
log_sum <- conf_summary(d, logRT, deviant, condition_label) %>% 
  rename(condition = condition_label)

# summary statistics (for plot overlay)
log_sum_stat <- log_sum %>% 
  group_by(condition) %>% 
  summarise(r = cor(deviant, mean) %>% round(2)) %>% 
  mutate(highlight = condition)

# no highlight
ggplot(log_sum) + aes(x = factor(deviant), y = mean, group = condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_text(aes(x = 7, y = 6.5, label = paste0("italic(r) == ", r)),
            parse = T, data = log_sum_stat) +
  facet_wrap(~ condition) +
  labs(x = deviant_lab, y = lrt_lab) +
  background_grid(major = 'y')

# ggsave_png(custom_name = "logrt_conditions", width = 9)

# with hihglight
log_sum %>% 
  crossing(highlight = unique(.$condition)) %>%
  ggplot() +
  aes(x = factor(deviant), y = mean, group = condition,
      alpha = highlight == condition) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  facet_wrap(~ condition) +
  labs(x = deviant_lab, y = lrt_lab) +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.1, 1)) +
  geom_text(aes(x = 7, y = 6.5, label = paste0("italic(r) == ", r)),
            parse = T, data = log_sum_stat) +
  background_grid(major = 'y') +
  theme(legend.position="none")

# ggsave_png(custom_name = "logrt_conditions_hi", width = 9)

#--------------------------------------------------------
# regression models (logRT_z)
#--------------------------------------------------------

# saturated model
mm_sat <- lmer(logRT_z ~ condition *
                 (deviant + preceding_stds + deviant_probability) +
                 (1|subjectID), d, REML = F)
summary(mm_sat)

# model selection
step(mm_sat)

# step-selected final model
mm_final <- d %>% 
  mutate(condition = condition_label) %>%
  lmer(logRT_z ~ deviant + preceding_stds + deviant_probability +
                   condition + condition:deviant_probability +
                   (1|subjectID), data = ., REML = F)
summary(mm_final)
# frontiers: Table 2
print_model(mm_final)

#--------------------------------------------------------
# plot model predictions
#--------------------------------------------------------

# plot effect of deviant probability
mm_pred_probability <- ggpredict(mm_final, terms = c('deviant_probability', 'condition'))
predprob <- ggplot(mm_pred_probability) +
  aes(x = x, y = predicted) +
  geom_line(aes(colour = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .1) +
  xlab('deviant probability') +
  ylab('reaction time\n(log-transformed, z-normalised)') +
  theme_dist_nolegend
  
# ggsave_png('predictions_probability', height = 5)

# plot effect of deviant distance
mm_pred_distance <- ggpredict(mm_final, terms = c('preceding_stds', 'condition'))
preddis <- ggplot(mm_pred_distance) +
  aes(x = x, y = predicted) +
  geom_line(aes(colour = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .1) +
  xlab('deviant distance') +
  ylab('reaction time\n(log-transformed, z-normalised)') +
  theme_dist_nolegend

# ggsave_png('predictions_distance', height = 5)

# plot effect of deviant position
mm_pred_position <- ggpredict(mm_final, terms = c('deviant', 'condition'))
predpos <- ggplot(mm_pred_position) +
  aes(x = x, y = predicted) +
  geom_line(aes(colour = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .1) +
  xlab('deviant position') +
  ylab('reaction time\n(log-transformed, z-normalised)') +
  theme_dist

# ggsave_png('predictions_position', height = 5)

plot_grid(predpos, preddis, predprob, ncol = 1, labels = 'AUTO')
# ggsave_png('model_predictions', height = 12)
# frontiers: Figure S2

#--------------------------------------------------------
# regression models (logRT)
#--------------------------------------------------------

# saturated model
mm_sat_log <- lmer(logRT ~ condition *
                 (deviant + preceding_stds + deviant_probability) +
                 (1|subjectID), d, REML = F)
summary(mm_sat_log)

# model selection
step(mm_sat_log)

# final model
mm_final_log <- lmer(logRT ~ deviant + preceding_stds + deviant_probability +
                       condition + condition:deviant_probability +
                   (1|subjectID), d, REML = F)
summary(mm_final_log)
print_model(mm_final_log)

#--------------------------------------------------------
# regression model: accuracy, only condition as predictor
#--------------------------------------------------------

# binary response accuracy as dependent variable (correct vs incorrect)
d_acc <- d_perform %>% 
  mutate(accurate = response_type == 'hit') %>% 
  left_join(dev_prob)

# logistic regression model
mm_sat_acc <- glmer(accurate ~ condition + (1|subjectID),
                    family = 'binomial', d_acc)
summary(mm_sat_acc)

# compare model with condition as predictor
# to null model without condition
mm_sat_acc_null <- glmer(accurate ~ (1|subjectID),
                    family = 'binomial', d_acc)
anova(mm_sat_acc, mm_sat_acc_null)

#--------------------------------------------------------
# regression model: accuracy
#--------------------------------------------------------

# logistic regression model.
mm_sat_acc_2 <- glmer(accurate ~ condition *
                      (deviant + preceding_stds + deviant_probability) +
                      (1|subjectID),
                    family = 'binomial', d_acc)

# accuracy_step <- cAIC4::stepcAIC(mm_sat_acc_2,
#                      direction = "both", trace = T, data = d_acc)

# stepwise selection (back+forward) of accuracy model
accuracy_step <- fitLMER.fnc(model = mm_sat_acc_2,
                 method = 'AIC',
                 backfit.on = 't',
                 ran.effects = list(ran.intercepts = "subjectID"))

summary(accuracy_step)
