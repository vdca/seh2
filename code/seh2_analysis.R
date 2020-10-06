#--------------------------------------------------------
# Globals
#--------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)
theme_set(theme_bw()) # ggplot theme
library(lmerTest)
library(xtable)

figdir <- '../plots/'

#------------------------------------------------------------
# helper functions
#------------------------------------------------------------

# summarySE()
# partially based on: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence 
## interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

summarySE <- function(data, measurevar, groupvars, conf.interval=.95) {
  
  # Mean by groups
  datac <- data %>% 
    mutate_(measurevar = measurevar) %>% 
    group_by_(.dots = groupvars) %>%
    summarise(mean = mean(measurevar, na.rm = T), sd = sd(measurevar, na.rm = T), N = n()) %>% 
    ungroup()
  datac[,measurevar] <- datac$mean
  
  # Calculate standard error of the mean
  datac$se <- datac$sd / sqrt(datac$N)
  
  # Confidence interval multiplier for standard error
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# save plots with cairo

ggsave.alt <- function(plot.lst = list(last_plot()), custom.name = "lastplot",
                       height = fh, width = fw, figscale = 1, ...) {
  if (is.null(names(plot.lst))) names(plot.lst) <- custom.name
  plot.names <- paste0(names(plot.lst), ".pdf")
  pmap(list(plot.names, plot.lst), ggsave,
       device=cairo_pdf, scale = figscale,
       height = height, width = width, ...)
}

#--------------------------------------------------------
# Load data
#--------------------------------------------------------

d <- read_csv('../data/seh2_response_data.csv') %>% 
  mutate(condition = as.factor(condition))

#--------------------------------------------------------
# Distributions
#--------------------------------------------------------

# Conditions were conducted on different machines:
# Conditions 1 and 2: Dell XPS M1330, Ubuntu 12.04.
# Condition 3: Lenovo T440s, Ubuntu 15.10.
# RTs are consistently shorter for condition 3.

ggplot(d) + aes(x = relRT, fill = condition, colour = condition) +
  geom_density(alpha = .5)

ggplot(d) + aes(x = zRT, fill = condition, colour = condition) +
  geom_density(alpha = .5)

summarySE(d, "relRT", groupvars = "condition")
summarySE(d, "zRT", "condition")

#--------------------------------------------------------
# By probe
#--------------------------------------------------------

## Summary by probe
dsum2 <- summarySE(d, "zRT", c("probe"))
ggplot(dsum2) + aes(x = factor(probe), y = mean, group = 1) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = "Position of probe within sequence",
       y = "Normalised reaction time (z)")
# ggsave("dsum_all.pdf", scale = 2)

## Summary by probe and condition
dsum <- summarySE(d, "zRT", c("probe", "condition"))
ggplot(dsum) + aes(x = factor(probe), y = mean, group = condition, color = condition) +
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
  theme(legend.position="none")
# ggsave(paste0(figdir, "dsum_conditions2.pdf"), width = 10, height = 6)

# summary statistics (for plot overlay)
dsum.stat <- dsum %>% 
  mutate(condition = paste("Condition", condition)) %>% 
  group_by(condition) %>% 
  summarise(r = cor(probe, zRT) %>% round(2)) %>% 
  mutate(highlight = condition)

# with hihglight
dsum %>% 
  mutate(condition = paste("Condition", condition)) %>% 
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
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  geom_text(aes(x = 7, y = 1.1, label = paste0("italic(r) == ", r)), parse = T, 
            data = dsum.stat)
# ggsave.alt(custom.name = "dsum_conditions3", width = 7.5, height = 5)

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

print(xtable(summary(mm3)$coefficients, digits = 4),
      scalebox = '.8', include.rownames = T, booktabs = T)

print(xtable(anova(mm3), digits = 4),
      scalebox = '.8', include.rownames = T, booktabs = T)
