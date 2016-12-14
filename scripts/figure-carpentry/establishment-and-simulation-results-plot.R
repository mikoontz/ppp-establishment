#### Title: experiment-and-simulation-results-plot

# Clear environment if necessary
# rm(list=ls())

library(lme4)
library(lsmeans)
library(multcompView)
library(tidyr)

#### Read and prepare data from the experiment ####
b <- read.csv("data/clean-establishment-data.csv")
b$number <- as.factor(b$number)
b$block <- as.factor(b$block)
b$gap <- as.factor(b$gap)
b.trim <- b[b$gap == "FALSE", ] # I think introduction gaps may be enough of an issue to just always exclude them

#### Read and prepare data from the simulations ####
sims_results <- read.csv("data/simulations/simulation_stats_tidy.csv")

sims_establish <- sims_results %>%
  filter(response == "extant_prop" & time_type == "absolute") %>%
  select(intro_regime, env, t_equals_7) %>%
  group_by(intro_regime) %>%
  summarize(establish_prop = mean(t_equals_7)) %>%
  arrange(intro_regime) %>%
  as.data.frame()

sims_establish$x_pos <- c(2, 1, 4, 3)
sims_establish <- sims_establish[order(sims_establish$x_pos), ]

#### Establishment ####
#### Establishment data from experiment ####
# Using establishment assessment at generation 7
prob_establishment <- glmer(extant7 ~ number + environment + number:environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) 

establishment_results <- lsmeans::lsmeans(prob_establishment, pairwise ~ number, adjust="none")
establishment_posthoc <- summary(establishment_results$lsmeans)

establishment_sig_letters <- lsmeans::cld(establishment_results, Letters = letters, sort = FALSE, adjust = "none")$.group
xvals <- 1:length(establishment_posthoc$lsmean)

# Minimum y value of the plot is the lowest possible value that can be plotted (lowest lower CI or 
# lowest result from simulations)
min_y <- min(c(plogis(establishment_posthoc$asymp.LCL)), sims_establish$establish_prop)

plot(x = xvals, y = plogis(establishment_posthoc$lsmean), 
     ylim=c(min_y, 1.02), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las=1, 
     pch=19, 
     xaxt="n", 
     xlab="Introduction regime", 
     ylab="Establishment probability", 
     bty="L")

axis(side=1, 
     at = xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE)

arrows(x0 = xvals, 
       y0 = plogis(establishment_posthoc$asymp.LCL), 
       y1 = plogis(establishment_posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2)

text(x = 1:4, y = 1.02, labels = establishment_sig_letters)

#### Add simulation results to the plot ####

segments(x0 = 1:4 - 0.25, x1 = 1:4 + 0.25, y0 = sims_establish$establish_prop, lwd = 4)


