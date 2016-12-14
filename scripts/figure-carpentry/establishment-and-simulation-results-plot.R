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
sims_results <- read.csv("data/simulations/simulation_stats_tidy.csv", stringsAsFactors = FALSE)

sims_establish <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "extant_prop" & time_type == "absolute") %>%
  select(propagule_number, env, t_equals_7) %>% # Here is where we get the simulation results from GEN 7
  group_by(propagule_number) %>%
  summarize(establish_prop = mean(t_equals_7)) %>% # Here is where we get the simulation results from GEN 7
  as.data.frame()

# Do not average over the environment treatments because there is a significant effect of environment on mean abundance
sims_abundance <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "mean_N" & time_type == "absolute") %>%
  select(propagule_number, env, mean_N = t_equals_7) %>% # Here is where we get the simulation results from GEN 7
  arrange(propagule_number, env)

#### Read samples data to put the equilibrium abundance value on the abundance plot ####



#### Establishment ####
#### Establishment data from experiment ####
# Using establishment assessment at generation 7
prob_establishment_model <- glmer(extant7 ~ number + environment + number:environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) 

establishment_results <- lsmeans::lsmeans(prob_establishment_model, pairwise ~ number, adjust="none")
establishment_posthoc <- summary(establishment_results$lsmeans)

# Get compact letter display for significant differences across levels
establishment_sig_letters <- lsmeans::cld(establishment_results, Letters = letters, sort = FALSE, adjust = "none")$.group
# Remove white space from significance letters object
establishment_sig_letters <- gsub(establishment_sig_letters, pattern = " ", replacement = "")

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


#### Abundance ####
#### Abundance data from experiment ####
abundance_model <- glmer(N6plus1 ~ number*environment + (1 | block), data=b.trim, family=poisson, control=glmerControl(optimizer="bobyqa"))

abundance_results <- lsmeans::lsmeans(abundance_model, pairwise ~ environment + number, adjust="none")
abundance_posthoc <- summary(abundance_results$lsmeans)

abundance_sig_letters <- lsmeans::cld(abundance_results, Letters = letters, sort = FALSE, adjust = "none")$.group
abundance_sig_letters <- gsub(abundance_sig_letters, pattern = " ", replacement = "")
xvals <- 1:length(abundance_posthoc$lsmean)

min_y <- min(c(exp(abundance_posthoc$asymp.LCL), sims_abundance$mean_N))
max_y <- max(exp(abundance_posthoc$asymp.UCL))

par(mar=c(4,4,3,1))

plot(x = xvals, y = exp(abundance_posthoc$lsmean), 
     ylim = c(min_y - 5, max_y + 4), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las = 1, 
     pch = c(1,19), 
     xaxt = "n", 
     xlab = "Introduction regime", 
     ylab = "Population abundance", 
     bty = "L")

axis(side = 1, 
     at = xvals[c(1,3,5,7)] + 0.5, 
     labels=c("20x1", "10x2", "5x4", "4x5"), 
     tick=FALSE)

arrows(x0 = xvals, 
       y0 = exp(abundance_posthoc$asymp.LCL), 
       y1 = exp(abundance_posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2)

legend("right", 
       legend = c("fluctuating", "stable"), 
       pch = c(1, 19), 
       bty = "n")

text(x = xvals, y = max_y + 2, labels = abundance_sig_letters)

#### Add simulation results to the plot ####
segments(x0 = 1:8 - 0.25, x1 = 1:8 + 0.25, y0 = sims_abundance$mean_N, lwd = 4)
