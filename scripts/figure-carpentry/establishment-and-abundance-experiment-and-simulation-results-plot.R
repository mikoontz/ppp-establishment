#### Title: experiment-and-simulation-results-plot

# Clear environment if necessary
# rm(list=ls())

library(lme4)
library(lsmeans)
library(multcompView)
library(tidyr)
library(dplyr)

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
sims_popSize <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "mean_N" & time_type == "absolute") %>%
  select(propagule_number, env, mean_N = t_equals_7) %>% # Here is where we get the simulation results from GEN 7
  arrange(propagule_number, env)

#### Read samples data to put the equilibrium abundance value on the abundance plot ####
samps <- read.csv("data/NBBg-samples/NBBg-samples-combined.csv")

equilibrium_popSize <- log(samps$R0) / samps$alpha

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
establishment_xvals <- 1:length(establishment_posthoc$lsmean)

# Minimum y value of the plot is the lowest possible value that can be plotted (lowest lower CI or 
# lowest result from simulations)
min_y <- min(c(plogis(establishment_posthoc$asymp.LCL)), sims_establish$establish_prop)
xlim <- range(establishment_xvals) + c(-0.5, 0.5)

pdf("figures/establishment-probability-experiment-and-simulations.pdf", height = 3.14961, width = 3.14961)
par(mar = c(3.5, 3.5, 0.5, 0.25), family = "Helvetica", mgp = c(2.25, 1, 0))

plot(x = establishment_xvals, y = plogis(establishment_posthoc$lsmean), 
     ylim=c(min_y - 0.03, 1.05), 
     xlim = xlim, 
     las=1, 
     pch=19, 
     xaxt = "n", 
     yaxt = "n",
     xlab = NA, 
     ylab = "Establishment probability")

mtext(side = 1,
      text = "Introduction regime",
      line = 2)

axis(side = 1, 
     at = establishment_xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE,
     line = -0.5)

axis(side = 2,
     at = c(0.6, 0.7, 0.8, 0.9, 1.0),
     las = 1,
     hadj = 0.75)

segments(x0 = establishment_xvals, 
         y0 = plogis(establishment_posthoc$asymp.LCL), 
         y1 = plogis(establishment_posthoc$asymp.UCL), 
         lwd = 2)

text(x = 1:4, y = 1.07, labels = establishment_sig_letters, pos = 1)

#### Add simulation results to the plot ####
points(x = establishment_xvals, y = sims_establish$establish_prop,
       pch = 17)

legend("bottomright",
       inset = c(0, -0.025),
       legend = c("microcosm", "simulation"),
       pch = c(19, 17),
       bty = "n")

dev.off()



#### Population size ####
#### Population size data from experiment ####
popSize_model <- glmer(N6plus1 ~ number*environment + (1 | block), data=b.trim, family=poisson, control=glmerControl(optimizer="bobyqa"))

popSize_results <- lsmeans::lsmeans(popSize_model, pairwise ~ environment + number, adjust="none")
popSize_posthoc <- summary(popSize_results$lsmeans)

popSize_sig_letters <- lsmeans::cld(popSize_results, Letters = letters, sort = FALSE, adjust = "none")$.group
popSize_sig_letters <- gsub(popSize_sig_letters, pattern = " ", replacement = "")

offset_xvals <- function(x, offset) {
  c(sapply(x, FUN = function(j) c(j - offset, j + offset)))
}
popSize_xvals <- offset_xvals(1:4, 0.1)

min_y <- min(c(exp(popSize_posthoc$asymp.LCL), sims_popSize$mean_N))
max_y <- max(exp(popSize_posthoc$asymp.UCL))
xlim <- range(popSize_xvals) + c(-0.4, 0.4)

pdf("figures/population-size-experiment-and-simulations.pdf", height = 3.14961, width = 3.14961)
par(mar = c(3.5, 3.25, 0.5, 0.25), family = "Helvetica", mgp = c(2.0, 1, 0))

plot(x = popSize_xvals, y = exp(popSize_posthoc$lsmean),
     type = "n", # Set the plot up, but do not print lines yet
     ylim = c(-10, max_y + 4), 
     xlim = xlim, 
     las = 1, 
     pch = 1, 
     xaxt = "n",
     yaxt = "n",
     xlab = NA, 
     ylab = "Population size")

mtext(side = 1,
      text = "Introduction regime",
      line = 2)

segments(x0 = popSize_xvals, 
         y0 = exp(popSize_posthoc$asymp.LCL), 
         y1 = exp(popSize_posthoc$asymp.UCL), 
         lwd = 2)

points(x = popSize_xvals[c(1, 3, 5, 7)], y = exp(popSize_posthoc$lsmean[c(1, 3, 5, 7)]),
       pch = 1,
       type = "b")

points(x = popSize_xvals[c(2, 4, 6, 8)], y = exp(popSize_posthoc$lsmean[c(2, 4, 6, 8)]),
     pch = 19,
     type = "b")

axis(side = 1, 
     at = 1:4, 
     labels = c("20x1", "10x2", "5x4", "4x5"), 
     tick = FALSE,
     line = -0.5)

axis(side = 2,
     at = seq(0, 60, by = 20),
     las = 1,
     hadj = 0.75)

legend("bottomleft",
       inset = c(0, -0.025),
       legend = c("fluctuating", "stable"), 
       pch = c(0, 15),
       bty = "n")

legend("bottom",
       inset = c(0, 0.15),
       legend = "equilibrium size", 
       lty = 2,
       lwd = 2,
       bty = "n")

legend("bottomright",
       inset = c(0, -0.025),
       legend = c("microcosm", "simulation"), 
       pch = c(19, 17),
       bty = "n")

text(x = popSize_xvals, y = max_y + 9, labels = popSize_sig_letters, pos = 1)

sims_xvals <- offset_xvals(popSize_xvals, offset = 0.05)
#### Add simulation results to the plot ####
# segments(x0 = sims_xvals[seq(1, length(sims_xvals), by = 2)], x1 = sims_xvals[seq(2, length(sims_xvals), by = 2)], y0 = sims_popSize$mean_N, lwd = 4)

points(x = popSize_xvals[c(1, 3, 5, 7)], y = sims_popSize$mean_N[c(1, 3, 5, 7)],
       pch = 2)

points(x = popSize_xvals[c(2, 4, 6, 8)], y = sims_popSize$mean_N[c(2, 4, 6, 8)],
       pch = 17)

segments(x0 = xlim[1], x1 = xlim[2], y0 = mean(equilibrium_popSize), y1 = mean(equilibrium_popSize), lty = "dashed", lwd = 2)

dev.off()
# Polygon of 95% credible interval around equilibrium population size Would need adjustment to y-axis.
# lwr <- quantile(equilibrium_popSize, probs = 0.025)
# upr <- quantile(equilibrium_popSize, probs = 0.975)
# polygon(x = c(0.25, 8.25, 8.25, 0.25), y = c(lwr, lwr, upr, upr), col = adjustcolor("black", alpha.f = 0.2) )
