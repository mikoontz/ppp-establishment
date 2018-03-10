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

#### Read samples data to put the equilibrium abundance value on the abundance plot ####
samps <- read.csv(here::here("data/NBBg-samples/NBBg-samples-combined.csv"))
equilibrium_popSize <- log(samps$R0) / samps$alpha

# Do not average over the environment treatments because there is a significant effect of environment on mean abundance
sims_popSize_absolute <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "mean_N" & time_type == "absolute") %>%
  select(propagule_number, env, mean_N = t_equals_9) %>% # Here is where we get the simulation results from GEN 9
  arrange(propagule_number, env)


popSize_model <- glmer(N8plus1 ~ number*environment + (1 | block), data=b.trim, family=poisson, control=glmerControl(optimizer="bobyqa"))

popSize_results <- lsmeans::lsmeans(popSize_model, pairwise ~ environment + number, adjust="none")
popSize_posthoc <- summary(popSize_results$lsmeans)

popSize_sig_letters <- lsmeans::cld(popSize_results, Letters = letters, sort = FALSE, adjust = "none")$.group
popSize_sig_letters <- gsub(popSize_sig_letters, pattern = " ", replacement = "")

offset_xvals <- function(x, offset) {
  c(sapply(x, FUN = function(j) c(j - offset, j + offset)))
}
popSize_xvals <- offset_xvals(1:4, 0.25)

min_y <- min(c(exp(popSize_posthoc$asymp.LCL), sims_popSize_absolute$mean_N))
max_y <- max(exp(popSize_posthoc$asymp.UCL))
xlim <- range(popSize_xvals) + c(-0.4, 0.4)

pdf("figures/population-size-absolute-9-and-relative-5-time-experiment-and-simulations.pdf", height = 3.14961, width = 2 * 3.14961)
# postscript("figures/population-size-absolute-9-and-relative-5-time-experiment-and-simulations.eps", height = 3.14961, width = 2 * 3.14961)
par(mar = c(1.5, 1, 0.5, 0), family = "Helvetica", mgp = c(2.25, 1, 0), mfrow = c(1, 2), oma = c(1.5, 2, 0, 2))

plot(x = popSize_xvals, y = exp(popSize_posthoc$lsmean),
     type = "n", # Set the plot up, but do not print lines yet
     ylim = c(-10, max_y + 10), 
     xlim = xlim, 
     las = 1, 
     pch = 1, 
     xaxt = "n",
     yaxt = "n",
     xlab = NA, 
     ylab = NA)

segments(x0 = popSize_xvals, 
         y0 = exp(popSize_posthoc$asymp.LCL), 
         y1 = exp(popSize_posthoc$asymp.UCL), 
         lwd = 2)

points(x = popSize_xvals[c(1, 3, 5, 7)], y = exp(popSize_posthoc$lsmean[c(1, 3, 5, 7)]),
       pch = 1)

points(x = popSize_xvals[c(2, 4, 6, 8)], y = exp(popSize_posthoc$lsmean[c(2, 4, 6, 8)]),
       pch = 19)

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

text(x = popSize_xvals, y = max_y + 13, labels = popSize_sig_letters, pos = 1)

sims_xvals <- offset_xvals(popSize_xvals, offset = 0.05)
#### Add simulation results to the plot ####

points(x = popSize_xvals[c(1, 3, 5, 7)], y = sims_popSize_absolute$mean_N[c(1, 3, 5, 7)],
       pch = 2)

points(x = popSize_xvals[c(2, 4, 6, 8)], y = sims_popSize_absolute$mean_N[c(2, 4, 6, 8)],
       pch = 17)

segments(x0 = xlim[1], x1 = xlim[2], y0 = mean(equilibrium_popSize), y1 = mean(equilibrium_popSize), lty = "dashed", lwd = 2)

#
# Now the population size at 5 generations after the final introduction event
#

sims_popSize_relative <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "mean_N" & time_type == "relative") %>%
  select(propagule_number, env, mean_N = t_equals_5) %>% # Here is where we get the simulation results from GEN 9
  arrange(propagule_number, env)

popSize_model <- glmer(N_5_after ~ number*environment + (1 | block), data=b.trim, family=poisson, control=glmerControl(optimizer="bobyqa"))

popSize_results <- lsmeans::lsmeans(popSize_model, pairwise ~ environment + number, adjust="none")
popSize_posthoc <- summary(popSize_results$lsmeans)

popSize_sig_letters <- lsmeans::cld(popSize_results, Letters = letters, sort = FALSE, adjust = "none")$.group
popSize_sig_letters <- gsub(popSize_sig_letters, pattern = " ", replacement = "")

offset_xvals <- function(x, offset) {
  c(sapply(x, FUN = function(j) c(j - offset, j + offset)))
}
popSize_xvals <- offset_xvals(1:4, 0.1)

min_y <- min(c(exp(popSize_posthoc$asymp.LCL), sims_popSize_relative$mean_N))
max_y <- max(exp(popSize_posthoc$asymp.UCL))
xlim <- range(popSize_xvals) + c(-0.4, 0.4)

plot(x = popSize_xvals, y = exp(popSize_posthoc$lsmean),
     type = "n", # Set the plot up, but do not print lines yet
     ylim = c(-10, max_y + 5), 
     xlim = xlim, 
     las = 1, 
     pch = 1, 
     xaxt = "n",
     yaxt = "n",
     xlab = NA, 
     ylab = NA)

mtext(side = 1,
      text = "Introduction regime",
      line = 0.5,
      outer = TRUE)

mtext(side = 2,
      text = "Population size",
      line = 1.15,
      outer = TRUE)

segments(x0 = popSize_xvals, 
         y0 = exp(popSize_posthoc$asymp.LCL), 
         y1 = exp(popSize_posthoc$asymp.UCL), 
         lwd = 2)

points(x = popSize_xvals[c(1, 3, 5, 7)], y = exp(popSize_posthoc$lsmean[c(1, 3, 5, 7)]),
       pch = 1)

points(x = popSize_xvals[c(2, 4, 6, 8)], y = exp(popSize_posthoc$lsmean[c(2, 4, 6, 8)]),
       pch = 19)

axis(side = 1, 
     at = 1:4, 
     labels = c("20x1", "10x2", "5x4", "4x5"), 
     tick = FALSE,
     line = -0.5)

axis(side = 4, las = 1)
# axis(side = 4,
#      at = seq(0, 60, by = 20),
#      las = 1,
#      hadj = 0.25)

# axis(side = 4,
#      at = seq(0, 60, by = 20),
#      las = 1,
#      hadj = 0.75)

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

points(x = popSize_xvals[c(1, 3, 5, 7)], y = sims_popSize_relative$mean_N[c(1, 3, 5, 7)],
       pch = 2)

points(x = popSize_xvals[c(2, 4, 6, 8)], y = sims_popSize_relative$mean_N[c(2, 4, 6, 8)],
       pch = 17)

segments(x0 = xlim[1], x1 = xlim[2], y0 = mean(equilibrium_popSize), y1 = mean(equilibrium_popSize), lty = "dashed", lwd = 2)

dev.off()
