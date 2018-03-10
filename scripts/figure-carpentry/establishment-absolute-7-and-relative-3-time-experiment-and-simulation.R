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

## Absolute time first
# Simulation data from generation 7
sims_establish_absolute <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "extant_prop" & time_type == "absolute") %>%
  select(propagule_number, env, t_equals_7) %>% # Here is where we get the simulation results from GEN 7
  group_by(propagule_number) %>%
  summarize(establish_prop = mean(t_equals_7)) %>% # Here is where we get the simulation results from GEN 7
  as.data.frame()


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
min_y <- min(c(plogis(establishment_posthoc$asymp.LCL)), sims_establish_absolute$establish_prop)
xlim <- range(establishment_xvals) + c(-0.5, 0.5)

pdf("figures/establishment-probability-absolute-7-and-relative-3-time-experiment-and-simulations.pdf", height = 3.14961, width = 2 * 3.14961)
# postscript("figures/establishment-probability-absolute-7-and-relative-3-time-experiment-and-simulations.eps", height = 3.14961, width = 2 * 3.14961)
par(mar = c(1.5, 1, 0.5, 0), family = "Helvetica", mgp = c(2.25, 1, 0), mfrow = c(1, 2), oma = c(1.5, 2, 0, 2))

plot(x = establishment_xvals, y = plogis(establishment_posthoc$lsmean), 
     ylim=c(min_y - 0.03, 1.05), 
     xlim = xlim, 
     las=1, 
     pch=19, 
     xaxt = "n", 
     yaxt = "n",
     xlab = NA, 
     ylab = NA)

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
points(x = establishment_xvals, y = sims_establish_absolute$establish_prop,
       pch = 17)

# legend("bottomright",
#        inset = c(0, -0.025),
#        legend = c("microcosm", "simulation"),
#        pch = c(19, 17),
#        bty = "n")


# Now relative time points
prob_establishment_model <- glmer(extant_3_after ~ number*environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa"))

# Simulation data from 3 generations after final introduction event:

#### Code to separate out environment treatments ####
# sims_establish_relative <- sims_results %>%
#   mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
#   filter(response == "extant_prop" & time_type == "relative") %>%
#   select(propagule_number, env, t_equals_3) %>% # Here is where we average establishment across environment treatment
#   rename(establish_prop = t_equals_3) %>% 
#   arrange(propagule_number, env) %>% 
#   as.data.frame()
# 
# 
# establishment_results <- lsmeans::lsmeans(prob_establishment_model, pairwise ~ environment*number, adjust="none")
# establishment_posthoc <- summary(establishment_results$lsmeans)
# 
# # Get compact letter display for significant differences across levels
# establishment_sig_letters <- lsmeans::cld(establishment_results, Letters = letters, sort = FALSE, adjust = "none")$.group
# # Remove white space from significance letters object
# establishment_sig_letters <- gsub(establishment_sig_letters, pattern = " ", replacement = "")
# 
# offset_xvals <- function(x, offset) {
#   c(sapply(x, FUN = function(j) c(j - offset, j + offset)))
# }
# 
# establishment_xvals <- offset_xvals(1:4, 0.2)
# 
# # Minimum y value of the plot is the lowest possible value that can be plotted (lowest lower CI or 
# # lowest result from simulations)
# xlim <- range(establishment_xvals) + c(-0.5, 0.5)
# 
# plot(x = establishment_xvals, y = plogis(establishment_posthoc$lsmean), 
#      ylim=c(min_y - 0.03, 1.05), 
#      xlim = xlim, 
#      las=1, 
#      pch=c(1, 19), 
#      xaxt = "n", 
#      yaxt = "n",
#      xlab = NA, 
#      ylab = NA)
# 

# # axis(side = 4,
# #      at = c(0.6, 0.7, 0.8, 0.9, 1.0),
# #      las = 1,
# #      hadj = 0.75)
# 
# segments(x0 = establishment_xvals, 
#          y0 = plogis(establishment_posthoc$asymp.LCL), 
#          y1 = plogis(establishment_posthoc$asymp.UCL), 
#          lwd = 2)
# 
# text(x = establishment_xvals, y = 1.07, labels = establishment_sig_letters, pos = 1)
# 
# Add simulation results to the plot
# 
# sims_xvals <- offset_xvals(establishment_xvals, offset = 0.05)
# Add simulation results to the plot
# # segments(x0 = sims_xvals[seq(1, length(sims_xvals), by = 2)], x1 = sims_xvals[seq(2, length(sims_xvals), by = 2)], y0 = sims_popSize$mean_N, lwd = 4)
# 
# points(x = establishment_xvals[c(1, 3, 5, 7)], y = sims_establish_relative$establish_prop[c(1, 3, 5, 7)],
#        pch = 2)
# 
# points(x = establishment_xvals[c(2, 4, 6, 8)], y = sims_establish_relative$establish_prop[c(2, 4, 6, 8)],
#        pch = 17)
# 
# 
# legend("bottomleft",
#        inset = c(0, -0.025),
#        legend = c("fluctuating", "stable"), 
#        pch = c(0, 15),
#        bty = "n")
# 
# legend("bottomright",
#        inset = c(0, -0.025),
#        legend = c("microcosm", "simulation"),
#        pch = c(19, 17),
#        bty = "n")

#### Code to average across environment treatments ####

sims_establish_relative <- sims_results %>%
  mutate(propagule_number = substr(intro_regime, start = nchar(intro_regime), stop = nchar(intro_regime))) %>%
  filter(response == "extant_prop" & time_type == "relative") %>%
  select(propagule_number, env, t_equals_3) %>% # Here is where we average establishment across environment treatment
  rename(establish_prop = t_equals_3) %>% 
  group_by(propagule_number) %>%
  summarize(establish_prop = mean(establish_prop)) %>% 
  as.data.frame()

establishment_results <- lsmeans::lsmeans(prob_establishment_model, pairwise ~ number, adjust="none")
establishment_posthoc <- summary(establishment_results$lsmeans)

# Get compact letter display for significant differences across levels
establishment_sig_letters <- lsmeans::cld(establishment_results, Letters = letters, sort = FALSE, adjust = "none")$.group
# Remove white space from significance letters object
establishment_sig_letters <- gsub(establishment_sig_letters, pattern = " ", replacement = "")

establishment_xvals <- 1:4

# Minimum y value of the plot is the lowest possible value that can be plotted (lowest lower CI or 
# lowest result from simulations)
xlim <- range(establishment_xvals) + c(-0.5, 0.5)

plot(x = establishment_xvals, y = plogis(establishment_posthoc$lsmean), 
     ylim=c(min_y - 0.03, 1.05), 
     xlim = xlim, 
     las = 1, 
     pch = 19, 
     xaxt = "n", 
     yaxt = "n",
     xlab = NA, 
     ylab = NA)

axis(side = 1, 
     at = 1:4, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE,
     line = -0.5)

segments(x0 = establishment_xvals, 
         y0 = plogis(establishment_posthoc$asymp.LCL), 
         y1 = plogis(establishment_posthoc$asymp.UCL), 
         lwd = 2)

text(x = establishment_xvals, y = 1.07, labels = establishment_sig_letters, pos = 1)

mtext(side = 1,
      text = "Introduction regime",
      line = 0.5,
      outer = TRUE)

mtext(side = 2,
      text = "Establishment probability",
      line = 1.15,
      outer = TRUE)

axis(side = 1,
     at = 1:4,
     labels = c("20x1","10x2","5x4","4x5"),
     tick = FALSE,
     line = -0.5)

axis(side = 4,
     at = c(0.6, 0.7, 0.8, 0.9, 1.0),
     las = 1,
     hadj = 0.25)

sims_xvals <- establishment_xvals

# Add simulation results to the plot
points(x = establishment_xvals, y = sims_establish_relative$establish_prop,
       pch = 17)

legend("bottomright",
       inset = c(0, -0.025),
       legend = c("microcosm", "simulation"),
       pch = c(19, 17),
       bty = "n")


dev.off()
