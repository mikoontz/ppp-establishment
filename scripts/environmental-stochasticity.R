### Title: Exploring environmental stochasticity in Tribolium system within Parsing Propagule Pressure experiment
###
### Author: Michael Koontz
###
### Date Created: 20140715
### Penultimate Update: 20160113
### Last Updated: 20161127 (all further updates reflected in version history)

### Purpose: Determine whether fluctuating environment treatment corresponded to a greater var(log(lambda))

# Clear environment if necessary
# rm(list=ls())

library(lme4)
library(lsmeans)

source("scripts/data-carpentry/generate-tidy-data.R")

beetles <- read.csv("data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("data/attributes.csv")

b <- tidy.beetles(beetles = beetles) # Convert long form data to short form (one generation per column)
data_col <- 2:(ncol(b$Nt) - 1) # The data columns (i.e. all time steps, but not including ID column for the tidy data.frames)

lambdas <- b$Ntp1[, data_col] / b$Nt[, data_col] # Calculate year to year population growth rates
stoch_tot <- data.frame(attributes, stoch_tot=apply(lambdas, 1, function(x) var(log(x)))) # Add total stochasticity (defined as var(log(lambdas))) to the attributes data frame

#------------
# Fitted models testing for effect of environment group (stable vs fluctuating) on actual total stochasticity
#------------
# Models are fit using only populations that never went temporarily or permanently extinct (n=667)

fm1 <- lmer(stoch_tot ~ environment + (environment | block), data=stoch_tot)
fm2 <- lmer(stoch_tot ~ environment + (1 | block), data=stoch_tot)

anova(fm1, fm2) # Drop random slope effect of environment

# Significantly different contrast in total stochasticity between stable and fluctuating environment (estimate=0.04707 p=0.0201)
results <- lsmeans::lsmeans(fm2, pairwise ~ environment)
(contr <- results$contrasts)

# 95% confidence interval on the contrast.
confint(contr)

#### Plot density of total stochasticity for each environment treatment ####
plot(density(stoch_tot$stoch_tot[stoch_tot$environment == "stable"], na.rm= TRUE))
lines(density(stoch_tot$stoch_tot[stoch_tot$environment == "fluctuating"], na.rm= TRUE), col = "red")
abline(v = summary(results)$lsmeans$lsmean[1], lwd = 2, col = "red")
abline(v = summary(results)$lsmeans$lsmean[2], lwd = 2, col = "black")

#### Compare to literature ####
# Values taken from Saether and Engin (2002)

saether.1 <- c(0.152, 0.026, 0.065, 0.375, 0.034, 0.015, 0.012, 0.090, 0.017, 0.294, 0.010, 0.120, 0.003, 0.024, 0.008, 0.012, 0.027, 0.069, 0.118, 0.011, 0.055, 0.134, 0.068, 0.004)
saether.2 <- c(0.041, 0.276, 0.014, 0.065, 0.411, 0.437, 0.114, 0.072, 0.005, 0.104, 0.014)

median(c(saether.1, saether.2))
