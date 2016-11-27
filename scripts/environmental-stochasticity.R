### Title: Exploring environmental stochasticity in Tribolium system within Parsing Propagule Pressure experiment
###
### Author: Michael Koontz
###
### Date Created: 20140715
### Last Updated: 20160113
###
### Purpose: Determine whether fluctuating environment treatment corresponded to a greater var(log(lambda))

rm(list=ls())

if (!require("lme4"))
{install.packages("lme4"); library(lme4)}
if (!require("lsmeans"))
{install.packages("lsmeans"); library(lsmeans)}

setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography")
source("Clean Data/generate-tidy-data.R")

beetles <- read.csv("Clean-Data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("Clean-Data/attributes.csv")

b <- tidy.beetles(beetles=beetles) # Convert long form data to short form (one generation per column)
data_col <- 2:(ncol(b$Nt) - 1) # The data columns (i.e. all time steps, but not including ID column for the tidy data.frames)

lambdas <- b$Ntp1[, data_col] / b$Nt[, data_col] # Calculate year to year population growth rates
stoch_tot <- data.frame(attributes, stoch_tot=apply(lambdas, 1, function(x) var(log(x)))) # Add total stochasticity (defined as var(log(lambdas))) to the attributes data frame

#------------
# Fitted models testing for effect of environment group (stable vs fluctuating) on actual total stochasticity
#------------

fm1 <- lmer(stoch_tot ~ environment + (1 + environment | block), data=stoch_tot)
fm2 <- lmer(stoch_tot ~ environment + (1 | block), data=stoch_tot)

anova(fm1, fm2) # Drop environment * block interaction

# Significantly different contrast in total stochasticity between stable and fluctuating environment (estimate=0.04707 p=0.0201)
lsmeans::lsmeans(fm2, pairwise ~ environment)$contrasts


