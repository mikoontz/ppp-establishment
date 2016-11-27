### Script to run Bayesian analysis with NBBG model of Tribolium population growth
###
### Author: Michael Koontz
###
### Date Created: 20141103
### Next to Last Update: 20141224
### Last Updated: 20150313
###
### The purpose of this code is to get data into a format that can be analyzed using a hierarchical Bayesian model, analyze those data using MCMC methods, and then present posterior distributions of the parameters of interest
###
###
###

# Clear working environment if necessary
# rm(list=ls())

source("scripts/NBBg-mcmc.R")
# source("NBBG population growth function.R")

data <- read.csv("data/initial-density-dependence.csv")

# Data collected for initial density dependence side experiment didn't use a unique identifier for each population, so set that now
data$ID <- 1:nrow(data)

# Set up data for analysis. Ensures column titles from actual data match the column titles that the MCMC code uses.
colnames(data)[colnames(data)=="census"] <- "Ntplus1"
colnames(data)[colnames(data)=="N0"] <- "migrants"
data$residents <- 0
data$Nt <- data$migrants + data$residents

# Setup MCMC parameters
n.mcmc <- 10000
inits <- list(c(R0=2, kE=22, kD=10, alpha=0.005),
              c(R0=5, kE=2, kD=25, alpha=0.01),
              c(R0=0.5, kE=10, kD=10, alpha=0.0005))

priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)
tune <- c(R0=0.1, kE=10, kD=1.5, alpha=0.001, RE=1)

#### Run MCMC algorithm ####
mcmc_output <- lapply(inits, FUN = function(x) NBBG.mcmc(data, priors.shape, priors.scale, x, tune, n.mcmc))
mcmc.list

# Percent of proposals accepted
num_accepted <- lapply(mcmc_output, FUN = function(x) x[["accept"]])
names(num_accepted) <- paste0("chain_", 1:length(num_accepted))

percent_accepted <- lapply(num_accepted, FUN = function(x) x / n.mcmc * 100)
percent_accepted

# Pull out samples of parameters
samples <- lapply(mcmc_output, FUN = function(x) x[["samps"]])
samples_list <- mcmc.list(samples)

# Pull out key parameters
key <- lapply(samples_list, function(x) x[, c("R0", "alpha", "kE", "kD")])

#### Trace plots of unburned samples for key parameters ####

plot(key)

#### Effective number of parameters for chains ####
effectiveSize(key)

#### Convergence diagnostics for MCMC chains ####

converge <- gelman.diag(x = key) # Should be at or *very* near 1 to indicate convergence.
gelman.plot(x = key) # Depicts shrink factor through time which can help ensure that convergence test above isn't a false positive. Must show decline through time, rather than a value of 1 (just by chance) throughout sampling

#### Burn in 20% of samples
burned <- lapply(samples_list, FUN = function(x) burn.in(x, n.mcmc*0.2))
burned_key <- lapply(burned, FUN = function(x) x[, c("R0", "alpha", "kE", "kD")])

#### Trace plots of burned in samples ####

plot(burned_key)

#### Plot the best fit curve for the estimated expectation of Ntplus1 given Nt
R0 <- mean(burned[[1]][, "R0"])
alpha <- mean(burned[[1]][, "alpha"])
x <- 1:200
mu_Ntp1 <- x * R0 * exp(-alpha * x)

plot(data$Nt, data$Ntplus1, pch=19)
lines(x = x, y = mu_Ntp1, col="red")

#### Write the samples data file for my data ####
# samples <- data.frame(R0=burned[["R0"]], kE=burned[["kE"]], kD=burned[["kD"]], alpha=burned[["alpha"]])
# write.csv(samples, 'NBBg samples.csv', row.names=FALSE)

#### Figuring out the modes and 95% credible intervals ####
report <- lapply(burned_key, FUN = summary)

# results using MJK Parsing Propagule Pressure data
# parameter         mode      lowerCI      upperCI
# 1        R0  1.134698685  1.004747623   1.26071576
# 2        kE 16.497861042 10.726579596 198.95984985
# 3        kD  1.517677349  0.830885042 100.68925330
# 4     alpha  0.008640162  0.007291477   0.01022171
# 5        RE  1.069107474  0.719188210   1.68939729

# results using data from Melbourne and Hastings 2008 Nature paper

# parameter         mean         mode      lowerCI      upperCI
# 1        R0  2.603156146  2.575061700  2.244819461  3.024520886
# 2        kE 29.488198170 21.381219238 10.870754031 75.050024319
# 3        kD  1.089430833  0.732032613  0.407613921  2.816244763
# 4     alpha  0.003717655  0.003707539  0.003398071  0.004049347
# 5        RE  2.602043286  2.457251521  1.644039932  3.887305169



###
###
### Randomly sample 10 population ID's and plot the predicted and known number of migrant females
# 
# get.mode <- function(x)
# {
# 	names(sort(-table(x)))[1]
# }
# 
# par(mfrow=c(2, 5))
# for (i in 1:length(ids))
# {
# 	plot(table(mcmc[["F.migrants"]][ids[i], ]), ylab="Frequency", main=paste("Expected Female Migrants:", test$Nt[ids[i]]/2, "\nEstimated Female Migrants:", get.mode(mcmc[["F.migrants"]][ids[i], ]), "\nActual Female Migrants:", test$F.migrants[ids[i]] ) )
# 	
# 	abline(v=test$F.migrants[ids[i]], col="blue", lwd=4)	
# 	abline(v=get.mode(mcmc[["F.migrants"]][ids[i], ]), col="red", lwd=4)	
# 
# }
