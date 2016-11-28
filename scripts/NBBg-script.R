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
n.mcmc <- 20000
inits <- list(c(R0=2, kE=22, kD=10, alpha=0.005),
              c(R0=5, kE=2, kD=25, alpha=0.015),
              c(R0=0.5, kE=35, kD=1, alpha=0.0005))

priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)
tune <- c(R0=0.1, kE=10, kD=1.5, alpha=0.001, RE=1)

#### Run MCMC algorithm ####
mcmc_output <- lapply(inits, FUN = function(x) NBBG.mcmc(data, priors.shape, priors.scale, x, tune, n.mcmc))

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

plot(key[[1]])

#### Effective number of parameters for chains ####
Neff <- lapply(key, effectiveSize)
names(Neff) <- paste0("chain_", 1:length(key))
Neff


#### Convergence diagnostics for MCMC chains ####

converge <- gelman.diag(x = key) # Should be at or *very* near 1 to indicate convergence.
converge
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# R0          1.01       1.04
# alpha       1.01       1.02
# kE          1.01       1.03
# kD          1.00       1.01
# 
# Multivariate psrf
# 
# 1.01

gelman.plot(x = key) # Depicts shrink factor through time which can help ensure that convergence test above isn't a false positive. Must show decline through time, rather than a value of 1 (just by chance) throughout sampling
# Shows clear decline and improvement through the sampling process.

#### Burn in 5% of samples
burned <- lapply(samples_list, FUN = function(x) burn.in(x, n.mcmc*0.05))
burned_key <- as.mcmc.list(lapply(burned, FUN = function(x) x[, c("R0", "alpha", "kE", "kD")]))

#### Trace plots of burned in samples ####

plot(burned_key[[1]])

#### Plot the best fit curve for the estimated expectation of Ntplus1 given Nt
R0 <- mean(burned[[1]][, "R0"])
alpha <- mean(burned[[1]][, "alpha"])
x <- 1:200
mu_Ntp1 <- x * R0 * exp(-alpha * x)

plot(data$Nt, data$Ntplus1, pch=19)
lines(x = x, y = mu_Ntp1, col="red")

#### Write the samples data file for my data ####
# Chains converge, so combine all samples into a single matrix-like object
combined_samples <- do.call(rbind, burned_key)

# Write the combined samples to a .csv
# write.csv(x = combined_samples,
#           file = 'data/NBBg-samples/NBBg-samples-combined.csv',
#           row.names = FALSE)

# Write each chain to its own .csv
# lapply(1:length(burned_key),
#        FUN = function(i)
#          write.csv(x = burned_key[[i]],
#                    file = paste0("data/NBBg-samples/NBBg-samples-chain", i),
#                    row.names = FALSE))

#### Figuring out the modes and 95% credible intervals ####
report <- lapply(burned_key, FUN = summary)
summary(burned_key)
# Iterations = 1:19000
# Thinning interval = 1 
# Number of chains = 3 
# Sample size per chain = 19000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean        SD  Naive SE Time-series SE
# R0     1.13110 0.0654174 2.740e-04      3.092e-03
# alpha  0.00876 0.0007555 3.165e-06      3.343e-05
# kE    19.66002 3.8511432 1.613e-02      9.229e-02
# kD     2.27827 0.9959276 4.171e-03      1.770e-02
# 
# 2. Quantiles for each variable:
#   
#   2.5%       25%       50%      75%    97.5%
# R0     1.009968  1.086325  1.128694  1.17402  1.26611
# alpha  0.007303  0.008243  0.008754  0.00927  0.01027
# kE    13.135438 16.911755 19.289537 22.02158 28.00732
# kD     1.015638  1.586667  2.070256  2.72190  4.81263

