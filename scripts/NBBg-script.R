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

source("scripts/NBBg-NIMBLE.R")
# source("NBBG population growth function.R")

#### Read and prepare data ####
data <- read.csv("data/initial-density-dependence.csv")
colnames(data) <- c("ID", "block", "env", "setup_order", "Nt", "Ntplus1", "check")
data <- subset(data, select = c(ID, Nt, Ntplus1))
data$residents <- 0
data$migrants <- data$Nt
head(data)

# Setup MCMC parameters
niter <- 50000
nchains <- 3
nburnin <- 10000

priors <- list(R0 = c(shape = 2.6, scale = 1),
               alpha = c(shape = 0.0037, scale = 1),
               kE = c(shape = 17.6, scale = 1),
               kD = c(shape = 1.07, scale = 1))


#### Run MCMC algorithm ####
key <- nbbgNIMBLE_run(data = data, priors = priors, nchains = nchains, niter = niter, nburnin = nburnin)

#### Trace plots for key parameters ####

plot(key[[1]])

#### Effective number of parameters for chains ####
Neff <- lapply(key, effectiveSize)
names(Neff) <- paste0("chain_", 1:length(key))
Neff

# Summed over all chains.
effectiveSize(key)


#### Convergence diagnostics for MCMC chains ####

converge <- gelman.diag(x = key) # Should be at or *very* near 1 to indicate convergence.
converge

gelman.plot(x = key) # Depicts shrink factor through time which can help ensure that convergence test above isn't a false positive. Must show decline through time, rather than a value of 1 (just by chance) throughout sampling
# Shows clear decline and improvement through the sampling process.



#### Write the samples data file for my data ####
# Chains converge, so combine all samples into a single matrix-like object
combined_samples <- do.call(rbind, key)

# Write the combined samples to a .csv
# write.csv(x = combined_samples,
#           file = 'data/NBBg-samples/NBBg-samples-combined.csv',
#           row.names = FALSE)
# 
# # Write each chain to its own .csv
# lapply(1:length(key),
#        FUN = function(i)
#          write.csv(x = key[[i]],
#                    file = paste0("data/NBBg-samples/NBBg-samples-chain", i),
#                    row.names = FALSE))

#### Figuring out the modes and 95% credible intervals ####
report <- lapply(key, FUN = summary)
report
