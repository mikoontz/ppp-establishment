### Script to run Bayesian analysis with NBBG model of Tribolium population growth
###
### Author: Michael Koontz
###
### Date Created: 20141103
### Last Updated: 20161122
###
### The purpose of this code is to validate the MCMC algorithm used to calculate NBBg parameters on the Parsing Propagule Pressure initial density dependence side experiment. We use known values of the parameters and simulate an NBBg model, then fit our hierarchical model to those simulated data to see if we can reobtain the known parameter values. We then fit our hierarchical model to data from Melbourne and Hastings (2008) to see if we can reobtain the parameter values they estimated from their analysis.

# Clear environment of defined variables if necessary
# rm(list = ls())

source("scripts/simulations/NBBg-population-dynamics-function.R")
source("scripts/NBBg-NIMBLE.R")

#### Simulated data to validate the MCMC algorithm ####
### Simulate data
# residents <- rep(0,2600)
# migrants <- data.frame(migrants0 = c(rep(4,500),
#                                      rep(5, 500), 
#                                      rep(10, 500), 
#                                      rep(20, 500), 
#                                      rep(50, 200), 
#                                      rep(100,200), 
#                                      rep(200, 100), 
#                                      rep(1000, 100)),
#                        migrants1 = 0)


past.residents <- rep(0,250)
migrants <- data.frame(migrants0 = c(rep(seq(2,100), times=2), 
                                     rep(seq(125, 500, by=25), times=2), 
                                     rep(seq(550, 1000, by=50), times=2)),
                       migrants1 = 0)


N <- data.frame(Nt = migrants$migrants0 + past.residents, Ntplus1 = 0)

# # Triple the amount of data simulated
# migrants <- rbind(migrants, migrants, migrants)
# N <- rbind(N, N, N)

# Known parameter values 
R0 <- 1.8
kE <- 11
kD <- 2
alpha <- 0.006

samps <- data.frame(R0, alpha, kE, kD)

pop_trajectory <- NBBg(Tf = 2, 
                  total.reps = nrow(N),
                  N = N,
                  past.residents = 0, 
                  migrants = migrants, 
                  p = 0.5, 
                  samps = samps)

test_data <- data.frame(ID = 1:nrow(pop_trajectory), 
                        pop_trajectory, 
                        migrants = migrants$migrants0, 
                        residents = past.residents)
head(test_data)

#### Plot the simulated data ####
x <- 1:1000
plot(test_data$Nt, test_data$Ntplus1, pch=16, cex = 0.5)
lines(x=x, y=x * R0 * exp(-alpha*x), col="red")

#### Set up MCMC parameters
niter <- 50000
nchains <- 4
nburnin <- 10000

#### Vague priors ####
priors <- list(R0 = c(shape = 0.001, scale = 1000),
               alpha = c(shape = 0.001, scale = 1000),
               kE = c(shape = 0.001, scale = 1000),
               kD = c(shape = 0.001, scale = 1000))

#### Weakly regularizing priors from Melbourne & Hastings (2008) ####
# These aid in model convergence

priors <- list(R0 = c(shape = 2.6, scale = 1),
               alpha = c(shape = 0.0037, scale = 1),
               kE = c(shape = 17.6, scale = 1),
               kD = c(shape = 1.07, scale = 1))

#### Run MCMC algorithm ####
key <- nbbgNIMBLE_run(data = test_data, priors = priors, nchains = nchains, niter = niter, nburnin = nburnin)

#### Trace plots of unburned samples for key parameters ####

plot(key[[1]])

#### Effective number of parameters for chains ####
Neff <- lapply(key, effectiveSize)
names(Neff) <- paste0("chain_", 1:length(key))
Neff

#### Summed over all chains. ####
effectiveSize(key)

#### Convergence diagnostics for MCMC chains ####

converge <- gelman.diag(x = key) # Should be at or *very* near 1 to indicate convergence.
converge


gelman.plot(x = key) # Depicts shrink factor through time which can help ensure that convergence test above isn't a false positive. Must show decline through time, rather than a value of 1 (just by chance) throughout sampling
# Shows clear decline and improvement through the sampling process.

lapply(key, FUN = summary)
# True values
samps



#### Data from Melbourne & Hastings (2008) ####

melbourne <- read.csv("data/melbourne_ricker_data.csv")
melbourne$At <- round(melbourne$At)
melbourne$Atp1 <- round(melbourne$Atp1)
melbourne$residents <- 0
melbourne$migrants <- melbourne$Nt
colnames(melbourne)[colnames(melbourne)=="Atp1"] <- "Ntplus1"
colnames(melbourne)[colnames(melbourne)=="At"] <- "migrants"

melbourne$residents <- 0
melbourne$Nt <- melbourne$migrants + melbourne$residents

plot(melbourne$Nt, melbourne$Ntplus1, pch=19)

### Define MCMC arguments
niter <- 50000
nchains <- 4
nburnin <- 10000

#### Weakly regularizing priors ####
# These aid in model convergence

priors <- list(R0 = c(shape = 1, scale = 1),
               alpha = c(shape = 0.005, scale = 1),
               kE = c(shape = 15, scale = 1),
               kD = c(shape = 15, scale = 1))

#### Run MCMC algorithm ####
key <- nbbgNIMBLE_run(data = melbourne, priors = priors, nchains = nchains, niter = niter, nburnin = nburnin)

#### Trace plots of unburned samples for key parameters ####

plot(key[[1]])

#### Effective number of parameters for chains ####
Neff <- lapply(key, effectiveSize)
names(Neff) <- paste0("chain_", 1:length(key))
Neff

#### Summed over all chains. ####
effectiveSize(key)

#### Convergence diagnostics for MCMC chains ####

converge <- gelman.diag(x = key) # Should be at or *very* near 1 to indicate convergence.
converge


gelman.plot(x = key) # Depicts shrink factor through time which can help ensure that convergence test above isn't a false positive. Must show decline through time, rather than a value of 1 (just by chance) throughout sampling
# Shows clear decline and improvement through the sampling process.

lapply(key, FUN = summary)
