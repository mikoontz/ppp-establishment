### Title: NBBg simulations script

# Clear environment if necessary
# rm(list=ls())

source("scripts/simulations/NBBg-simulation-functions.R")
source("scripts/data-carpentry/generate-establishment-responses-for-analysis.R")

#### Set up parameters ####

# Samples from posterior distribution of R0, alpha, kE, and kD parameters
samps <- read.csv('data/NBBg-samples/NBBg-samples-combined.csv')
reps <- 10 # Number of reps per introduction scenario
Tf <- 10 # Number of time steps, INCLUDING initial introduction
total.to.introduce <- 20
p <- 0.5 # Probability of being female
pp <- get.propagule.pressure(total.to.introduce, Tf, excludeMaxNum = FALSE)

#### Simulate dynamics in stable environment ####
# Baseline dynamics under 'stable' environment with density dependence. All NBBg parameters are the estimated values from the Bayesian data analysis. Model error (i.e. uncertainty in the parameter values) are propagated through the simulation since each modeled population in each time step draws a parameter set (maintaining correlation structure) from the MCMC samples

N <- simNBBg(samps, reps, Tf, total.to.introduce, p = 0.5, save_objects = FALSE, excludeMaxNum = FALSE)

#### Simulate dynamics in more variable environment ####
# Multiply kE samples by 100/121 to increase variance of Re by 20% and standard deviation of Re by 10%
var_samps <- samps
var_samps$kE <- samps$kE * (100/121) # Add ~10% to standard deviation of R
# HIvar_samps <- samps
# HIvar_samps$kE <- samps$kE * (1/4) # Add ~100% to standard deviation of R

N_var <- simNBBg(samps=var_samps, reps, Tf, total.to.introduce, p=0.5, save_objects=FALSE, data_descriptor="_var", excludeMaxNum = FALSE)
# N_HIvar <- simNBBg(samps=HIvar_samps, reps, Tf, total.to.introduce, p=0.5, save_objects=TRUE, data_descriptor="_HIvar")

#### Simulation summaries ####
intro.regime <- get.intro.regime(N = N, propagule.pressure = pp)

#### Establshment proportions and population abundance at fixed time points ####
# I.e. is the population still present at generation x?
# We can look starting at generation (max number of introductions) until the final generation in the experiment
time_points <- max(pp$number):(Tf - 1)

#### Extant proportion at fixed time points for stable environment ####
# Population still present at generations 5, 6, 7, 8, and 9 for stable environment treatment
extant_prop_after_x <- simulation_stats(time_points = time_points,
                                        FUN = function(x, ...) !extinct.after.x(x, ...),
                                        N = N[, -1], # Just pass the equivalent of the Ntp1 dataframe (i.e. no initial introduction column)
                                        intro.regime = intro.regime,
                                        gap = rep(FALSE, nrow(pp) * reps),
                                        gens.censused = rep(Tf, nrow(pp) * reps),
                                        col_names = paste0("extant_prop_after_", time_points))
#### Extant proportion at fixed time points for fluctuating environment ####
# Population still present at generations 5, 6, 7, 8, and 9 for fluctuating environment treatment
extant_prop_after_x_var <- simulation_stats(time_points = time_points,
                                            FUN = function(x, ...) !extinct.after.x(x, ...),
                                            N = N_var[, -1], # Just pass the equivalent of the Ntp1 dataframe (i.e. no initial introduction column)
                                            intro.regime = intro.regime,
                                            gap = rep(FALSE, nrow(pp) * reps),
                                            gens.censused = rep(Tf, nrow(pp) * reps),
                                            col_names = paste0("extant_prop_after_", time_points))
#### Population sizes at fixed time points for stable populations ####
# First, convert all Ntp1 censuses of 0 to NA so they don't get included in averages
N_for_abundance <- N
N_for_abundance[, -1] <- sapply(N_for_abundance[, -1], function(j) replace(j, which(j == 0), NA))

mean_N_after_x <- apply(X = N_for_abundance[, time_points + 1], 
                        MARGIN = 2, 
                        FUN = function(x) tapply(x, INDEX = intro.regime, FUN = mean, na.rm = TRUE))
colnames(mean_N_after_x) <- paste0("N_after", time_points)
#### Population sizes at fixed time points for fluctuating populations ####
N_for_abundance_var <- N_var
N_for_abundance_var[, -1] <- sapply(N_for_abundance_var[, -1], function(j) replace(j, which(j == 0), NA))

mean_N_after_x_var <- apply(X = N_for_abundance_var[, time_points + 1], 
                        MARGIN = 2, 
                        FUN = function(x) tapply(x, INDEX = intro.regime, FUN = mean, na.rm = TRUE))
colnames(mean_N_after_x_var) <- paste0("N_after", time_points)

##### Establishment proportions and population abundance at relative time points #####
# (x number of generations after final introduction event)#### Establshment proportions at relative time points
# We can assess this for all populations from the 1st generation after introductions end to the (Tf - maximum propagule number treatment)th generation after introductions end
time_points <- 1:(Tf - pp$number[nrow(pp)])
#### Establishment proportion at relative time points for stable environment ####
# Population still present 1, 2, 3, 4, and 5 generations after final introduction event for stable environment treatment
extant_prop_x_after <- simulation_stats(time_points = time_points,
                                        FUN = function(x, ...) !extinct.x.after.intro(x, ...),
                                        N = N[, -1], # Just pass the equivalent of the Ntp1 dataframe (i.e. no initial introduction column)
                                        intro.regime = intro.regime,
                                        gap = rep(FALSE, nrow(pp) * reps),
                                        gens.censused = rep(Tf, nrow(pp) * reps),
                                        col_names = paste0("extant_prop_", time_points, "after"))
#### Establishment proportion at relative time points for fluctuating environment ####
# Population still present 1, 2, 3, 4, and 5 generations after final introduction event for fluctuating environment treatment
extant_prop_x_after_var <- simulation_stats(time_points = time_points,
                                            FUN = function(x, ...) !extinct.x.after.intro(x, ...),
                                            N = N_var[, -1], # Just pass the equivalent of the Ntp1 dataframe (i.e. no initial introduction column)
                                            intro.regime = intro.regime,
                                            gap = rep(FALSE, nrow(pp) * reps),
                                            gens.censused = rep(Tf, nrow(pp) * reps),
                                            col_names = paste0("extant_prop_", time_points, "after"))
#### Population sizes at relative time points for stable populations ####
mean_N_x_after <- simulation_stats(time_points = time_points,
                                   FUN = N_x.after.intro,
                                   Ntp1 = N[, -1], 
                                   intro.regime = intro.regime,
                                   col_names = paste0("N_", time_points, "after"))
#### Population sizes at relative time points for fluctuating populations ####
mean_N_x_after_var <- simulation_stats(time_points = time_points,
                                       FUN = N_x.after.intro,
                                       Ntp1 = N_var[, -1], 
                                       intro.regime = intro.regime,
                                       col_names = paste0("N_", time_points, "after"))
