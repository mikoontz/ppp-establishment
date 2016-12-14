### Title: NBBg simulations script

# Clear environment if necessary
# rm(list=ls())

source("scripts/simulations/NBBg-simulation-functions.R")
source("scripts/data-carpentry/generate-establishment-responses-for-analysis.R")

#### Set up parameters ####

# Samples from posterior distribution of R0, alpha, kE, and kD parameters
samps <- read.csv('data/NBBg-samples/NBBg-samples-combined.csv')
reps <- 500000 # Number of reps per introduction scenario
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
                                        gens.censused = rep(Tf - 1, nrow(pp) * reps), # 1 fewer census than total time steps
                                        col_names = paste0("extant_prop_after_", time_points))
#### Extant proportion at fixed time points for fluctuating environment ####
# Population still present at generations 5, 6, 7, 8, and 9 for fluctuating environment treatment
extant_prop_after_x_var <- simulation_stats(time_points = time_points,
                                            FUN = function(x, ...) !extinct.after.x(x, ...),
                                            N = N_var[, -1], # Just pass the equivalent of the Ntp1 dataframe (i.e. no initial introduction column)
                                            intro.regime = intro.regime,
                                            gap = rep(FALSE, nrow(pp) * reps),
                                            gens.censused = rep(Tf - 1, nrow(pp) * reps), # 1 fewer census than total time steps
                                            col_names = paste0("extant_prop_after_", time_points, "_var"))
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
colnames(mean_N_after_x_var) <- paste0("N_after", time_points, "_var")

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
                                        gens.censused = rep(Tf - 1, nrow(pp) * reps), # 1 fewer census than total time steps
                                        col_names = paste0("extant_prop_", time_points, "after"))
#### Establishment proportion at relative time points for fluctuating environment ####
# Population still present 1, 2, 3, 4, and 5 generations after final introduction event for fluctuating environment treatment
extant_prop_x_after_var <- simulation_stats(time_points = time_points,
                                            FUN = function(x, ...) !extinct.x.after.intro(x, ...),
                                            N = N_var[, -1], # Just pass the equivalent of the Ntp1 dataframe (i.e. no initial introduction column)
                                            intro.regime = intro.regime,
                                            gap = rep(FALSE, nrow(pp) * reps),
                                            gens.censused = rep(Tf - 1, nrow(pp) * reps), # 1 fewer census than total time steps
                                            col_names = paste0("extant_prop_", time_points, "after_var"))
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
                                       col_names = paste0("N_", time_points, "after_var"))

results <- data.frame(
  extant_prop_after_x,
  extant_prop_after_x_var,
  mean_N_after_x,
  mean_N_after_x_var,
  extant_prop_x_after,
  extant_prop_x_after_var,
  mean_N_x_after,
  mean_N_x_after_var)

results

#### Tidy the results ####
sims_results_tidy <- function(sims_results, response, time_type, env) {
  
  # Defense. Make sure column names are appropriate.
  if (!(response %in% c("extant_prop", "mean_N"))) {
    stop("Response must be extant_prop or mean_N")
  }
  
  if (!(time_type %in% c("absolute", "relative"))) {
    stop("Time type must be absolute or relative")
  } 
  
  if (!(env %in% c("stable", "fluctuating"))) {
    stop("Environment must be stable or fluctuating")
  }
  
  # Find the numbered columns and rename them something generic.
  char_start <- regexpr(pattern = "[0-9]+", colnames(sims_results))
  char_stop <- attributes(char_start)$match.length + char_start - 1
  generation <- substr(colnames(sims_results), start = char_start, stop = char_stop)
  colnames(sims_results) <- paste0("t_equals_", generation)
  
  # Add some NA columns to flesh out the whole time series. Absolute time types will need missing values for 1st through 4th time points
  if (time_type == "absolute") {
    tmp <- as.data.frame(matrix(NA, nrow = nrow(sims_results), ncol = max(pp$number) - 1)) # 4-row data frame of NAs
    colnames(tmp) <- paste0("t_equals_", (1:(max(pp$number) - 1))) # Named t_equals_1 through t_equals_4
    sims_results <- data.frame(tmp, sims_results) # These new columns go first
  }
  
  # Relative time types will need missing values for 46th through 50th time points
  if (time_type == "relative") {
    tmp <- as.data.frame(matrix(NA, nrow = nrow(sims_results), ncol = max(pp$number) - 1)) # 4-row data frame of NAs
    # Named t_equals_47 through t_equals_50
    colnames(tmp) <- paste0("t_equals_", (max(as.numeric(generation) + 1):(max(as.numeric(generation) + max(pp$number) - 1))))
    sims_results <- data.frame(sims_results, tmp) # These columns go at the end of the time series.
  }
  
  # Fill in variable columns with their appropriate data
  sims_results$intro_regime <- apply(pp, 1, function(x) paste(rev(x), collapse = "x"))
  sims_results$response <- response
  sims_results$time_type <- time_type
  sims_results$env <- env
  
  return(sims_results)
}

#### Tidy the individual results data frames ####
extant_prop_after_x_tidy <- sims_results_tidy(extant_prop_after_x, response = "extant_prop", time_type = "absolute", env = "stable")
extant_prop_after_x_var_tidy <- sims_results_tidy(extant_prop_after_x_var, response = "extant_prop", time_type = "absolute", env = "fluctuating")

mean_N_after_x_tidy <- sims_results_tidy(mean_N_after_x, response = "mean_N", time_type = "absolute", env = "stable")
mean_N_after_x_var_tidy <- sims_results_tidy(mean_N_after_x_var, response = "mean_N", time_type = "absolute", env = "fluctuating")

extant_prop_x_after_tidy <- sims_results_tidy(extant_prop_x_after, response = "extant_prop", time_type = "relative", env = "stable")
extant_prop_x_after_var_tidy <- sims_results_tidy(extant_prop_x_after_var, response = "extant_prop", time_type = "relative", env = "fluctuating")

mean_N_x_after_tidy <- sims_results_tidy(mean_N_x_after, response = "mean_N", time_type = "relative", env = "stable")
mean_N_x_after_var_tidy <- sims_results_tidy(mean_N_x_after_var, response = "mean_N", time_type = "relative", env = "fluctuating")

results_tidy <-
  rbind(extant_prop_after_x_tidy,
        extant_prop_after_x_var_tidy,
        mean_N_after_x_tidy,
        mean_N_after_x_var_tidy,
        extant_prop_x_after_tidy,
        extant_prop_x_after_var_tidy,
        mean_N_x_after_tidy,
        mean_N_x_after_var_tidy)

#### Rearrange columns so that variable columns are first ####
results_tidy <- results_tidy[, c((ncol(results_tidy) - 4) + 1:4, 1:(ncol(results_tidy) - 4))]

#### Write data to files ####
# write.csv(x = results, file = "data/simulations/simulation_stats_raw.csv", row.names = FALSE)
# write.csv(x = results_tidy, file = "data/simulations/simulation_stats_tidy.csv", row.names = FALSE)

# for (i in unique(intro.regime)) {
#   file_tag <- paste(rev(pp[pp$number == i, ]), collapse = "x")
#   write.csv(N[intro.regime == i, ], file = paste0("data/simulations/N_", file_tag, "_regime.csv"), row.names = FALSE)
# }
# 
# for (i in unique(intro.regime)) {
#   file_tag <- paste(rev(pp[pp$number == i, ]), collapse = "x")
#   write.csv(N_var[intro.regime == i, ], file = paste0("data/simulations/N_var_", file_tag, "_regime.csv"), row.names = FALSE)
# }
# 
# for (i in unique(intro.regime)) {
#   file_tag <- paste(rev(pp[pp$number == i, ]), collapse = "x")
#   write.csv(migrants[intro.regime == i, ], file = paste0("data/simulations/migrants_", file_tag, "_regime.csv"), row.names = FALSE)
# }

