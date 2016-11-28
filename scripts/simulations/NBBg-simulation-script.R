### Title: NBBg simulations script

# Clear environment if necessary
# rm(list=ls())

source("scripts/simulations/NBBg-simulation-functions.R")

#-----------
# Set up parameters
#-----------

# Samples from posterior distribution of R0, alpha, kE, and kD parameters
samps <- read.csv('data/NBBg-samples.csv')
reps <- 500 # Number of reps per introduction scenario
Tf <- 10 # Number of time steps, INCLUDING initial introduction
total.to.introduce <- 20
p <- 0.5 # Probability of being female
pp <- get.propagule.pressure(total.to.introduce, Tf)

#-------------
# Baseline dynamics under 'stable' environment with density dependence. All NBBg parameters are the estimated values from the Bayesian data analysis. Model error (i.e. uncertainty in the parameter values) are propagated through the simulation since each modeled population in each time step draws a parameter set (maintaining correlation structure) from the MCMC samples
#-------------

N <- simNBBg(samps, reps, Tf, total.to.introduce, p=0.5, save_objects=FALSE)

#------
# more variable environment
#------

# Multiply kE samples by 100/121 to increase variance of Re by 20% and standard deviation of Re by 10%
var_samps <- samps
var_samps$kE <- samps$kE * (100/121) # Add ~10% to standard deviation of R
# HIvar_samps <- samps
# HIvar_samps$kE <- samps$kE * (1/4) # Add ~100% to standard deviation of R

N_var <- simNBBg(samps=var_samps, reps, Tf, total.to.introduce, p=0.5, save_objects=FALSE, data_descriptor="_var")
# N_HIvar <- simNBBg(samps=HIvar_samps, reps, Tf, total.to.introduce, p=0.5, save_objects=TRUE, data_descriptor="_HIvar")

#--------
# Get establshment proportions through time (is the population still present, x generations after the final introduction event)
#--------
# Proportion of extant populations through time after the final introduction event
# We can assess this for all populations from the 1st generation after introductions end to the final recorded generation (Tf - maximum propagule number treatment)
time_points <- 1:(Tf-pp$number[nrow(pp)])

extant <- sapply(time_points, function(x) 1-extinct.x.after.intro(x=x, N=N, propagule.pressure=pp)$percent.extinct)

e_var <- sapply(time_points, function(x) 1-extinct.x.after.intro(x=x, N=N_var, propagule.pressure=pp)$percent.extinct)

#-------------
# Population sizes at different time steps
#-------------

make_table <- function(e1=5, e2=50, a1=5)
{
  intro.regime <- get.intro.regime(N=N, propagule.pressure = pp)
  
  # Get all population abundances a1 generations after final intro
  N_abund <- N_x.after.intro(x=a1, N=N, intro.regime=intro.regime)
  N_var_abund <- N_x.after.intro(x=a1, N=N_var, intro.regime=intro.regime)
  N_noDD_abund <- N_x.after.intro(x=a1, N=N_noDD, intro.regime=intro.regime)
  N_noDDvar_abund <- N_x.after.intro(x=a1, N=N_noDDvar, intro.regime=intro.regime)
  
  # Get mean of population abundances a1 generations after final intro
  N_mu <- aggregate(abundance ~ intro.regime, data=N_abund, FUN=mean)
  N_var_mu <- aggregate(abundance ~ intro.regime, data=N_var_abund, FUN=mean)
  N_noDD_mu <- aggregate(abundance ~ intro.regime, data=N_noDD_abund, FUN=mean)
  N_noDDvar_mu <- aggregate(abundance ~ intro.regime, data=N_noDDvar_abund, FUN=mean)
  
  DDdf <- 
    data.frame(extant[, e1], 
                     N_mu$abundance, 
                     extant[, e2], 
                     e_var[, e1], 
                     N_var_mu$abundance, 
                     e_var[, e2] )
  noDDdf <- 
    data.frame(e_noDD[, e1], 
                     N_noDD_mu$abundance, 
                     e_noDD[, e2], 
                     e_noDDvar[, e1], 
                     N_noDDvar_mu$abundance, 
                     e_noDDvar[, e2] )
 
  return(list(DD=DDdf, noDD=noDDdf))
}

establishment_abundance_table_5_50 <- make_table()
establishment_abundance_table_1_50 <- make_table(e1=1, a1=1)

# write.csv(establishment_abundance_table_5_50, "simulations/establishment_abundance_table_5_50.csv", row.names=FALSE)
# write.csv(establishment_abundance_table_1_50, "simulations/establishment_abundance_table_1_50.csv", row.names=FALSE)

mean_when_extinct <- function(x=5)
{
  intro.regime <- get.intro.regime(N=N, propagule.pressure = pp)
  
  cols <- 1:(x + max(intro.regime))
  N_when <- determine.when.extinct(N[, cols], intro.regime)
  N_var_when <- determine.when.extinct(N_var[, cols], intro.regime)
  N_noDD_when <- determine.when.extinct(N_noDD[, cols], intro.regime)
  N_noDDvar_when <- determine.when.extinct(N_noDDvar[, cols], intro.regime)
  TTEdf <- data.frame()

  for (i in unique(intro.regime))
  {
    TTE <- mean(N_when[intro.regime==i & N_when <= i+x] - i, na.rm=TRUE)
    TTE_var <- mean(N_var_when[intro.regime==i & N_var_when <= i+x] - i, na.rm=TRUE)
    TTE_noDD <- mean(N_noDD_when[intro.regime==i & N_noDD_when <= i+x] - i, na.rm=TRUE)
    TTE_noDDvar <- mean(N_noDDvar_when[intro.regime==i & N_noDDvar_when <= i+x] - i, na.rm=TRUE)
    TTEdf <- rbind(TTEdf, data.frame(TTE, TTE_var, TTE_noDD, TTE_noDDvar))
  }
  
  return(TTEdf)
}

TTEdf5 <- mean_when_extinct(x=5)
# write.csv(TTEdf5, "simulations/TTEdf5.csv", row.names=FALSE)
TTEdf50 <- mean_when_extinct(x=50)
# write.csv(TTEdf50, "simulations/TTEdf50.csv", row.names=FALSE)


# Very low amount of density dependence to prevent overflow errors
# e_loDDvar <- sapply(time_points, function(x) 1-extinct.x.after.intro(x=x, N=N_loDDvar, propagule.pressure=pp)$percent.extinct)

# High variance
# e_HIvar <- sapply(time_points, function(x) 1-extinct.x.after.intro(x=x, N=N_HIvar, propagule.pressure=pp)$percent.extinct)
# 
# e_noDDHIvar <- sapply(time_points, function(x) 1-extinct.x.after.intro(x=x, N=N_noDDHIvar, propagule.pressure=pp)$percent.extinct)

# Very low amount of density dependence to prevent overflow errors
# e_loDDHIvar <- sapply(time_points, function(x) 1-extinct.x.after.intro(x=x, N=N_loDDHIvar, propagule.pressure=pp)$percent.extinct)

# write.csv(extant, "simulations/extant.csv", row.names=FALSE)
# write.csv(e_noDD, "simulations/e_noDD.csv", row.names=FALSE)
# write.csv(e_var, "simulations/e_var.csv", row.names=FALSE)
# write.csv(e_noDDvar, "simulations/e_noDDvar.csv", row.names=FALSE)
# write.csv(e_HIvar, "simulations/e_var.csv", row.names=FALSE)
# write.csv(e_noDDHIvar, "simulations/e_noDDvar.csv", row.names=FALSE)


#-----------
# Plot results
#-----------

# Long term: density dependence versus not
# colors <- 1:4
# # colors <- c("dodgerblue", "gold", "darkred", "darkgreen")
# par(mfrow=c(1,2), oma=c(3,3,0,0), mar=c(3,3,2,1))
# matplot(time_points, t(extant), type="l", lty=1, col=colors, lwd=2, las=1, xlab=NA, ylab=NA, main="With negative density dependence")
# matplot(time_points, t(e_noDD), type="l", lty=1, col=colors, lwd=2, las=1, xlab=NA, ylab=NA, main="No density dependence")
# mtext(side=1, text="Generations after final introduction event", outer=TRUE)
# mtext(side=2, text="Proportion of populations extant", outer=TRUE)
# 
# # Short term: variation vs not in same plot; density dependence vs. not in diff panels
# interval <- 1:5
# matplot(time_points[interval], t(extant[ , interval]), type="l", lty=1, col=colors, lwd=2, las=1, xlab=NA, ylab=NA, main="With negative density dependence")
# matplot(time_points[interval], t(e_var[ , interval]), type="l", lty=2, col=colors, lwd=2, add=TRUE)
# 
# matplot(time_points[interval], t(e_noDD[ , interval]), type="l", lty=1, col=colors, lwd=2, las=1, xlab=NA, ylab=NA, main="No density dependence")
# matplot(time_points[interval], t(e_noDDvar[ , interval]), type="l", lty=2, col=colors, lwd=2, add=TRUE)
# 
# # Long term: variation vs not in same plot; density dependence vs. not in diff panels
# interval <- 45:50
# matplot(time_points[interval], t(extant[ , interval]), type="l", lty=1, col=colors, lwd=2, las=1, xlab=NA, ylab=NA, main="With negative density dependence")
# matplot(time_points[interval], t(e_var[ , interval]), type="l", lty=2, col=colors, lwd=2, add=TRUE)
# 
# matplot(time_points[interval], t(e_noDD[ , interval]), type="l", lty=1, col=colors, lwd=2, las=1, xlab=NA, ylab=NA, main="No density dependence")
# matplot(time_points[interval], t(e_noDDvar[ , interval]), type="l", lty=2, col=colors, lwd=2, add=TRUE)