### Script to run Bayesian analysis with NBBG model of Tribolium population growth
###
### Author: Michael Koontz
###
### Date Created: 20141103
### Last Updated: 20161122
###
### The purpose of this code is to validate the MCMC algorithm used to calculate NBBg parameters on the Parsing Propagule Pressure initial density dependence side experiment. We use known values of the parameters and simulate an NBBg model, then fit our hierarchical model to those simulated data to see if we can reobtain the known parameter values. We then fit our hierarchical model to data from Melbourne and Hastings (2008) to see if we can reobtain the parameter values they estimated from their analysis.

source("scripts/simulations/NBBg-population-dynamics-function.R")
source("scripts/NBBg-mcmc.R")

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

# Known parameter values 
R0 <- 2.6
kE <- 17.6
kD <- 1.07
alpha <- 0.0037
samps <- data.frame(R0, alpha, kE, kD)


test_data <- NBBg(Tf = 2, 
                  total.reps = nrow(df),
                  N = N,
                  past.residents = 0, 
                  migrants = migrants, 
                  p = 0.5, 
                  samps = samps)

test_data <- data.frame(ID = 1:nrow(test_data), 
                        test_data, 
                        migrants = migrants$migrants0, 
                        residents = past.residents)
head(test_data)

x <- 1:1000

# par(mfrow=c(3,2))
plot(test_data$Nt, test_data$Ntplus1, pch=16, cex = 0.5)
lines(x=x, y=x * R0 * exp(-alpha*x), col="red")

n.mcmc <- 1000
priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)
tune <- c(R0=0.2, kE=5, kD=0.75, alpha=0.0001, RE=1.5)

inits <- list(c(R0=2, kE=22, kD=10, alpha=0.005),
              c(R0=5, kE=2, kD=25, alpha=0.015),
              c(R0=0.5, kE=35, kD=1, alpha=0.0005))

mcmc_test_data <- NBBG.mcmc(data = test_data, 
                            priors.shape = priors.shape, 
                            priors.scale = priors.scale, 
                            inits = inits[[1]], 
                            tune = tune,
                            n.mcmc = n.mcmc)

mcmc_test_data[["accept"]]

key_params <- mcmc_test_data$samps[, c("R0", "kE", "kD", "alpha")]
summary(key_params)

#### Data from Melbourne & Hastings (2008) ####

# melbourne <- read.csv("melbourne_ricker_data.csv")
# head(melbourne)
# melbourne$At <- round(melbourne$At)
# melbourne$Atp1 <- round(melbourne$Atp1)
# melbourne$residents <- 0
# melbourne$migrants <- melbourne$Nt
# 
# ### Define MCMC arguments
# 
# n.mcmc <- 5000
# inits <- c(R0=2, kE=22, kD=10, alpha=0.005)

# colnames(melbourne)[colnames(melbourne)=="Atp1"] <- "Ntplus1"
# colnames(melbourne)[colnames(melbourne)=="At"] <- "migrants"
# melbourne$residents <- 0
# melbourne$Nt <- melbourne$migrants + melbourne$residents
# melbourne
# plot(melbourne$Nt, melbourne$Ntplus1, pch=19)
#--------------
# For Brett's data
# priors.shape <- c(R0=0.001, kE=0.001, kD=0.001, alpha=0.001)
# priors.scale <- c(R0=1000, kE=1000, kD=1000, alpha=1000)
# tune <- c(R0=0.2, kE=15, kD=1, alpha=0.0002, RE=1.2)
# mcmc.melbourne <- NBBG.mcmc(melbourne, priors.shape, priors.scale, inits, tune, n.mcmc)
# mcmc.melbourne[["accept"]]

# par(mfrow=c(3,2))
# plot(melbourne$Nt, melbourne$Ntplus1, pch=19)
# lines(x=x, y=x * results[results$parameter=="R0", "mean"] * exp(-results[results$parameter=="alpha", "mean"]*x), col="red")
# plot(mcmc.melbourne[["R0"]], type="l", main="R0")
# plot(mcmc.melbourne[["kE"]], type="l", main="kE")
# plot(mcmc.melbourne[["kD"]], type="l", main="kD")
# plot(mcmc.melbourne[["alpha"]], type="l", main=expression(alpha))
# plot(mcmc.melbourne[["RE"]][1, ], type="l", main="RE")
# burned <- burn.in(mcmc.melbourne[c("R0", "kE", "kD", "alpha")], 500)
# burned.df <- burn.in.df(mcmc.melbourne[c("RE", "F.migrants", "F.residents")], 150)
# dev.off()

#### Write the samples file for Brett's data ####
# samples <- data.frame(R0=burned[["R0"]], kE=burned[["kE"]], kD=burned[["kD"]], alpha=burned[["alpha"]])
# 
# write.csv(samples, 'NBBg samples melbourne.csv', row.names=FALSE)