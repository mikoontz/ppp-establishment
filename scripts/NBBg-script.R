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

colnames(data)[colnames(data)=="census"] <- "Ntplus1"
colnames(data)[colnames(data)=="N0"] <- "migrants"
data$residents <- 0
data$Nt <- data$migrants + data$residents


#--------------
# For my data
#--------------
n.mcmc <- 1000
inits <- c(R0=2, kE=22, kD=10, alpha=0.005)
priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)
tune <- c(R0=0.1, kE=10, kD=1, alpha=0.0004, RE=1)

head(data)

mcmc <- NBBG.mcmc(data, priors.shape, priors.scale, inits, tune, n.mcmc)

burned <- burn.in(mcmc[c("R0", "kE", "kD", "alpha")], 2000)
burned.df <- burn.in.df(mcmc[c("RE", "F.migrants", "F.residents")], 2000)

# mcmc[["accept"]]

#### Write the samples data file for my data ####
# samples <- data.frame(R0=burned[["R0"]], kE=burned[["kE"]], kD=burned[["kD"]], alpha=burned[["alpha"]])
# write.csv(samples, 'NBBg samples.csv', row.names=FALSE)

#### Trace plots and a sample trace plot for RE (1st population) ####


# pdf('test.pdf')
# par(mfrow=c(3,2))
# plot(data$Nt, data$Ntplus1, pch=19)
# # lines(x=x, y=x * R0 * exp(-alpha*x), col="red")
# plot(mcmc[["R0"]], type="l", main="R0")
# plot(mcmc[["kE"]], type="l", main="kE")
# plot(mcmc[["kD"]], type="l", main="kD")
# plot(mcmc[["alpha"]], type="l", main=expression(alpha))
# plot(mcmc[["RE"]][2, ], type="l", main="RE")
# burned <- burn.in(mcmc[c("R0", "kE", "kD", "alpha")], 2000)
# burned.df <- burn.in.df(mcmc[c("RE", "F.migrants", "F.residents")], 2000)
# dev.off()
#### Figuring out the modes and 95% credible intervals ####
# parameters <- c("R0", "kE", "kD", "alpha")
# results <- data.frame(parameter=c(parameters, "RE"), mean=0, mode=0, lowerCI=0, upperCI=0)
# 
# for (k in 1:length(parameters) )
# {
#   results$mean[results$parameter==parameters[k]] <- mean(burned[[parameters[k]]])
#   results$mode[results$parameter==parameters[k]] <- density(burned[[parameters[k]]])$x[which.max(density(burned[[parameters[k]]])$y)]
#   results[results$parameter==parameters[k], 4:5] <- quantile(burned[[parameters[k]]], prob=c(0.025, 0.975))
# }
# 
# results$mean[results$parameter=="RE"] <- mean(burned.df[["RE"]])
# results$mode[results$parameter=="RE"] <- density(burned.df[["RE"]])$x[which.max(density(burned.df[["RE"]])$y)]
# results[results$parameter=="RE", 4:5] <- quantile(burned.df[["RE"]], prob=c(0.025, 0.975))

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

#### Simulated data to validate the MCMC algorithm ####
### Simulate data

# df <- data.frame(residents=rep(0,2600), migrants=c(rep(4,500), rep(5, 500), rep(10, 500), rep(20, 500), rep(50, 200), rep(100,200), rep(200, 100), rep(1000, 100)) )
# df <- data.frame(residents=rep(0,250), migrants=c(rep(seq(2,100), times=2), rep(seq(125, 500, by=25), times=2), rep(seq(550, 1000, by=50), times=2)) )
# 
# R0 <- 2.6
# kE <- 17.6
# kD <- 1.07
# alpha <- 0.0037
# 
# test <- Tribolium.NBBG(data=df, R0=R0, kE=kE, kD=kD, alpha=alpha)
# x <- 1:1000
# 
# # par(mfrow=c(3,2))
# plot(test$Nt, test$Ntplus1, pch=".")
# lines(x=x, y=x * R0 * exp(-alpha*x), col="red")
# head(test)

# priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
# priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)
# tune <- c(R0=0.2, kE=5, kD=0.75, alpha=0.0001, RE=1.5)
# 
# mcmc.test <- NBBG.mcmc(test, priors.shape, priors.scale, inits, tune, n.mcmc)
# mcmc <- mcmc.test
# mcmc[["accept"]]


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

#### Testing environmental stochasticity effect ####

# We expect the fluctuating environment to have a greater value for kE than the stable environment

# d <- read.csv("Census data through generation 10.csv")
# d <- subset(d, subset=Generation==1)
# d <- subset(d, select=c(ID, Stochastic, N0, Census))
# colnames(d)[2:4] <- c("env", "migrants", "Ntplus1")
# d$residents <- 0
# d$Nt <- d$residents + d$migrants
# head(d)
# 
# fluc <- subset(d, subset=env=="fluctuating")
# stab <- subset(d, subset=env=="stable")
# 
# n.mcmc <- 10000
# inits <- c(R0=2, kE=22, kD=10, alpha=0.005)
# priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
# priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)

# fluc.tune <- c(R0=0.08, kE=8, kD=1, alpha=0.008, RE=1)
# stab.tune <- c(R0=0.08, kE=8, kD=1, alpha=0.008, RE=1)
# 
# mcmc.fluc <- NBBG.mcmc(fluc, priors.shape, priors.scale, inits, fluc.tune, n.mcmc)
# mcmc.stab <- NBBG.mcmc(stab, priors.shape, priors.scale, inits, stab.tune, n.mcmc)
# 
# mcmc.fluc[["accept"]]
# mcmc.stab[["accept"]]

# Fluctuating environment first
# x <- 1:20
# par(mfrow=c(3,2))
# plot(fluc$Nt, fluc$Ntplus1, pch=19)
# # lines(x=x, y=x * R0 * exp(-alpha*x), col="red")
# plot(mcmc.fluc[["R0"]], type="l", main="R0")
# plot(mcmc.fluc[["kE"]], type="l", main="kE")
# plot(mcmc.fluc[["kD"]], type="l", main="kD")
# plot(mcmc.fluc[["alpha"]], type="l", main=expression(alpha))
# plot(mcmc.fluc[["RE"]][2, ], type="l", main="RE")
# burned.fluc <- burn.in(mcmc.fluc[c("R0", "kE", "kD", "alpha")], 2000)
# burned.df.fluc <- burn.in.df(mcmc.fluc[c("RE", "F.migrants", "F.residents")], 2000)
# 
# plot(stab$Nt, stab$Ntplus1, pch=19)
# plot(mcmc.stab[["R0"]], type="l", main="R0")
# plot(mcmc.stab[["kE"]], type="l", main="kE")
# plot(mcmc.stab[["kD"]], type="l", main="kD")
# plot(mcmc.stab[["alpha"]], type="l", main=expression(alpha))
# plot(mcmc.stab[["RE"]][2, ], type="l", main="RE")
# burned.stab <- burn.in(mcmc.stab[c("R0", "kE", "kD", "alpha")], 2000)
# burned.df.stab <- burn.in.df(mcmc.stab[c("RE", "F.migrants", "F.residents")], 2000)
# 
# stack <- c(density(burned.fluc[["kE"]])$y, density(burned.stab[["kE"]])$y)
# max.y.idx <- which.max(stack)
# max.y <- stack[max.y.idx]
# 
# par(mfrow=c(1,1))
# plot(density(burned.fluc[["kE"]]), lwd=2, col="red", main="Parameter estimate for environmental stochasticity", ylim=c(0,max.y))
# lines(density(burned.stab[["kE"]]), lwd=2, col="blue")
# legend("topright", legend=c("Stable", "Fluctuating"), col=c("blue", "red"), lwd=3)

# mu is R0
# var is R0^2 / kE

# fRE.mu <- burned.fluc[["R0"]]
# fRE.s2 <- fRE.mu^2 / burned.fluc[["kE"]]
# 
# sRE.mu <- burned.stab[["R0"]]
# sRE.s2 <- sRE.mu^2 / burned.stab[["kE"]]
#   
# t.test(burned.fluc[["kE"]], burned.stab[["kE"]])
# 
# stack <- c(density(burned.fluc[["R0"]])$y, density(burned.stab[["R0"]])$y)
# max.y.idx <- which.max(stack)
# max.y <- stack[max.y.idx]
# 
# par(mfrow=c(1,1))
# plot(density(burned.stab[["R0"]]), lwd=2, col="blue", main="Parameter estimate for R0", ylim=c(0,max.y))
# lines(density(burned.fluc[["R0"]]), lwd=2, col="red")
# legend("topright", legend=c("Stable", "Fluctuating"), col=c("blue", "red"), lwd=3)
# 
# stack <- c(density(sRE.s2)$y, density(fRE.s2)$y)
# max.y.idx <- which.max(stack)
# max.y <- stack[max.y.idx]
# 
# par(mfrow=c(1,1))
# plot(density(sRE.s2), lwd=2, col="blue", main="Estimate for derived parameter: variance of RE", ylim=c(0,max.y))
# lines(density(fRE.s2), lwd=2, col="red")
# legend("topright", legend=c("Stable", "Fluctuating"), col=c("blue", "red"), lwd=3)



# 
# 
# 

#### Trace plots after removing burn in ####

# par(mfrow=c(3,2))
# plot(test$Nt, test$Ntplus1, pch=19)
# lines(x=x, y=x * R0 * exp(-alpha*x), col="red")
# plot(burned[["R0"]], type="l", main="R0")
# plot(burned[["kE"]], type="l", main="kE")
# plot(burned[["kD"]], type="l", main="kD")
# plot(burned[["alpha"]], type="l", main=expression(alpha))
# plot(burned.df[["RE"]][1, ], type="l", main="RE")
# 
# length(burned.df)

####  Plots of priors and posteriors ####
# par(mfrow=c(3,2))
# plot(density(burned[["R0"]]), main="R0")
# # abline(v=R0, col="red")
# R0.x <- seq(0, 12, 0.5)
# R0.y <- dgamma(R0.x, shape=priors.shape["R0"], scale=priors.scale["R0"])
# lines(R0.x, R0.y, col="blue")
# 
# plot(density(burned[["kE"]]), main="kE")
# # abline(v=kE, col="red")
# kE.x <- seq(5,35,0.5)
# kE.y <- dgamma(kE.x, shape=priors.shape["kE"], scale=priors.scale["kE"])
# lines(kE.x, kE.y, col="blue")
# 
# plot(density(burned[["kD"]]), main="kD")
# # abline(v=kD, col="red")
# kD.x <- seq(0,20,0.05)
# kD.y <- dgamma(kD.x, shape=priors.shape["kD"], scale=priors.scale["kD"])
# lines(kD.x, kD.y, col="blue")
# 
# plot(density(burned[["alpha"]]), main=expression(alpha))
# # abline(v=alpha, col="red")
# alpha.x <- seq(0, 0.4, 0.001)
# alpha.y <- dgamma(alpha.x, shape=priors.shape["alpha"], scale=priors.scale["alpha"])
# lines(alpha.x, alpha.y, col="blue")

# RE.x <- seq(0, 20, by=0.05)
# RE.y <- dgamma(RE.x, shape=kE, scale=R0/kE)
# plot(RE.x, RE.y, col="red", main="RE", type="l")
# lines(density(burned.df[["RE"]]))
# plot(density(burned.df[["RE"]]))




###
###
### Randomly sample 10 population ID's and plot the predicted and known number of migrant females

get.mode <- function(x)
{
	names(sort(-table(x)))[1]
}

par(mfrow=c(2, 5))
for (i in 1:length(ids))
{
	plot(table(mcmc[["F.migrants"]][ids[i], ]), ylab="Frequency", main=paste("Expected Female Migrants:", test$Nt[ids[i]]/2, "\nEstimated Female Migrants:", get.mode(mcmc[["F.migrants"]][ids[i], ]), "\nActual Female Migrants:", test$F.migrants[ids[i]] ) )
	
	abline(v=test$F.migrants[ids[i]], col="blue", lwd=4)	
	abline(v=get.mode(mcmc[["F.migrants"]][ids[i], ]), col="red", lwd=4)	

}
