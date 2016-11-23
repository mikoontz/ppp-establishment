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

# Set up data for analysis. Ensures column titles from actual data match the column titles that the MCMC code uses.
colnames(data)[colnames(data)=="census"] <- "Ntplus1"
colnames(data)[colnames(data)=="N0"] <- "migrants"
data$residents <- 0
data$Nt <- data$migrants + data$residents

# Setup MCMC parameters
n.mcmc <- 100000
inits <- c(R0=2, kE=22, kD=10, alpha=0.005)
priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)
tune <- c(R0=0.1, kE=10, kD=1.5, alpha=0.001, RE=1)

head(data)

mcmc <- NBBG.mcmc(data, priors.shape, priors.scale, inits, tune, n.mcmc)

burned <- burn.in(mcmc[c("R0", "kE", "kD", "alpha")], 2000)
burned.df <- burn.in.df(mcmc[c("RE", "F.migrants", "F.residents")], 2000)

# Percent of proposals accepted
(mcmc[["accept"]] / n.mcmc) * 100

R0 <- mean(burned[["R0"]])
alpha <- mean(burned[["alpha"]])
x <- 1:200
mu_Ntp1 <- x * R0 * exp(-alpha * x)

par(mfrow=c(3,2))
plot(data$Nt, data$Ntplus1, pch=19)
lines(x = x, y = mu_Ntp1, col="red")
plot(burned[["R0"]], type="l", main="R0")
plot(burned[["kE"]], type="l", main="kE")
plot(burned[["kD"]], type="l", main="kD")
plot(burned[["alpha"]], type="l", main=expression(alpha))
plot(burned.df[["RE"]][2, ], type="l", main="RE")

# Effective number of parameters for a single chain using method in coda package
burned_matrix <- as.matrix(as.data.frame(burned))
spec <- spectrum0.ar(burned_matrix)$spec
Neff <- ifelse(spec == 0, 0, nrow(burned_matrix) * apply(burned_matrix, 2, var) / spec)
Neff


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
