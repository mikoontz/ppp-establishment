### Script to run Bayesian analysis with NBBG model of Tribolium population growth
###
### Author: Michael Koontz
###
### Date Created: 20141103
### Last Updated: 20161122
###
### Purpose is to estimate parameters of the NBBg model separately for the populations in stable environment versus fluctuating environments
#### Testing environmental stochasticity effect ####

# We expect the fluctuating environment to have a lower value for kE than the stable environment which means MORE variability
source("scripts/data-carpentry/generate-tidy-data.R")
source("scripts/NBBg-mcmc.R")

beetles <- read.csv("data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("data/attributes.csv")

d <- tidy.beetles(beetles)
N0 <- d$Nt[, 1:2]
N1 <- d$Ntp1[, 1:2]
b <- merge(N0, N1, by = "ID")
b <- merge(b, attributes, "ID")
b <- subset(b, select = c("ID", "0", "1", "environment"))
colnames(b) <- c("ID", "Nt", "Ntplus1", "env")
b$residents <- 0
b$migrants <- b$Nt

fluc <- subset(b, subset=env=="fluctuating")
stab <- subset(b, subset=env=="stable")

n.mcmc <- 1000
inits <- c(R0=2, kE=22, kD=10, alpha=0.005)
priors.shape <- c(R0=2.6, kE=17.6, kD=1.07, alpha=0.0037)
priors.scale <- c(R0=1, kE=1, kD=1, alpha=1)

fluc.tune <- c(R0=0.08, kE=8, kD=1, alpha=0.008, RE=1)
stab.tune <- c(R0=0.08, kE=8, kD=1, alpha=0.008, RE=1)

mcmc.fluc <- NBBG.mcmc(fluc, priors.shape, priors.scale, inits, fluc.tune, n.mcmc)
mcmc.stab <- NBBG.mcmc(stab, priors.shape, priors.scale, inits, stab.tune, n.mcmc)

mcmc.fluc[["accept"]]
mcmc.stab[["accept"]]

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


