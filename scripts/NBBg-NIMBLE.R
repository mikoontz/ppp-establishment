library(nimble)

#### Read and prepare data ####
data <- read.csv("data/initial-density-dependence.csv")
colnames(data) <- c("ID", "block", "env", "setup_order", "Nt", "Ntplus1", "check")
data <- subset(data, select = c(ID, Nt, Ntplus1))
data$residents <- 0
data$migrants <- data$Nt
head(data)
#### Small sample for NIMBLE Developers Team ####
# dput(data[sample(nrow(data), 20), ])
# data <- structure(list(ID = c(32L, 14L, 8L, 13L, 16L, 9L, 23L, 12L, 25L, 
#                               29L, 9L, 19L, 11L, 12L, 3L, 21L, 13L, 28L, 12L, 1L), 
#                        Nt = c(200L, 20L, 10L, 20L, 30L, 10L, 50L, 20L, 50L, 75L, 10L, 30L, 20L, 20L, 5L, 50L, 20L, 75L, 20L, 5L),
#                        Ntplus1 = c(36L, 12L, 15L, 24L, 27L, 6L, 16L, 18L, 26L, 35L, 17L, 27L, 19L, 21L, 12L, 16L, 8L, 47L, 18L, 4L), 
#                        residents = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                        migrants = c(200L, 20L, 10L, 20L, 30L, 10L, 50L, 20L, 50L, 75L, 10L, 30L, 20L, 20L, 5L, 50L, 20L, 75L, 20L, 5L)), 
#                   .Names = c("ID", "Nt", "Ntplus1", "residents", "migrants"), 
#                   row.names = c(92L, 74L, 38L, 13L, 16L, 69L, 83L, 42L, 55L, 59L, 102L, 112L, 104L, 12L, 33L, 81L, 43L, 58L, 72L, 94L), 
#                   class = "data.frame")

#### F_mated_adjust nimbleFunction ####
# This nimbleFunction should assign the number of mated females. Default is all females in the population, 
# but this needs adjustment if no males are present. Then the number of mated females is the number
# of migrant females, because they arrive from a large external source population with males

F_mated_adjust <- nimbleFunction(run = function(Nt = integer(), F_migrants = integer(), F_residents = integer()) {
  returnType(integer())
  
  F_mated <- F_migrants + F_residents
  
  if (F_mated == Nt) {
    F_mated <- F_migrants
    }
  
  return(F_mated)
})

#### nimbleCode ####
nbbgCode <- nimbleCode({ 
  
  for (i in 1:reps){

    F_residents[i] ~ dbinom(size = residents[i], prob = p) # Number of resident females from binomial distribution
    F_migrants[i] ~ dbinom(size = migrants[i], prob = p) # Number of migrant females from binomial distribution
    
    # Number of mated females equal to all of the females in the population except when population is all females
    # Then only the migrants are considered mated (since they come from a large external source population that
    # includes males)

    F_mated[i] <- F_mated_adjust(Nt = Nt[i],
                                 F_migrants = F_migrants[i],
                                 F_residents = F_residents[i])

    # The particular intrinsic growth rate for a population comes from a gamma distribution to incorporate 
    # environmental stochasticity
    RE[i] ~ dgamma(shape = kE, scale = R0 / kE)
    
    # This is the Ricker skeleton describing expectation of population size including density dependent 
    # process (cannibalism, whose rate is described directly as alpha)
    mu[i] <- 1 / p * F_mated[i] * RE[i] * exp(-alpha * Nt[i])

    # Negative binomial distribution to account for demographic heterogeneity is implemented as a gamma mixture
    # of Poissons. Notably, the overdispersion parameter (size parameter in negative binomial, equivalent to shape
    # parameter of the gamma distribution of Poisson's intensity parameter) is density dependent.
    lambda[i] ~ dgamma(shape = kD * F_mated[i], scale = mu[i] / (kD * F_mated[i]))
    
    # Demographic stochasticity is accounted for using a Poisson distribution
    Ntplus1[i] ~ dpois(lambda = lambda[i])

  }
  
  # Priors from Melbourne & Hastings (2008)
  R0 ~ dgamma(shape = 2.6, scale = 1)
  alpha ~ dgamma(shape = 0.0037, scale = 1)
  kD ~ dgamma(shape = 1.07, scale = 1)
  kE ~ dgamma(shape = 17.6, scale = 1)
  
})

#### constants ####
nbbgConsts <- list(reps = nrow(data),
                   p = 0.5)

#### data ####
nbbgData <- list(Nt = data$Nt,
                 Ntplus1 = data$Ntplus1,
                 residents = data$residents,
                 migrants = data$migrants)


#### initial values ####
nbbgInits <- list(R0 = 1,
                  alpha = 0.005,
                  kD = 2,
                  kE = 15,
                  RE = rep(1, times = nrow(data)),
                  lambda = rep(1, times = nrow(data)),
                  F_migrants = floor(data$migrants / 2),
                  F_residents = floor(data$residents / 2),
                  F_mated = floor(data$Nt / 2))

#### nimbleModel ####
nbbg <- nimbleModel(code = nbbgCode, name = 'NBBg', constants = nbbgConsts,
                    data = nbbgData, inits = nbbgInits)

#### build the nimble model ####
nbbgMCMC <- buildMCMC(nbbg)

#### function to generate initial values ####
inits_generate <- function(data) {
  R0 <- runif(n = 1, min = 0.5, max = 3.0)
  alpha <- runif(n = 1, min = 0.001, max = 0.01)
  kE <- runif(n = 1, min = 0.5, max = 25)
  kD <- runif(n = 1, min = 0.5, max = 25)
  RE <- rep(R0, times = nrow(data))
  lambda <- RE
  F_migrants <- floor(data$migrants / 2)
  F_residents <- floor(data$residents / 2)
  F_mated <- floor(data$Nt / 2)
  
  initsList <-  list(R0 = R0,
                     alpha = alpha,
                     kD = kD,
                     kE = kE,
                     RE = RE,
                     lambda = lambda,
                     F_migrants = F_migrants,
                     F_residents = F_residents,
                     F_mated = F_mated)
}

Cnbbg <- compileNimble(nbbg)
CnbbgMCMC <- compileNimble(nbbgMCMC, project = nbbg)

#### Fit the model; multiple chains ####

sampsList <- runMCMC(CnbbgMCMC, niter = 30000, nchains = 4, inits = inits_generate(data), returnCodaMCMC = TRUE)






# #### Fit the model ####
# CnbbgMCMC$run(20000)
# 
# #### Extract the samples ####
# samps <- as.data.frame(as.matrix(CnbbgMCMC$mvSamples))
# key <- samps[, c("R0", "alpha", "kD", "kE")]
# 
# #### Plot the samples ####
# samps_plot <- function(x) {
#   plot(key[[x]], type = "l", xlab = "Iteration", ylab = "Value")
#   mtext(side = 2, names(key[x]), line = 5)
#   
#   plot(density(key[[x]]), main = NA)
#   abline(v = mean(key[[x]]), col = "red", lwd = 2)
#   return(invisible())
# }
# 
# par(mfrow = c(ncol(key), 2), oma = c(0, 3, 2, 0))
# lapply(1:length(key), samps_plot)
# mtext(side = 3, outer = TRUE, text = "Trace Plot(s)", at = 0.25)
# mtext(side= 3, outer = TRUE, text = "Density Plot(s)", at = 0.75)

