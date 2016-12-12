library(nimble)
library(coda)

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

F_mated_adjust <- nimbleFunction(run = function(Nt = integer(), F_migrants = integer(), F_residents = integer()) {
  returnType(integer())
  
  F_mated <- F_migrants + F_residents
  
  if (F_mated == Nt) {
    F_mated <- F_migrants
  }
  
  return(F_mated)
})

nbbgNIMBLE_run <- function(data, priors, nbbgInits = NA, nchains = 1, niter = 10000, nburnin = 1000) {
  #### F_mated_adjust nimbleFunction ####
  # This nimbleFunction should assign the number of mated females. Default is all females in the population, 
  # but this needs adjustment if no males are present. Then the number of mated females is the number
  # of migrant females, because they arrive from a large external source population with males
  
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
    R0 ~ dgamma(shape = R0shape, scale = R0scale)
    alpha ~ dgamma(shape = alphaShape, scale = alphaScale)
    kD ~ dgamma(shape = kDshape, scale = kDscale)
    kE ~ dgamma(shape = kEshape, scale = kEscale)
    
  })
  
  #### constants ####
  nbbgConsts <- list(reps = nrow(data),
                     p = 0.5,
                     R0shape = priors$R0["shape"],
                     R0scale = priors$R0["scale"],
                     alphaShape = priors$alpha["shape"],
                     alphaScale = priors$alpha["scale"],
                     kDshape = priors$kD["shape"],
                     kDscale = priors$kD["scale"],
                     kEshape = priors$kE["shape"],
                     kEscale = priors$kE["scale"])
  
  #### data ####
  nbbgData <- list(Nt = data$Nt,
                   Ntplus1 = data$Ntplus1,
                   residents = data$residents,
                   migrants = data$migrants)
  
  
  #### initial values, arbitrary; these will be updated when using multiple chains ####
  if (is.na(nbbgInits)) {
    nbbgInits <- list(R0 = 1,
                      alpha = 0.005,
                      kD = 2,
                      kE = 15,
                      RE = rep(1, times = nrow(data)),
                      lambda = rep(1, times = nrow(data)),
                      F_migrants = floor(data$migrants / 2),
                      F_residents = floor(data$residents / 2),
                      F_mated = floor(data$Nt / 2))
  }
  
  #### nimbleModel ####
  nbbg <- nimbleModel(code = nbbgCode, name = 'NBBg', constants = nbbgConsts,
                      data = nbbgData, inits = nbbgInits)
  
  #### nimbleModel configuration ####
  
  nbbgConfig <- configureMCMC(nbbg)
  # nbbgConfig$printSamplers()
  
  #### build the nimble model ####
  nbbgMCMC <- buildMCMC(nbbgConfig)
  
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
  # Return coda object so we can use effective size and convergence diagnostic criteria
  
  sampsList <- runMCMC(CnbbgMCMC, niter = niter, nchains = nchains, inits = inits_generate(data), nburnin = nburnin, returnCodaMCMC = TRUE)
  
  #### Pull out samples for key parameters ####
  if (is.mcmc.list(sampsList)) {
    key <- lapply(sampsList, function(x) x[, c("R0", "alpha", "kE", "kD")])
  } else {
    key <- sampsList[, c("R0", "alpha", "kE", "kD")]
  }
  
  return(key)
}