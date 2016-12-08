### MCMC Code for NBBG model of Tribolium population growth
###
### Author: Michael Koontz
###
### Date Created: 20141103
### Antepenultimate Update: 20141114
### Penultimate Update: 20150604 (changed proposals of discrete latent variables (number of resident females, number of migrant females) to sampling from a discrete uniform rather than rounding a sample from a continuous uniform)
### Last Update: 20161123 and now using version control (see result of diff call to see changes)

### The purpose of this code is to generate a chain of MCMC samples to estimate the posterior distributions for parameters in a population dynamics model. By estimating what variables might be at play in the system, we can hope to isolate the role of environmental stochasticity and then manipulate it during simulations.

###
### Parameters to estimate:
###
### R0: per capita population growth rate in an average environment
### alpha: adult search rate for eggs (the density dependent parameter)
### kD: shape parameter for gamma distribution for demographic heterogeneity (higher number means it matters less)
### kE: shape parameter for gamma distribution for environmental stochasticity (higher number means it matters less)

###
### Parameters:
###
### RE: per capita population growth rate across all environments (incorporates kE value)
### F.migrants: the number of female migrants
### F.residents: the number of female residents


### Data:
###
### Ntplus1: Population size at time t+1
### Nt: Population size at time t
### 
###

library(coda)

NBBG.mcmc <- function(data, priors.shape, priors.scale, inits, tune, n.mcmc, p=0.5)
{
	reps <- nrow(data)

	### Set up storage variables for parameters and hyperparameters
	
	R0.save <- rep(0, n.mcmc)
	R0.save[1] <- inits["R0"]
	R0 <- inits["R0"]
	
	kE.save <- rep(0, n.mcmc)
	kE.save[1] <- inits["kE"]
	kE <- inits["kE"]
	
	kD.save <- rep(0, n.mcmc)
	kD.save[1] <- inits["kD"]
	kD <- inits["kD"]
	
	alpha.save <- rep(0, n.mcmc)
	alpha.save[1] <- inits["alpha"]
	alpha <- inits["alpha"]
	
	RE.save <- matrix(0, nrow = n.mcmc, ncol = reps, dimnames = list(NULL, paste0("RE_Pop", data$ID)))
	RE.save[1, ] <- inits["R0"]
	RE <- RE.save[1, ]

	names(RE) <- "RE"
	
	F.mated.save <- matrix(0, nrow = n.mcmc, ncol = reps, dimnames = list(NULL, paste0("F_mated_Pop", data$ID)))
	F.migrants.save <- matrix(0, nrow = n.mcmc, ncol = reps, dimnames = list(NULL, paste0("F_migrants_Pop", data$ID)))
	F.residents.save <- matrix(0, nrow = n.mcmc, ncol = reps, dimnames = list(NULL, paste0("F_residents_Pop", data$ID)))
	
	# Set initial conditions for Female migrants (half of all migrants)
	F.migrants <- ceiling(data$migrants * p)
	F.migrants.save[1, ] <- F.migrants

	# Set initial conditions for Female residents (half of all residents)
	F.residents <- ceiling(data$residents * p)
	F.residents.save[1, ] <- F.residents

	F.mated <- F.migrants + F.residents # Default number of mated females is all females
	F.mated.idx <- which(F.mated == data$Nt) # Indices where population is all females
	F.mated[F.mated.idx] <- F.migrants[F.mated.idx] # Populations with all females use only the migrant females for the number of mated females (because we assume they arrive mated from the large external source)
	F.mated.save[1, ] <- F.mated
	
	mu <- 1/p * F.mated * RE * exp(-alpha * data$Nt)

	accept <- c(R0=0, alpha=0, kE=0, kD=0, RE=0, F.migrants=0, F.residents=0)
	
	### Begin MCMC loop
		
	for (i in 2:n.mcmc)
	{
		### Print progress
		
		if(i%%100==0) cat(i," "); flush.console()
		
		###
		### Update alpha
		###
	
		# Proposal
		
		alpha.star <- rnorm(1, mean=alpha, sd=tune["alpha"])
		
		if (alpha >= 0)
		{
		# Because alpha.star changes mu, we need to recalculate
		mu.star <- 1/p * F.mated * RE * exp(-alpha.star * data$Nt)

		# Calculate mh ratio
		
		mh1 <- sum(dnbinom(data$Ntplus1, mu=mu.star, size=kD * F.mated + (F.mated == 0), log=TRUE), na.rm=TRUE) + dgamma(alpha.star, shape=priors.shape["alpha"], scale=priors.scale["alpha"], log=TRUE)
		
		mh2 <- sum(dnbinom(data$Ntplus1, mu=mu, size=kD * F.mated + (F.mated == 0), log=TRUE), na.rm=TRUE) + dgamma(alpha, shape=priors.shape["alpha"], scale=priors.scale["alpha"], log=TRUE)
					
		mh <- exp(mh1 - mh2)
		
		if (mh > runif(1))
		{
			alpha <- alpha.star
			mu <- mu.star
			accept["alpha"] <- accept["alpha"] + 1
		}
		} # End the 'if alpha is non-negative' statement

		###
		### Update kD
		###
		
		# Proposal
		
		kD.star <- rnorm(1, mean=kD, sd=tune["kD"])
		# kD.star <- rgamma(1, shape = tune["kD"], scale = kD / tune["kD"])
		  
		if (kD.star >= 0)
		{
		
		# Calculate mh ratio
		# correction1 <- dgamma(kD, shape = tune["kD"], scale = kD.star / tune["kD"], log = TRUE)
		# correction2 <- dgamma(kD.star, shape = tune["kD"], scale = kD / tune["kD"], log = TRUE)
		  
		mh1 <- sum(dnbinom(data$Ntplus1, mu=mu, size=kD.star * F.mated + (F.mated == 0), log=TRUE), na.rm=TRUE) + dgamma(kD.star, shape=priors.shape["kD"], scale=priors.scale["kD"], log=TRUE)
			
		mh2 <- sum(dnbinom(data$Ntplus1, mu=mu, size=kD * F.mated + (F.mated == 0), log=TRUE), na.rm=TRUE) + dgamma(kD, shape=priors.shape["kD"], scale=priors.scale["kD"], log=TRUE)
		
		# mh1 <- mh1 + correction1
		# mh2 <- mh2 + correction2
		
		mh <- exp(mh1 - mh2)
		
		if (mh > runif(1))
		{
			kD <- kD.star
			accept["kD"] <- accept["kD"] + 1
		}
		} # End if kD is non-negative
		
		###
		### Update RE
		###
	
		# Proposal
		
		RE.star <- rnorm(reps, mean=RE, sd=tune["RE"])
		RE.star[which(RE.star < 0)] <- NA # Any proposals of RE that are less than 0 (impossible values) are instead NA. This will propagate the nothingness through the mu.star calculations, the likelihoods, and the mh ratios such that none of them throw errors or warnings. We also need to ensure that none of those proposals get accepted by chance, which is done by making sure to use the which() function when determining for which of the populations to accept or reject the proposal RE. Because the RE proposals that were negative are now NA, they will be excluded from the candidate set of proposals to accept or reject. Thus, the only updates will be to proposals that meet 2 criteria: They're mh.ratio > runif(1) AND their mh.ratio is NOT NA.
		
		# Because RE changes mu, we need to recalculate
		
		mu.star <- 1/p * F.mated * RE.star * exp(-alpha * data$Nt)

###### Calculate mh ratio #####
		
		mh1 <- dnbinom(data$Ntplus1, mu=mu.star, size=kD * F.mated + (F.mated == 0), log=TRUE) + dgamma(RE.star, shape=kE, scale=R0/kE, log=TRUE)
		
		mh2 <- dnbinom(data$Ntplus1, mu=mu, size=kD * F.mated + (F.mated == 0), log=TRUE) + dgamma(RE, shape=kE, scale=R0/kE, log=TRUE)
		
		mh <- exp(mh1 - mh2)
		mh.idx <- which(mh > runif(n=reps) )
		
		RE[mh.idx] <- RE.star[mh.idx]
		mu[mh.idx] <- mu.star[mh.idx]
		accept["RE"] <- accept["RE"] + (length(mh.idx) / reps)


		###
		### Update R0
		###
		
		# Proposal

		R0.star <- rnorm(1, mean=R0, sd=tune["R0"])
		
		if (R0.star >= 0)
		{
		
		mh1 <- sum(dgamma(RE, shape=kE, scale=R0.star/kE, log=TRUE)) + dgamma(R0.star, shape=priors.shape["R0"], scale=priors.scale["R0"], log=TRUE) 
		
		mh2 <- sum(dgamma(RE, shape=kE, scale=R0/kE, log=TRUE)) + dgamma(R0, shape=priors.shape["R0"], scale=priors.scale["R0"], log=TRUE)
		
		mh <- exp(mh1 - mh2)
	
		# Decide whether to accept or reject proposal
		if (mh > runif(1))
		{
			R0 <- R0.star
			accept["R0"] <- accept["R0"] + 1
		}
		} # End if R0 proposal is non-negative statement
		
		###
		### Update kE
		###
		
		# Proposal
	
		kE.star <- rnorm(1, mean=kE, sd=tune["kE"])
				
		if (kE.star >= 0)
		{

		# Calculate mh ratio
	
		mh1 <- sum(dgamma(RE, shape=kE.star, scale=R0/kE.star, log=TRUE)) + dgamma(kE.star, shape=priors.shape["kE"], scale=priors.scale["kE"], log=TRUE)
		
		mh2 <- sum(dgamma(RE, shape=kE, scale=R0/kE, log=TRUE)) + dgamma(kE, shape=priors.shape["kE"], scale=priors.scale["kE"], log=TRUE) 
		
		mh <- exp(mh1 - mh2)
		
		# Decide whether to accept or reject proposal
		if (mh > runif(1))
		{
			kE <- kE.star
			accept["kE"] <- accept["kE"] + 1
		}
		} # End if kE proposal is non-negative
		
		###
		### Update F.migrants
		###

		# Discrete uniform proposal distribution
		
    F.migrants.star <- sapply(data$migrants, FUN=function(x) sample(x, size=1))
    
    # Because F.migrants changes F.mated, we need to recalculate
    F.mated.star <- F.migrants.star + F.residents

    F.mated.idx.star <- which((F.migrants.star + F.residents) == data$Nt) # Indices where proposal number of migrant females makes the population be all females
    
    F.mated.star[F.mated.idx.star] <- F.migrants.star[F.mated.idx.star]

    # Because F.mated changes mu, we need to recalculate		
		mu.star <- 1/p * F.mated.star * RE * exp(-alpha * data$Nt)

		# Calculate mh ratio
		
		mh1 <- dnbinom(data$Ntplus1, mu=mu.star, size= kD * F.mated.star + (F.mated.star == 0), log=TRUE) + dbinom(F.migrants.star, size=data$migrants, prob=p, log=TRUE) 
			
		mh2 <- dnbinom(data$Ntplus1, mu=mu, size=kD * F.mated + (F.mated == 0), log=TRUE) + dbinom(F.migrants, size=data$migrants, prob=p, log=TRUE)
	
		mh <- exp(mh1 - mh2)
		mh.idx <- which(mh > runif(reps))
		
		F.migrants[mh.idx] <- F.migrants.star[mh.idx]
		F.mated[mh.idx] <- F.mated.star[mh.idx]
		mu[mh.idx] <- mu.star[mh.idx]

		accept["F.migrants"] <- accept["F.migrants"] + length(mh.idx)/reps
		
		###
		### Update F.residents
		###
	
    F.residents.star <- sapply(data$residents, FUN=function(x) sample(x, size=1))
    
    # Because F.migrants changes F.mated, we need to recalculate
    F.mated.star <- F.migrants + F.residents.star
		
    F.mated.idx.star <- which((F.migrants + F.residents.star) == data$Nt) # Indices where proposal number of resident females makes the population be all females
    
    # Because F.mated changes mu, we need to recalculate		
    mu.star <- 1/p * F.mated.star * RE * exp(-alpha * data$Nt)

		mh1 <- dnbinom(data$Ntplus1, mu=mu.star, size=kD * F.mated.star + (F.mated.star == 0), log=TRUE) + dbinom(F.residents.star, size=data$residents, prob=p, log=TRUE)
		
		mh2 <- dnbinom(data$Ntplus1, mu=mu, size=kD * F.mated + (F.mated == 0), log=TRUE) + dbinom(F.residents, size=data$residents, prob=p, log=TRUE)
	
		mh <- exp(mh1 - mh2)
		mh.idx <- which(mh > runif(reps))
		
		F.residents[mh.idx] <- F.residents.star[mh.idx]
		F.mated[mh.idx] <- F.mated.star[mh.idx]
		mu[mh.idx] <- mu.star[mh.idx]
	
		accept["F.residents"] <- accept["F.residents"] + length(mh.idx)/reps
	
					
		###
		### Save samples
		###
		
		R0.save[i] <- R0
		kE.save[i] <- kE
		kD.save[i] <- kD
		alpha.save[i] <- alpha
		RE.save[i, ] <- RE
		F.migrants.save[i, ] <- F.migrants
		F.residents.save[i, ] <- F.residents
		F.mated.save[i, ] <- F.mated
		
	}
	
	samps <- cbind(R0 = R0.save, kE = kE.save, kD = kD.save, alpha = alpha.save, F.migrants.save, F.residents.save, F.mated.save)
	
	list(samps = mcmc(samps), accept = accept)
	
}

burn.in <- function(chain, burn.in)
{
  return(mcmc(chain[-(1:burn.in), ]))
}










