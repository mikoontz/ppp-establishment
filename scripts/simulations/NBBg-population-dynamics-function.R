### Function to simulate Tribolium population dynamics
###
### Author: Michael Koontz
###
### Date Created: 20141107
### Last Updated: 20141114
###
### The purpose of this code is to simulate population growth or decline given parameter values for different kinds of stochasticity and starting population size (including number of individuals that are resident and number that are migrants)

###
###
###

### The data argument should be a dataframe with 2 columns: Number of residents at time t and number of migrants at time t

Tribolium.NBBG <- function(data, R0, kE, kD, alpha)
{
	reps <- nrow(data)
	p <- 0.5
	data$Nt <- data$residents + data$migrants
	
	# Create new column and assign all values to 0 for now
	data$Ntplus1 <- 0
	
	# Randomly draw resident females and migrant females separately
	data$F.residents <- rbinom(reps, size=data$residents, prob=p)
	data$F.migrants <- rbinom(reps, size=data$migrants, prob=p)
	
	# Total females is the sum of resident and migrant females
	data$F <- data$F.residents + data$F.migrants
	
	# Two cases where population will be extinct due to skewed sex ratio: (1) all males or (2) all resident females
	idx <- (data$F.residents == data$Nt) | (data$F == 0)
	idx.count <- reps - sum(idx)
	
	RE <- rgamma(n=reps, shape=kE, scale=R0/kE)
	
	# Special case where there are all females, but a mixture of migrants and residents. Only the migrant females would be able to produce eggs.
	
	idx2 <- (data$F == data$Nt) & (data$F.migrants != data$Nt)
	
	eggs <- 1/p * data$F * RE
	eggs[idx2] <- 1/p * data$F.migrants[idx2] * RE[idx2]
	mu <- eggs*exp(-alpha * data$Nt)
		
	data$Ntplus1[!idx] <- rnbinom(idx.count, mu=mu[!idx], size=kD * data$F[!idx])
	
	### And for the special case with all Females but a mixture of residents and migrants
	data$Ntplus1[idx2] <- rnbinom(sum(idx2), mu=mu[idx2], size=kD * data$F.migrants[idx2])

	data$Ntplus1[idx] <- 0
	data
}