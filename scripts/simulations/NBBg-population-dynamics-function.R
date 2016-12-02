### Function to simulate Tribolium population dynamics
###
### Author: Michael Koontz
###
### Date Created: 20141107
### Last Updated: 20141114
###
### The purpose of this code is to simulate population growth or decline over Tf generations given parameter values for different kinds of stochasticity and starting population size (including number of individuals that are resident and number that are migrants)

# Simulate NBBg dynamics
NBBg <- function(Tf, total.reps, past.residents, p, N, migrants, samps, overflowGuard=FALSE)
{
  for (t in 2:Tf)
  {
    # Define number of females for mixed population of new migrants and residents from previous generation
    
    F.residents <- rbinom(total.reps, size = past.residents, prob = p)
    F.migrants <- rbinom(total.reps, size = migrants[,t-1], prob = p)
    F.mated <- F.residents + F.migrants
    
    # Two cases where sex ratio affects how simulations will run
    # If population is all females, it isn't necessarily doomed to extinction because migrant females arrived mated from the large external source population. Thus, if a population is all females then the number of mated females is equal to the number of migrant females. This has two consequences:
    # (1) A smaller mean number of offspring because 1/p * F.mated is smaller than 1/p * F.all
    # (2) A larger effect of demographic heterogeneity, because the kD * F.mated is smaller than kD * F.all in the overdispersion parameter in the negative binomial distribution for demographic heterogeneity
    
    idx.F.mated.change <- which((F.residents + F.migrants) == N[, t-1])
    F.mated[idx.F.mated.change] <- F.migrants[idx.F.mated.change]
    
    idx <- sample(1:nrow(samps), size = total.reps, replace = TRUE)
    Re.shape <- samps[idx, "kE"]
    Re.scale <- samps[idx, "R0"] / Re.shape
    Re <- rgamma(total.reps, shape = Re.shape, scale = Re.scale)
    
    # Note that mu will be 0 when there are all males in the population
    mu <- 1/p * F.mated * Re * exp(-samps[idx, "alpha"] * N[ , t-1])
    
    N[, t] <- rnbinom(n=total.reps, mu=mu, size=samps[idx, "kD"] * F.mated + (F.mated == 0)) # Ensures we'll get a number if F.mated == 0 rather than an NaN. N[, t] will still always be 0 if F.mated is 0.
    
    
    # In scenarios with no density dependence and high environment fluctuation, populations can explode to sizes that are beyond what can be handled by the computer-- that is, there are overflow errors where very (very very) large numbers just become NaN and then the other functions think those populations went extinct. If overflowGuard is TRUE, we'll convert population sizes greater than some absurdly large population size to something that is still absurdly large, but that can be handled by the computer. Keeping populations below about 1 billion seems to prevent overflow.
    if (overflowGuard)
    {
      N[(N[,t] > 1e9), t] <- 1e9
    }
    # Now all new individuals are considered "past.residents" from the perspective of the next generation 
    past.residents <- N[,t]
  }
  
  return(N)
}