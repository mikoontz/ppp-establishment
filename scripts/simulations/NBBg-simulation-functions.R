# Title: NBBg simulate functions
# 
# Functions for simulating NBBg dynamics including all forms of demographic stochasticity AND model undertainty in parameter distributions

#----------
# Functions
#----------

source("scripts/simulations/NBBg-population-dynamics-function.R")

# Returns a data.frame with the possible combinations of propagule number and propagule size to reach a fixed total number of individuals to introduce.

get.propagule.pressure <- function(total.to.introduce, Tf, excludeMaxNum=TRUE)
{
  # Determine the possible combinations of propagule number and propagule size. Only evenly divided combinations for now.
  # Propagule number and size are the index number where the modulus of total.to.introduce and potential number/size is 0
  propagule.number <- which(total.to.introduce %% 1:total.to.introduce == 0)
  propagule.size <- rev(propagule.number)
  
  propagule.pressure <- data.frame(number=propagule.number, size=propagule.size)
  
  # Take away the possibility of 1 individual being introduced in each time step
  propagule.pressure <- propagule.pressure[propagule.pressure$size!=1, ]
  
  # Take away the possibility of more introduction events than time steps
  propagule.pressure <- subset(propagule.pressure, subset=number<=(Tf-1))
  
  # Exclude maximum number of introductions if you like (aside from 1xN, which is automatically excluded since we model a sexual species. For instance, no 2 individuals by 10 introductions when the total.to.introduce is 20
  if (excludeMaxNum)
    propagule.pressure <- propagule.pressure[-nrow(propagule.pressure), ]
  
  return(propagule.pressure)
}

# This function sets up storage matrices, determines the migrants matrix, and initiates the first generation of the population. Returns a number of set up variables.
set.up <- function(total.to.introduce, Tf, reps, excludeMaxNum=TRUE)
{
  # First get different propagule number/size combinations
  propagule.pressure <- get.propagule.pressure(total.to.introduce, Tf, excludeMaxNum)
  
  # Set up the number of migrants per generation based on combinations of propagule number and size
  total.reps <- dim(propagule.pressure)[1]*reps
  migrants <- matrix(0, nrow=total.reps, ncol=Tf)
  
  # Initialize the first time step's worth of migrants
  migrants[, 1] <- rep(propagule.pressure$size, each=reps)
  
  # Fill in the rest of the migrants matrix
  for (i in 2:dim(propagule.pressure)[1])
  {
    # Only fill in rows where introductions are still forthcoming
    # Possible rows are 1:total number of reps
    all.rows <- 1:total.reps
    # Rows we want to fill in are total rows minus the sets of 'reps' populations who finished receiving introductions
    rows.of.interest <- all.rows[-(1:((i-1)*reps))]
    migrants[rows.of.interest, 1:propagule.pressure$number[i]] <- rep(propagule.pressure$size[i], times=reps)
  }
  
  # Start with no resident beetles. That is, no unmated beetles. All initial introductions come from external population.
  past.residents <- numeric(total.reps)
  
  #Initial population size is the number of migrants in the first generation
  N <- matrix(NA, nrow=total.reps, ncol=Tf)
  N[,1] <- migrants[,1] + past.residents
  
  return(list(N=N, migrants=migrants, past.residents=past.residents, total.reps=total.reps))
}

# This function extracts the N0 values for each generation
get_Nt <- function(N, migrants)
{
  Nt <- N
  Nt[, 2:ncol(Nt)] <- N[, 2:ncol(Nt)] + migrants[, 2:ncol(migrants)]
  return(Nt[, -ncol(Nt)])
}

get_Ntp1 <- function(N)
{
  return(N[, -1])
}


# This function will determine the number of introductions for each population. Used internally in other functions, mostly.
# Inputs: matrix of population sizes (only the first time step needs to be filled in). Could also pass in the migrants matrix.
# Output: vector of propagule numbers-- number of introductions for each population
# Use match to find the row in propagule.pressure that matches each of the initial generation numbers, then access those rows and the "number" column of propagule.pressure to find out how many generations the introductions occurred. 
get.intro.regime <- function(N, propagule.pressure)
{
  intro.regime.idx <- match(N[, 1], propagule.pressure$size)
  intro.regime <- propagule.pressure[intro.regime.idx, "number"]
  
  return(intro.regime)
}

### Convenience wrapper that calls required functions above with desired parameters and returns a short form data frame of population sizes. Also includes ability to save objects for later retrieval.
### Takes in parameters required to run a NBBg (see Melbourne & Hastings 2008) model for Tf generations while introducing a fixed total number of individuals in different combinations of introduction regime.

simNBBg <- function(samps, reps, Tf, total.to.introduce = 20, p = 0.5, excludeMaxNum = TRUE, save_objects = FALSE, data_descriptor = "", dirpath = "data/simulations/", overflowGuard = FALSE)
{
  # Derive setup variables and storage objects using basic arguments passed to simNBBg
  setup <- set.up(total.to.introduce, Tf, reps, excludeMaxNum = excludeMaxNum)
  
  # Run the actual simulation
  N <- NBBg(Tf=Tf, total.reps = setup$total.reps, past.residents = setup$past.residents, p = p, N = setup$N, migrants = setup$migrants, samps = samps, overflowGuard = overflowGuard)
  
  # Save N and migrants data frame or not
  if (save_objects) {
    if (!dir.exists(dirpath)) {
      dir.create(dirpath)
    }
    if (!dir.exists(paste0(dirpath, "/N"))) {
      dir.create(paste0(dirpath, "/N"))
    }
    if (!dir.exists(paste0(dirpath, "/migrants"))) {
      dir.create(paste0(dirpath, "/migrants"))
    }
    
    write.csv(N, paste0(dirpath, "/N", data_descriptor, ".csv"), row.names=FALSE)
    write.csv(setup$migrants, paste0(dirpath, "/migrants", data_descriptor, ".csv"), row.names=FALSE)
  }
  return(N)
  
}