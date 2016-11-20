### Title: Generate extinction response variables for analysis
###
### Author: Michael Koontz
### Email: mikoontz@gmail.com
###
### Date Created: 20150501
### Last Updated: 20161120
###
### Purpose: Use long form data and tidy beetle functions to determine when populations went extinct, whether they were extinct 5 generations after the final introduction event, and whether they were extinct by the 9th generation. For use in logistic regression analysis.

# Clear your environment or restart R first
# rm(list=ls())

source('scripts/generate-tidy-data.R')

#-----------------------
#-----------------------

# Returns the number of generations before each population went extinct. (First time a census showed 0 individuals, AND no further introductions coming). This function returns an NA for an element corresponding to a population that never went extinct

# Arguments
#   N is the matrix of population sizes; the censuses values
#   intro.regime is a vector of propagule numbers. I.e. the number of introduction events for each population
#   gap is a logical vector specifying whether there was a gap in the introductions in generation 2
#   gens.censused is a numeric vector representing the total number of generations censused for each population. ID 1:480 were censused for 10 generations, ID 481:945 were censused for 9 generations

determine.when.extinct <- function(N, intro.regime, gap, gens.censused=c(rep(10, 480), rep(9, 437)) )
{
  Tf <- ncol(N)
  # Determine the first generation that each population went permanently extinct. 
  # Offset by 2 because there are really only Tf-1 generations and minus 1 more so that we get the first EXTINCT generation, not the last EXTANT generation 
  when.extinct <- Tf - (apply(N[,Tf:1], MARGIN=1, FUN=function(x) match(FALSE,x==0) ) - 2)
  
  # Note that if a population went to 0 each time you tried to introduce new individuals, the when.extinct value would say extinction happened too soon. In fact, the extinction isn't permanent until no more introductions will occur either. So the maximum of the when.extinct value and the number of introductions (plus 1 to look at the first GENERATION where no more introductions came).
  
  when.extinct <- pmax(when.extinct, (intro.regime + 1 + as.numeric(gap)))

# Determine populations that never went extinct. The when.extinct vector will be equal to gens.censused+1 when the population has a positive number of individuals at time (gens.censused)
# This is where it was important to know that some populations weren't censused as many times as others (9 vs. 10 generations)

  never.extinct <- which(when.extinct == (gens.censused+1))
  
  # Make the populations that never went extinct into NAs
  when.extinct[never.extinct] <- NA
  
  return(when.extinct)
}

#--------------------
#--------------------

# Function determines whether population was extinct x generations after its FINAL introduction event. Importantly, the timing of the final introduction event depends on how many introduction events there were, as well as whether there were any gaps in the introduction effort. Returns a logical vector with TRUE when population was extinct within x generations after final introduction and FALSE otherwise

# Arguments:
#   x: the number of generations after the final introduction event
#   N: the population census data (Ntp1)
#   intro.regime: the number of introduction events
#   gap: whether there was a gap in the introduction events in generation 2

extinct.x.after.intro <- function(x=5, N, intro.regime, gap)
{
  Tf <- ncol(N)
  # When did introductions finish for each rep?
  # Add x value to when introductions finished.
  # Essentially, the evaluation of (intro.regime + as.numeric(gap)) represents the first generation in the Ntp1 data frame where we can look and know that there won't be further introductions. For instance, for a 20x1 introduction (which can't have a gap), that expression evaluates to (1 + 0 = 1). We can treat the generation 1 column in the Ntp1 dataframe as the first generation where those individuals won't be augmented by future introductions. For a 5x4 introduction scenario with a gap in introduction in the 2nd generation, that expression evaluates to (4 + 1 = 5). We can treat the generation 5 column in the Ntp1 dataframe as the first generation where no further introductions will augment the population.
  # Since the evaluation of (intro.regime + as.numeric(gap)) represents the first generation to look at for post-introduction extinction, we want to look (x-1) generations later to ask the question, was the population extinct x generations after introductions were complete?
  gen.of.interest <- (intro.regime + x + as.numeric(gap)) - 1
  
  when.extinct <- determine.when.extinct(N, intro.regime, as.numeric(gap))
  
  # Determine populations that went extinct by the generation of interest
  extinct.x.after <- (when.extinct <= gen.of.interest)
  extinct.x.after[is.na(extinct.x.after)] <- FALSE
  
#   # Subset the N so that we only look at reps where enough time (Tf) had passed after the final introduction
#   analysis.extinct <- extinct[gen.of.interest <= Tf]
#   analysis.intro.regime <- intro.regime[gen.of.interest <= Tf]
#   analysis.gen <- gen.of.interest[gen.of.interest <= Tf]
  
  return(extinct.x.after)
}

#--------------------
#--------------------

# Function determines whether population was extinct x generations after the FIRST introduction event (i.e. the start of the experiment). Note that the timing of this question is independent of the introduction type or whether gaps were present in the introduction effort. Returns a logical vector with TRUE when population was extinct within x generations after final introduction and FALSE otherwise

#   N is the matrix of population sizes; the censuses values
#   gens.censused is a numeric vector representing the total number of generations censused for each population. ID 1:480 were censused for 10 generations, ID 481:945 were censused for 9 generations

extinct.after.x <- function(x=9, N, intro.regime, gap, gens.censused=c(rep(10, 480), rep(9, 437)))
{
  # Check if arguments imply the same number of populations
  if (length( unique(c(nrow(N), length(gens.censused), length(gap), length(intro.regime)) )) !=1)
    stop("Number of rows of N must be the same as the lengths of intro.regime, gap, and gens.censused. Non-matching lengths imply different numbers of populations.")
  
  when.extinct <- determine.when.extinct(N, intro.regime, gap, gens.censused)

  extinct <- rep(FALSE, times=nrow(N)) 
  
  # All populations that have a number designation for when they became extinct get a TRUE
  extinct[!is.na(when.extinct) & when.extinct <= x] <- TRUE
  
  return(extinct)
}

#--------------------
#--------------------
# Function returns summary statistics for extinction probability

# Arguments
#   extinct: a vector of TRUE or FALSE representing whether a population was extinct or not
#   INDEX: a list of grouping indices to aggregate by. Each element in the list must have the same length as the extinct vector

extinction.stats <- function(extinct, INDEX)
{
  df <- aggregate(x=list(percent.extinct=extinct), by=INDEX, FUN=function(x) c(length(which(x)), length(x), length(which(x))/length(x)) )
  
  groups <- length(INDEX)  
  labels <- names(df)[1:groups]
              
  df <- cbind(df[,1:groups], as.data.frame(df[,groups+1]))
  
  names(df)[1:groups] <- labels
  names(df)[(groups+1):(groups+3)] <- c("n[extinct]", "n[total]", "percent.extinct")
  return(df)
}

#--------------------
# Function returns a dataframe containing the unique ID, the introduction regime of the population, and the population size x generations after the final introduction.
# Function will break and return an error if x is too big (i.e. there weren't x generations of data past the final introdution event)
# This function will work on simulated data or on microcosm data as long as the IDs and gap (1 or 0) are specified.

N_x.after.intro <- function(x=1, N, intro.regime, ID=1:nrow(N), gap=numeric(nrow(N)), zeros=TRUE)
{
  if (x > (ncol(N) - max(intro.regime)))
    stop("You are looking too far after the final introduction scenario. Can't look this far for all treatments. Try decreasing x.")
  
  numbers <- unique(intro.regime)
  mat <- N*0 # Data frame of 0's the same dimensions as N
  
  # Iterate through the unique values of propagule number
  for (i in numbers)
  {
    # For population time series where the introduction regime is the current propagule number AND there's no introduction gap, the correct time point is i+x
    mat[(intro.regime == i & gap == 0), i+x] <- 1
    # For population time series where the introduction regime is the current propagule number AND there WAS an introduction gap, the correct time point is offset by 1 (i.e., i+x+1)
    if (any(intro.regime == i & gap == 1))
      mat[(intro.regime == i & gap == 1), i+x+1] <- 1
  }
  
  
  abundance <- apply(mat * N, 1, sum)
  df <- data.frame(ID, intro.regime, abundance)
  idx <- df$abundance == 0 | is.na(df$abundance)
  
  if (zeros)
    return(df)
  else
    return(df[!idx, ])
}

#--------------------


# b <- merge(df, abund[, c("ID", "abundance")], by="ID")
# b$number <- as.factor(b$number)
# b$block <- as.factor(b$block)
# b$gap <- as.factor(b$gap)

attributes <- read.csv('data/attributes.csv')
beetles <- read.csv('data/Tribolium-propagule-pressure-data.csv')

tidyb <- tidy.beetles(beetles = beetles)

N <- tidyb$Ntp1
Nt <- tidyb$Nt
b <- merge(N, attributes, by="ID")
census.columns <- which(colnames(b) %in% as.character(1:10))

col <- as.numeric(b$gap) + b$number
b$temp.extinctions <- 0
b$latest.extinct <- 0
b$loss <- 0

for (i in which(tidyb$migrants[, 2] != 20))
{
  temp <- sum(N[i, 2:col[i]] == 0, na.rm=TRUE)
  b$temp.extinctions[i] <- temp
  
  if (temp > 0) {
    b$latest.extinct[i] <- max(which(N[i, 2:col[i]] == 0))
  
    b$loss[i] <- (b$latest.extinct[i] - as.numeric(b$gap[i])) * b$size[i]
  }
}

#------------------------------------
# Alternate way to calculate temporary extinctions and lost inputs
#------------------------------------
# b$n.temp.extinctions <- 0
# b$loss <- 0
# 
# for (i in 1:nrow(tidyb$Ntp1))
# {
#   b$n.temp.extinctions[i] <- sum(tidyb$Ntp1[i, 2:col[i]] == 0)
#   b$loss[i] <- sum(tidyb$migrants[i,2:ncol(tidyb$migrants)]) - (match(0, tidyb$Ntp1[i, col[i]:2]) - 1) * b$size[i]
# }
# 
# b[is.na(b$loss), "loss"] <- 0
#-------------------------------------

#-------------------------------------
# Original code to get the abundance and establishment at a single time point, 5 generations after the final introduction event for each introduction regime
#--------------------------------------

# Example function run to determine which populations went extinct within 5 generations after their final introduction
# extinct5 <- extinct.x.after.intro(x=5, N=b[,census.columns], intro.regime=b$number, gap=b$gap)

# Example function run to determine which populations went extinct before the 9th generation of the experiment (this is the last generation that all populations were censused)
# extinct.after9 <- extinct.after.x(x=9, N=b[,census.columns], intro.regime=b$number, gap=b$gap)

# Example function run to get population abundance 5 generations after final introduction event.
# First convert NAs to 0's in the population size data.frame
# N_for_abundance <- b[, c(which(colnames(b) == "size"), census.columns)]
# N_for_abundance[is.na(N_for_abundance)] <- 0

# abund <- N_x.after.intro(x=5, N=N_for_abundance, intro.regime=b$number, ID=b$ID, gap=as.numeric(b$gap), zeros=TRUE)
# abund <- subset(abund, select=c(ID, abundance))

# b$extinct5after <- extinct5
# b$extinctafter9 <- extinct.after9
# b <- merge(b, abund, by="ID")

# extinction.stats(extinct5, INDEX=list(number=b$number))
# extinction.stats(extinct5, INDEX=list(env=b$env))
# extinction.stats(extinct5, INDEX=list(number=b$number, env=b$env))

#-----------------------------------------

# Example function run to determine when extinctions occurred
when.extinct <- determine.when.extinct(N=b[,census.columns], intro.regime=b$number, gap=b$gap)

b$when.extinct <- when.extinct
established <- as.data.frame(!sapply(1:9, FUN = extinct.after.x, N=b[,census.columns], intro.regime=b$number, gap=b$gap))
names(b)[2:11] <- paste0("N", names(b[2:11]))
names(established) <- paste0("established", 1:9)

b <- head(data.frame(b, established))


bb <- subset(b, select=c(ID, when.extinct, extinct5after, extinctafter9, temp.extinctions, latest.extinct, loss, abundance))

log.regression <- merge(attributes, bb, by="ID")

# write.csv(log.regression, 'figshare/extinctions.csv', row.names=FALSE)
