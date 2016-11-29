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

source('scripts/data-carpentry/generate-tidy-data.R')

##### Helper functions #####

#### determine.when.extinct() ####
#### function description ####
# Returns the number of generations before each population went extinct. (First time a census showed 0 individuals, AND no further introductions coming). This function returns an NA for an element corresponding to a population that never went extinct

# Arguments
#   N is the matrix of population sizes in the NEXT generation; the censuses values representing Ntp1 (and not including the initial introduction)
#   intro.regime is a vector of propagule numbers. I.e. the number of introduction events for each population
#   gap is a logical vector specifying whether there was a gap in the introductions in generation 2
#   gens.censused is a numeric vector representing the total number of generations censused for each population. ID 1:480 were censused for 10 generations, ID 481:945 were censused for 9 generations
#### function definition ####
determine.when.extinct <- function(N, intro.regime, gap, gens.censused = c(rep(10, 480), rep(9, 437)), Tf = ncol(N))
{
  # Determine the first generation that each population went permanently extinct. 
  # Offset by 2 because there are really only Tf-1 generations of Ntp1 censuses and minus 1 more so that we get the first EXTINCT generation, not the last EXTANT generation 
  when.extinct <- (Tf - apply(N[,Tf:1], MARGIN=1, FUN=function(x) match(FALSE, x == 0)) + 2) 
  when.extinct[is.na(when.extinct)] <- 1
    # Note that if a population went to 0 each time you tried to introduce new individuals, the when.extinct value would say extinction happened too soon. In fact, the extinction isn't permanent until no more introductions will occur either. So the maximum of the when.extinct value and the number of introductions (plus 1 to look at the first GENERATION where no more introductions came).
  
  when.extinct <- pmax(when.extinct, (intro.regime + as.numeric(gap)))

# Determine populations that never went extinct. The when.extinct vector will be equal to gens.censused+1 when the population has a positive number of individuals at time (gens.censused)
# This is where it was important to know that some populations weren't censused as many times as others (9 vs. 10 generations)
  never.extinct <- which(when.extinct == (gens.censused + 1))
  
  # Make the populations that never went extinct into NAs
  when.extinct[never.extinct] <- NA
  
  return(when.extinct)
}

#### simulation_stats() ####
#### function description ####
# Convenience wrapper for getting establishment and population size response variables from experiment data or simulation
#### function definition ####
simulation_stats <- function(time_points, FUN, col_names, ...) {
  values <- sapply(X = time_points, 
                    FUN = FUN,
                    ...)
  
  mean_by_regime <- apply(X = values, 
                          MARGIN = 2, 
                          FUN = function(x) tapply(X = x, 
                                                   INDEX = intro.regime, 
                                                   FUN = mean,
                                                   na.rm = TRUE))
  
  mean_by_regime <- as.data.frame(mean_by_regime)
  colnames(mean_by_regime) <- col_names
  
  return(mean_by_regime)
}

# Function returns summary statistics for extinction probability

# Arguments
#   extinct: a vector of TRUE or FALSE representing whether a population was extinct or not
#   INDEX: a list of grouping indices to aggregate by. Each element in the list must have the same length as the extinct vector
# extinction.stats <- function(extinct, INDEX)
# {
#   df <- aggregate(x=list(percent.extinct=extinct), by=INDEX, FUN=function(x) c(length(which(x)), length(x), length(which(x))/length(x)) )
#   
#   groups <- length(INDEX)  
#   labels <- names(df)[1:groups]
#               
#   df <- cbind(df[,1:groups], as.data.frame(df[,groups+1]))
#   
#   names(df)[1:groups] <- labels
#   names(df)[(groups+1):(groups+3)] <- c("n[extinct]", "n[total]", "percent.extinct")
#   return(df)
# }

##### Response variable functions #####

#### extinct.x.after.intro() ####
#### function description ####
# Function determines whether population was extinct x generations after its FINAL introduction event. Importantly, the timing of the final introduction event depends on how many introduction events there were, as well as whether there were any gaps in the introduction effort. Returns a logical vector with TRUE when population was extinct within x generations after final introduction and FALSE otherwise

# Arguments:
#   x: the number of generations after the final introduction event
#   N: the population census data (Ntp1)
#   intro.regime: the number of introduction events
#   gap: whether there was a gap in the introduction events in generation 2
#   gens.censused: The final generation of census for each population
#### function definition ####
extinct.x.after.intro <- function(x=5, N, intro.regime, gap, gens.censused=c(rep(10, 480), rep(9, 437)))
{
  if (x > (ncol(N) - max(intro.regime) + 1))
    stop("You are looking too far after the final introduction scenario. Can't look this far for all treatments. Try decreasing x.")
  
  Tf <- ncol(N)
  # When did introductions finish for each rep?
  # Add x value to when introductions finished.
  # Essentially, the evaluation of (intro.regime + as.numeric(gap)) represents the first generation in the Ntp1 data frame where we can look and know that there won't be further introductions. For instance, for a 20x1 introduction (which can't have a gap), that expression evaluates to (1 + 0 = 1). We can treat the generation 1 column in the Ntp1 dataframe as the first generation where those individuals won't be augmented by future introductions. For a 5x4 introduction scenario with a gap in introduction in the 2nd generation, that expression evaluates to (4 + 1 = 5). We can treat the generation 5 column in the Ntp1 dataframe as the first generation where no further introductions will augment the population.
  # Since the evaluation of (intro.regime + as.numeric(gap)) represents the first generation to look at for post-introduction extinction, we want to look (x-1) generations later to ask the question, was the population extinct x generations after introductions were complete?
  gen.of.interest <- (intro.regime + x + as.numeric(gap)) - 1
  
  when.extinct <- determine.when.extinct(N, intro.regime, as.numeric(gap), gens.censused = gens.censused)
  
  # Determine populations that went extinct by the generation of interest
  extinct.x.after <- (when.extinct <= gen.of.interest)
  extinct.x.after[is.na(extinct.x.after)] <- FALSE
  
  #   # Subset the N so that we only look at reps where enough time (Tf) had passed after the final introduction
  #   analysis.extinct <- extinct[gen.of.interest <= Tf]
  #   analysis.intro.regime <- intro.regime[gen.of.interest <= Tf]
  #   analysis.gen <- gen.of.interest[gen.of.interest <= Tf]
  
  return(extinct.x.after)
}

#### extinct.after.x() ####
#### function description ####
# Function determines whether population was extinct x generations after the FIRST introduction event (i.e. the start of the experiment). Note that the timing of this question is independent of the introduction type or whether gaps were present in the introduction effort. Returns a logical vector with TRUE when population was extinct within x generations after final introduction and FALSE otherwise

#   N is the matrix of population sizes; the censuses values
#   gens.censused is a numeric vector representing the total number of generations censused for each population. ID 1:480 were censused for 10 generations, ID 481:945 were censused for 9 generations
#### function definition ####
extinct.after.x <- function(x=9, N, intro.regime, gap, gens.censused=c(rep(10, 480), rep(9, 437)))
{
  if (x > (ncol(N) + 1))
    stop("You are looking too far after the final introduction scenario. Can't look this far for all treatments. Try decreasing x.")
  
  # Check if arguments imply the same number of populations
  if (length( unique(c(nrow(N), length(gens.censused), length(gap), length(intro.regime)) )) !=1)
    stop("Number of rows of N must be the same as the lengths of intro.regime, gap, and gens.censused. Non-matching lengths imply different numbers of populations.")
  
  when.extinct <- determine.when.extinct(N, intro.regime, gap, gens.censused)
  
  extinct <- rep(FALSE, times=nrow(N)) 
  
  # All populations that have a number designation for when they became extinct get a TRUE
  extinct[!is.na(when.extinct) & when.extinct <= x] <- TRUE
  
  return(extinct)
}

#### N_x.after.intro() ####
#### function description ####
# Function returns a dataframe containing the unique ID, the introduction regime of the population, and the population size x generations after the final introduction.
# Function will break and return an error if x is too big (i.e. there weren't x generations of data past the final introdution event)
# This function will work on simulated data or on microcosm data as long as the IDs and gap (1 or 0) are specified.
# The 'zeroes' argument lets 0's be included instead of NAs. Default is FALSE so they don't affect population mean calculations
#### function definition ####
N_x.after.intro <- function(x = 1, Ntp1, intro.regime, ID = 1:nrow(N), gap = numeric(nrow(N)), zeros = FALSE)
{
  if (x > (ncol(Ntp1) - max(intro.regime) + 1))
    stop("You are looking too far after the final introduction scenario. Can't look this far for all treatments. Try decreasing x.")
  
  numbers <- unique(intro.regime)
  mat <- Ntp1 * 0 # Data frame of 0's the same dimensions as N
  # Iterate through the unique values of propagule number
  
  for (i in numbers)
  {
    # For population time series where the introduction regime is the current propagule number AND there's no introduction gap, the correct time point is i+x-1 since we are using the Ntp1 data frame
    mat[(intro.regime == i), i + x - 1 ] <- 1
    # For population time series where the introduction regime is the current propagule number AND there WAS an introduction gap, the correct time point is offset by 1 (i.e., i+x+1)
    if (any(intro.regime == i & gap == 1))
      mat[(intro.regime == i & gap == 1), i + x] <- 1
  }
  
  abundance <- apply(mat * Ntp1, 1, sum, na.rm = TRUE)
  idx <- abundance == 0
  abundance_no_zeroes <- abundance
  abundance_no_zeroes[idx] <- NA 
  
  if (zeros)
      return(abundance)
    else
      return(abundance_no_zeroes)
}

#### Calculate temporary extiction response variables ####
#### temp_extinction() ####
#### function description ####
# Calculates number of times a population went extinct before being demographically rescued by additional inputs
#### function definition ####
temp_extinction <- function(Ntp1, intro.regime, gap = rep(FALSE, nrow(Ntp1))) {
  final_intro_column <- as.numeric(gap) + intro.regime
  num_temp_extinctions <- 
    sapply(X = 1:length(final_intro_column), 
           FUN = function(j) sum(Ntp1[j, 2:final_intro_column[j]] == 0, na.rm = TRUE))
  return(num_temp_extinctions)
}

#### loss_from_temp_extinction() ####
#### function description ####
# Calculates number of lost inputs when a population goes temporarily extinct. Populations that went temporarily extinct when there were fewer additional inputs scheduled to arrive "lost" more inputs than those populations that went temporarily extinct and were revived earlier.
#### function definition ####
loss_from_temp_extinction <- function(Ntp1, migrants, intro.regime, gap = rep(FALSE, nrow(Ntp1))) {
  # Column representing final introduction
  final_intro_column <- as.numeric(gap) + intro.regime
  # Initialize latest temporary extinction generation at 0
  latest_temp_extinct <- rep(0, nrow(Ntp1))
  # Use temp_extinction() function to get whether or not extinctions occurred
  temp <- temp_extinction(Ntp1 = Ntp1, intro.regime = intro.regime, gap = gap)
  # Population indicies with temporary extinctions
  idx <- which(temp > 0)
  
  # Calculate the latest temporary extinction by looking at each row of Ntp1 and the census columns during which introductions were still occurring, counting how many generations had a Ntp1 == 0
  latest_temp_extinct[idx] <- sapply(X = idx, 
                                     FUN = function(j) max(which(Ntp1[j, 2:final_intro_column[j]] == 0)))
  
  # Need to ignore the ID column, so the column representing the latest temporary extinction is 1 + the generation of the latest temporary extinction
  latest_temp_extinct_column <- latest_temp_extinct + 1
  
  # Initialize amount of loss at 0 
  loss <- rep(0, nrow(Ntp1))
  
  # Calculate loss as the sum of migrants up to the latest generation with a temporary extinction
  loss[idx] <- sapply(X = idx,
                      FUN = function(j) sum(migrants[j, 2:latest_temp_extinct_column[j]])) 
  
  return(loss)
}
#### Main script ####

attributes <- read.csv('data/attributes.csv')
beetles <- read.csv('data/Tribolium-propagule-pressure-data.csv')

tidyb <- tidy.beetles(beetles = beetles)

Ntp1 <- tidyb$Ntp1 # Represents population abundance tp1 generations after start of experiment
Nt <- tidyb$Nt
migrants <- tidyb$migrants

b <- merge(Ntp1, attributes, by="ID")
census.columns <- grep(pattern = "[0-9]", x = colnames(b)) # All columns with a number in them are census columns (representing Ntp1)




#### Add response variables to data frame ####
b$temp.extinctions <- temp_extinction(Ntp1 = Ntp1, intro.regime = b$number, gap = b$gap)

b$loss <- loss_from_temp_extinction(Ntp1 = Ntp1, migrants = migrants, intro.regime = b$number, gap = b$gap)

b$when.extinct <- determine.when.extinct(N = b[ ,census.columns], intro.regime = b$number, gap = b$gap)

extant <- as.data.frame(!sapply(1:9, 
                                FUN = extinct.after.x, 
                                N = b[, census.columns], 
                                intro.regime = b$number, 
                                gap = b$gap))
names(extant) <- paste0("extant", 1:9)

b <- data.frame(b, extant)

# Also include the relative extinction calculation in the final data set. These columns represent whether the population was extant or not 1, 2, 3, 4, or 5 generations after the introduction regime completed.
extant_x_after <- as.data.frame(!sapply(1:5, 
                                        FUN = extinct.x.after.intro, 
                                        N = b[,census.columns], 
                                        intro.regime = b$number, 
                                        gap = b$gap))

names(extant_x_after) <- paste0("extant_", 1:5, "_after")

b <- data.frame(b, extant_x_after)

# Also include the relative abundance calculations in the final data set. These columns represent the abundance of the population 1, 2, 3, 4, or 5 generations after the final introduction event for a population's partiuclar introduction regime.

abund_x_after <- as.data.frame(sapply(1:5, 
                                      FUN = N_x.after.intro, 
                                      N = b[, census.columns], 
                                      intro.regime = b$number, 
                                      ID = b$ID, 
                                      gap = as.numeric(b$gap)))

names(abund_x_after) <- paste0("N_", 1:5, "_after")

b <- data.frame(b, abund_x_after)

#### Finally, convert all Ntp1 censuses of 0 to NA for the purposes of the response variable ####
b[, census.columns] <- sapply(b[, census.columns], function(j) replace(j, which(j == 0), NA))


bb <- subset(b, select= -c(block, color, number, size, environment, special, gap, notes, drought))



clean_establishment_data <- merge(attributes, bb, by="ID")
head(clean_establishment_data)
tail(clean_establishment_data, 10)

# write.csv(clean_establishment_data, 'data/clean-establishment-data.csv', row.names=FALSE)
