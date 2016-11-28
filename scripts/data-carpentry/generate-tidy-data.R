# Title: generate tidy data
# 
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150330
# Last Updated: 20150331

# This function takes the entered Tribolium flour beetle data from the "Eco-evolutionary consequences of multiple introductions" main experiment (which is in long form), and puts it in a 2 dimensional form to more easily represent the time series. Each population is defined by a unique ID, which is included in each of the produced dataframes here. This makes for easy merging (using ID as a key) with the 'attributes.csv' file (which includes block number, treatment types, whether a gap in the introduction occurred, whether the populations experienced the reduced incubator moisture during generation 2, etc.)
# Note this code also corrects 12 population trajectories that had their population drop to 1 individual, but rebounded due to egg contamination (eggs on a lab surface getting into the microcosm). We coerced those populations to be extinct after they dropped in size to 1 individual.
# Affected populations: 13, 45, 87, 98, 303, 362, 500, 523, 640, 758, 777, 825

# Requires the tidyr package for reshaping the data.
# Input is the long-form entered data in a dataframe object type.
# Returns a list of dataframes representing the different values that are unique to each ID/Generation combination.
# Further manipulations on these data frames are possible using other functions. For instance, cbind(Ntp1$ID, Ntp1[,2:11]/Nt[,2:11] would be equivalent to a dataframe of lambda values for each time step.

# Load tidyr library
library(tidyr)

tidy.beetles <- function(beetles, deal_with_loners=TRUE)
{
#---------- 
# N[t+1] dataframe
# ---------

# The final column has NA for ID 481:945 because only blocks 1 and 2 were set up at the end of generation 9 to yield N[t+1] data for generation 10

# Subset to relevant parts
b <- subset(beetles, select=c(ID, Generation, Census))
# Spread it
Ntp1 <- spread(b, Generation, Census)
names(Ntp1) <- c("ID", paste0("N", 0:9, "plus1"))
# Check it
head(Ntp1)
tail(Ntp1)

#---------- 
# N[t] dataframe
# ---------

# The final column has NA for ID 481:945 because only blocks 1 and 2 were set up at the end of generation 9 to yield N[t+1] data for generation 10

b <- subset(beetles, select=c(ID, Generation, N0))
Nt <- spread(b, Generation, N0)
names(Nt) <- c("ID", paste0("N", 0:9))
head(Nt)
tail(Nt)

#---------- 
# Migration dataframe
# ---------

b <- subset(beetles, select=c(ID, Generation, Addition))
migrants <- spread(b, Generation, Addition)
# Add the initial introduction onto the dataframe
migrants$'0' <- Nt$N0
# Reorder dataframe columns and leave off final generation addition, since those populations weren't even set up for another generation
total.columns <- ncol(migrants)
migrants <- migrants[, c(1, total.columns, 2:(total.columns-2))]
names(migrants) <- c("ID", paste0("migrants", 0:9))
head(migrants, 120)
tail(migrants)
      
#---------- 
# Environment % standard media mixture dataframe
# ---------

b <- subset(beetles, select=c(ID, Generation, Environment))
environment <- spread(b, Generation, Environment)
names(environment) <- c("ID", paste0("env", 0:9))

head(environment)

#---------- 
# Census taker dataframe
# ---------

b <- subset(beetles, select=c(ID, Generation, Person))
person <- spread(b, Generation, Person)
names(person) <- c("ID", paste0("person", 0:9))
head(person)
tail(person)

#---------- 
# Setup order dataframe
# ---------

b <- subset(beetles, select=c(ID, Generation, Setup.Order))
setup.order <- spread(b, Generation, Setup.Order)
names(setup.order) <- c("ID", paste0("setup.order", 0:9))
head(setup.order)
tail(setup.order)

#---------
# Make populations of size 1 go extinct
#---------
loners_df <- data.frame()
if (deal_with_loners)
{
# Due to some egg contamination, populations of size 1 sometimes persisted. Make these extinct with this code.
# Two conditions must be met: Population had 1 individual, no more introductions were coming
# Make next generation size 0, then make the rest of the columns NA
# Affected population IDs: 13, 45, 87, 98, 303, 500, 523, 640, 758

# First determine when introductions stopped
# Total columns in migrants data frame
columns <- ncol(migrants)

# By rows, look for the index of the first non-zero number of migrants.
idx <- apply(migrants[, (columns:2)], MARGIN=1, FUN=function(x) match(FALSE, x==0))

# The total columns minus this index value represents the last generation whose Ntp1 census data was influenced by migrants
last.migrants <- columns - idx

loner.col <- rep(NA, nrow(Ntp1))

# For each row, determine the index for the first generation where the population has only 1 individual, but only look at the columns AFTER no more migrants arrive (that is, the (last.migrants + 1) column until the final column)
# We have to add back in the columns that we skipped (last.migrants columns) to this index value
# Also, adding 1 to the last.migrants value accounts for the "ID" column
for (i in 1:nrow(Ntp1))
{
  loner.col[i] <- match(1, Ntp1[i, (last.migrants[i] + 1):columns]) + last.migrants[i]
}

# Indices of populations that had 1 individual after introductions were finished (includes populations that actually did go extinct as well as populations that should have gone extinct)
loner.idx <- which(!is.na(loner.col))

# Column of Ntp1 dataframe for each of the above indices that should be 0
extinct.col <- loner.col[loner.idx] + 1
# Note that the equivalent column in Nt is extinct.col+1

# Column of Ntp1 dataframe for each of the loner indicies that should be NA. All columns after this should be NA as well.
NA.col <- extinct.col + 1
counter <- 0
running_sum <- 0
culprits <- numeric(length(loner.idx))

for (j in 1:length(loner.idx))
{
  # Only proceed if loner population wasn't in the final generation.
  if (extinct.col[j] <= columns)
  {
    # Only proceed if population didn't go extinct when it should have
    if (as.numeric(Ntp1[loner.idx[j], extinct.col[j]] > 0))
      {
        # Update counter of number of culprit populations
        counter <- counter + as.numeric(Ntp1[loner.idx[j], extinct.col[j]] > 0)
        # Add total Ntp1 to running sum
        running_sum <- running_sum + Ntp1[loner.idx[j], extinct.col[j]]
        # Flag this population as a culprit
        culprits[j] <- 1
        
        #  Force the appropriate Ntp1 time point to extinction
        Ntp1[loner.idx[j], extinct.col[j]] <- 0
        
        # Only proceed if there are still generations after extinction that need to be forced to NAs
        # This also covers proceeding if there are columns in the Nt dataframe that need to be forced to 0
        if (NA.col[j] <= columns)
        {
          # Force appropriate Nt columns to 0
          Nt[loner.idx[j], (extinct.col[j]+1):columns] <- 0
          # Force appropriate Ntp1 columns to NA
          Ntp1[loner.idx[j], NA.col[j]:columns] <- NA
        }
      } # End check on whether loners didn't go extinct
  } # End check on whether loners occured before end of experiment
} # End for loop iterating through each loner population
loners_df <- data.frame(loner.idx, culprits, extinct.col)

} # End if statement about making loners go extinct

return(list(Nt=Nt, Ntp1=Ntp1, migrants=migrants, environment=environment, person=person, setup.order=setup.order, loners_df=loners_df))
}
