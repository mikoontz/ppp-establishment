### Title: measurement error
###
### Author: Michael Koontz
### 
### Date Created: 20140604
### Last Updated: 20170127
###
# Purpose: Evaluates the measurement error in censusing using repeat measures 
# of 15 populations of unknown size

# Fifteen populations were set up with unknown, but variable population sizes 
# and 8 different observers censused each population. Five populations were
# "small", 5 populations were "medium", and 5 populations were "large"

# Read file
error <- read.csv("data/measurement-error.csv")
error$size_cat <- c(rep("s", 5), rep("m", 5), rep("l", 5))
head(error)

# Plot the range of census values
plot(x = rep(1, times = nrow(error)), y = error$Census, pch = 19)
boxplot(error$Census, add = TRUE, col = NA)

# Mean and variance of censused population sizes. If variance is 0, then all 
# observers recorded the same population size 
(group.mn <- tapply(error$Census, error$Box, FUN=function(x) mean(x, na.rm=TRUE)))
(variance <- tapply(error$Census, error$Box, FUN=function(x) var(x, na.rm=TRUE)))

# The squared difference between each observer's census and the group mean for
# each population
(ss_resids <- tapply(X = error$Census, 
                  INDEX = error$Person, 
                  FUN=function(x) ( ((x - group.mn)^2) ) ))

ss_resids_df <- do.call(rbind, ss_resids)

# Coefficient of variation of observer error per population (standard deviation of
# sum of squares residuals per popoulation divided by group mean per population)
CoV <- apply(ss_resids_df, MARGIN = 2, FUN = function(x) mean(sqrt(x), na.rm = TRUE)) / group.mn
mean(CoV)
  
fm1 <- lm(Census ~ Person + size_cat, data = error)
summary(fm1)
