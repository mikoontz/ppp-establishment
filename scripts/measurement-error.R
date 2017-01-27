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
library(dplyr)
library(tidyr)

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

# Coefficient of variation of observer error per population (standard deviation
# of sum of squares residuals per popoulation per observer divided by group mean
# per population)
CoV <- apply(ss_resids_df, MARGIN = 2, FUN = function(x) sqrt(x)) %>%
  apply(MARGIN = 1, FUN = function(y) y / group.mn) %>%
  t()

pct_err_long <- CoV %>%
  as.data.frame() %>%
  mutate(observer = rownames(.)) %>%
  gather(key = "box", value = "pct_err", -observer) %>%
  mutate(size_cat = "s") %>%
  mutate(size_cat = ifelse(box %in% 1:5, yes = "s", no = ifelse(box %in% 6:10, yes = "m", no = "l"))) %>%
  mutate(group_mn = group.mn[box]) # Add the group means for each box

# Treat the population size as a categorical predictor with 3 levels (small, medium, large)
fm1 <- lm(pct_err ~ observer + size_cat, data = pct_err_long)
summary(fm1)

# Treat the group means as a continuous predictor
fm2 <- lm(pct_err ~ observer + group_mn, data = pct_err_long)
summary(fm2)

# Conclusion: No difference in measurement error amongst observers, no effect of
# population size on measurement error %

# Average percent difference for any one observer and any one population was 0.26%
fm3 <- lm(pct_err ~ 1, data = pct_err_long)
summary(fm3)
