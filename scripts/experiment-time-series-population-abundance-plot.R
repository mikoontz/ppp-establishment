### Experiment time series plot
###
### Plots population abundance of extant populations and 1 standard error margin for each generation in the parsing propagule pressure experiment

# Clear environment if necessary
# rm(list=ls())

source("scripts/generate-tidy-data.R")

beetles <- read.csv("data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("data/attributes.csv")  

# b <- read.csv("Clean-Data/extinctions.csv")

library(viridis)

tidy.b <- tidy.beetles(beetles=beetles)

Nt <- tidy.b$Nt
Ntp1 <- tidy.b$Ntp1
migrants <- tidy.b$migrants
env <- tidy.b$environment
