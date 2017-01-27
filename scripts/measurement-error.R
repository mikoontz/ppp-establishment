### Title: measurement error
###
### Author: Michael Koontz
### 
### Date Created: 20140604
### Last Updated: 20140604
###
### Purpose: Evaluates the measurement error in censusing using repeat measures of 15 populations of unknown size

setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography")
error <- read.csv("Entered Data/csv/Measurement Error/20140602 measurement error.csv")

head(error)

par(mar=c(5,5,4,3)+0.1)
plot(error$Box, error$Census, pch=19)

variance <- tapply(error$Census, error$Box, FUN=function(x) var(x, na.rm=TRUE))

group.mn <- tapply(error$Census, error$Box, FUN=function(x) mean(x, na.rm=TRUE))
group.mn

resids <- tapply(error$Census, error$Person, FUN=function(x) ( ((x-group.mn)^2) ) )
resids

df <- do.call(rbind, resids)
mean(apply(df, 2, function(x) (mean(x, na.rm=TRUE))^(1/2) ) / group.mn)
apply(df, 1, function(x) sd(x/group.mn, na.rm=TRUE))

lapply(resids, FUN=function(x) sum(x^2, na.rm=TRUE))

var(variance[1:5])
var(variance[6:10])
var(variance[10:15])