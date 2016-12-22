### Experiment time series plot
###
### Plots population abundance of extant populations and 1 standard error margin for each generation in the parsing propagule pressure experiment

# Clear environment if necessary
# rm(list=ls())

source("scripts/data-carpentry/generate-tidy-data.R")

beetles <- read.csv("data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("data/attributes.csv")  

# b <- read.csv("Clean-Data/extinctions.csv")

library(viridis)

tidy.b <- tidy.beetles(beetles=beetles)

Nt <- tidy.b$Nt
Ntp1 <- tidy.b$Ntp1
migrants <- tidy.b$migrants
env <- tidy.b$environment

series <- merge(Nt, Ntp1)

# The vector of the order of columns to make it true: N0, N0p1, N1, N1p1, N2, N2p1, etc.
col_order <- numeric(length = ncol(series))
col_order[1:2] <- 1:2
col_order[seq(3, length(col_order), by=2)] <- 12:21
col_order[seq(4, length(col_order), by=2)] <- 3:11  

# Rearrange the time series data frame to put columns in proper order for plotting
series <- series[, col_order]
series <- series[, -ncol(series)] # Remove the final column, since we don't need to plot what the 10th generation would have began as

# Look at all rows in series and replace everything after an NA with all NAs row-wise. This will ensure no plotting occurs after extinction
extinctToNA <- function(x) {
  if (any(is.na(x))) { # If any of the time series is NA (indicating extinction)
    col <- min(which(is.na(x))) # Find the first index value where it went extinct
    x[col:length(x)] <- NA # Change all of the generations after that point to NA
  }
  return(t(x)) # Return the transposed original or modified x to keep same format with each row being a different population
}
# Apply the function rowwise
time_series <- t(apply(series, MARGIN = 1, FUN= extinctToNA))
colnames(time_series) <- colnames(series)

# Add attribute data to time series
data <- merge(time_series, attributes)

# Subset data to eliminate populations experiencing an introduction gap
data <- data[data$gap == FALSE, ]

# We need to calculate the mean and standard error of population size by generation and introduction regime. So we need mean abundance, standard deviation of abundance, and sample size. We'll use calculations that exclude populations that have been extinct (represented as NAs in the data) and populations that JUST went extinct (represented as zeroes before becoming NAs in the next generation).

means_no_zeroes <- apply(data[, 2:20], 
                         MARGIN= 2, 
                         FUN= function(x) tapply(x, 
                                                 INDEX= data$number, 
                                                 FUN= function(y) mean(y[y>0],  na.rm= TRUE)))

counts_no_zeroes <- apply(data[, 2:20], 
                          MARGIN= 2, 
                          FUN= function(x) tapply(x, 
                                                  INDEX= data$number, 
                                                  FUN= function(y) length(which(!is.na(y) & y > 0))))

sigma_no_zeroes <- apply(data[, 2:20], 
                         MARGIN= 2, 
                         FUN= function(x) tapply(x, 
                                                 INDEX= data$number, 
                                                 FUN= function(y) sd(y, na.rm= TRUE)))

##### Calculate the standard error and time steps ####
se_no_zeroes <- sigma_no_zeroes / sqrt(counts_no_zeroes)

time <- c(0, rep(1:9, each=2))

#### Make time series plot ####
# pdf("figures/experiment-time-series-population-abundance.pdf", height=5, width=5)
# postscript("figures/experiment-time-series-population-abundance.eps", family = "Times", height=5, width=5)

tiff("figures/experiment-time-series-population-abundance.tif", units= "in", res= 300, height=5, width=6)
par(mar = c(4.1, 4.1, 1, 2), family = "Helvetica")
matplot(x = time, y = t(time_series[, 2:ncol(time_series)]), 
        type = "l", 
        lty = "solid", 
        col = adjustcolor("grey", alpha.f= 0.1), 
        las = 1, 
        xlab = "Generation", 
        ylab = "Population size", 
        xaxt = "n", 
        ylim = c(0,80))

legend("topleft", 
       legend = c("20 individuals introduced in 1 event", "10x2", "4x5", "5x4"), 
       bty = "n", 
       col = viridis(4), 
       lty = 1,
       lwd = 3)

axis(side = 1, 
     at = 0:9, 
     labels = c("P", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9") )

# ---------
# Plot a 1 standard error envelope around the mean population abundances
# ---------
for (j in 1:nrow(means_no_zeroes))
{
  
  polygon(x= c(time, rev(time)), y= c((means_no_zeroes[j, ] - se_no_zeroes[j, ]), rev(means_no_zeroes[j, ] + se_no_zeroes[j, ])), col=adjustcolor(viridis(4)[j], alpha.f= 0.5), border = NA)
}
matplot(time, t(means_no_zeroes), add= TRUE, col=viridis(4), type="l", lwd=1.5, lty="solid")

dev.off()

