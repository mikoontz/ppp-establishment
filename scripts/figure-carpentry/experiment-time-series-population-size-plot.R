### Experiment time series plot
###
### Plots population abundance of extant populations and 1 standard error margin for each generation in the parsing propagule pressure experiment

# Clear environment if necessary
# rm(list=ls())

#### Setup plotting variables ####
library(viridis)
source("scripts/data-carpentry/generate-tidy-data.R")

beetles <- read.csv("data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("data/attributes.csv")  

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

ts_getSummary <- function(data, INDEX = data$number) {
  
  means_no_zeroes <- apply(data[, 2:20], 
                           MARGIN= 2, 
                           FUN= function(x) tapply(x, 
                                                   INDEX= INDEX, 
                                                   FUN= function(y) mean(y[y>0],  na.rm= TRUE)))
  
  counts_no_zeroes <- apply(data[, 2:20], 
                            MARGIN= 2, 
                            FUN= function(x) tapply(x, 
                                                    INDEX= INDEX, 
                                                    FUN= function(y) length(which(!is.na(y) & y > 0))))
  
  sigma_no_zeroes <- apply(data[, 2:20], 
                           MARGIN= 2, 
                           FUN= function(x) tapply(x, 
                                                   INDEX= INDEX, 
                                                   FUN= function(y) sd(y, na.rm= TRUE)))
  
  ##### Calculate the standard error and time steps ####
  se_no_zeroes <- sigma_no_zeroes / sqrt(counts_no_zeroes - 1)
  
  return(list(means_no_zeroes = means_no_zeroes, se_no_zeroes = se_no_zeroes))
}

# Passing both the number and environment columns as a list to the INDEX argument will get the summary statistics for each introduction regime/environment treatment combination
ts_summary <- ts_getSummary(data, INDEX = list(data$number, data$env))
se_no_zeroes <- ts_summary$se_no_zeroes
means_no_zeroes <- ts_summary$means_no_zeroes

time <- c(0, rep(1:9, each=2))

#### Make a 2-panel plot of population size through generation 7 where each panel has a different environment treatment ####

# pdf("figures/experiment-time-series-population-size-seven-generations-environment-facet.pdf", height = 3.14961, width = 3.14961 * 2)
# 
# par(mfrow = c(1, 2), mar = c(3.5, 1.5, 2, 0.5), family = "Helvetica", mgp = c(3.3, 1, 0), oma = c(0, 1.5, 0, 0))
# 
# 
# # Plot populations from the stable environment first
# matplot(x = time[1:15], y = t(data[data$env == "stable", 2:16]), 
#         type = "l", 
#         lty = "solid", 
#         col = adjustcolor("grey", alpha.f= 0.2), 
#         las = 1, 
#         xlab = NA, 
#         ylab = NA, 
#         xaxt = "n", 
#         yaxt = "n",
#         ylim = c(0, 90),
#         cex.axis = 1)
# 
# mtext(side = 3, text = "Stable environment", cex = 1, line = 0.5)
# 
# legend("topleft", 
#        legend = c("20x1", "10x2", "4x5", "5x4"), 
#        bty = "n", 
#        col = viridis(4), 
#        pch = 15)
# 
# axis(side = 1, 
#      at = 0:7,
#      cex.axis = 1)
# 
# axis(side = 2,
#      las = 1,
#      at = seq(0, 90, by = 30),
#      hadj = 0.6)
# 
# 
# # Plot a 1 standard error envelope around the mean population abundances
# 
# for (j in 5:8)
# {
#   
#   polygon(x= c(time[1:15], time[15:1]), y= c((means_no_zeroes[j, 1:15] - se_no_zeroes[j, 1:15]), means_no_zeroes[j, 15:1] + se_no_zeroes[j, 15:1]), col=adjustcolor(viridis(4)[j - 4], alpha.f= 0.5), border = NA)
# }
# 
# matplot(time[1:15], t(means_no_zeroes[5:8, 1:15]), add= TRUE, col=viridis(4), type="l", lwd=1.5, lty="solid")
# 
# # Second panel is all of the populations in fluctuating environments
# 
# matplot(x = time[1:15], y = t(data[data$env == "fluctuating", 2:16]), 
#         type = "l", 
#         lty = "solid", 
#         col = adjustcolor("grey", alpha.f= 0.2), 
#         las = 1, 
#         xlab = NA, 
#         ylab = NA, 
#         xaxt = "n", 
#         yaxt = "n",
#         ylim = c(0, 90),
#         cex.axis = 1)
# 
# mtext(side = 3, text = "Fluctuating environment", cex = 1, line = 0.5)
# 
# axis(side = 1, 
#      at = 0:7,
#      cex.axis = 1)
# 
# axis(side = 2,
#      las = 1,
#      at = seq(0, 90, by = 30),
#      hadj = 0.6)
# 
# mtext(side = 1, text = "Generation", outer = TRUE, line = -1.25, cex = 1)
# mtext(side = 2, text = "Population size", outer = TRUE, line = 0.5, cex = 1)
# 
# # Plot a 1 standard error envelope around the mean population abundances
# 
# for (j in 1:4)
# {
#   
#   polygon(x= c(time[1:15], time[15:1]), y= c((means_no_zeroes[j, 1:15] - se_no_zeroes[j, 1:15]), means_no_zeroes[j, 15:1] + se_no_zeroes[j, 15:1]), col=adjustcolor(viridis(4)[j], alpha.f= 0.5), border = NA)
# }
# 
# matplot(time[1:15], t(means_no_zeroes[1:4, 1:15]), add= TRUE, col=viridis(4), type="l", lwd=1.5, lty="solid")
# 
# dev.off()

#### Make a 2-panel plot of population size through generation 7 where each panel has a different environment treatment and populations that eventually went extinct are highlighted. ####

pdf("figures/experiment-time-series-population-size-seven-generations-environment-facet-with-extinctions.pdf", height = 3.14961, width = 3.14961 * 2)

par(mfrow = c(1, 2), mar = c(3.5, 1.5, 2, 0.5), family = "Helvetica", mgp = c(3.3, 1, 0), oma = c(0, 1.5, 0, 0))

plotColors <- viridis(4)

# Plot populations from the stable environment first
matplot(x = time[1:15], y = t(data[data$env == "stable", 2:16]), 
        type = "l", 
        lty = "solid", 
        col = adjustcolor("grey", alpha.f= 0.2), 
        las = 1, 
        xlab = NA, 
        ylab = NA, 
        xaxt = "n",
        yaxt = "n",
        ylim = c(0, 90),
        cex.axis = 1)

matplot(x = time[1:15], y = t(data[is.na(data$N8) & data$environment == "stable", 2:16]),
        type = "l",
        col = adjustcolor(magma(5)[4], alpha.f= 0.2),
        lty = "solid",
        add = TRUE)

mtext(side = 3, text = "Stable environment", cex = 1, line = 0.5)

# legend("topleft", 
#        legend = c("20x1", "10x2", "4x5", "5x4", "extinct"), 
#        bty = "n", 
#        col = c(plotColors, magma(5)[4]), 
#        pch = 15)

legend("topleft", 
       legend = c("20x1", "10x2", "4x5", "5x4"), 
       bty = "n", 
       col = plotColors, 
       pch = 15)

legend("topright", 
       legend = "extinct by gen. 7", 
       bty = "n", 
       col = magma(5)[4], 
       pch = 15)

axis(side = 1, 
     at = 0:7,
     cex.axis = 1)

axis(side = 2,
     las = 1,
     at = seq(0, 90, by = 30),
     hadj = 0.6)


# Plot a 1 standard error envelope around the mean population abundances

for (j in 5:8)
{
  
  polygon(x= c(time[1:15], time[15:1]), y= c((means_no_zeroes[j, 1:15] - se_no_zeroes[j, 1:15]), means_no_zeroes[j, 15:1] + se_no_zeroes[j, 15:1]), col=adjustcolor(plotColors[j - 4], alpha.f= 0.5), border = NA)
}

matplot(time[1:15], t(means_no_zeroes[5:8, 1:15]), add= TRUE, col=plotColors[1:4], type="l", lwd=1.5, lty="solid")

# Second panel is all of the populations in fluctuating environments

par(mar = c(3.5, 1.5, 2, 1))

matplot(x = time[1:15], y = t(data[data$env == "fluctuating", 2:16]), 
        type = "l", 
        lty = "solid", 
        col = adjustcolor("grey", alpha.f= 0.2), 
        las = 1, 
        xlab = NA, 
        ylab = NA, 
        xaxt = "n",
        yaxt = "n",
        ylim = c(0, 90))

matplot(x = time[1:15], y = t(data[is.na(data$N8) & data$environment == "fluctuating", 2:16]),
        type = "l",
        col = adjustcolor(magma(5)[4], alpha.f= 0.2),
        lty = "solid",
        add = TRUE)

mtext(side = 3, text = "Fluctuating environment", cex = 1, line = 0.5)

axis(side = 1, 
     at = 0:7,
     cex.axis = 1)

axis(side = 2,
     las = 1,
     at = seq(0, 90, by = 30),
     hadj = 0.6)

mtext(side = 1, text = "Generation", outer = TRUE, line = -1.25, cex = 1)
mtext(side = 2, text = "Population size", outer = TRUE, line = 0.5, cex = 1)

# Plot a 1 standard error envelope around the mean population abundances

for (j in 1:4)
{
  
  polygon(x= c(time[1:15], time[15:1]), y= c((means_no_zeroes[j, 1:15] - se_no_zeroes[j, 1:15]), means_no_zeroes[j, 15:1] + se_no_zeroes[j, 15:1]), col=adjustcolor(plotColors[j], alpha.f= 0.5), border = NA)
}

matplot(time[1:15], t(means_no_zeroes[1:4, 1:15]), add= TRUE, col=plotColors[1:4], type="l", lwd=1.5, lty="solid")

dev.off()

#### Make time series plot through Generation 7 ####
# pdf("figures/experiment-time-series-population-size-seven-generations.pdf", height=5, width=6)
# par(mar = c(4.7, 4.7, 1, 1), family = "Helvetica", mgp = c(3.3, 1, 0))
# 
# matplot(x = time[1:15], y = t(time_series[, 2:16]), 
#         type = "l", 
#         lty = "solid", 
#         col = adjustcolor("grey", alpha.f= 0.1), 
#         las = 1, 
#         xlab = "Generation", 
#         ylab = "Population size", 
#         xaxt = "n", 
#         ylim = c(0,80),
#         cex.lab = 1.5,
#         cex.axis = 1.5)
# 
# legend("topleft", 
#        legend = c("20 individuals introduced in 1 event", "10x2", "4x5", "5x4"), 
#        bty = "n", 
#        col = viridis(4), 
#        lty = 1,
#        lwd = 3)
# 
# axis(side = 1, 
#      at = 0:7)
# 
# 
# # Plot a 1 standard error envelope around the mean population abundances
# 
# for (j in 1:nrow(means_no_zeroes))
# {
#   
#   polygon(x= c(time[1:15], time[15:1]), y= c((means_no_zeroes[j, 1:15] - se_no_zeroes[j, 1:15]), means_no_zeroes[j, 15:1] + se_no_zeroes[j, 15:1]), col=adjustcolor(viridis(4)[j], alpha.f= 0.5), border = NA)
# }
# matplot(time[1:15], t(means_no_zeroes[, 1:15]), add= TRUE, col=viridis(4), type="l", lwd=1.5, lty="solid")
# 
# dev.off()


# #### Make time series plot through Generation 9 ####
# pdf("figures/experiment-time-series-population-size.pdf", height=5, width=6)
# par(mar = c(4.7, 4.7, 1, 1), family = "Helvetica", mgp = c(3.3, 1, 0))
# 
# matplot(x = time, y = t(time_series[, 2:ncol(time_series)]), 
#         type = "l", 
#         lty = "solid", 
#         col = adjustcolor("grey", alpha.f= 0.1), 
#         las = 1, 
#         xlab = "Generation", 
#         ylab = "Population size", 
#         xaxt = "n", 
#         ylim = c(0,80),
#         cex.lab = 1.5,
#         cex.axis = 1.5)
# 
# legend("topleft", 
#        legend = c("20 individuals introduced in 1 event", "10x2", "4x5", "5x4"), 
#        bty = "n", 
#        col = viridis(4), 
#        lty = 1,
#        lwd = 3)
# 
# axis(side = 1, 
#      at = 0:9)
# 
# 
# # Plot a 1 standard error envelope around the mean population abundances
# 
# for (j in 1:nrow(means_no_zeroes))
# {
#   
#   polygon(x= c(time, rev(time)), y= c((means_no_zeroes[j, ] - se_no_zeroes[j, ]), rev(means_no_zeroes[j, ] + se_no_zeroes[j, ])), col=adjustcolor(viridis(4)[j], alpha.f= 0.5), border = NA)
# }
# matplot(time, t(means_no_zeroes), add= TRUE, col=viridis(4), type="l", lwd=1.5, lty="solid")
# 
# dev.off()

