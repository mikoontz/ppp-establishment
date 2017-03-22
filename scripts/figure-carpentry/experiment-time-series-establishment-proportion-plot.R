### Title: establishment-through-time

#### Load libraries ####
# Clear environment if necessary
# rm(list=ls())

library(tidyr)
library(dplyr)
library(viridis)

#### Read experiment data ####
b <- read.csv("data/clean-establishment-data.csv", stringsAsFactors = FALSE) %>%
  filter(gap == FALSE) %>%
  select(ID, size, number, extant1, extant2, extant3, extant4, extant5, extant6, extant7, extant8, extant9) %>%
  gather(key = generation, value = extant, -c(1, 2, 3)) %>%
  group_by(size, number, generation) %>%
  summarize(establish_prop = mean(extant)) %>%
  spread(key = generation, value = establish_prop) %>%
  mutate(intro_regime = paste(size, number, sep = "x")) %>%
  ungroup() %>%
  arrange(number) %>%
  select(12, 3:11) %>%
  as.data.frame()
b

# Set the time points where populations were sustained by additinoal introductions to NA so they don't plot
b[2, 2] <- NA
b[3, 2:4] <- NA
b[4, 2:5] <- NA

#### Through 7 generations ####
pdf("figures/experiment-time-series-establishment-proportion-absolute-time-type-seven-generations.pdf", height = 3.14961, width = 3.14961)
par(mar = c(3.5, 3.5, 0.5, 0.5), family = "Helvetica", mgp = c(2.25, 1, 0))

matplot(x = 1:7, y = 100 * t(b[, c(-1, -9, -10)]), 
        lty = 1, 
        pch = 19, 
        lwd = 2, 
        col = viridis(4), 
        type = "b", 
        xlab = "Generation", 
        ylab = "Percent established", 
        xaxt = "n", 
        yaxt = "n", 
        xlim = c(1, 7), 
        ylim = c(75, 102),
        cex = 0.75)

legend("bottomleft", 
       legend = c("20x1", "10x2", "4x5", "5x4"), 
       bty = "n", 
       col = viridis(4), 
       pch = 15)

axis(side = 1, at = 1:7, labels = 1:7)
axis(side = 2, at = c(80, 90, 100), las = 1)

# No need for the introductionr regime labels if we can have color
# text(x = 7.25, y = 100 * b[, 8], labels = b[, 1], pos = 4)

dev.off()

#### Through all 9 generations ####
# pdf("figures/experiment-time-series-establishment-proportion-absolute-time-type.pdf", height=5, width=6)
# par(mar = c(4.7, 4.7, 1, 1), family = "Helvetica", mgp = c(3.3, 1, 0))
# 
# matplot(x = 1:9, y = 100*t(b[, -1]), 
#         lty = 1, 
#         pch = 19, 
#         lwd = 2, 
#         col = 1, 
#         type = "b", 
#         xlab = "Generation", 
#         ylab = "Percent established", 
#         xaxt = "n", 
#         yaxt = "n", 
#         xlim = c(1, 10.5), 
#         ylim = c(69, 102),
#         bty = "L",
#         cex.lab = 1.5)
# 
# axis(side = 1, at = 1:9, labels = c(1, NA, 3, NA, 5, NA, 7, NA, 9), cex.axis = 1.5)
# axis(side = 2, at = c(70, 80, 90, 100), las = 1, cex.axis = 1.5)
# 
# text(x = 10, y = 100 * b[, ncol(b)], labels = b[, 1], cex = 1.5)
# 
# dev.off()


#### Establishment when assessed at relative time points ####

b <- read.csv("data/clean-establishment-data.csv", stringsAsFactors = FALSE) %>%
  filter(gap == FALSE) %>%
  select(ID, size, number, extant_1_after, extant_2_after, extant_3_after, extant_4_after, extant_5_after) %>%
  gather(key = generation, value = extant, -c(1, 2, 3)) %>%
  group_by(size, number, generation) %>%
  summarize(establish_prop = mean(extant)) %>%
  spread(key = generation, value = establish_prop) %>%
  mutate(intro_regime = paste(size, number, sep = "x")) %>%
  ungroup() %>%
  arrange(number) %>%
  select(ncol(.), 3:(ncol(.)-1)) %>%
  as.data.frame()
b

pdf("figures/experiment-time-series-establishment-proportion-relative-time-type.pdf", height=5, width=6)
par(mar = c(5, 5.5, 1, 1), family = "Helvetica", mgp = c(4, 1, 0))

matplot(x = 1:5, y = 100*t(b[, -1]), 
        lty = 1, 
        pch = 19, 
        lwd = 2, 
        col = 1, 
        type = "b", 
        xlab = "Generations since\nfinal introduction event", 
        ylab = "Percent established", 
        xaxt = "n", 
        yaxt = "n", 
        xlim = c(1, 6), 
        bty = "L",
        cex.lab = 1.5)

axis(side = 1, at = 1:5, cex.axis = 1.5)
axis(side = 2, las = 1, cex.axis = 1.5)
text(x = 5.5, y = 100 * b[, ncol(b)], labels = b[, 1], cex = 1.5)
dev.off()
