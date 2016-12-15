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

tiff("figures/experiment-time-series-establishment-proportion-absolute-time-type.tif", units= "in", res= 300, height=5, width=6)
par(mar=c(4,4,1,1), family = "Helvetica")
matplot(x = 1:9, y = 100*t(b[, -1]), 
        lty = 1, 
        pch = 19, 
        lwd = 2, 
        col = 1, 
        type = "b", 
        xlab = "Generation", 
        ylab = "Percent established", 
        xaxt = "n", 
        yaxt = "n", 
        xlim = c(1, 10), 
        bty = "L")

axis(side = 1, at = 1:9, labels = c("F1", "", "F3", "", "F5", "", "F7", "", "F9"))
axis(side = 2, las = 1)
text(x = 10, y = 100 * b[, ncol(b)], labels = b[, 1])
dev.off()


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

tiff("figures/experiment-time-series-establishment-proportion-relative-time-type.tif", units= "in", res= 300, height=5, width=6)
par(mar=c(4,4,1,1), family = "Helvetica")
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
        bty = "L")

axis(side = 1, at = 1:5)
axis(side = 2, las = 1)
text(x = 6, y = 100 * b[, ncol(b)], labels = b[, 1])
dev.off()
