### Title: Population abundance 5 generations after introductions end
### 
### Author: Michael Koontz
###
### Date Created: 20160124
### Last Updated: 20161121 
###
#### Purpose: Executes a generalized linear mixed model to assess whether treatments have different times to  extinction. First determines if it is appropriate to pool the two environment treatments (stable vs. fluctuating) and if we can ignore the drought. ####

# There are 3 analyses in this script that parallel the 3 analyses from the establishment-probability.R script that are most likely to end up in the final paper.

### Analysis 1 assesses population abundance at the end of the experiment, in generation F9

### Analysis 2 assesses population abundance at generation 6, a couple of generations after all introduction events have completed and prior to suspected major influences of other biological processes

### Analysis 3 assesses population abundance 5 generations after the final introduction event for each introduction regime (generation F5 for 20x1, F6 for 10x2, F8 for 5x4, and F9 for 4x5)

### Possible analysis 4 woiuld assess population abundance at generation 7, which allows 1 more generation of accumulated biological processes, but still appears to be before the major decline in population sizes which could be hinting that other biological processes begin to donimate the dynamics.

#### Load packages, set up data for analysis ####
# Clear environment if necessary
# rm(list=ls())

library(lme4)
library(lsmeans)
library(multcompView)

b <- read.csv("data/clean-establishment-data.csv")

b$number <- as.factor(b$number)
b$block <- as.factor(b$block)
b$gap <- as.factor(b$gap)


##### Analysis 1: Population abundance assessment at each generation, include generation as a fixed effect, population ID as random effect #####
#### Analysis 1: Assessing effect of the introduction gap ####
bb <- subset(b, select = c("ID", "block", "number", "environment", "gap", paste0("N", 0:8, "plus1")))

# Gather the extant/extinct assessments currently in short form across columns into a 2-column pair of keys and values.
long_b <- gather(bb, key = generation, value = abundance, -c(ID, block, number, environment, gap))
# Convert the gathered "generation" column to a numeric value representing which filial generation the row refers to
start <- regexpr(pattern = "[0-9]+", long_b$generation)
stop <- attr(x = start, which = "match.length") + start - 1
long_b$generation <- as.numeric(substr(long_b$generation, start = start, stop = stop)) + 1 # add the one to account for the "plus1"

# We are only interested in generations 5 through 9, which represent time points with all introductions completed
long_b <- subset(long_b, generation >= 5)

gap.potential <- subset(long_b, subset=(block!=3)&(number%in%c(2,4)))

fm1 <- glmer(abundance ~ number * environment * as.factor(generation) + gap + (1 | ID) + (1 | block), data = gap.potential, family = "poisson", control=glmerControl(optimizer="bobyqa"))

fm2 <- update(fm1, . ~ . - gap)
anova(fm1, fm2) # Looks like model with gap is NOT more likely than the model without, so we could proceed. I think it's probably best to still drop the populations that experienced an introduction gap, since that's what we did for the establishment probability analysis. The difference between these 2 analyses could itself be interesting...

b_trim <- long_b[long_b$gap == "FALSE", ]
#### Analysis 1: Determining the random effects structure ####
# Even the simplest possible model that we'd be interested in won't fit using this framework. Let's abandon it in favor of fitting models for each generation that we care about.
fm1 <- glmer(abundance ~ (number + environment) * as.factor(generation) + number:environment + (1 | ID) + (1 | block), data = b_trim, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(fm1)


#### Analysis 2: Population abundance at generation F9 ####

#### Analysis 2: Effect of introduction gap ####
# First test the effect of introduction gap on the time to extinction
# We proceed the exact same way as we do with the extinction probability analysis, except we use a Poisson likelihood instead of a binomial likelihood.
#-------------------

# Subset dataset to only populations that COULD have experienced a gap
gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# We use a simple random effects structure here, because there is far less data in this subsetted dataset.
m1 <- glmer(N9 ~ number*environment*gap + (1 | block), data=gap.potential, family=poisson, control=glmerControl(optimizer="bobyqa"))

# Update the model fit to remove the 3-way interaction
m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1,m2) # LRT suggests model with 3-way interaction is WAY more likely than model with additive effect of gap. 

#-------
# Remove gap populations 
#-------

bb <- subset(b, subset = (gap == FALSE))
#### Analysis 2: Influence of fixed effects ####
# We'll use the same random effects structure as we did for the establishment probability analysis
# There are not as many data points in this dataset, so we keep the random effects structure simple with just a random intercept of temporal block
#---------------
# Look at the effects of the fixed effects to guide interpretation, but leave all covariates in the model in the end.
#---------------

m6 <- glmer(N9 ~ number*environment + (1 | block), data=bb, family=poisson, control=glmerControl(optimizer="bobyqa"))

m7 <- update(m6, formula= .~. - number:environment)

anova(m6, m7) # Significant interaction of number of introductions and environmental stability

# The final model which includes all fixed effects
final <- m6
#### Analysis 2: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ environment + number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters)$.group[order(cld(results)$number, cld(results)$environment)]
xvals <- 1:length(posthoc$lsmean)
min_y <- min(exp(posthoc$asymp.LCL))

par(mar=c(4,4,3,1))

plot(x = xvals, y = exp(posthoc$lsmean), 
     ylim = c(0, 55), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las = 1, 
     pch = c(1,19), 
     xaxt = "n", 
     xlab = "Introduction regime", 
     ylab = "Population abundance", 
     bty = "L")

axis(side = 1, 
     at = xvals[c(1,3,5,7)] + 0.5, 
     labels=c("20x1", "10x2", "5x4", "4x5"), 
     tick=FALSE)

arrows(x0 = xvals, 
       y0 = exp(posthoc$asymp.LCL), 
       y1 = exp(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2, 
       col = rep(cols, each = 2))

legend("bottomright", 
       legend = c("fluctuating", "stable"), 
       pch = c(1, 19), 
       bty = "n")

text(x = xvals, y = 53, labels = sig_letters)

#### Analysis 3: Population abundance 5 generations after final intro event ####

#### Analysis 3: Effect of introduction gap ####
# First test the effect of introduction gap on the time to extinction
# We proceed the exact same way as we do with the extinction probability analysis, except we use a Poisson likelihood instead of a binomial likelihood.
#-------------------

# Subset dataset to only populations that COULD have experienced a gap
gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# We use a simple random effects structure here, because there is far less data in this subsetted dataset.
m1 <- glmer(N_5_after ~ number*environment*gap + (1 | block), data=gap.potential, family=poisson, control=glmerControl(optimizer="bobyqa"))

# Update the model fit to remove the 3-way interaction
m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1,m2) # LRT suggests model with 3-way interaction is WAY more likely than model with additive effect of gap. 

#-------
# Remove gap populations 
#-------

bb <- subset(b, subset=gap==FALSE)
#### Analysis 3: Influence of fixed effects ####
# We'll use the same random effects structure as we did for the establishment probability analysis
# There are not as many data points in this dataset, so we keep the random effects structure simple with just a random intercept of temporal block
#---------------
# Look at the effects of the fixed effects to guide interpretation, but leave all covariates in the model in the end.
#---------------

m6 <- glmer(N_5_after ~ number*environment + (1 | block), data=bb, family=poisson, control=glmerControl(optimizer="bobyqa"))

m7 <- update(m6, formula= .~. - number:environment)

anova(m6, m7) # Significant interaction of number of introductions and environmental stability

# The final model which includes all fixed effects
final <- m6
#### Analysis 3: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ environment + number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters)$.group[order(cld(results)$number, cld(results)$environment)]
xvals <- 1:length(posthoc$lsmean)
min_y <- min(exp(posthoc$asymp.LCL))

par(mar=c(4,4,3,1))

plot(x = xvals, y = exp(posthoc$lsmean), 
     ylim = c(0, 55), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las = 1, 
     pch = c(1,19), 
     xaxt = "n", 
     xlab = "Introduction regime", 
     ylab = "Population abundance", 
     bty = "L")

axis(side = 1, 
     at = xvals[c(1,3,5,7)] + 0.5, 
     labels=c("20x1", "10x2", "5x4", "4x5"), 
     tick=FALSE)

arrows(x0 = xvals, 
       y0 = exp(posthoc$asymp.LCL), 
       y1 = exp(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2, 
       col = rep(cols, each = 2))

legend("bottomright", 
       legend = c("fluctuating", "stable"), 
       pch = c(1, 19), 
       bty = "n")

text(x = xvals, y = 53, labels = sig_letters)

#### Add in simulation results for plot ####
# simResults <- read.csv("Clean-Analyses/simulations/establishment_abundance_table_5_50.csv")
# abundance <- numeric(8)
# abundance[seq(1,7,2)] <- simResults[, 5]
# abundance[seq(2,8,2)] <- simResults[, 2]
# 
# # Replot with simulation expectations.
# # pdf("Clean-Plots/modeled-population-abundance-with-simulations-dot-whisker.pdf", height=3.75, width=3.75)
# par(mar=c(3,4,2,1), xpd=NA, oma=c(2,0,0,0))
# 
# plot(x=xvals, y=exp(posthoc$lsmean), ylim=c(0, 52), xlim=c(0.5,4.5), las=1, pch=c(1,19), xaxt="n", xlab=NA, ylab="Population abundance", col=rep(cols, each=2), bty="L")
# axis(side=1, at=1:4, labels=c("20x1", "10x2", "5x4", "4x5"), tick=FALSE)
# arrows(x0=xvals, y0=exp(posthoc$asymp.LCL), y1=exp(posthoc$asymp.UCL), code=3, length=0.025, angle=90, lwd=2, col=rep(cols, each=2))
# legend("topright", legend=c("fluctuating", "stable"), pch=c(1, 19), bty="n")
# 
# segments(x0=(xvals)-0.025, x1=(xvals)+0.025, y0=abundance, lwd=4)
# text(x=xvals, y=exp(posthoc$asymp.UCL)+2, labels=sig_labs, cex=0.5)
# 
# mtext(side=1, text="Introduction regime", line=3)

# dev.off()
