#### Title: Extinction probability analysis ####
###
### Author: Michael Koontz
### Email: mikoontz@gmail.com
###
### Date Created: 20150502
### Penultimate Updated: 20160124 (Changed frame of reference from extinction to establishment)
### Last Updated: 20161120
###
#### Purpose: Use tidy extinction data to assess the probability of extinction depending on different covariates. ####
###
### There are 6 different analysis presented here. 

### The first uses generation 9 (the end of the experiment) as the fixed point where we assess whether or not a population has established.

### The second analysis uses generation 6 as the fixed point where we assess whether or not a population has established. Generation 6 is after all introductions have completed, but before we see the rapid decline in population abundance that may indicate an additional biological process is at work.

### The third analysis uses all generations after introductions are complete (5 through 9) and includes a random effect of population ID to account for the repeated measures. This analysis was halted early on due to an inability to fit the models to the data.

### The fourth analysis uses generation 5 as the fixed point where we assess whether or not a population has established. Generation 5 is the first generation after all introductions have completed for all populations.

### The fifth analysis uses a relative timepoint to measure extant/extinct. Here, we assess establishment 5 generations after the final introduction event for each introduction regime. Since the 20x1 introduction regime took 1 generation to implement (new individuals to a patch in generation P), we look at generation 5 (F5) for establishment. Since the 4x5 introduction regime took 5 generations to implement (new individuals to a patch in generations P, F1, F2, F3, F4), we look at generation 9 (F9) for establishment. This analysis preempts the counterargument that the fewer, larger introductions had more time to go extinct, thus we see they have a lower estimated establishment probability at any given absolute time point.

### Sixth analysis assess establishment at generation 7

#### Load libraries, read data, and set it up for analysis ####
# Clear environment if necessary
# rm(list=ls())

library(lme4)
library(lsmeans)
library(multcompView)
library(tidyr)

b <- read.csv("data/clean-establishment-data.csv")
b$number <- as.factor(b$number)
b$block <- as.factor(b$block)
b$gap <- as.factor(b$gap)


##### Analysis 1: Extant/extinct assessment at generation 9 #####

#### Analysis 1: Assessing effect of the introduction gap ####

# Subset data to represent the only blocks and introduction treatments that COULD have experienced introduction gaps
gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# We use a simple random effects structure here, because there is far less data in this subsetted dataset. Use extant assessment at generation 10 (9 full censuses after initial introduction).
m1 <- glmer(extant9 ~ number*environment*gap + (1 | block), data=gap.potential, family=binomial, control=glmerControl(optimizer="bobyqa"))

# Update the model fit to remove the 3-way interaction
m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1,m2) # LRT suggests model with 3-way interaction is not more likely than model with additive effect of gap. Drop 3-way interaction with gap.

m3 <- update(m2, formula= .~. -gap:environment)

anova(m2, m3) # LRT suggests model with interaction between gap and environment is not more likely than model without, so we drop the gap:environment interaction.

m4 <- update(m3, formula= .~. -gap:number)

anova(m3, m4) # LRT suggests model with interaction between gap and propagule number is not more likely than model without, so we drop the gap:number interaction.

m5 <- update(m4, formula= .~. -gap)

anova(m4, m5) # LRT suggests that the model with an additive effect of gap is NOT more likely than the model without, so we drop the effect of gap. Note this is DIFFERENT than the analysis where extant/extinct is assessed relative to when introductions finish.

#--------------------
#### Analysis 1: Plot effect (if any) of introduction gap ####
# Quick plot to show marginal effect of introduction gap (very overlapping CIs when looking at generation 9)
#--------------------
# Recall that plogis() has the effect of being an inverse logit function and qlogis() has the effect of being a logit function.

ls <- summary(lsmeans::lsmeans(m4, pairwise ~ as.factor(gap), adjust="none")$lsmeans)

cont <- lsmeans::lsmeans(m4, pairwise ~ as.factor(gap), adjust="none")$contrasts
cont.mu <- plogis(summary(cont)$estimate)
cont.SE <- c(plogis(summary(cont)$estimate - summary(cont)$SE), plogis(summary(cont)$estimate + summary(cont)$SE))
cont.SE

mu <- plogis(ls$lsmean)
CI <- rbind(plogis(ls$asymp.LCL), plogis(ls$asymp.UCL))

# pdf("Clean Plots/introduction gap effect.pdf", height=3.75, width=3.75)
par(mar=c(3,4,3,1))
plot(x=c(1,2), y=mu, ylim=c(0.45, 1), xlim=c(0.5,2.5), las=1, pch=19, xaxt="n", xlab=NA, ylab="Probability of Establishment")
axis(side=1, at=c(1,2), labels=c("No gap", "Gap"), tick=FALSE)
arrows(x0=c(1,2), y0=CI[1,], y1=CI[2,], code=3, length=0.1, angle=90, lwd=2)
# dev.off()

mtext(side=3, text="Effect of introduction gap on establishment probability\nafter 9 filial generations")

#-----------------------
# No need to trim the data to exclude those populations that experienced an introduction gap.
#-----------------------

# Subset to remove populations with an introduction gap
b.trim <- b

#-------------------
#### Analysis 1: Determining the random effects structure ####
### Test whether there are interactions between the covariates and the random effect of block
#-------------------

# The 'keep it maximal' random effects structure; Possibly too many parameters to justify this approach.
m5a <- glmer(extant9 ~ number*environment + (number*environment | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) # Doesn't converge

m5b <- update(m5a, formula= .~ number*environment + (number + environment | block))

anova(m5a, m5b) # Doesn't appear to be a three way interaction bw number, environment, and block.

m5c <- update(m5b, formula= . ~ number*environment + (number | block))

anova(m5b, m5c) # No interaction between block and environment

m5d <- update(m5c, formula= . ~ number*environment + (environment | block))

anova(m5c, m5d) # No interaction between block and number

# There appear to be no interactions between block and the fixed covariates, so we use a simple random intercept structure for the random effect of temporal block.

#----------------
#### Analysis 1: Influence of fixed effects ####
#----------------

# Use LRT tests to guide interpretation, but all fixed effects will remain in the model in the end
m6 <- glmer(extant9 ~ number*environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa"))

m7 <- update(m6, formula= .~. - number:environment)

anova(m6, m7) # Doesn't appear to be an interaction between propagule number and environmental stability

m8 <- update(m7, formula= .~. -environment)

anova(m7, m8) # Environment doesn't seem to affect extinction probability

m9 <- update(m8, formula= .~. -number)

anova(m8, m9) # Propagule number seems to be very important for extinction probability

# The final model which includes all fixed effects and uses the trimmed dataset without populations that experienced a gap in the introduction period
final <- m6

#### Analysis 1: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters)$.group[order(cld(results)$number)]
xvals <- 1:length(posthoc$lsmean)
min_y <- min(plogis(posthoc$asymp.LCL))

plot(x = xvals, y=plogis(posthoc$lsmean), 
     ylim=c(min_y, 1.0), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las=1, 
     pch=19, 
     xaxt="n", 
     xlab="Introduction regime", 
     ylab="Establishment probability", 
     bty="L")

axis(side=1, 
     at = xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE)

arrows(x0 = xvals, 
       y0 = plogis(posthoc$asymp.LCL), 
       y1 = plogis(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2)

text(x = 1:4, y = 1.02, labels = sig_letters)

mtext(side=3, 
      text="Effect of propagule number on extinction probability\nafter 9 filial generations", 
      line = 2)
#### Analysis 1: Effect of temporary extinctions ####

bb <- subset(b.trim, subset=number %in% c(2,4,5))
bb$temp <- ifelse(bb$temp.extinctions > 0, yes=TRUE, no=FALSE)
bb$temp <- as.factor(bb$temp)

(nrow(b.trim) - sum(b.trim$extant9)) / nrow(b.trim)
# Total extinct through 9 generations: 150/917 = 0.164
(nrow(b.trim) - sum(b.trim$extant_5_after)) / nrow(b.trim)
# Total extinct 5 generations after final intro: 108/917 = 0.1178

sum(b.trim$temp.extinctions > 0) / nrow(bb)
# Temporary extinction before new introduction: 115/677 = 0.1699
sum(b.trim$temp.extinctions > 1) / nrow(bb)
# Temporarily extinct twice before new introduction: 13/677 = 0.0192
sum(b.trim$temp.extinctions > 2) / nrow(bb)
# Temporarily extinct thrice before new introduction: 1/677 = 0.00148

# Proportion extant at end of experiment
aggregate(extant_5_after ~ temp + number, FUN=mean, data=bb) 
aggregate(loss ~ number, FUN=mean, data=bb)
aggregate(temp.extinctions ~ number, FUN=sum, data=bb)

tempA <- glmer(extant9 ~ number*environment + temp + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

reducedModel <- glmer(extant9 ~ number*environment + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

anova(tempA, reducedModel) # Significant effect of temporary extinction versus not

lossA <- glmer(extant9 ~ number*environment + loss + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

anova(lossA, reducedModel) # Significant effect of AMOUNT of loss. More genetic inputs lost meant a higher extinction probability
summary(lossA) # One more potential beetle lost means a decrease in 0.11 log odds of establishment.
# plogis(0.92) # Mean probability of establishment
# plogis(0.92 - 0.11) # Approximate probability of establishment with 1 fewer individual
# plogis(0.92 - .11*5) # Approximate probability of establishment with 5 fewer individuals



##### Analysis 2: Extant/extinct assessment at generation 6 #####

#### Analysis 2: Assessing effect of the introduction gap ####
# Subset data to represent the only blocks and introduction treatments that COULD have experienced introduction gaps
gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# We use a simple random effects structure here, because there is far less data in this subsetted dataset. Use extant assessment at generation 6.
m1 <- glmer(extant6 ~ number*environment*gap + (1 | block), data=gap.potential, family=binomial, control=glmerControl(optimizer="bobyqa"))
# Model fails to converge...

# Update the model fit to remove the 3-way interaction
m2 <- update(m1, formula= .~. - number:environment:gap)
# This model is nearly unidentifiable.

anova(m1,m2) # LRT suggests model with 3-way interaction is not more likely than model with additive effect of gap. Drop 3-way interaction with gap. Also, the previous model failed to converge so really we'll "start" with the model that doesn't have the 3-way interaction

m3 <- update(m2, formula= .~. -gap:environment)

anova(m2, m3) # LRT suggests model with interaction between gap and environment is not more likely than model without, so we drop the gap:environment interaction. Also, m2 and m3 didn't fit without a Warning about degenerate Hessians.

m4 <- update(m3, formula= .~. -gap:number) # This is the first model that fits without throwing a Warning message. Let's use this as the fullest model.

anova(m3, m4) # LRT suggests model with interaction between gap and propagule number IS more likely than the model without, so we have to deal with this.

#--------------------
#### Analysis 2: Plot effect (if any) of introduction gap ####
# Quick plot to show marginal effect of introduction gap (very overlapping CIs when looking at generation 9)
#--------------------
# Recall that plogis() has the effect of being an inverse logit function and qlogis() has the effect of being a logit function.

ls <- summary(lsmeans::lsmeans(m4, pairwise ~ as.factor(gap) + number, adjust="none")$lsmeans)

cont <- lsmeans::lsmeans(m4, pairwise ~ as.factor(gap) + number, adjust="none")$contrasts
cont.mu <- plogis(summary(cont)$estimate)
cont.SE <- c(plogis(summary(cont)$estimate - summary(cont)$SE), plogis(summary(cont)$estimate + summary(cont)$SE))
cont.SE

mu <- plogis(ls$lsmean)
CI <- rbind(plogis(ls$asymp.LCL), plogis(ls$asymp.UCL))

# pdf("Clean Plots/introduction gap effect.pdf", height=3.75, width=3.75)
par(mar=c(5,4,2,1))
plot(x=1:4, y=mu, ylim=c(0.45, 1), xlim=c(0.5,4.5), las=1, pch=19, xaxt="n", xlab=NA, ylab="Probability of Establishment")
axis(side=1, at=1:4, labels=rep(c("No gap", "Gap"), times = 2), tick=FALSE)
mtext(side = 1, at = c(1.5, 3.5), text = c("10x2", "5x4"), line = 3)

arrows(x0=1:4, y0=CI[1,], y1=CI[2,], code=3, length=0.1, angle=90, lwd=2)
# dev.off()

mtext(side=3, text="Effect of introduction gap on establishment probability\nafter 6 generations")

#-----------------------
# We need to trim the data to exclude those populations that experienced an introduction gap.
#-----------------------

# Subset to remove populations with an introduction gap
b.trim <- b[b$gap == "FALSE", ]
# dim(b.trim) # 842 populations didn't experience an introduction gap.

#-------------------
#### Analysis 2: Determining the random effects structure ####
# Test whether there are interactions between the covariates and the random effect of block
#-------------------

# The 'keep it maximal' random effects structure; Possibly too many parameters to justify this approach.
m5a <- glmer(extant6 ~ number*environment + (number*environment | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa"))

m5b <- update(m5a, formula= .~ number*environment + (number + environment | block)) # This is the most complicated random effects structure we can fit.

anova(m5a, m5b) # Doesn't appear to be a three way interaction bw number, environment, and block based on LRT, but we had trouble fitting the m5a model.

m5c <- update(m5b, formula= . ~ number*environment + (number | block))

anova(m5b, m5c) # No interaction between block and environment

m5d <- update(m5c, formula= . ~ number*environment + (environment | block))

anova(m5c, m5d) # No interaction between block and number

# There appear to be no interactions between block and the fixed covariates, so we use a simple random intercept structure for the random effect of temporal block.

#----------------
#### Analysis 2: Influence of fixed effects ####

# Use LRT tests to guide interpretation, but all fixed effects will remain in the model in the end
m6 <- glmer(extant6 ~ number + environment + number:environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa"))

m7 <- update(m6, formula= .~. - number:environment)

anova(m6, m7) # Doesn't appear to be an interaction between propagule number and environmental stability

m8 <- update(m7, formula= .~. -environment)

anova(m7, m8) # Environment doesn't seem to affect extinction probability

m9 <- update(m8, formula= .~. -number)

anova(m8, m9) # Propagule number seems to be very important for extinction probability

# The final model which includes all fixed effects and uses the trimmed dataset without populations that experienced a gap in the introduction period
final <- m6
#### Analysis 2: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters)$.group[order(cld(results)$number)]
xvals <- 1:length(posthoc$lsmean)

plot(x = xvals, y=plogis(posthoc$lsmean), 
     ylim=c(0.7, 1.0), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las=1, 
     pch=19, 
     xaxt="n", 
     xlab="Introduction regime", 
     ylab="Establishment probability", 
     col=cols, 
     bty="L")

# cols <- c("dodgerblue", "brown1", "gold", "green")
cols <- "black"

axis(side=1, 
     at = xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE)

arrows(x0 = xvals, 
       y0 = plogis(posthoc$asymp.LCL), 
       y1 = plogis(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2,
       col = cols)

text(x = 1:4, y = 1.02, labels = sig_letters)

mtext(side=3, 
      text="Effect of propagule number on extinction probability\nafter 6 generations", 
      line = 2)

##### Analysis 3: Extant/extinct assessment at each generation, include generation as a fixed effect, population ID as random effect #####

#### Analysis 3: Assessing effect of the introduction gap ####
bb <- subset(b, select = c("ID", "block", "number", "environment", "gap", paste0("extant", 1:9)))

# Gather the extant/extinct assessments currently in short form across columns into a 2-column pair of keys and values.
long_b <- gather(bb, key = generation, value = extant, -c(ID, block, number, environment, gap))
# Convert the gathered "generation" column to a numeric value representing which filial generation the row refers to
long_b$generation <- as.numeric(substr(long_b$generation, start = 7, stop = nchar(long_b$generation)))

# We are only interested in generations 5 through 9, which represent time points with all introductions completed
long_b <- subset(long_b, generation >= 5)

gap.potential <- subset(long_b, subset=(block!=3)&(number%in%c(2,4)))

str(gap.potential)
fm1 <- glmer(extant ~ number * environment * generation + gap + (1 | ID) + (1 | block), data = gap.potential, family = "binomial", control=glmerControl(optimizer="bobyqa"))

fm2 <- update(fm1, . ~ . - gap)
anova(fm1, fm2) # Looks like model with gap is more likely than the model without, so we'll drop those populations that experienced an introduction gap.

b_trim <- long_b[long_b$gap == "FALSE", ]
#### Analysis 3: Determining the random effects structure ####
# Even the simplest possible model that we'd be interested in won't fit using this framework. Let's abandon it in favor of fitting models for each generation that we care about.
fm1 <- glmer(extant ~ number * environment + generation + (1 | ID) + (1 | block), data = b_trim, family = "binomial", control=glmerControl(optimizer="bobyqa"))

# What about ensuring ID is a factor?
b_trim$ID <- as.factor(b_trim$ID)
fm2 <- glmer(extant ~ number * environment + generation + (1 | ID) + (1 | block), data = b_trim, family = "binomial", control=glmerControl(optimizer="bobyqa"))
#No convergence.

# What about making generation a factor?
fm3 <- glmer(extant ~ number * environment + as.factor(generation) + (1 | ID) + (1 | block), data = b_trim, family = "binomial", control=glmerControl(optimizer="bobyqa"))
# Still no convergence.
# Additional check if no convergence should be worrisome
relgrad <- with(fm3@optinfo$derivs, solve(Hessian, gradient))
max(abs(relgrad))


# What kind of variation is in the data?
aggregate(extant ~ number + environment + generation, FUN = mean, data = b_trim)
# Some, but not a lot when broken down by generation like this. Perhaps that's why the model won't converge?

##### Analysis 4: Extant/extinct assessment at generation 5 #####
## This analysis may not prove useful, as no populations in the 4x5 introduction regime went extinct by this time.
#### Analysis 4: Assessing effect of the introduction gap ####
# Subset data to represent the only blocks and introduction treatments that COULD have experienced introduction gaps
gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# We use a simple random effects structure here, because there is far less data in this subsetted dataset. Use extant assessment at generation 10 (9 full censuses after initial introduction).
m1 <- glmer(extant5 ~ number*environment*gap + (1 | block), data=gap.potential, family=binomial, control=glmerControl(optimizer="bobyqa"))

# Update the model fit to remove the 3-way interaction
m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1,m2) # LRT suggests model with 3-way interaction is not more likely than model with additive effect of gap. Also, model with 3-way interaction failed to converge. Drop 3-way interaction with gap.

m3 <- update(m2, formula= .~. -gap:environment)

anova(m2, m3) # LRT suggests model with interaction between gap and environment is not more likely than model without, so we drop the gap:environment interaction.

m4 <- update(m3, formula= .~. -gap:number)

anova(m3, m4) # LRT suggests model with interaction between gap and propagule number is not more likely than model without, so we drop the gap:number interaction.

m5 <- update(m4, formula= .~. -gap)

anova(m4, m5) # LRT suggests that the model with an additive effect of gap is NOT more likely than the model without, so we drop the effect of gap.
#### Analysis 4: Plot effect (if any) of introduction gap ####
# Quick plot to show marginal effect of introduction gap (very overlapping CIs when looking at generation 9)
#--------------------
# Recall that plogis() has the effect of being an inverse logit function and qlogis() has the effect of being a logit function.

ls <- summary(lsmeans::lsmeans(m4, pairwise ~ as.factor(gap), adjust="none")$lsmeans)

cont <- lsmeans::lsmeans(m4, pairwise ~ as.factor(gap), adjust="none")$contrasts
cont.mu <- plogis(summary(cont)$estimate)
cont.SE <- c(plogis(summary(cont)$estimate - summary(cont)$SE), plogis(summary(cont)$estimate + summary(cont)$SE))
cont.SE

mu <- plogis(ls$lsmean)
CI <- rbind(plogis(ls$asymp.LCL), plogis(ls$asymp.UCL))

# pdf("Clean Plots/introduction gap effect.pdf", height=3.75, width=3.75)
par(mar=c(5,4,2,1))
plot(x=1:2, y=mu, ylim=c(0.45, 1), xlim=c(0.5,2.5), las=1, pch=19, xaxt="n", xlab=NA, ylab="Probability of Establishment")
axis(side=1, at=1:4, labels=rep(c("No gap", "Gap"), times = 2), tick=FALSE)

arrows(x0=1:4, y0=CI[1,], y1=CI[2,], code=3, length=0.1, angle=90, lwd=2)
# dev.off()

mtext(side=3, text="Effect of introduction gap on establishment probability\nafter 5 generations")

#-----------------------
# No need to trim the data to exclude those populations that experienced an introduction gap.
#-----------------------

# Subset to remove populations with an introduction gap
b.trim <- b
#### Analysis 4: Determining the random effects structure ####
# Test whether there are interactions between the covariates and the random effect of block
#-------------------

# The 'keep it maximal' random effects structure; Possibly too many parameters to justify this approach.
m5a <- glmer(extant5 ~ number*environment + (number*environment | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) # This model doesn't converge

m5b <- update(m5a, formula= .~ number*environment + (number + environment | block)) # This model doesn't converge

anova(m5a, m5b) # Doesn't appear to be a three way interaction bw number, environment, and block based on LRT, but we had trouble fitting the m5a model.

m5c <- update(m5b, formula= . ~ number*environment + (number | block)) # This model doesn't converge

anova(m5b, m5c) # No interaction between block and environment

m5d <- update(m5c, formula= . ~ number*environment + (environment | block)) # This is the most complicated random effects structure whose model converges

anova(m5c, m5d) # No interaction between block and number
# There appear to be no interactions between block and the fixed covariates, so we use a simple random intercept structure for the random effect of temporal block.
# m6 <- update(m5d, formula = . ~ number * environment + (1 | block))
#### Analysis 4: Influence of fixed effects ####
# Use LRT tests to guide interpretation, but all fixed effects will remain in the model in the end
m6 <- glmer(extant5 ~ number + environment + number:environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) # Model fails to converge. Problem is that the 4x5 introduction regime has NO extinctions yet at generation 5

m7 <- update(m6, formula= .~. - number:environment) # Model fails to converge

anova(m6, m7) # Doesn't appear to be an interaction between propagule number and environmental stability

m8 <- update(m7, formula= .~. -environment) # This model does converge

anova(m7, m8) # Environment doesn't seem to affect extinction probability

m9 <- update(m8, formula= .~. -number) # This model also converges

anova(m8, m9) # Propagule number seems to be very important for extinction probability

# The final model which includes all fixed effects and uses the trimmed dataset without populations that experienced a gap in the introduction period. Note though that this model didn't converge.

final <- m6
#### Analysis 4: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters)$.group[order(cld(results)$number)]
xvals <- 1:length(posthoc$lsmean)
min_y <- min(plogis(posthoc$asymp.LCL))

plot(x = xvals, y=plogis(posthoc$lsmean), 
     ylim=c(min_y, 1.0), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las=1, 
     pch=19, 
     xaxt="n", 
     xlab="Introduction regime", 
     ylab="Establishment probability", 
     col=cols, 
     bty="L")

# cols <- c("dodgerblue", "brown1", "gold", "green")
cols <- "black"

axis(side=1, 
     at = xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE)

arrows(x0 = xvals, 
       y0 = plogis(posthoc$asymp.LCL), 
       y1 = plogis(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2,
       col = cols)

text(x = 1:4, y = 1.02, labels = sig_letters)

mtext(side=3, 
      text="Effect of propagule number on extinction probability\nafter 5 filial generations", 
      line = 2)

##### Analysis 5: Extant/extinct assessment 5 generations after final introduction event #####

#### Analysis 5: Assessing effect of the introduction gap ####

# Subset data to represent the only blocks and introduction treatments that COULD have experienced introduction gaps
gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# We use a simple random effects structure here, because there is far less data in this subsetted dataset. Use extant assessment at generation 10 (9 full censuses after initial introduction).

m1 <- glmer(extant_5_after ~ number*environment*gap + (1 | block), data=gap.potential, family=binomial, control=glmerControl(optimizer="bobyqa"))

# Update the model fit to remove the 3-way interaction
m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1,m2) # LRT suggests model with 3-way interaction is not more likely than model with additive effect of gap. Drop 3-way interaction with gap.

m3 <- update(m2, formula= .~. -gap:environment)

anova(m2, m3) # LRT suggests model with interaction between gap and environment is not more likely than model without, so we drop the gap:environment interaction.

m4 <- update(m3, formula= .~. -gap:number)

anova(m3, m4) # LRT suggests model with interaction between gap and propagule number is not more likely than model without, so we drop the gap:number interaction.

m5 <- update(m4, formula= .~. -gap)

anova(m4, m5) # LRT suggests that the model with an additive effect of gap IS more likely than model without, so we deal with this by removing populations that experienced an introduction gap.

#--------------------
#### Analysis 5: Plot effect (if any) of introduction gap ####
# Quick plot to show marginal effect of introduction gap (very overlapping CIs when looking at generation 9)
#--------------------
# Recall that plogis() has the effect of being an inverse logit function and qlogis() has the effect of being a logit function.

ls <- summary(lsmeans::lsmeans(m4, pairwise ~ as.factor(gap), adjust="none")$lsmeans)

cont <- lsmeans::lsmeans(m4, pairwise ~ as.factor(gap), adjust="none")$contrasts
cont.mu <- plogis(summary(cont)$estimate)
cont.SE <- c(plogis(summary(cont)$estimate - summary(cont)$SE), plogis(summary(cont)$estimate + summary(cont)$SE))
cont.SE

mu <- plogis(ls$lsmean)
CI <- rbind(plogis(ls$asymp.LCL), plogis(ls$asymp.UCL))

# pdf("Clean Plots/introduction gap effect.pdf", height=3.75, width=3.75)
par(mar=c(3,4,3,1))
plot(x=c(1,2), y=mu, ylim=c(0.45, 1), xlim=c(0.5,2.5), las=1, pch=19, xaxt="n", xlab=NA, ylab="Probability of Establishment")
axis(side=1, at=c(1,2), labels=c("No gap", "Gap"), tick=FALSE)
arrows(x0=c(1,2), y0=CI[1,], y1=CI[2,], code=3, length=0.1, angle=90, lwd=2)
# dev.off()

mtext(side=3, text="Effect of introduction gap on establishment probability\n5 generations after the final introductione event")

### Introduction gap significantly decreases the probability of establishment.

#-----------------------
# We need to trim the data to exclude those populations that experienced an introduction gap.
#-----------------------

# Subset to remove populations with an introduction gap
b.trim <- b[b$gap == "FALSE", ]

#-------------------
#### Analysis 5: Determining the random effects structure ####
### Test whether there are interactions between the covariates and the random effect of block
#-------------------

# The 'keep it maximal' random effects structure; Possibly too many parameters to justify this approach.
m5a <- glmer(extant_5_after ~ number*environment + (number*environment | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) # as we might expect, the model doesn't fit.

m5b <- update(m5a, formula= .~ number*environment + (number + environment | block)) # This one does converge and has the most complex random effects structure we could use.

anova(m5a, m5b) # Doesn't appear to be a three way interaction bw number, environment, and block.

m5c <- update(m5b, formula= . ~ number*environment + (number | block))

anova(m5b, m5c) # No interaction between block and environment

m5d <- update(m5c, formula= . ~ number*environment + (environment | block))

anova(m5c, m5d) # No interaction between block and number

# There appear to be no interactions between block and the fixed covariates, so we use a simple random intercept structure for the random effect of temporal block.
#----------------
#### Analysis 5: Influence of fixed effects ####
#----------------

# Use LRT tests to guide interpretation, but all fixed effects will remain in the model in the end
m6 <- glmer(extant_5_after ~ number*environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa"))

m7 <- update(m6, formula= .~. - number:environment)

anova(m6, m7) # Doesn't appear to be an interaction between propagule number and environmental stability

m8 <- update(m7, formula= .~. -environment)

anova(m7, m8) # Environment doesn't seem to affect extinction probability

m9 <- update(m8, formula= .~. -number)

anova(m8, m9) # Propagule number seems to be very important for extinction probability

# The final model which includes all fixed effects and uses the trimmed dataset without populations that experienced a gap in the introduction period
final <- m6
#### Analysis 5: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters)$.group[order(cld(results)$number)]
xvals <- 1:length(posthoc$lsmean)
min_y <- min(plogis(posthoc$asymp.LCL))

plot(x = xvals, y=plogis(posthoc$lsmean), 
     ylim=c(min_y, 1.0), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las=1, 
     pch=19, 
     xaxt="n", 
     xlab="Introduction regime", 
     ylab="Establishment probability",
     bty="L")

axis(side=1, 
     at = xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE)

arrows(x0 = xvals, 
       y0 = plogis(posthoc$asymp.LCL), 
       y1 = plogis(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2)

text(x = 1:4, y = 1.0, labels = sig_letters)

mtext(side=3, 
      text="Effect of propagule number on extinction probability\nafter 9 filial generations", 
      line = 2)
#### Analysis 5: Effect of temporary extinctions ####

bb <- subset(b.trim, subset=number %in% c(2,4,5))
bb$temp <- ifelse(bb$temp.extinctions > 0, yes=TRUE, no=FALSE)
bb$temp <- as.factor(bb$temp)

# Note for all summaries of 5 generations after last intro event, the gap populations have been removed from the dataset
(nrow(b.trim) - sum(b.trim$extant9) ) / nrow(b.trim)
# Total extinct through 9 generations: 130/842 = 0.154 
(nrow(b.trim) - sum(b.trim$extant_5_after)) / nrow(b.trim)
# Total extinct 5 generations after final intro: 90/842 = 0.107

sum(b.trim$temp.extinctions > 0) / nrow(bb)
# Temporary extinction before new introduction: 104/602 = 0.1699
sum(b.trim$temp.extinctions > 1) / nrow(bb)
# Temporarily extinct twice before new introduction: 12/602 = 0.0199
sum(b.trim$temp.extinctions > 2) / nrow(bb)
# Temporarily extinct thrice before new introduction: 1/602 = 0.00148

dim(bb)

# Proportion extant at end of experiment
aggregate(extant_5_after ~ temp + number, FUN=mean, data=bb) 
aggregate(loss ~ number, FUN=mean, data=bb)
aggregate(temp.extinctions ~ number, FUN=sum, data=bb)

tempA <- glmer(extant_5_after ~ number*environment + temp + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

reducedModel <- glmer(extant_5_after ~ number*environment + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

anova(tempA, reducedModel)
# No significant effect of temporary extinction versus not

lossA <- glmer(extant_5_after ~ number*environment + loss + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

anova(lossA, reducedModel)
# No significant effect of AMOUNT of loss. Check this against the amount of loss at different time points (up to generation 6, for instance)


##### Analysis 6: Extant/extinct assessment at generation 7 #####
# Go directly to excluding populations that experienced an introduction gap
b.trim <- b[b$gap == "FALSE", ]
# dim(b.trim) # 842 populations didn't experience an introduction gap.

m1 <- glmer(extant7 ~ number + environment + number:environment + (1 | block), data=b.trim, family=binomial, control=glmerControl(optimizer="bobyqa")) 

m2 <- update(m1, formula= .~. - number:environment)
anova(m1, m2) # Doesn't appear to be an interaction between propagule number and environmental stability

m3 <- update(m2, formula= .~. -environment)
anova(m2, m3) # Environment doesn't seem to affect extinction probability

m4 <- update(m3, formula= .~. -number)
anova(m3, m4) # Propagule number seems to be very important for extinction probability

# The final model which includes all fixed effects and uses the trimmed dataset without populations that experienced a gap in the introduction period. Note though that this model didn't converge.

final <- m1
#### Analysis 6: Interpretation and contrasts ####

results <- lsmeans::lsmeans(final, pairwise ~ number, adjust="none")
results
posthoc <- summary(results$lsmeans)

sig_letters <- lsmeans::cld(results, Letters = letters, sort = FALSE, adjust = "none")$.group
xvals <- 1:length(posthoc$lsmean)
min_y <- min(plogis(posthoc$asymp.LCL))

plot(x = xvals, y=plogis(posthoc$lsmean), 
     ylim=c(min_y, 1.02), 
     xlim = range(xvals) + c(-0.5, 0.5), 
     las=1, 
     pch=19, 
     xaxt="n", 
     xlab="Introduction regime", 
     ylab="Establishment probability", 
     bty="L")

axis(side=1, 
     at = xvals, 
     labels = c("20x1","10x2","5x4","4x5"), 
     tick = FALSE)

arrows(x0 = xvals, 
       y0 = plogis(posthoc$asymp.LCL), 
       y1 = plogis(posthoc$asymp.UCL), 
       code = 3, 
       length = 0.1, 
       angle = 90, 
       lwd = 2)

text(x = 1:4, y = 1.02, labels = sig_letters)

mtext(side=3, 
      text="Effect of propagule number on extinction probability\nafter 7 filial generations", 
      line = 2)

(nrow(b.trim) - sum(b.trim$extant7)) / nrow(b.trim)
# Total extinct through 7 generations: 101/842 = 12.0 %

#### Analysis 6: Effect of temporary extinctions ####

bb <- subset(b.trim, subset=number %in% c(2,4,5))
bb$temp <- ifelse(bb$temp.extinctions > 0, yes=TRUE, no=FALSE)
bb$temp <- as.factor(bb$temp)
bb$effective_colonists <- 20 - bb$loss

tempA <- glmer(extant7 ~ number*environment + temp + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

reducedModel <- glmer(extant7 ~ number*environment + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

anova(tempA, reducedModel)
# Significant effect of temporary extinction versus not
summary(tempA)

effective_colonists_model <- glmer(extant7 ~ number*environment + effective_colonists + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(effective_colonists_model)

anova(effective_colonists_model, reducedModel)
# No significant effect of AMOUNT of loss.