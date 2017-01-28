# Title: temporary extinctions
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150606
# Last Updated: 20170127
#
# Test 'Eco-evolutionary consequences of multiple introductions for colonizing individuals' main experiment data for an effect of temporary extinctions on extinction probability
#



#-----------
# Load libraries
#-----------

library(lme4)
library(lsmeans)

#-----------
# Set up objects
#-----------

#### Effect on establishment probability ####
b <- read.csv("data/clean-establishment-data.csv")
b$number <- as.factor(b$number)
b$block <- as.factor(b$block)
b$gap <- as.factor(b$gap)

bb <- subset(b, subset = number %in% c(2,4,5) & gap == "FALSE")
bb$temp <- ifelse(bb$temp.extinctions > 0, yes=TRUE, no=FALSE)
bb$temp <- as.factor(bb$temp)
bb$effective_colonists <- 20 - bb$loss

tempA <- glmer(extant7 ~ number*environment + temp + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))
reducedModel <- glmer(extant7 ~ number*environment + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))

anova(tempA, reducedModel) ### Significant effect of temporary extinction versus not
summary(tempA)

### Describe the effect further: confidence interval around estimate of difference
extinct_effect <- lsmeans(tempA, pairwise ~ temp, adjust="none") 
confint(extinct_effect$contrasts)

(lsmeans(tempA, pairwise ~ temp, adjust="none", transform = "response"))

effective_colonists_model <- glmer(extant7 ~ number*environment + effective_colonists + (1 | block), data=bb, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(effective_colonists_model)

anova(effective_colonists_model, reducedModel)
# No significant effect of AMOUNT of loss.

#### Effect on size of extant populations ####

tempA_popSize <- glmer(N6plus1 ~ number*environment + temp + (1 | block), data=bb, family=poisson, control=glmerControl(optimizer="bobyqa"))
reducedModel_popSize <- glmer(N6plus1 ~ number*environment + (1 | block), data=bb, family=poisson, control=glmerControl(optimizer="bobyqa"))

anova(tempA_popSize, reducedModel_popSize)
# Significant effect of temporary extinction versus not
summary(tempA_popSize)

### Describe the effect further: confidence interval around estimate of difference
extinct_effect_popSize <- lsmeans(tempA_popSize, pairwise ~ temp, adjust="none") 
confint(extinct_effect_popSize$contrasts)

(lsmeans(tempA_popSize, pairwise ~ temp, adjust="none", transform = "response"))


effective_colonists_model_popSize <- glmer(N6plus1 ~ number*environment + effective_colonists + (1 | block), data=bb, family=poisson, control=glmerControl(optimizer="bobyqa"))
summary(effective_colonists_model_popSize)

anova(effective_colonists_model_popSize, reducedModel_popSize)
confint(effective_colonists_model_popSize)
# Significant effect of AMOUNT of loss.
