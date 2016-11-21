# Title: temporary extinctions
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150606
# Last Updated: 20160626
#
# Test 'Eco-evolutionary consequences of multiple introductions for colonizing individuals' main experiment data for an effect of temporary extinctions on extinction probability
#

source("scripts/generate-tidy-data.R")

beetles <- read.csv("data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("data/attributes.csv")  

#-----------
# Load libraries
#-----------

library(lme4)

#-----------
# Set up objects
#-----------

# We need the number of individuals in each generation, the number of migrants, and the number of individuals in the next generation (to get the final census value)

tidy.b <- tidy.beetles(beetles=beetles)

b <- read.csv("data/clean-establishment-data.csv")

Nt <- tidy.b$Nt
Ntp1 <- tidy.b$Ntp1
colnames(Ntp1) <- colnames(Nt)
migrants <- tidy.b$migrants
b$temp.extinct <- as.factor(b$temp.extinctions > 0)

b <- merge(b, Ntp1, by="ID")

head(b)

# Convert to factors for analysis
b$number <- as.factor(b$number)
b$block <- as.factor(b$block)
b$gap <- as.factor(b$gap)

gap.potential <- subset(b, subset=(block!=3)&(number%in%c(2,4)))

# Temporary extinction or not? (Boolean covariate)
m1 <- glmer(extant_5_after ~ temp.extinct + number * environment * gap + (1 | block), data = gap.potential, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1)

m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1, m2) # Model with 3-way interaction with gap not more likely, so drop it.

m3 <- update(m2, formula= .~. - number:gap)

anova(m2, m3) # Model with number:gap not more likely than model without, so drop it.

m4 <- update(m3, formula= .~. - environment:gap)

anova(m3, m4) # Model with environment:gap not more likely than model without, so drop the term

m5 <- update(m4, formula= .~. - gap)

anova(m4, m5) # Model with additive effect of gap is more likely than model without, so leave it in.

summary(m4)

gap_model <- glmer(extant_5_after ~ temp.extinct + number * environment + gap + (1 | block), data = b, subset = (number !=1 ), family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(gap_model)

# Gap not included in the model, but data subsetted such that no populations that experienced a gap are part of the analysis
no_gap_model <- glmer(extant_5_after ~ temp.extinct + number * environment + (1 | block), data = b, subset = (gap == "FALSE" & number != 1), family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(no_gap_model)

#-----------
#-----------

# Amount of lost input (rather than just lost input vs no lost input) as a covariate. That is, a 10x2 that went extinct after the first input "lost" 10 beetles. A 4x5 introduction that went extinct in the 3rd generation, but still had one more introduction coming "lost" 15 beetles.
# # Subset to not include the 20x1 introduction scenario, since it wasn't possible for them to have a temporary extinction.
m1 <- glmer(extant_5_after ~ loss + number * environment * gap + (1 | block), data = gap.potential, family = binomial, control = glmerControl(optimizer = "bobyqa"))

m2 <- update(m1, formula= .~. - number:environment:gap)

anova(m1, m2) # Model with 3-way interaction with gap not more likely, so drop it.

m3 <- update(m2, formula= .~. - number:gap)

anova(m2, m3) # Model with number:gap not more likely than model without, so drop it.

m4 <- update(m3, formula= .~. - environment:gap)

anova(m3, m4) # Model with environment:gap not more likely than model without, so drop the term

m5 <- update(m4, formula= .~. - gap)

anova(m4, m5) # Model with additive effect of gap is more likely than model without, so leave it in.

loss_gap_model <- glmer(extant_5_after ~ loss + number * environment + gap + (1 | block), data = b, subset = (number != 1), family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(loss_gap_model)

loss_no_gap_model <- glmer(extant_5_after ~ loss * number * environment + (1 | block), data = b, subset = (gap == "FALSE" & number != 1), family = binomial, control = glmerControl(optimizer="bobyqa"))
summary(loss_no_gap_model)


# # Gap seems important so I'm going to again remove those populations that had an introduction gap AND don't include the 20x1 introduction scenario, since it wasn't possible for them to have a temporary extinction.
# m2 <- glmer(extant5after ~ temp.extinct*gap + environment*number + (1 | block), subset=(number!=1), data=b, family=binomial)
# summary(m2)
# # Neither temporary extinctions nor its interactions seem all that important
# 
# 
# m3 <- glmer(extant5after ~ loss + gap + number*environment + (1 | block), data=b, family=binomial, subset=(number!=1), control=glmerControl(optimizer="bobyqa"))
# summary(m3)
# 
# m3 <- glmer(extant_5_after ~ loss + gap + number*environment + (1 | block), data=b, family=binomial, subset=(number!=1), control=glmerControl(optimizer="bobyqa"))
# summary(m3)

bb <- subset(b, subset= (number!=1))

aggregate(extant_5_after ~ temp.extinct + number, FUN=mean, data=b)
aggregate(loss ~ number, FUN=mean, data=bb)
aggregate(temp.extinctions ~ number, FUN=sum, data=bb)

# Conclusion: Greater amount of "lost" inputs yields decreased probability of population establishment