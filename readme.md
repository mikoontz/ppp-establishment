Table of Contents for `ppp-establishment` Repository
================
Michael Koontz
March 23, 2017

-   [Abstract](#abstract)
-   [Keywords](#keywords)
-   [General Methods Framework](#general-methods-framework)
-   [Files in Repository](#files-in-repository)
    -   [scripts](#scripts)
        -   [data-carpentry](#data-carpentry)
            -   [generate-establishment-responses-for-analysis](#generate-establishment-responses-for-analysis)
            -   [generate-tidy-data](#generate-tidy-data)
        -   [figure-carpentry](#figure-carpentry)
            -   [experiment-time-series-population-size-plot](#experiment-time-series-population-size-plot)
            -   [experiment-time-series-establishment-proportion-plot](#experiment-time-series-establishment-proportion-plot)
            -   [establishment-and-size-experiment-and-simulation-results-plot](#establishment-and-size-experiment-and-simulation-results-plot)
        -   [simulations](#simulations)
            -   [NBBg-population-dynamics-function](#nbbg-population-dynamics-function)
            -   [NBBg-simulation-functions](#nbbg-simulation-functions)
            -   [NBBg-simulation-script](#nbbg-simulation-script)
        -   [environmental-stochasticity](#environmental-stochasticity)
        -   [establishment-probability](#establishment-probability)
        -   [measurement-error](#measurement-error)
        -   [NBBg-NIMBLE](#nbbg-nimble)
        -   [NBBg-script-validation](#nbbg-script-validation)
        -   [NBBg-script](#nbbg-script)
        -   [population-size](#population-size)
        -   [temporary-extinctions](#temporary-extinctions)
    -   [data](#data)
        -   [simulations](#simulations-1)
            -   [simulation\_stats\_tidy](#simulation_stats_tidy)
            -   [simulation\_stats\_raw](#simulation_stats_raw)
            -   [N\_20x1\_regime](#n_20x1_regime)
            -   [N\_var\_20x1\_regime](#n_var_20x1_regime)
            -   [N\_10x2\_regime](#n_10x2_regime)
            -   [N\_var\_10x2\_regime](#n_var_10x2_regime)
            -   [N\_5x4\_regime](#n_5x4_regime)
            -   [N\_var\_5x4\_regime](#n_var_5x4_regime)
            -   [N\_4x5\_regime](#n_4x5_regime)
            -   [N\_var\_4x5\_regime](#n_var_4x5_regime)
            -   [migrants\_20x1\_regime](#migrants_20x1_regime)
            -   [migrants\_10x2\_regime](#migrants_10x2_regime)
            -   [migrants\_5x4\_regime](#migrants_5x4_regime)
            -   [migrants\_4x5\_regime](#migrants_4x5_regime)
        -   [NBBg-samples](#nbbg-samples)
            -   [NBBg-samples-combined](#nbbg-samples-combined)
            -   [NBBg-samples-chain1](#nbbg-samples-chain1)
            -   [NBBg-samples-chain2](#nbbg-samples-chain2)
            -   [NBBg-samples-chain3](#nbbg-samples-chain3)
        -   [attributes](#attributes)
        -   [clean-establishment-data](#clean-establishment-data)
        -   [initial-density-dependence](#initial-density-dependence)
        -   [measurement-error](#measurement-error-1)
        -   [melbourne-ricker-data](#melbourne-ricker-data)
        -   [Tribolium-propagule-pressure-data](#tribolium-propagule-pressure-data)
        -   [Tribolium-source-population-data](#tribolium-source-population-data)
    -   [figures](#figures)
        -   [establishment-probability-experiment-and-simulations](#establishment-probability-experiment-and-simulations)
        -   [experiment-time-series-establishment-proportion-absolute-time-type-seven-generations](#experiment-time-series-establishment-proportion-absolute-time-type-seven-generations)
        -   [experiment-time-series-establishment-proportion-absolute-time-type](#experiment-time-series-establishment-proportion-absolute-time-type)
        -   [experiment-time-series-establishment-proportion-relative-time-type](#experiment-time-series-establishment-proportion-relative-time-type)
        -   [experiment-time-series-population-size-seven-generations-environment-facet-with-extinctions](#experiment-time-series-population-size-seven-generations-environment-facet-with-extinctions)
        -   [experiment-time-series-population-size-seven-generations-environment-facet](#experiment-time-series-population-size-seven-generations-environment-facet)
        -   [experiment-time-series-population-size-seven-generations](#experiment-time-series-population-size-seven-generations)
        -   [experiment-time-series-population-size](#experiment-time-series-population-size)
        -   [population-size-experiment-and-simulations](#population-size-experiment-and-simulations)
    -   [written-notes](#written-notes)
        -   [NBBg-model-specification](#nbbg-model-specification)

Abstract
========

Predicting whether individuals will colonize a novel habitat is of fundamental ecological interest and is crucial to both conservation efforts and invasive species management. The only consistent predictor of colonization success is the number of individuals introduced, also called propagule pressure. Propagule pressure increases with the number of introductions and the number of individuals per introduction (the size of the introduction), but it is unresolved which process is a stronger driver of colonization success. Furthermore their relative importance may depend upon the environment, with multiple introductions potentially enhancing colonization of fluctuating environments. To evaluate the relative importance of the number and size of introductions and its dependence upon environmental variability, we paired demographic simulations with a microcosm experiment. Using Tribolium flour beetles as a model system, we introduced a fixed number of individuals into replicated novel habitats of stable or fluctuating quality, varying the number of introductions through time and size of each introduction. We evaluated establishment probability and the size of extant populations after 7 generations. In the simulations and microcosms, we found that establishment probability increased with more, smaller introductions, but was not affected by biologically realistic fluctuations in environmental quality. Population size was not significantly affected by environmental variability in the simulations, but populations in the microcosms grew larger in a stable environment, especially with more introduction events. In general, the microcosm experiment yielded higher establishment probability and larger populations than the demographic simulations. We suggest that genetic mechanisms likely underlie these differences and thus deserve more attention in efforts to parse propagule pressure. Our results highlight the importance of preventing further introductions of undesirable species to invaded sites, and suggest conservation efforts should focus on increasing the number of introductions or re-introductions of desirable species rather than increasing the size of those introduction events.

Keywords
========

invasion, reintroduction, propagule pressure, biocontrol, conservation, microcosm, population dynamics, stochasticity, simulation

General Methods Framework
=========================

In simulations and a microcosm experiment, we evaluated the outcome of introducing 20 total individuals to one of two environmental contexts (a stable or fluctuating novel environment), varying the number of introduction events used to distribute those individuals through time. This total was low enough to allow some population extinction within a reasonable timeframe, but high enough to be representative of documented introductions in the literature (Simberloff 1989, Simberloff 2009; Grevstad 1999; Berggren 2001; Taylor et al. 2005; Drake et al. 2005). The introduction regimes were: 20 individuals introduced in the first generation, 10 individuals introduced in each of the first 2 generations, 5 individuals introduced in each of the first 4 generations, and 4 individuals introduced in each of the first 5 generations. To create the fluctuating environment, we imposed a magnitude of variability corresponding to environmental stochasticity in nature, which leads to frequent, mild to moderate perturbations in population growth rate due to external forces (sensu Lande 1993). Populations were tracked for seven generations following the initial introduction.

Establishment probability and the size of established populations were used as measures of colonization success. Populations were deemed ‘established’ for any time step in which they were extant and population size was noted for all extant populations in every time step. For analysis, establishment probability and mean population size were assessed in generation 7.

Files in Repository
===================

This repository (`mikoontz\ppp-establishment`) represents the data and analysis from the simulations and the experiment.

scripts
-------

### data-carpentry

These scripts clean and fortify the data to make it ready for analyses.

#### generate-establishment-responses-for-analysis

File type: .R

Takes short form data and appends population attributes (e.g. treatments, notes, presence of introduction gap) as well as various response variables of interest (established/extinct for each fixed generation, established/extinct for each relative generation after introductions completed) for further analysis.

#### generate-tidy-data

File type: .R

Converts long form census data (one row represents the census for each population/generation combination) into short form (one row represents censuses for each population with different generations in different columns) for more intuitive access to data.

### figure-carpentry

Scripts that generate the figures in the paper

#### experiment-time-series-population-size-plot

File type: .R

Script to generate the <a href="#experiment-time-series-population-size">experiment-time-series-population-size</a> figure.

#### experiment-time-series-establishment-proportion-plot

File type: .R

Generates the <a href="#experiment-time-series-establishment-proportion-absolute-time-type">experiment-time-series-establishment-proportion-absolute-time-type</a> and <a href="#experiment-time-series-establishment-proportion-relative-time-type.tif">experiment-time-series-establishment-proportion-relative-time-type</a> figures.

#### establishment-and-size-experiment-and-simulation-results-plot

File type: .R

Generates the <a href="#establishment-probability-experiment-and-simulations.tif">establishment-probability-experiment-and-simulations.tif</a> and <a href="#population-size-experiment-and-simulations.tif">population-size-experiment-and-simulations.tif</a> figures.

### simulations

Scripts related to the simulation component to the project.

#### NBBg-population-dynamics-function

File type: .R

Simulates N<sub>t+1</sub> given N<sub>t</sub> and the values of 4 key parameters of the NBBg model (R<sub>0</sub>, *α*, kE, and kD).

#### NBBg-simulation-functions

File type: .R

Helper functions to setting up the simulation of the NBBg model using estimated parameter values (and their distributions) from a separate fitted model.

#### NBBg-simulation-script

File type: .R

Script runs the simulation functions and summarizes output.

### environmental-stochasticity

File type: .R

Calculates total stochasticity through time using lambda values for each population that remained extant throughout entire experiment. Uses a linear mixed effects model to assess the effect of environmental treatment (stable or fluctuating) on this total stochasticity value. We'd expect our fluctuating treatment to impose more variation in lambda values.

### establishment-probability

File type: .R

Mixed effects logistic regression analyses with establishment as response (1 or 0) with introduction regime and environmental stability as covariates. The script includes several analyses, all with a check to see whether the introduction gap in generation F2 was important.

Assessments of establishment or extinction were made at fixed time points throughout the experiment (generations F5, F6, F7, F8, and F9), at a relative time point (5 generations after the final introduction event for each introduction regime), and at each time point with the population ID as a random effect in a repeated measures framework.

### measurement-error

File type: .R

Analysis of measurement error as a function of observer and size of population being censused.

### NBBg-NIMBLE

File type: .R

Metropolis-Hastings algorithm to estimate parameters in NBBg model from Melbourne and Hastings (2008) using the `nimble` package.

### NBBg-script-validation

File type: .R

Test of the NBBg `nimble` code on simulated dataset with known parameters and on a published dataset with previously estimated parameters (Melbourne and Hastings, 2008).

### NBBg-script

File type: .R

Script to run `nimble` algorithm, assess MCMC diagnostics, and write samples to a file for use in simulations.

### population-size

File type: .R

Mixed effects Poisson regression analyses with population size as response with introduction regime and environmental stability as covariates.The script includes several analyses, all with a check to see whether the introduction gap in generation F2 was important.

Assessments of population size were made at fixed time points throughout the experiment (generations F5, F6, F7, F8, and F9) and at a relative time point (5 generations after the final introduction event for each introduction regime).

### temporary-extinctions

File type: .R

Assessment of whether temporary extinctions in the multiply-introduced populations affected establishment probability or population size. Tests effect of temporary extinction versus no temporary extinction (categorical variable with 2 levels) as well as effect of "amount of loss" (continuous variable) representing how many of the 20 possible individuals introduced did not have a genetic contribution to the population when response was assessed. That is, the "amount of loss" represents how much of a bottleneck the population passed through during introduction.

data
----

### simulations

These files are broken down by introduction regime to keep them small enough for easy upload/download from GitHub.

#### simulation\_stats\_tidy

File type: .csv

Tidy summary output from simulation runs of 500,000 replications per treatment (4 introduction regimes \* 2 environment stabilities)

#### simulation\_stats\_raw

File type: .csv

Raw summary output from simulation runs of 500,000 replications per treatment (4 introduction regimes \* 2 environment stabilities)

#### N\_20x1\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 20x1 introduction regime into a relatively stable novel environment

#### N\_var\_20x1\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 20x1 introduction regime into a fluctuating novel environment

#### N\_10x2\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 10x2 introduction regime into a relatively stable novel environment

#### N\_var\_10x2\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 10x2 introduction regime into a fluctuating novel environment

#### N\_5x4\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 5x4 introduction regime into a relatively stable novel environment

#### N\_var\_5x4\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 5x4 introduction regime into a fluctuating novel environment

#### N\_4x5\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 4x5 introduction regime into a relatively stable novel environment

#### N\_var\_4x5\_regime

File type: .csv

Population trajectories for 500,000 simulated *Tribolium* populations introduced via 4x5 introduction regime into a fluctuating novel environment

#### migrants\_20x1\_regime

File type: .csv

Matrix representing the migrant *Tribolium* individuals into the populations introduced via 20x1 introduction regime.

#### migrants\_10x2\_regime

File type: .csv

Matrix representing the migrant *Tribolium* individuals into the populations introduced via 10x2 introduction regime.

#### migrants\_5x4\_regime

File type: .csv

Matrix representing the migrant *Tribolium* individuals into the populations introduced via 5x4 introduction regime.

#### migrants\_4x5\_regime

File type: .csv

Matrix representing the migrant *Tribolium* individuals into the populations introduced via 4x5 introduction regime.

### NBBg-samples

#### NBBg-samples-combined

File type: .csv

Samples from parameter posterior distribution from NBBg model fit. All converged chains combined into a single matrix-like object.

#### NBBg-samples-chain1

File type: .csv

Chain 1 samples from parameter posterior distribution from NBBg model fit. Separate chains could be retrieved and diagnostics run.

#### NBBg-samples-chain2

File type: .csv

Chain 2 samples from parameter posterior distribution from NBBg model fit.

#### NBBg-samples-chain3

File type: .csv

Chain 3 samples from parameter posterior distribution from NBBg model fit.

### attributes

File type: .csv

These are the population attributes that describe each population. They can be easily merged (aka joined) with other population information (*e.g.* response variables) by the unique identifier "ID".

Column descriptions:

ID: unique population identifier

block: which of the 4 temporal blocks the population belonged to (1, 2, 3, or 4)

color: the color code of the block (blue, green, pink, or yellow)

number: the number of introductions in the introduction regime of the population (1, 2, 4, or 5)

size: the size of each introduction in the introduction regime of the population (20, 10, 5, or 4)

environment: the environment context (stable or fluctuating) of the population

special: was there anything idiosyncratic about this population?

gap: was there an introduction gap during the multiple introductions (these populations were excluded from analysis)

notes: more detailed notes on idiosyncracies

drought: was there a dramatic loss of relative humidity between generations 2 and 3 due to total evaporation of the water trays? (turned out to not have any detectable effect by any measure)

### clean-establishment-data

File type: .csv

The data used for establishment and population size analysis, which was generated with the <a href="#generate-establishment-responses-for-analysis.r">generate-establishment-responses-for-analysis.R</a> script.

### initial-density-dependence

File type: .csv

Data for estimating parameters in the NBBg stochastic population dynamics model. Represents 125 populations with a starting population size ranging from 5 to 200 individuals, and a full census at the next time step.

### measurement-error

File type: .csv

Data from experiment to quantify measurement error during census as a function of observer and size of the population being censused.

### melbourne-ricker-data

File type: .csv

Census data from Melbourne and Hastings (2008) used to validate our implementation of the NBBg model using the NIMBLE language.

### Tribolium-propagule-pressure-data

File type: .csv

Raw, long-form census data from parsing propagule pressure experiment. Each row represents the census value (and other attributes) for a single population in a single generation.

Column descriptions:

ID: unique population identifier

Week: During which experiment week did the census take place? There were 2 blocks of surveys per week for 2 weeks in a row, then no censuses for 3 weeks, then the next generation's censuses start

Environment: The percent of the growth medium that was the natal medium (95% wheat flour, 5% corn flour)

Generation: To which generation did the censused adults belong? (e.g. initial introduction would be 0, offspring from initial introduction are generation 1, etc.)

N0: The number of adulte beetles from the previous generation

Census: The number of adult beetles in the current generation

Person: Initials of the person doing the census

Setup.Order: Subsets of populations (consisting of no more than 10 populations each) were censused in a random order in each block; this number represents the position of the subset of populations in that census order to which the population belonged

Addition: How many beetles will be added to the population this generation from the external source population

Last.addition: The number of beetles introduced to the population in the previous generation

Proportion.new: What proportion of the total beetles from the previous generation were introduced from the external source population

### Tribolium-source-population-data

File type: .csv

Raw, long-form census data from the source population maintenance. Note that not all patches were censused each generation, just a sample to monitor population growth rates.

Column descriptions:

ID: convenience identifier; these are NOT unique and did NOT correspond to the same population throughout the course of the maintenance of this population. All IDs within a block were mixed together after census, before selecting colonists, and before refounding the source population for incubation

Week: During which experiment week did the census take place? There were 2 blocks of surveys per week for 2 weeks in a row, then no censuses for 3 weeks, then the next generation's censuses start

Block: which of the 4 temporal blocks the population was in

Color: the color code for each block

Environment: The percent of the growth medium that was the natal medium (95% wheat flour, 5% corn flour)

Generation: To which generation did the censused adults belong? (e.g. initial introduction would be 0, offspring from initial introduction are generation 1, etc.)

N0: The number of adulte beetles from the previous generation

Census: The number of adult beetles in the current generation

Check: Data entry check on the Census column

Person: Initials of the person doing the census

notes: Anything out of ordinary for the census

figures
-------

### establishment-probability-experiment-and-simulations

File type: .pdf

Plot of establishment probability results from experiment and simulations assessed at generation 7.

### experiment-time-series-establishment-proportion-absolute-time-type-seven-generations

File type: .pdf

Establishment proportion of experimental populations as a function of the total number of generations elapsed in the experiment through generation 7.

### experiment-time-series-establishment-proportion-absolute-time-type

File type: .pdf

Establishment proportion of experimental populations as a function of the total number of generations elapsed in the experiment through generation 9.

### experiment-time-series-establishment-proportion-relative-time-type

File type: .pdf

Establishment proportion of experimental populations as a function of number of generations since the final introduction event.

### experiment-time-series-population-size-seven-generations-environment-facet-with-extinctions

File type: .pdf

Population size through time for all 842 populations faceted by environment treatment. Mean population size and standard error are shown in each facet for each introduction regime. The population size trajectories for populations that go extinct by generation 7 are highlighted in orange.

### experiment-time-series-population-size-seven-generations-environment-facet

File type: .pdf

Population size through time for all 842 populations faceted by environment treatment. Mean population size and standard error are shown in each facet for each introduction regime.

### experiment-time-series-population-size-seven-generations

File type: .pdf

The population trajectories for the 842 populations that didn't experience a gap in the introduction regime between generations 1 and 2. (That is, no population augmentation of generation 1 adults). Also plots mean population size and a 1 standard error envelope around that mean for each introduction regime. Plot through generation 7.

### experiment-time-series-population-size

File type: .pdf

The population trajectories for the 842 populations that didn't experience a gap in the introduction regime between generations 1 and 2. (That is, no population augmentation of generation 1 adults). Also plots mean population size and a 1 standard error envelope around that mean for each introduction regime. Plot through generation 9.

### population-size-experiment-and-simulations

File type: .pdf

Plot of population size results from experiment and simulations assessed at generation 7.

written-notes
-------------

### NBBg-model-specification

File type: .Rmd

The specification for the NBBg stochastic hierarchical population dynamics model as well as the specification of the full conditional distributions for each estimated parameter. Other files by the same name with different extensions are part of the .pdf rendering process.
