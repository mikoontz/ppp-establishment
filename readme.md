Table of Contents for `ppp-establishment` Repository
================
Michael Koontz
January 27, 2017

-   [Introduction](#introduction)
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
        -   [melbourne\_ricker\_data](#melbourne_ricker_data)
        -   [Tribolium-propagule-pressure-data](#tribolium-propagule-pressure-data)
    -   [figures](#figures)
        -   [experiment-time-series-population-size](#experiment-time-series-population-size)
        -   [experiment-time-series-population-size-seven-generations](#experiment-time-series-population-size-seven-generations)
        -   [experiment-time-series-establishment-proportion-absolute-time-type](#experiment-time-series-establishment-proportion-absolute-time-type)
        -   [experiment-time-series-establishment-proportion-absolute-time-type-seven-generations](#experiment-time-series-establishment-proportion-absolute-time-type-seven-generations)
        -   [experiment-time-series-establishment-proportion-relative-time-type](#experiment-time-series-establishment-proportion-relative-time-type)
        -   [establishment-probability-experiment-and-simulations](#establishment-probability-experiment-and-simulations)
        -   [population-size-experiment-and-simulations](#population-size-experiment-and-simulations)
    -   [written-notes](#written-notes)
        -   [NBBg-model-specification](#nbbg-model-specification)

Introduction
============

We introduced 917 populations of *Tribolium* flour beetles comprising 20 individuals apiece to novel environments that were either stable or randomly fluctuating through time, varying the number of introduction events used to distribute them (1, 2, 4, or 5 events). Thus, there were 4 levels of "introduction regime" (20 individuals introduced in the 1st generation, 10 individuals introduced in the first 2 generations, 5 individuals introduced in the first 4 generations, and 4 individuals introduced in the first 5 generations), and 2 levels of "environmental stability" (stable or randomly fluctuating through time).

Our goal was to investigate how the introduction regime and environmental stability of the recipient environment might affect establishment probability and the population size of the established populations.

This repository (`mikoontz\ppp-establishment`) represents the data and analysis from the experiment.

Files in Repository
===================

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

Simulates N<sub>t+1</sub> given N<sub>t</sub> and the values of 4 key parameters of the NBBg model (R<sub>0</sub>, \(\alpha\), kE, and kD).

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

### clean-establishment-data

File type: .csv

The data used for establishment and population size analysis, which was generated with the <a href="#generate-establishment-responses-for-analysis.r">generate-establishment-responses-for-analysis.R</a> script.

### initial-density-dependence

File type: .csv

Data for estimating parameters in the NBBg stochastic population dynamics model. Represents 125 populations with a starting population size ranging from 5 to 200 individuals, and a full census at the next time step.

### measurement-error

File type: .csv

Data from experiment to quantify measurement error during census as a function of observer and size of the population being censused.

### melbourne\_ricker\_data

File type: .csv

Census data from Melbourne and Hastings (2008) used to validate our implementation of the NBBg model using the NIMBLE language.

### Tribolium-propagule-pressure-data

File type: .csv

Raw, long-form census data from parsing propagule pressure experiment. Each row represents the census value (and other attributes) for a single population in a single generation.

figures
-------

### experiment-time-series-population-size

File type: .pdf

The population trajectories for the 842 populations that didn't experience a gap in the introduction regime between generations 1 and 2. (That is, no population augmentation of generation 1 adults). Also plots mean population size and a 1 standard error envelope around that mean for each introduction regime. Plot through generation 9.

### experiment-time-series-population-size-seven-generations

File type: .pdf

The population trajectories for the 842 populations that didn't experience a gap in the introduction regime between generations 1 and 2. (That is, no population augmentation of generation 1 adults). Also plots mean population size and a 1 standard error envelope around that mean for each introduction regime. Plot through generation 7.

### experiment-time-series-establishment-proportion-absolute-time-type

File type: .pdf

Establishment proportion of experimental populations as a function of the total number of generations elapsed in the experiment through generation 9.

### experiment-time-series-establishment-proportion-absolute-time-type-seven-generations

File type: .pdf

Establishment proportion of experimental populations as a function of the total number of generations elapsed in the experiment through generation 7.

### experiment-time-series-establishment-proportion-relative-time-type

File type: .pdf

Establishment proportion of experimental populations as a function of number of generations since the final introduction event.

### establishment-probability-experiment-and-simulations

File type: .pdf

Plot of establishment probability results from experiment and simulations assessed at generation 7.

### population-size-experiment-and-simulations

File type: .pdf

Plot of population size results from experiment and simulations assessed at generation 7.

written-notes
-------------

### NBBg-model-specification

File type: .Rmd

The specification for the NBBg stochastic hierarchical population dynamics model as well as the specification of the full conditional distributions for each estimated parameter. Other files by the same name with different extensions are part of the .pdf rendering process.
