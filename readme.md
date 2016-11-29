Table of Contents for ppp-establishment Repository
================
Michael Koontz
11/27/2016

-   [Introduction](#introduction)
-   [Files in Repository](#files-in-repository)
    -   [scripts](#scripts)
        -   [data-carpentry](#data-carpentry)
            -   [generate-establishment-responses-for-analysis.R](#generate-establishment-responses-for-analysis.r)
            -   [generate-tidy-data.R](#generate-tidy-data.r)
        -   [simulations](#simulations)
            -   [NBBg-population-dynamics-function.R](#nbbg-population-dynamics-function.r)
            -   [NBBg-simulation-functions.R](#nbbg-simulation-functions.r)
            -   [NBBg-simulation-script.R](#nbbg-simulation-script.r)
        -   [environmental-stochasticity.R](#environmental-stochasticity.r)
        -   [establishment-probability.R](#establishment-probability.r)
        -   [experiment-time-series-population-abundance-plot.R](#experiment-time-series-population-abundance-plot.r)
        -   [population-abundance.R](#population-abundance.r)
        -   [temporary-extinctions.R](#temporary-extinctions.r)
        -   [NBBg-mcmc.R](#nbbg-mcmc.r)
        -   [NBBg-script.R](#nbbg-script.r)
        -   [NBBg-script-validation.R](#nbbg-script-validation.r)
        -   [NBBg-environment-stability-effect.R](#nbbg-environment-stability-effect.r)
    -   [data](#data)
        -   [simulations](#simulations-1)
            -   [simulation\_stats.csv](#simulation_stats.csv)
            -   [N.csv](#n.csv)
            -   [N\_var.csv](#n_var.csv)
        -   [NBBg-samples](#nbbg-samples)
            -   [NBBg-samples-combined.csv](#nbbg-samples-combined.csv)
            -   [NBBg-samples-chain1.csv](#nbbg-samples-chain1.csv)
            -   [NBBg-samples-chain2.csv](#nbbg-samples-chain2.csv)
            -   [NBBg-samples-chain3.csv](#nbbg-samples-chain3.csv)
        -   [attributes.csv](#attributes.csv)
        -   [clean-establishment-data.csv](#clean-establishment-data.csv)
        -   [initial-density-dependence.csv](#initial-density-dependence.csv)
        -   [Tribolium-propagule-pressure-data.csv](#tribolium-propagule-pressure-data.csv)
    -   [figures](#figures)
        -   [experiment-time-series-population-abundance.tif](#experiment-time-series-population-abundance.tif)
    -   [written-notes](#written-notes)
        -   [NBBg-model-specification.tex](#nbbg-model-specification.tex)

Introduction
============

We introduced 917 populations of *Tribolium* flour beetles comprising 20 individuals apiece to novel environments that were either stable or randomly fluctuating through time, varying the number of introduction events used to distribute them (1, 2, 4, or 5 events). Thus, there were 4 levels of "introduction regime" (20 individuals introduced in the 1st generation, 10 individuals introduced in the first 2 generations, 5 individuals introduced in the first 4 generations, and 4 individuals introduced in the first 5 generations), and 2 levels of "environmental stability" (stable or randomly fluctuating through time).

Our goal was to investigate how the introduction regime and environmental stability of the recipient environment might affect establishment probability and the population abundance of the established populations.

This repository (`mikoontz\ppp-establishment`) represents the data and analysis from the experiment.

Files in Repository
===================

scripts
-------

### data-carpentry

These scripts clean and fortify the data to make it ready for analyses.

#### generate-establishment-responses-for-analysis.R

Takes short form data and appends population attributes (e.g. treatments, notes, presence of introduction gap) as well as various response variables of interest (established/extinct for each fixed generation, established/extinct for each relative generation after introductions completed) for further analysis.

#### generate-tidy-data.R

Converts long form census data (one row represents the census for each population/generation combination) into short form (one row represents censuses for each population with different generations in different columns) for more intuitive access to data.

### simulations

#### NBBg-population-dynamics-function.R

Simulates N<sub>t+1</sub> given N<sub>t</sub> and the values of 4 key parameters of the NBBg model (R<sub>0</sub>, \(\alpha\), kE, and kD).

#### NBBg-simulation-functions.R

Functions to simulate N<sub>t+1</sub> given N<sub>t</sub> and the values of 4 key parameters of the NBBg model (R<sub>0</sub>, \(\alpha\), kE, and kD) and also allows incorporation of model uncertainty.

#### NBBg-simulation-script.R

Script runs the simulation functions and summarizes output.

### environmental-stochasticity.R

Calculates total stochasticity through time using lambda values for each population that remained extant throughout entire experiment. Uses a linear mixed effects model to assess the effect of environmental treatment (stable or fluctuating) on this total stochasticity value. We'd expect our fluctuating treatment to impose more variation in lambda values.

### establishment-probability.R

Mixed effects logistic regression analyses with establishment as response (1 or 0) with introduction regime and environmental stability as covariates. The script includes several analyses, all with a check to see whether the introduction gap in generation F2 was important.

Assessments of establishment or extinction were made at fixed time points throughout the experiment (generations F5, F6, F7, F8, and F9), at a relative time point (5 generations after the final introduction event for each introduction regime), and at each time point with the population ID as a random effect in a repeated measures framework.

### experiment-time-series-population-abundance-plot.R

Script to generate the <a href="#experiment-time-series-population-abundance.tif">experiment-time-series-population-abundance.tif</a> figure.

### population-abundance.R

Mixed effects Poisson regression analyses with population abundance as response with introduction regime and environmental stability as covariates.The script includes several analyses, all with a check to see whether the introduction gap in generation F2 was important.

Assessments of population size were made at fixed time points throughout the experiment (generations F5, F6, F7, F8, and F9) and at a relative time point (5 generations after the final introduction event for each introduction regime).

### temporary-extinctions.R

Assessment of whether temporary extinctions in the multiply-introduced populations affected establishment probability or population size. Tests effect of temporary extinction versus no temporary extinction (categorical variable with 2 levels) as well as effect of "amount of loss" (continuous variable) representing how many of the 20 possible individuals introduced did not have a genetic contribution to the population when response was assessed. That is, the "amount of loss" represents how much of a bottleneck the population passed through during introduction.

### NBBg-mcmc.R

Metropolis-Hastings algorithm to estimate parameters in NBBg model from Melbourne and Hastings (2008).

### NBBg-script.R

Script to run Metropolis-Hastings algorithm, assess MCMC diagnostics, and write samples to a file for use in simulations.

### NBBg-script-validation.R

Test of the NBBg MCMC code on simulated dataset with known parameters and on a published dataset with previously estimated parameters (Melbourne and Hastings, 2008).

### NBBg-environment-stability-effect.R

Uses the NBBg MCMC code to estimate parameters of the NBBg model fit to the raw experimental data from Generation P to Generation F1 in the fluctuating environment versus the stable environment. We would expect the kE parameter to be larger in the stable environment, indicating less variability in N<sub>t+1</sub> arising from environmental stochasticity.

data
----

### simulations

#### simulation\_stats.csv

Summary output from simulation runs of 500,000 replications per treatment (4 introduction regimes \* 2 environment stabilities)

#### N.csv

Population trajectories for 2,000,000 simulated *Tribolium* populations introduced via 1 of 4 different introduction regimes into a stable environment.

#### N\_var.csv

Population trajectories for 2,000,000 simulated *Tribolium* populations introduced via 1 of 4 different introduction regimes into a fluctuating environment.

### NBBg-samples

#### NBBg-samples-combined.csv

Samples from parameter posterior distribution from NBBg model fit. All converged chains combined into a single matrix-like object.

#### NBBg-samples-chain1.csv

Chain 1 samples from parameter posterior distribution from NBBg model fit. Separate chains could be retrieved and diagnostics run.

#### NBBg-samples-chain2.csv

Chain 2 samples from parameter posterior distribution from NBBg model fit.

#### NBBg-samples-chain3.csv

Chain 3 samples from parameter posterior distribution from NBBg model fit.

### attributes.csv

These are the population attributes that describe each population. They can be easily merged (aka joined) with other population information (*e.g.* response variables) by the unique identifier "ID".

### clean-establishment-data.csv

The data used for establishment and population abundance analysis, which was generated with the <a href="#generate-establishment-responses-for-analysis.r">generate-establishment-responses-for-analysis.R</a> script.

### initial-density-dependence.csv

Data for estimating parameters in the NBBg stochastic population dynamics model. Represents 125 populations with a starting population size ranging from 5 to 200 individuals, and a full census at the next time step.

### Tribolium-propagule-pressure-data.csv

Raw, long-form census data from parsing propagule pressure experiment. Each row represents the census value (and other attributes) for a single population in a single generation.

figures
-------

### experiment-time-series-population-abundance.tif

The population trajectories for the 842 populations that didn't experience a gap in the introduction regime between generations F1 and F2. (That is, no population augmentation of F1 adults). Also plots mean population abundance and a 1 standard error envelope around that mean for each introduction regime.

written-notes
-------------

### NBBg-model-specification.tex

The specification for the NBBg stochastic hierarchical population dynamics model as well as the specification of the full conditional distributions for each estimated parameter. Other files by the same name with different extensions are part of the .pdf rendering process.
