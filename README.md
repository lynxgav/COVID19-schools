# COVID19-schools

## DOI for this GitHub repository
doi: https://doi.org/10.5281/zenodo.4541431

## Overview

We investigate the impact of reducing school and other (non-school-related) contacts on COVID-19 control using an age-structured model fitted to age-specific seroprevalence and hospital admission data from the Netherlands. 

The analyses were published as

> Rozhnova G, van Dorp CH, Bruijning-Verhagen P, Bootsma MCJ, van de Wijgert JHHM, Bonten MJM, Kretzschmar ME. Model-based evaluation of school- and non-school-related measures to control the COVID-19 pandemic. Nature Communications. 2021;12(1):1614. https://doi.org/10.1038/s41467-021-21899-6.

## Data

The data are added to the `data` folder for convenience.

### Seroprevalence data

We digitized data from the publication

> Vos ERA, den Hartog G, Schepp RM, *et al*
  Nationwide seroprevalence of SARS-CoV-2 and identification of risk factors in the general population of the Netherlands during the first epidemic wave
  *J Epidemiol Community Health* Published Online First: 28 November 2020. doi: https://doi.org/10.1136/jech-2020-215678

### Contact matrices

We used contact patterns from the preprint

> Jantien A. Backer, Liesbeth Mollema, R.A. Eric Vos, Don Klinkenberg, Fiona R.M. van der Klis, Hester E. de Melker, Susan van den Hof, Jacco Wallinga
  The impact of physical distancing measures against COVID-19 transmission on contacts and mixing patterns in the Netherlands: 
  repeated cross-sectional surveys in 2016/2017, April 2020 and June 2020
  *medRxiv* 2020.05.18.20101501; doi: [https://doi.org/10.1101/2020.05.18.20101501](https://doi.org/10.1101/2020.05.18.20101501)
  
The school-specific contact matrix was taken from the publication 

> Prem K, Cook AR, Jit M. Projecting social contact matrices in 152 countries using contact surveys and demographic data. *PLOS Computational Biology* 2017;13(9):1-21. doi: https://doi.org/10.1371/journal.pcbi.1005697

### Demographic data

We used publicly available data from the Statistics Netherlands (CBS): https://www.cbs.nl.

### Hospitalization data

The hospital data included COVID-19 hospitalizations by date of admission and stratified by age in the Netherlands (OSIRIS database). 

## Model 

### Inference

Model inference was done with R Version 3.6.0 using R Studio Version 1.3.959 (Interface to R) and Stan using rstan R package Version 2.19.3 (R interface to Stan) and cmdstanr R package Version 0.1.3 on Windows 10 Home Version 1903. The scripts can be found in the `scripts` directory. The R and Stan scripts are based on scripts used for the publication

> van Boven M, Teirlinck AC, Meijer A, Hooiveld M, van Dorp CH, Reeves RM, Campbell H, van der Hoek W; RESCEU Investigators. 
  Estimating Transmission Parameters for Respiratory Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination. 
  *J Infect Dis.* 2020 Oct 7;222(Supplement_7):S688-S694. doi: https://doi.org/10.1093/infdis/jiaa424

### Analysis

Analysis of the model was performed using Mathematica Version Number 10.0.2.0 on Platform Mac OS X El Capitan Version 10.11.5. The notebook SchoolAnalyses.nb can be found in the `notebooks` directory.

### Figures

Figures for the manuscript were produced in the notebook SchoolAnalyses.nb. Figures can be found in the `figures` directory.
Notebook JointPosteriorPlot.ipynb produces Figure S5 of the manuscript (posterior correlations between the parameters). 

### Output

Output files produced in R or Mathematica can be found in the `output` directory.

## Other

### OS requirements

Mac OS X El Capitan Version 10.11.5 for Mathematica notebook<br />
Windows 10 Home Version 1903 for R and Stant scripts

### Hardware requirements

Our study requires only a standard computer with enough RAM to support the in-memory operations.

### Installation guide

R Version 3.6.0 https://www.r-project.org/<br />
R Studio Version 1.3.959 (Interface to R) https://rstudio.com/<br />
rstan R package Version 2.19.3 (R interface to Stan) https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html<br />
cmdstanr R package Version 0.1.3 on Windows 10 Home Version 1903 https://mc-stan.org/cmdstanr/<br />
Mathematica 10.0.2.0 https://www.wolfram.com/mathematica/<br />
