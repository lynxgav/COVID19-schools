# COVID19-schools

We used an age-structured model fitted to age-specific seroprevalence and hospital admission data from the Netherlands to investigate the impact of reducing school and other (non-school-related) contacts on controlling the COVID-19 pandemic.

All details are described in the preprint

> Rozhnova G, van Dorp CH, Bruijning-Verhagen, P, Bootsma MCJ, van de Wijgert JHHM, Bonten MJM, Kretzschmar ME,
  Model-based evaluation of school- and non-school-related measures to control the COVID-19 pandemic


## Data

The data are added to the `data` folder for convenience.

### Seroprevalence data

We digitized data from the publication

> Vos ERA, den Hartog G, Schepp RM, *et al*
  Nationwide seroprevalence of SARS-CoV-2 and identification of risk factors in the general population of the Netherlands during the first epidemic wave
  *J Epidemiol Community Health* Published Online First: 28 November 2020. doi: [10.1136/jech-2020-215678](10.1136/jech-2020-215678)

### Contact matrices

We used contact patterns from the preprint:

> Jantien A. Backer, Liesbeth Mollema, R.A. Eric Vos, Don Klinkenberg, Fiona R.M. van der Klis, Hester E. de Melker, Susan van den Hof, Jacco Wallinga
  The impact of physical distancing measures against COVID-19 transmission on contacts and mixing patterns in the Netherlands: 
  repeated cross-sectional surveys in 2016/2017, April 2020 and June 2020
  *medRxiv* 2020.05.18.20101501; doi: [https://doi.org/10.1101/2020.05.18.20101501](https://doi.org/10.1101/2020.05.18.20101501)
  
### Hospitalization data

The hospital data included COVID-19 hospitalizations by date of admission and stratified by age during the period of 69 days following the first official case in the Netherlands (OSIRIS database). The applications for this dataset should be forwarded to the National Institute for Public Health and the Environment (COVID-19EPI@rivm.nl).

## Model 

### Inference

Model inference is done with R and Stan. The scripts can be found in the `scripts` directory.
The R and Stan scripts are based on scripts used for the publication

> van Boven M, Teirlinck AC, Meijer A, Hooiveld M, van Dorp CH, Reeves RM, Campbell H, van der Hoek W; RESCEU Investigators. 
  Estimating Transmission Parameters for Respiratory Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination. 
  *J Infect Dis.* 2020 Oct 7;222(Supplement_7):S688-S694. doi: [10.1093/infdis/jiaa424](10.1093/infdis/jiaa424)

### Analysis

Analysis of the model is performed using Mathematica Version Number 10.0.2.0 on Platform Mac OS X El Capitan Version 10.11.5. The notebook SchoolAnalyses.nb can be found in the `notebooks` directory.

