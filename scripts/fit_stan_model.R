## Corona evidence synthesis for Dutch Hospitalization and Seroprevalence data.
##
## This R script is adapted from scripts developed for the publication
## 
## Michiel van Boven, Anne C Teirlinck, Adam Meijer, Mariëtte Hooiveld, 
## Christiaan H van Dorp, Rachel M Reeves, Harry Campbell, Wim van der Hoek, 
## RESCEU Investigators, Estimating Transmission Parameters for Respiratory 
## Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination, 
## The Journal of Infectious Diseases, Volume 222, Issue Supplement_7, 
## 1 November 2020, Pages S688–S694, https://doi.org/10.1093/infdis/jiaa424
##


# load packages
library(tidyverse)
library(cowplot)
library(lubridate)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(readxl)
library(ggplot2)
library(loo) ## model comparison


## install cmdstanr to get access to newer version of Stan
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)
## manually set the path to the cmdstan install
#set_cmdstan_path("~/cmdstan-2.24.1/")
## or install cmdstan
#install_cmdstan()


## Set working directory. NB: not required if we open the R project file!
## in Rstudio: File -> Open Project (or recent projects)
getwd()

# some preliminaries
numberofageclasses <- 10
t0 <- 0

# Load hospitalisations by date of admission
hospitalisations <- read.csv("data/hospitalizations_15062020.csv")
daystofirstcase <- which(hospitalisations$datum_zhopname == "2020-02-27")
daysunperturbed <- which(hospitalisations$datum_zhopname == "2020-03-20") 
daystoend <- which(hospitalisations$datum_zhopname == "2020-04-30") 

date("2020-03-15") - date("2020-02-21")

hospitalisations_selected <- hospitalisations %>%
  slice(daystofirstcase : daystoend) %>%
  select(-c(X, datum_zhopname)) %>%
  as.matrix
dim(hospitalisations_selected)

# Pad with extra days of zero cases to enforce integrator starting low 
# and to enable ODEs approach stable class distribution (I is small)
numberofextradays <- 5
extradays = matrix(0, nrow = numberofextradays, ncol = numberofageclasses)
hospitalisations_selected = rbind(extradays, hospitalisations_selected)
numberofdays <- dim(hospitalisations_selected)[1]

write.csv(hospitalisations_selected, file = "output/hospitalisations_selected.csv", row.names = FALSE)

# Demographic data from the Netherlands
demography <- pull(read_excel("data/DEMOGRAPHY NL NEW.xlsx", skip = 0, sheet = "Population", 
                              col_names = FALSE, col_types = "numeric"), var = 1) # June 2020
sum(demography)

# Contact matrices
ContactData <- read_tsv("data/contact_matrices.tsv")
ContactDataBaselineAll <- ContactData %>% filter(survey == "baseline", contact_type == "all")
ContactDataPhysicalDistancingAll <- ContactData %>% filter(survey == "physical distancing", contact_type == "all")
ContactMatrixBaselineAll <- matrix(ContactDataBaselineAll$m_est, nrow = 10, ncol = 10)
ContactMatrixPhysicalDistancingAll <- matrix(ContactDataPhysicalDistancingAll$m_est, nrow = 10, ncol = 10)
contactmatrix_unperturbed <- ContactMatrixBaselineAll
contactmatrix_distancing <- ContactMatrixPhysicalDistancingAll

## split distancing matrix into household and community contacts

ContactDataPhysicalDistancingHouse <- ContactData %>% filter(survey == "physical distancing", contact_type == "household")
contactmatrix_distancing_house <- matrix(ContactDataPhysicalDistancingHouse$m_est, nrow = 10, ncol = 10)
ContactDataPhysicalDistancingComty <- ContactData %>% filter(survey == "physical distancing", contact_type == "community")
contactmatrix_distancing_comty <- matrix(ContactDataPhysicalDistancingComty$m_est, nrow = 10, ncol = 10)

## In the ODEs, we need C_pc with p the age class of the participant
## and c the age class of the contact:
##
## dS_p / dt = -epsilon * S_p * \sum_c C_pc I_c
##
## R default matrix method uses column-major order and so does Stan.
## The contact matrix file is in row-major-order (contact age groups, then participant age groups)
## Hence, the raw contact matrix has to be transposed

contactmatrix_unperturbed <- t(contactmatrix_unperturbed)
contactmatrix_distancing <- t(contactmatrix_distancing)
## same for split matrix
contactmatrix_distancing_house <- t(contactmatrix_distancing_house)
contactmatrix_distancing_comty <- t(contactmatrix_distancing_comty)

write.csv(contactmatrix_unperturbed, file = "output/contactmatrix_unperturbed.csv", row.names = FALSE)
write.csv(contactmatrix_distancing, file = "output/contactmatrix_distancing.csv", row.names = FALSE)

## import serological data
SeroData <- read_tsv("data/digitized_sero_data.tsv")

## hyper parameters (for the priors)

## use 3 susceptability classes
#<20 : OR=0.23 (0.11,0.46)
#20-59 : OR=0.64 (0.43,0.97)
#60+ : OR=1 (reference class)

beta_short_hyp <- c(0.23, 0.64, 1.0) # 3 classes
#beta_short_hyp <- c(0.23, 0.23, 0.64, 1.0) # 4 classes (with identical prior for [0,10) and [10,20))

susc_classes = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3)
#susc_classes = c(1, 1, 2, 3, 3, 3, 3, 4, 4, 4)

hosp_classes <- c(1, 1, 1, 2, 3, 4, 5, 6, 7, 8)

## hyper parameters for gamma and alpha

## 99% T_I between 5 and 15 days
a_gamma <- 22.6; b_gamma <- 2.44

## 95% T_I between 5.3 and 8.2 days
#a_gamma <- 80; b_gamma <- 12

## 95% T_I between 4.9 and 9.5 days
#a_gamma <- 35; b_gamma <- 5


## 95% T_E between 2.0 and 4.2 days
#a_alpha <- 30; b_alpha <- 10

## 99% T_E between 2 and 5 days
a_alpha <- 32.25; b_alpha <- 9.75


datalist <- list(
  t0 = t0, ## should be smaller than the smallest integration time
  ts_hosp = 1:numberofdays, 
  numdayshosp = numberofdays,
  A = numberofageclasses,
  J = 3, ## Erlang shape parameter for infectious phase
  F = 1, ## Erlang shape parameter for exposed phase
  demography = demography, 
  Cunp = contactmatrix_unperturbed, 
  Cdis = contactmatrix_distancing,
  hospitalisations = hospitalisations_selected, 
  Ahosp = 8,
  hosp_classes = hosp_classes,
  Asusc = 3,
  ref_class = 3,
  susc_classes = susc_classes,
  numdayssero = 1,
  ts_sero = as.array(30), ## TODO: read from file
  sero_num_sampled = matrix(SeroData$SampleSize, nrow=1),
  sero_num_pos = matrix(SeroData$Positive, nrow=1),
  beta_short_hyp = beta_short_hyp,
  log_beta_sd = 0.1,
  a_gamma = a_gamma,
  b_gamma = b_gamma,
  a_alpha = a_alpha,
  b_alpha = b_alpha,
  m_zeta = 1.0, ## we think zeta should be around 1
  s_zeta = 0.1,
  rel_tol = 1e-5,
  abs_tol = 1e-7, 
  max_num_steps = 1e5,
  mode = 0
)


nu_short_scale_guess <- 0.08
nu_short_simplex_guess <- c(0.001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04) ## 8 classes

## make sure that sum(nu_short_simplex_guess) = 1
nu_short_simplex_guess <- nu_short_simplex_guess / sum(nu_short_simplex_guess)

beta_short_raw_guess <- c(0.23, 0.64) ## 3 age classes (one reference class)


## check if the initial parameter guess leads to an R0 that is reasonable
## as an approximation for R0, take rho(diag(S) * C) * epsilon/gamma
## where rho is the spectral radius

calc_R0 <- function(epsilon, beta, gamma, alpha, C) {
  ## FIXME: compute R0 properly!!
  SC <- diag(beta) %*% as.matrix(C)
  eval <- eigen(SC)$values[1]
  R0 <- eval * epsilon / gamma
  return(R0)
}

beta_short_guess <- c(beta_short_raw_guess, 1)
beta_guess <- beta_short_guess[susc_classes]
epsilon_guess <- 0.07
gamma_guess <- 1/4.0
R0_guess <- calc_R0(epsilon_guess, beta_guess, gamma_guess, 
                    alpha_guess, contactmatrix_unperturbed)
print("approx R0 for initial values:")
print(R0_guess)

# initial values 
initials <- function(){
  return(list(
    epsilon = epsilon_guess,
    gamma = gamma_guess,
    zeta = 0.15,
    alpha = 0.33,
    inoculum = 0.0002,
    x0 = 30,
    k = 0.5,
    r = 20,
    nu_scale = nu_short_scale_guess,
    nu_simplex = nu_short_simplex_guess,
    beta_short_raw = beta_short_raw_guess
  ) 
  )
}



## model with Erlang-distributed infectious period
#stan_model_file <- "scripts/stan-models/corona_erlang_hosp_sero_susc.stan"

## model with Erlang-distributed infectious period and neg-binom hosp
stan_model_file <- "scripts/stan-models/corona_erlang_od-hosp_sero_susc.stan"

## model with Erlang-distributed infectious and exposed period and poisson hosp
#stan_model_file <- "scripts/stan-models/corona_erlangEI_hosp_sero_susc.stan"

## model with Erlang-distributed infectious and exposed period and neg-binom hosp
#stan_model_file <- "scripts/stan-models/corona_erlangEI_od-hosp_sero_susc.stan"


# run Stan model

iter <- 1500
warmup <- 1000
thin <- 1

sm <- cmdstan_model(stan_model_file)

fit <- sm$sample(
  data = datalist,
  init = initials,
  chains = 4,
  iter_sampling = iter - warmup,
  iter_warmup = warmup,
  thin = thin,
  refresh = 1,
  adapt_delta = 0.9,
  #fixed_param = TRUE,
  max_treedepth = 15,
  output_dir = "cmdstan-cache", ## make sure this directory exists
  step_size = 0.01
)

## diagnose!

fit$cmdstan_diagnose()


## convert to rstan object to not break the code....

fit.rstan <- rstan::read_stan_csv(fit$output_files())

# print fit
print(fit.rstan, pars=c("epsilon", "gamma", "zeta", "x0", "k", "nu_short", "inoculum", "p_short", "beta_short"), digits = 4)

# traces
rstan::traceplot(fit.rstan, pars=c("epsilon", "gamma", "zeta", "nu_scale", "inoculum", "x0", "k", "beta_short", "r"))
#rstan::traceplot(fit.rstan, pars=c("epsilon", "gamma", "zeta", "nu_scale", "inoculum", "x0", "k", "beta_short"))

# pair plots
pairs(fit.rstan, pars=c("inoculum", "nu_scale", "beta_short_raw", "epsilon", "gamma", "alpha"))

pairs(fit.rstan, pars=c("nu_scale", "nu_simplex"))


# calculate WBIC; note set mode = 1
print(fit.rstan, pars="log_lik", digits=5)

# save and load fit; note: large files
save(fit.rstan, file = "output/fit_J3_wide-prior-alpha-gamma_DSD.rda")

# extract parameters
params = rstan::extract(fit.rstan)

## In case a chain fails, we want to exclude it from the 
## further analysis. This function extracts unwanted chains
## from a parameter list returned by rstan::extract.
## FIXME: there must be a better way!
filter.chains <- function(traces, exclude, iter, warmup, thin) {
  stride <- (iter - warmup) / thin
  ## make a vector of indices we want to exclude
  excl.idxs <- as.vector(sapply(exclude, function(idx) {
    ((idx-1)*stride+1) : (idx*stride)
  }))
  ## for each parameter, remove the excluded samples
  filt.traces <- lapply(traces, function(trace) {
    shape <- dim(trace)
    if ( length(shape) == 1 ) {
      return(trace[-excl.idxs])      
    } else if ( length(shape) == 2 ) {
      return(trace[-excl.idxs,])
    } else if ( length(shape) == 3 ) {
      return(trace[-excl.idxs,,])
    } else if ( length(shape) == 4 ) {
      return(trace[-excl.idxs,,,])
    } else {
      ## FIXME: How do we properly do this in R??
      print("warning: n-dim array with n > 4 not implemented")  
    }
  })
  filt.traces
}

## use filter.chains to remove some chains
#params <- filter.chains(params, c(1,3), iter, warmup, thin) ## exclude chains 1 and 3

# write selected output
output = as.data.frame(fit.rstan, pars=c("epsilon", "zeta", "gamma", "nu_short" , "inoculum", "p_short", "x0", "k", "alpha", "beta_short", "r")) 
write.csv(output, file = "output/ganna_J3_wide-prior-alpha-gamma_DSD.csv", row.names = FALSE)

# select sample closest to MAP
hist(params$lp__, nclass = 30)
max(params$lp__)
posmax <- as.numeric(which(params$lp__ == max(params$lp__)))
MAP <- as.data.frame(fit.rstan)[posmax, c("inoculum", 
                                    "epsilon", 
                                    "gamma")] 
MAP

# color-blind friendly scheme
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# make a figure of hospitalisations and model fit
rm(df_Hospitalizations)
df_Hospitalizations <- params$simulated_hospitalisations
yHosp_sim025 <- matrix(NA, nrow = numberofdays, ncol = numberofageclasses)
yHosp_sim50 <- matrix(NA, nrow = numberofdays, ncol = numberofageclasses)
yHosp_sim975 <- matrix(NA, nrow = numberofdays, ncol = numberofageclasses)
for (t in 1 : numberofdays) {
  for (c in 1 : numberofageclasses) {
    yHosp_sim025[t,c] <- quantile(df_Hospitalizations[,t,c], 0.025)
    yHosp_sim50[t,c] <- quantile(df_Hospitalizations[,t,c], 0.5)
    yHosp_sim975[t,c] <- quantile(df_Hospitalizations[,t,c], 0.975)
  }
}

df_Expect_Hosps <- params$expected_hospitalisations
yHosp_hat025 <- matrix(NA, nrow = numberofdays, ncol = numberofageclasses)
yHosp_hat50 <- matrix(NA, nrow = numberofdays, ncol = numberofageclasses)
yHosp_hat975 <- matrix(NA, nrow = numberofdays, ncol = numberofageclasses)
for (t in 1 : numberofdays) {
  for (c in 1 : numberofageclasses) {
    yHosp_hat025[t,c] <- quantile(df_Expect_Hosps[,t,c], 0.025)
    yHosp_hat50[t,c] <- quantile(df_Expect_Hosps[,t,c], 0.5)
    yHosp_hat975[t,c] <- quantile(df_Expect_Hosps[,t,c], 0.975)
  }
}


ts = 1:numberofdays

dfHosp <- list()

for (j in 1:numberofageclasses) {
  dfHosp[[j]] = data.frame(list(ts = ts,
                                ## simulated values
                                y_sim50=yHosp_sim50[,j],
                                y_sim025=yHosp_sim025[,j],
                                y_sim975=yHosp_sim975[,j],
                                ## expected values
                                y_hat50=yHosp_hat50[,j],
                                y_hat025=yHosp_hat025[,j],
                                y_hat975=yHosp_hat975[,j],
                                ## data
                                hospitalizations = hospitalisations_selected[,j]
  ))
}

## Manually combine age classes

dfHosp[[numberofageclasses+1]] = data.frame(list(ts = ts, 
                                                 ## simulations
                                                 y_sim50=yHosp_sim50[,1]+yHosp_sim50[,2]+yHosp_sim50[,3],
                                                 y_sim025=yHosp_sim025[,1]+yHosp_sim025[,2]+yHosp_sim025[,3],
                                                 y_sim975=yHosp_sim975[,1]+yHosp_sim975[,2]+yHosp_sim975[,3],
                                                 ## expected values
                                                 y_hat50=yHosp_hat50[,1]+yHosp_hat50[,2]+yHosp_hat50[,3],
                                                 y_hat025=yHosp_hat025[,1]+yHosp_hat025[,2]+yHosp_hat025[,3],
                                                 y_hat975=yHosp_hat975[,1]+yHosp_hat975[,2]+yHosp_hat975[,3],
                                                 ## data
                                                 hospitalizations = hospitalisations_selected[,1]+hospitalisations_selected[,2]+
                                                   hospitalisations_selected[,3]
))


# some plot preliminaries
plothosp <- list()
margs <- c(b=1,l=0,t=0,r=0)
size <- 12
maxy <- 250
date.start <- as.Date('2020-02-21')

for (j in 1 : (numberofageclasses+1)) {
  plothosp[[j]] <- ggplot(data = dfHosp[[j]], fill=cbPalette[[j]], mapping = aes(x = ts+date.start)) +
    geom_ribbon(aes(ymin = y_hat025, ymax = y_hat975),  alpha = 0.3) +
    geom_ribbon(aes(ymin = y_sim025, ymax = y_sim975),  alpha = 0.2) +
    geom_line(aes(y=y_hat50),colour = cbPalette[[j]]) +
    geom_point(aes(y = hospitalizations), colour = cbPalette[[j]]) +  
    labs(x="Day", y = "Hospitalisations") + 
    #coord_cartesian(xlim = c(0, numberofdaysunperturbed), ylim = c(0, 100)) +
    theme_bw(base_size = size) + 
    theme(plot.margin = unit(margs, "pt"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())  
}


plot_hosp <- 
  plot_grid(
    plothosp[[1]] + ggtitle("0-4 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[2]] + ggtitle("5-9 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[3]] + ggtitle("10-19 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[4]] + ggtitle("20-29 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[5]] + ggtitle("30-39 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[6]] + ggtitle("40-49 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[7]] + ggtitle("50-59 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[8]] + ggtitle("60-69 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[9]] + ggtitle("70-79 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[10]] + ggtitle("80+ yr") + theme(plot.title = element_text(hjust = 0.5)),
    nrow = 2, 
    ncol= 5
  )

## look at all individual age classes (is pclasses OK?)
plot_hosp

plot_hosp_comb <- 
  plot_grid(
    plothosp[[11]] + ggtitle("0-19 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[4]] + ggtitle("20-29 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[5]] + ggtitle("30-39 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[6]] + ggtitle("40-49 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[7]] + ggtitle("50-59 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[8]] + ggtitle("60-69 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[9]] + ggtitle("70-79 yr") + theme(plot.title = element_text(hjust = 0.5)),
    plothosp[[10]] + ggtitle("80+ yr") + theme(plot.title = element_text(hjust = 0.5)),
    nrow = 2, 
    ncol= 4
  )

plot_hosp_comb

png(filename = "x.png",
    width = 1100, height = 1620, units = "px", pointsize = 12,
    bg = "white",  res = NA, type = c("cairo"))
plot_hosp
dev.off()

###

# logistic function for probability of infection
probinfection = function(t, x0, k){
  return(1 / (1 + exp((-k) * (t - x0))))
} 

# tibble for plotting
n_samples <- length(params$lp__)

probinfection_curves <- tibble(
  sample = 1:n_samples,
  x0 = params$x0,
  k = params$k
)

plot_data <-
  pmap_df(probinfection_curves,
          function(sample, x0, k) {
            tibble(sample = sample,
                   x = seq(10, 40, by = 1),
                   y = probinfection(x, x0, k))
          }) 


# plot estimated curves of prob(infection)
prob_infection_plot <- ggplot(data = plot_data) +
  geom_line(aes(group = sample, x = x, y = y), alpha = .02, color = cbPalette[[1]]) +
  geom_vline(xintercept=23, linetype="dashed", color = "#D55E00", alpha = 1.0) +
  labs(x = "Date", y = "Lockdown weight") + 
  scale_x_continuous(limits = c(10, 40), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

prob_infection_plot


## Make a violin plot of the serological data and model predictions

age.class.names <- SeroData$AgeClass ## [0,5), [5,10), ...

df_sero <- as.data.frame(params$expected_serodata[,1,])
names(df_sero) <- age.class.names
df_sero_long <- gather(df_sero, age.class)

sero_plot <- 
  ggplot() +
  geom_point(data=SeroData, aes(x=AgeClass, y=Positive), size=3) +
  geom_violin(data=df_sero_long, aes(x=age.class, y=value)) +
  scale_x_discrete(limits=age.class.names) ## avoid sorting strings

sero_plot


## Are the simulations of the sero data from Stan compatible with our data?
## We can compute confidence intervals for the fractions of non-susceptible
## individuals based on the sero data alone (using Jeffreys method).
## This will be displayed as age-specific error bars.
## Then, we can compare this with the point estimates from each simulated 
## value. These are shown as a violin plot.
df_sero_frac <- as.data.frame(sweep(as.matrix(df_sero), 2, SeroData$SampleSize, "/"))
df_sero_frac_long <- gather(df_sero_frac, age.class)

SeroData$Frac <- SeroData$Positive / SeroData$SampleSize
SeroData$Lower <-qbeta(0.025, SeroData$Positive + 0.5, SeroData$SampleSize - SeroData$Positive + 0.5)
SeroData$Upper <-qbeta(0.975, SeroData$Positive + 0.5, SeroData$SampleSize - SeroData$Positive + 0.5)


sero_plot2 <- 
  ggplot(SeroData) +
  geom_point(aes(x=AgeClass, y=Frac)) +
  geom_errorbar(aes(x=AgeClass, ymin=Lower, ymax=Upper), width=0.3) +
  geom_violin(data=df_sero_frac_long, aes(x=age.class, y=value), color='red', alpha=0.5) +
  scale_x_discrete(limits=age.class.names) ## avoid sorting strings

sero_plot2

## look at zeta for household and community transmission
## Only applicable to the split model!

ggplot(as.data.frame(params)) + geom_density(aes(x=zeta_house), fill='blue') +
  geom_density(aes(x=zeta_comty), fill='red') + xlim(0,2)


## plot the estimated susc ORs

df_susc <- as.data.frame(params$beta)
names(df_susc) <- age.class.names
df_susc_long <- gather(df_susc, age.class)

susc_plot <- 
  ggplot() +
  geom_violin(data=df_susc_long, aes(x=age.class, y=value)) +
  scale_x_discrete(limits=age.class.names) ## avoid sorting strings

susc_plot




## look at hospitalization rates 

df_hosp <- as.data.frame(log10(params$nu))
names(df_hosp) <- age.class.names
df_hosp_long <- gather(df_hosp, age.class)

hosp_rate_plot <- ggplot() +
  geom_violin(data=df_hosp_long, aes(x=age.class, y=value)) +
  scale_x_discrete(limits=age.class.names) ## avoid sorting strings

hosp_rate_plot



## compute posterior density of R0

beta_est <- params$beta
gamma_est <- params$gamma
alpha_est <- params$alpha
epsilon_est <- params$epsilon

R0_est <- c()

for ( i in 1:length(epsilon_est) ) {
  R0 <- calc_R0(epsilon_est[i], beta_est[i,], gamma_est[i], alpha_est[i], 
                contactmatrix_unperturbed)
  R0_est <- c(R0_est, R0)
}

hist(R0_est)


## plot prior distributions vs marginal posteriors

## 1/gamma
df.prior.posterior <- gather(data.frame(posterior=1/params$gamma, 
                                        prior=1/params$prior_sample_gamma))
pp.plot.gamma <- ggplot(data=df.prior.posterior) +
  geom_density(aes(x=value, fill=key), alpha=0.3)

## 1/alpha
df.prior.posterior <- gather(data.frame(posterior=1/params$alpha, 
                                        prior=1/params$prior_sample_alpha))
pp.plot.alpha <- ggplot(data=df.prior.posterior) +
  geom_density(aes(x=value, fill=key), alpha=0.3)

## zeta 
df.prior.posterior <- gather(data.frame(posterior=params$zeta, 
                                        prior=params$prior_sample_zeta))
pp.plot.zeta <- ggplot(data=df.prior.posterior) +
  geom_density(aes(x=value, fill=key), alpha=0.3)


pp.plots <- plot_grid(
  pp.plot.gamma + ggtitle("1/gamma"),
  pp.plot.alpha + ggtitle("1/alpha"), 
  pp.plot.zeta + ggtitle("zeta"), 
  nrow=1,
  ncol=3
)

pp.plots

################### COMPARE MODELS ##################

## compute WAIC to compare models
LLarray <- extract_log_lik(fit, parameter_name='log_lik_vec', merge_chains=F)

## I've been using WAIC to compare models a lot
waic(LLarray)

## But according the the Elders of Stan, the PSIS-LOO method is better
r_eff <- relative_eff(exp(LLarray), cores = 4)

loo(LLarray, r_eff=r_eff)

## If you want to compare models, use the loo_compare() method 
## to get an accurate standard error estimate

## TODO: WBIC
