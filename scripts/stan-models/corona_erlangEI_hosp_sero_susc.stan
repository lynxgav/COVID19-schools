/* Corona evidence synthesis for Dutch Hospitalization and Seroprevalence data.
 * This Stan script is adapted from scripts developed for the publication
 * 
 * Michiel van Boven, Anne C Teirlinck, Adam Meijer, Mariëtte Hooiveld, 
 * Christiaan H van Dorp, Rachel M Reeves, Harry Campbell, Wim van der Hoek, 
 * RESCEU Investigators, Estimating Transmission Parameters for Respiratory 
 * Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination, 
 * The Journal of Infectious Diseases, Volume 222, Issue Supplement_7, 
 * 1 November 2020, Pages S688–S694, https://doi.org/10.1093/infdis/jiaa424
 *
 * This model uses some features from Stan 2.4 to make life easier.
 * Use a poisson distribution for hospitalizations
 * The infectious period and the exposed period are Erlang-distributed.
 */

functions {
	/* ODE model */
  vector corona_model(real t, vector y, real epsilon, real alpha, real gamma, 
                      real k, real x0, real zeta, vector beta_short, vector nu_short,
                      data matrix Cunp, data matrix Cdis,
                      int[] susc_classes, int[] hosp_classes, int A, int J, int F) {
    vector[(1+F+J)*A] dydt;
    real q = inv_logit(k * (t-x0));
    matrix[A,A] C = (1-q) * Cunp + q * zeta * Cdis;
    vector[A] incidence;
    vector[A] Itot = y[(F+1)*A+1:(F+2)*A];
    for ( j in 2:J ) { // add other I-stages to Itot if J > 1
      Itot += y[(F+j)*A+1:(F+j+1)*A];
    }
    incidence = (C * Itot * epsilon) .* beta_short[susc_classes] .* y[1:A];
    dydt[1:A] = -incidence;                                  // Susceptible
    // first exposed phase
    dydt[A+1:2*A] = incidence - alpha*F * y[A+1:2*A];       // Exposed
    // other E-stages (if F > 1)
    for ( j in 2:F ) {
      dydt[j*A+1:(j+1)*A] = alpha*F * (y[(j-1)*A+1:j*A] - y[j*A+1:(j+1)*A]);
    }
    // first infectious stage
    dydt[(1+F)*A+1:(2+F)*A] = alpha*F * y[F*A+1:(1+F)*A] 
        - (gamma*J + nu_short[hosp_classes]) .* y[(1+F)*A+1:(2+F)*A];
    // other infectous stages (if J > 1)
    for ( j in 2:J ) {
      dydt[(F+j)*A+1:(1+F+j)*A] = gamma*J * y[(F+j-1)*A+1:(F+j)*A] 
        - (gamma*J + nu_short[hosp_classes]) .* y[(F+j)*A+1:(1+F+j)*A];
    }
    return dydt;
  }
}

data {
  /* preliminaries  */ 
  real t0; // start integration at t0=0 (now 2020-02-21)
  int <lower = 1> A;                       
  int <lower = 1> numdayshosp;              
  
  /* for grouping of classes wrt reporting in the hospital */
  int <lower = 1> Ahosp;
  int <lower = 1, upper = Ahosp> hosp_classes[A];
  
  /* for grouping of susceptibility classes */
  int <lower = 1> Asusc;
  int <lower = 1, upper = Asusc> susc_classes[A];
  int <lower = 1, upper = Asusc> ref_class; // has OR of 1

  /* contact matrices in the pre-lockdown and lockdown */
  matrix[A, A] Cunp; // unperturbed contact matrix
  matrix[A, A] Cdis; // distancing contact matrix
  
  /* demography */
  real demography[A]; // demographic composition of the Netherlands
  
  /* observation times */
  real ts_hosp[numdayshosp];                                      
  
  int<lower=1> F; // number of compartments for Erlang-distributed exposed period
  int<lower=1> J; // number of compartments for Erlang-distributed infectious period
  
  /* fix latent period for faster sampling */
  //real <lower = 0> alpha;                                    
  
  /* sampling temperature */
  int <lower = 0, upper = 1> mode; // 0 = estimation, 1 = WBIC calculation
  
  /* hospitalisation data by admission date */
  int<lower=0> hospitalisations[numdayshosp, A];    

  /* serological data from the Netherlands */
  int<lower=0> numdayssero;
  real ts_sero[numdayssero]; // times of sero measurements
  int<lower=0> sero_num_sampled[numdayssero, A];
  int<lower=0> sero_num_pos[numdayssero, A];

  /* ODE integrator settings */
  real<lower = 0> rel_tol;
  real<lower = 0> abs_tol;
  int max_num_steps;
  
  /* Hyper parameters */
  real<lower=0> beta_short_hyp[Asusc]; // typical ORs
  real<lower=0> log_beta_sd; // sd for log-normal prior for ORs
  
  real<lower=0> a_gamma; // shape parameter for prior of gamma
  real<lower=0> b_gamma; // scale parameter for prior of gamma
  
  real<lower=0> a_alpha; // shape parameter for prior of alpha
  real<lower=0> b_alpha; // scale parameter for prior of alpha
  
  real<lower=0> m_zeta; // location parameter for prior of zeta
  real<lower=0> s_zeta; // scale parameter for prior of zeta
}

transformed data {
  /* sampling temperature */  
  real<lower = 0, upper = 1> watanabe_beta; // 0 = normal; 1 = WBIC

  int numdays_all = numdayshosp + numdayssero;
  
  real ts_all[numdays_all]; // merged times for integrating
  int sort_indices[numdays_all]; // used for sorting ts_all for ODE solver
  int time_indices_hosp[numdayshosp]; // for finding predictions for hosp observations
  int time_indices_sero[numdayssero]; // for finding predictions for sero observations
  
  /* sampling mode. 0 = normal; 1 = WBIC */
  if ( mode == 0 ) {
    watanabe_beta = 1.0;
  }
  else { // mode == 1
    /* NB: if any of the age classes for the sero data has a 0 sample
     * size, this should be exluded from the number of observations 
     */
    watanabe_beta = 1.0/log(A * numdays_all);
  }
  
  ts_all = append_array(ts_hosp, ts_sero);
  // sort_indices is used to sort the concatenated times ts_all (required for ODE integrator)
  sort_indices = sort_indices_asc(ts_all);
  // we now have to "invert" the permutation so that we can find predictions for hosp data and sero data
  for ( i in 1:numdays_all ) {
    int idx = sort_indices[i];
    if ( idx <= numdayshosp ) {
      time_indices_hosp[idx] = i;  
    } else {
      time_indices_sero[idx-numdayshosp] = i;  
    }
  }
}

parameters {
  real<lower=0> alpha;                                             // rate of becoming infectious
  real<lower=0> gamma;                                             // recovery rate
  real<lower=0, upper=1> epsilon;                                  // probability of transmission per contact    
  real<lower=0> zeta;                                              // proportionality parameter for the lockdown contact matrix
  real<lower=0> x0;                                                // x0 of logist from pre-lockdown to lockdown
  real<lower=0> k;                                                 // k of logist from pre-lockdown to lockdown
  real<lower=0> nu_scale;                                          // average hosp rate
  simplex[Ahosp] nu_simplex;                                       // relative hosp rates
  vector<lower=0, upper=1>[Asusc-1] beta_short_raw;                // age-dependent susceptibility 
  real<lower=1e-7, upper=5e-4> inoculum;                           // approx 1-10k initial infections
}

transformed parameters {
  vector<lower=0, upper=1>[(1+F+J)*A] y0;                      // initial conditions 
  vector[(1+F+J)*A] y_hat[numdays_all];                        // prevalences at time t of age a (3 classes)
  vector<lower=0>[Ahosp] p_short;                              // prob of hospitalisation
  vector<lower=0>[A] nu;                                       // hospitalisation rates              
  vector<lower=0>[Asusc] beta_short;                           // age-dependent susceptibility ORs
  vector<lower=0>[A] beta;                                     // susceptibilities for all age classes

  /* for easy reference */
  matrix[numdays_all, A] Susceptible; 
  matrix[numdays_all, A] Exposed = rep_matrix(0.0, numdays_all, A); 
  matrix[numdays_all, A] Infected = rep_matrix(0.0, numdays_all, A); 
  matrix[numdayshosp, A] log_likes_hosp;  // log-likelihood contributions of hospitalisations  
  matrix[numdayssero, A] log_likes_sero;  // log-likelihood contributions of serological data 

  vector[Ahosp] nu_short = nu_scale * nu_simplex;
  /* reduced number age-specific parameters for severe disease */
  nu = nu_short[hosp_classes]; // elongate nu_short
  
  // add a OR of 1 to the beta_short_raw vector at the reference class
  beta_short = append_row(beta_short_raw[1:ref_class-1], 
    append_row(rep_vector(1.0, 1),beta_short_raw[ref_class:]));
  
  beta = beta_short[susc_classes]; // full vector of age-dependent susc

  /* initial conditions */
  for (a in 1:A) { // TODO: extend/include long run-up period to approach stable class distribution
    y0[a] = 1.0 - inoculum;                          // Susceptible
    for ( j in 1:F ) {
      y0[j*A+a] = 0.5/F * inoculum;                  // Exposed 
    }
    for ( j in 1:J ) {
      y0[(F+j)*A+a] = 0.5/J * inoculum;              // Infected
    }
  }
  
  /* integrate ODEs and take result */
  y_hat = ode_adams_tol(corona_model, y0, t0, ts_all[sort_indices], 
      rel_tol, abs_tol, max_num_steps,
      epsilon, alpha, gamma, k, x0, zeta, beta_short, nu_short, Cunp, Cdis, 
      susc_classes, hosp_classes, A, J, F);
  
  /* extract trajectories */
  for ( i in 1:numdays_all ) {
      Susceptible[i,:] = y_hat[i][1:A]';
      for ( j in 1:F ) {
        Exposed[i,:] += y_hat[i][j*A+1:(j+1)*A]'; 
      }
      for ( j in 1:J ) {
        Infected[i,:] += y_hat[i][(F+j)*A+1:(F+j+1)*A]';
      }
  }
  
  /* probability of hospitalisation for rate-based reporting model */
  p_short = nu_short ./ (nu_short + gamma);

  /* likelihood contributions */
  for ( i in 1:numdayshosp ) {
    int idx = time_indices_hosp[i];
    for ( a in 1:A ) {
      /* prevalence-based reporting */
      real x = nu[a] * Infected[idx, a] * demography[a];
      log_likes_hosp[i, a] = poisson_lpmf(hospitalisations[i, a] | x);
    }
  }
    
  /* likelihood of sero data */
  for ( i in 1:numdayssero ) {
    int idx = time_indices_sero[i];
    for ( a in 1:A ) {
      log_likes_sero[i,a] = 
          binomial_lpmf(sero_num_pos[i,a] | sero_num_sampled[i,a], 1-Susceptible[idx,a]);
    }
  }
}

model {
  /* prior distributions */
  alpha ~ inv_gamma(a_alpha, b_alpha);  // 1/latent period; 95% prior coverage 2.2-4.4 days
  gamma ~ inv_gamma(a_gamma, b_gamma);  // 1/infectious period; 95% prior coverage 4.2-15 days

  nu_short ~ normal(0, 5);
  // manual jacobian correction: log|dnu_short/dnu_simplex| = Ahosp * nu_scale
  target += Ahosp * nu_scale;

  
  for ( a in 1:Asusc ) {
    beta_short[a] ~ lognormal(log(beta_short_hyp[a]), log_beta_sd);   // the log-OR has a normal distribution
  }

  zeta ~ normal(m_zeta, s_zeta);                                      // a priori expect no additional reduction, max 2se=0.5 at most
  x0 ~ normal(23, 7);                                                 // x0 of logistic transition from pre-lockdown to lockdown
                                                                      // day23 is currently 15 March (lockdown declared) - generalise
  k ~ exponential(1);                                                 // k of logistic transition from pre-lockdown to lockdown
  target += watanabe_beta * sum(log_likes_hosp);
  target += watanabe_beta * sum(log_likes_sero);
}

generated quantities {
  real log_lik = sum(log_likes_hosp) + sum(log_likes_sero);           // estimates WBIC when mode = 1
  real log_lik_vec[A*numdays_all] = 
      append_array(to_array_1d(log_likes_hosp), to_array_1d(log_likes_sero));
  real expected_hospitalisations[numdayshosp, A];   // for credible intervals
  int simulated_hospitalisations[numdayshosp, A];   // for prediction intervals of hospitalisations
  real expected_serodata[numdayssero, A];            // for credible intervals of serological data
  int simulated_serodata[numdayssero, A];            // for prediction intervals of serological data
  
  /* sample from priors to compare marginal posteriors with priors on parameters */
  real prior_sample_gamma = inv_gamma_rng(a_gamma, b_gamma); // samples from the prior of gamma
  real prior_sample_alpha = inv_gamma_rng(a_alpha, b_alpha); // samples from the prior of alpha
  real prior_sample_zeta = normal_rng(m_zeta, s_zeta); // samples from the prior of zeta

  /* hospitalisations */
  for ( i in 1:numdayshosp ) { 
    int idx = time_indices_hosp[i];
    for ( a in 1:A ) {
      real x = nu[a] * Infected[idx, a] * demography[a];
      expected_hospitalisations[i, a] = x;
      simulated_hospitalisations[i, a] = poisson_rng(x);
    } 
  }  
  
  /* serological data */
  for ( i in 1:numdayssero ) {
    int idx = time_indices_sero[i];
    for ( a in 1:A ) {
      int num_sam = sero_num_sampled[i, a];
      real prob_pos = (1-Susceptible[idx, a]);
      expected_serodata[i, a] = num_sam * prob_pos;
      simulated_serodata[i, a] = binomial_rng(num_sam, prob_pos);
    }
  }
} 
