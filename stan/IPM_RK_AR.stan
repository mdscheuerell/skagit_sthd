
functions {
  // spawner-recruit functions
  real ricker(real log_a, real b, real S) {
    real R;
    R = S * exp(log_a - S / b);
    return(R);
  }
}

data {
  // number of years of spawner observations
  int<lower=1> n_yrs;                       
  // number of adult age classes
  int<lower=1> AA;                       
  // minimum adult age
  int<lower=1> age_min;          
  // maximum adult age
  int<lower=age_min> age_max;
  // early years to skip
  int<lower=0> age_skip;
  // number of years to forecast
  int<lower=0> n_fore;          
  // number of missing recruits to impute
  int<lower=1> n_rec_missing;
  // number of years with missing spawner abundance obs
  int<lower=0,upper=n_yrs> n_S_NA;          
  // rows with missing spawner abundance obs
  int<lower=1,upper=n_yrs> which_S_NA[max(n_S_NA,1)]; 
  // observed annual total spawner abundance (not density)
  vector<lower=0>[n_yrs] dat_esc;
  // observed annual spawner age distributions
  matrix<lower=0>[n_yrs-(age_min+age_skip),AA] dat_age;
  // total catch of adults (no harvest on jacks)
  vector[n_yrs+n_fore] dat_harv;                         
}

parameters {
  // log intrinsic productivity (units: R/S)
  real<lower=-2,upper=10> log_a;
  // asymptotic recruitment (units: S)
  // real<lower=0,upper=1> b;       
  real<lower=1> b;       
  // process error SD
  real<lower=1e-4,upper=3> sigma_proc;
  // process error z-scores
  vector[n_yrs-age_min+n_fore] proc_err_z;
  // AR(1) coef for proc errors
  real<lower=0,upper=1> phi;
  // observation error SD
  real<lower=1e-4,upper=3> sigma_obs;
  // hyper-mean age distribution by cohort
  simplex[AA] mu_p;                     
  // process error precision of age distribution by cohort
  real<lower=0.1,upper=50> prec_p;          
  // age distributions by cohort (i.e., brood year)
  simplex[AA] pp[n_yrs-age_min+n_fore];                    
  // true total spawner abundance in years 1:(age_min+age_skip)
  vector<lower=0>[age_min+age_skip] S_tot_init;       
  // unknown (missing) spawner abundance obs
  vector<lower=0>[max(n_S_NA,1)] dat_esc_NA;    
  // true recruits for early missing broods
  vector<lower=0>[n_rec_missing] R_early;           
}

transformed parameters {
  // carrying capacity
  // real<lower=0> K;                     
  // true total spawner abundance
  vector<lower=0>[n_yrs+n_fore] spawners;           
  // true total run size
  vector<lower=0>[n_yrs+n_fore] tot_run;           
  // age distribution for calendar year returns
  matrix[n_yrs-(age_min+age_skip),AA] age_props;                      
  // expected recruits by brood year
  vector<lower=0>[n_yrs-age_min+n_fore] R_tot_hat;
  // true recruits by brood year
  vector<lower=0>[n_yrs-age_min+n_fore] R_tot;           
  // process errors
  vector[n_yrs-age_min+n_fore] proc_err;
  // temp variable: true spawners by age
  vector<lower=0>[AA] run;
  // recruits by age
  vector<lower=0>[AA] Rec[n_yrs-age_min+n_fore];            
  
  // if(log_a>0) {
  //   K = log_a / b;
  // }
  // else {
  //   K = 0;
  // }
  
  {
    // init dummy counter
    int cnt;
    cnt = 0;                     
    
    for(t in 1:(n_yrs+n_fore)) {
      if(t <= age_min+age_skip) {
        // 1: early years with no projected recruits
        tot_run[t] = S_tot_init[t];
      }
      else {
        if(t <= age_max) {
          // 2: middle years with incomplete recruits
          // modeled recruits
          for(a in 1:(t-age_min)) {
            run[a] = Rec[t-age_min,a];
          }
          // imputed recruits
          for(a in (t-age_min+1):AA) {
            cnt = cnt + 1;
            run[a] = R_early[cnt];
          }
        }
        else {
          // 3: later years with complete recruits 
          for(a in 1:AA) {
            run[a] = Rec[t-age_min,a];
          }
        }
        // total spawners in years with projected recruits
        tot_run[t] = sum(run);
        // projected calendar-yr age dist
        if(t <= n_yrs) {
          for(a in 1:AA) {
            age_props[t-age_min-age_skip,a] = log(run[a]/tot_run[t]);
          } 
        }
      }  // end else for middle & late years
      // estimated calendar-yr spawners
      spawners[t] = tot_run[t] - dat_harv[t];
      // get projected recruits by age
      if(t <= (n_yrs-age_min+n_fore)) {
        // process errors are AR(1)
        if(t == 1) {
          proc_err[t] = sigma_proc * proc_err_z[t] / sqrt(1 - phi^2);
        }
        else {
          proc_err[t] = phi * proc_err[t-1] + sigma_proc * proc_err_z[t];
        }
        // total recruits
        R_tot_hat[t] = ricker(log_a, b, spawners[t]);
        R_tot[t] = R_tot_hat[t] * exp(proc_err[t]);
        // brood-yr recruits by age
        Rec[t] = R_tot[t] * pp[t];
      }
      
    }  // end loop over t
    
  }  // end group for dummy counter
  
}  // end transformed params

model {
  // observed total spawner abundance with missing obs filled in
  vector[n_yrs] dat_esc_aug;     

  // Priors
  log_a ~ normal(0,5);
  // b ~ normal(0,2);
  // phi ~ normal(0,5);
  S_tot_init ~ lognormal(0,10);
  R_early ~ lognormal(0,10);
  if(n_S_NA == 0) {
    dat_esc_NA ~ lognormal(0,10);
  }
  
  // Process models
  // hierarchical model of recruit age distn
  for(t in 1:(n_yrs-age_min+n_fore)) {
    pp[t] ~ dirichlet(mu_p*prec_p);
  }
  // total recruitment by brood year
  // R_tot ~ lognormal(log(R_tot_hat), sigma_proc);
  proc_err_z ~ normal(0,1);
  
  // Observation models
  dat_esc_aug = dat_esc;
  if(n_S_NA > 0) {
    for(t in 1:n_S_NA) {
      dat_esc_aug[which_S_NA[t]] = dat_esc_NA[t];
    }
  }
  // observed total spawners
  dat_esc_aug ~ lognormal(log(spawners), sigma_obs);
  // observed age distribution
  target += sum(dat_age .* age_props);
}

