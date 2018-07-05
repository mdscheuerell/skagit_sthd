functions {
  # spawner-recruit functions
  real ricker(real log_a, real b, real S) {
    real R;
    R <- S*exp(log_a-b*S);
    return(R);
  }
}

data {
  # number of years of spawner observations
  int<lower=1> TT;                       
  # number of adult age classes
  int<lower=1> AA;                       
  # minimum adult age
  int<lower=1> age_min;          
  # maximum adult age
  int<lower=age_min> age_max;
  # early years to skip
  int<lower=0> age_skip;
  # number of years to forecast
  int<lower=0> n_fore;          
  # number of missing recruits to impute
  int<lower=1> n_rec_missing;
  # number of years with missing spawner abundance obs
  int<lower=0,upper=TT> n_S_NA;          
  # rows with missing spawner abundance obs
  int<lower=1,upper=TT> which_S_NA[max(n_S_NA,1)]; 
  # observed annual total spawner abundance (not density)
  vector<lower=0>[TT] dat_esc;
  # observed annual spawner age distributions
  matrix<lower=0>[TT-(age_min+age_skip),AA] dat_age;
  # total catch of adults (no harvest on jacks)
  vector[TT+n_fore] dat_harv;                         
}

parameters {
  # log intrinsic productivity (units: R/S)
  real<lower=-1,upper=5> log_a;
  # asymptotic recruitment (units: S)
  real<lower=1e-6,upper=5e-3> b;       
  # process error SD
  real<lower=1e-4,upper=3> sigma_proc;
  # observation error SD
  real<lower=1e-4,upper=3> sigma_obs;    
  # hyper-mean age distribution by cohort
  simplex[AA] mu_p;                     
  # process error precision of age distribution by cohort
  real<lower=1,upper=50> prec_p;          
  # age distributions by cohort (i.e., brood year)
  simplex[AA] pp[TT-age_min+n_fore];                    
  # true total spawner abundance in years 1:(age_min+age_skip)
  vector<lower=0>[age_min+age_skip] S_tot_init;       
  # unknown (missing) spawner abundance obs
  vector<lower=0>[max(n_S_NA,1)] dat_esc_NA;    
  # true recruit abundance by brood year
  vector<lower=0>[TT-age_min+n_fore] R_tot;           
  # true recruits for early missing broods
  vector<lower=0>[n_rec_missing] R_early;           
}


transformed parameters {
  # carrying capacity
  real<lower=0> K;                     
  # true total spawner abundance
  vector<lower=0>[TT+n_fore] spawners;           
  # true total run size
  vector<lower=0>[TT+n_fore] tot_run;           
  # age distribution for calendar year returns
  matrix[TT-(age_min+age_skip),AA] age_props;                      
  # expected recruit abundance (not density) by brood year
  vector<lower=0>[TT-age_min+n_fore] R_tot_hat;
  # temp variable: true spawners by age
  vector<lower=0>[AA] run;
  # recruits by age
  vector<lower=0>[AA] Rec[TT-age_min+n_fore];            

  if(log_a>0)
  {
  	K <- log_a/b;
  }
  else
  {
  	K <- 0;
  }

  {
  # init dummy counter
  int cnt;
  cnt <- 0;                     
  
  for(t in 1:(TT+n_fore))
  {
    if(t <= age_min+age_skip)
    # 1: early years with no projected recruits
    {
      tot_run[t] <- S_tot_init[t];
    }
    else
    {
      if(t <= age_max)
      # 2: middle years with incomplete recruits
      {
	    # modeled recruits
	    for(a in 1:(t-age_min))
	    {
	      run[a] <- Rec[t-age_min,a];
	    }
	    # imputed recruits
	    for(a in (t-age_min+1):AA)
	    {
	      cnt <- cnt + 1;
	      run[a] <- R_early[cnt];
	    }
      }
      else
      # 3: later years with complete recruits 
      {
	    for(a in 1:AA)
	    {
	      run[a] <- Rec[t-age_min,a];
	    }
      }
    # total spawners in years with projected recruits
    tot_run[t] <- sum(run);
    # projected calendar-yr age dist
    if(t <= TT)
    {
      for(a in 1:AA)
      {
        age_props[t-age_min-age_skip,a] <- log(run[a]/tot_run[t]);
      } 
    }
    }  # end else for late years
    # estimated calendar-yr spawners 
      spawners[t] <- tot_run[t] - dat_harv[t];
    # get projected recruits by age
    if(t <= (TT-age_min+n_fore))
    {
	  # total recruits
	  R_tot_hat[t] <- ricker(log_a, b, spawners[t]);
	  # brood-yr recruits by age
	  Rec[t] <- R_tot[t] * pp[t];
    }

  }  # end loop over t

  }  # end group for dummy counter

}  # end transformed params

model {
  # observed total spawner abundance with missing obs filled in
  vector[TT] dat_esc_aug;     
  vector[TT] l_spawners;
  
  # Priors
  S_tot_init ~ lognormal(0,10);
  R_early ~ lognormal(0,10);
  if(n_S_NA == 0) dat_esc_NA ~ lognormal(0,10);
  
  # Process models
  # hierarchical model of recruit age distn
  for(t in 1:(TT-age_min+n_fore))
  {
    pp[t] ~ dirichlet(mu_p*prec_p);
  }
  # total recruitment by brood year
  R_tot ~ lognormal(log(R_tot_hat), sigma_proc);
  
  # Observation models
  dat_esc_aug <- dat_esc;
  if(n_S_NA > 0)
  {
    for(t in 1:n_S_NA)
    {
      dat_esc_aug[which_S_NA[t]] <- dat_esc_NA[t];
    }
  }
  for(t in 1:TT)
  {
  	l_spawners[t] <- spawners[t];
  }
  # observed total spawners
  dat_esc_aug ~ lognormal(log(l_spawners), sigma_obs);
  # observed age distribution
  increment_log_prob(sum(dat_age .* age_props));

}








