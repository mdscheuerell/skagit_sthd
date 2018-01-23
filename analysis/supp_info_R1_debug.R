## ----set_options, echo = FALSE, message = FALSE--------------------------
options(width = 100)
knitr::opts_chunk$set(message = FALSE)
set.seed(123)
if(file.exists("cnt_time.txt")) {
  file.remove("cnt_time.txt")
}

## ----load_pkgs, message = FALSE, warning = FALSE-------------------------
if(!require("R2jags")) {
  install.packages("R2jags")
  library("R2jags")
}
if(!require("readr")) {
  install.packages("readr")
  library("readr")
}
if(!require("gsl")) {
  install.packages("gsl")
  library("gsl")
}
if(!require("loo")) {
  install.packages("loo")
  library("loo")
}

## ----define_Re2prec------------------------------------------------------
## better round
Re2prec <- function(x, map = "round", prec = 1) {
  ## 'map' can be round, floor, or ceiling
  ## 'prec' is nearest value (eg, 0.1 means to nearest tenth); default 1 gives normal behavior
  if(prec<=0) { stop("\"prec\" cannot be less than or equal to 0") }
  do.call(map,list(x/prec))*prec
}

## return spring transition index
get_STI <- function(x, day_max = 200) {
	return(min(which(x==min(x[1:day_max]))))
}

## colVars; from Gelman
## returns the column-wise variance of a matrix
colVars <- function(a) {
	n <- dim(a)[[1]]
	c <- dim(a)[[2]]
	mm <- matrix(.colMeans(a, n, c), n, c, byrow = TRUE)
	return(.colMeans(((a - mm) ^ 2), n, c) * n / (n - 1))
}

## waic; from Gelman
## computes WAIC based on pointwise log-like
waic <- function(log_lik) {
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw /
    matrix(colMeans(loo_weights_raw),nrow = S,ncol = n,byrow = TRUE)
  loo_weights_regularized <- pmin(loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized) /
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic = total["waic"],
              elpd_waic = total["elpd_waic"],
              p_waic = total["p_waic"],
              elpd_loo = total["elpd_loo"],
              p_loo = total["p_loo"],
              pointwise = pointwise,
              total = total,
              se = se))
}

## ----get_user_inputs-----------------------------------------------------
## 1. file with escapement data
## [n_yrs x 2] matrix of obs counts; 1st col is calendar yr
fn_esc <- "skagit_sthd_esc.csv"

## 2. file with age comp data
## [n_yrs x (1+A)]; 1st col is calendar yr
fn_age <- "skagit_sthd_age.csv"
## min & max ages
age_min <- 3
age_max <- 8
## years, if any, of age-comp to skip; see below
age_skip <- 0

## 3. file with harvest data
## [n_yrs x 2] matrix of obs catch; 1st col is calendar yr
fn_harv <- "skagit_sthd_catch.csv"

## 4. file with hatchery release data
## [n_yrs x (1+MM)]; 1st col is calendar yr
fn_hrel <- "skagit_sthd_hrel.csv"

## time lags (years) for covariates
flow_lag <- 1
marine_lag <- 2
hrel_lag <- 2

## number of years of forecasts
n_fore <- 0

## upper threshold for Gelman & Rubin's (1992) potential scale reduction factor (Rhat).
Rhat_thresh <- 1.1

## URL for example data files
## set to NULL if using a local folder/directory
ex_url <- "https://raw.githubusercontent.com/mdscheuerell/skagit_sthd/master/data/"

## ----get_escapement_data-------------------------------------------------
## escapement
dat_esc <- read_csv(paste0(ex_url,fn_esc))
## years of data
dat_yrs <- dat_esc$year
## number of years of data
n_yrs <- length(dat_yrs)
## get first & last years
yr_frst <- min(dat_yrs)
yr_last <- max(dat_yrs)
## log of escapement
ln_dat_esc <- c(log(dat_esc$escapement),rep(NA,n_fore))

## ----get_age_data--------------------------------------------------------
## age comp data
dat_age <- read_csv(paste0(ex_url,fn_age))
## drop year col & first age_min+age_skip rows
dat_age <- dat_age[-(1:(age_min+age_skip)),-1]
## num of age classes
A <- age_max-age_min+1
## add row(s) of NA's for forecast years
if(n_fore > 0) {
  dat_age <- rbind(dat_age,
                   matrix(0, n_fore, A,
                          dimnames = list(n_yrs+seq(n_fore),colnames(dat_age))))
}
## total num of age obs by cal yr
dat_age[,"sum"] <- apply(dat_age, 1, sum)
## row indices for any years with no obs age comp
idx_NA_yrs <- which(dat_age$sum<A, TRUE)
## replace 0's in yrs w/o any obs with NA's
dat_age[idx_NA_yrs,(1:A)] <- NA
## change total in yrs w/o any obs from 0 to A to help dmulti()
dat_age[idx_NA_yrs,"sum"] <- A
## convert class
dat_age <- as.matrix(dat_age)

## ----get_harvest---------------------------------------------------------
## harvest
dat_harv <- read_csv(paste0(ex_url, fn_harv))
## drop year col & first age_max rows
dat_harv <- c(dat_harv$catch, rep(0,n_fore))

## ----get_flow_url--------------------------------------------------------
## flow site
flow_site <- 12178100
## get URL for flow data from USGS
flow_url <- paste0("https://waterdata.usgs.gov/nwis/dv",
                   "?cb_00060=on",
                   "&format=rdb",
                   "&site_no=",flow_site,
                   "&begin_date=",yr_frst,"-01-01",
                   "&end_date=",yr_last,"-12-31")

## ----get_flow_metadata---------------------------------------------------
## raw flow data from USGS
flow_raw <- read_lines(flow_url)
## lines with metadata
hdr_flow <- which(lapply(flow_raw, grep, pattern = "\\#")==1, arr.ind = TRUE)
## print flow metadata
print(flow_raw[hdr_flow], quote = FALSE)

## ----get_flows-----------------------------------------------------------
## flow data for years of interest
dat_flow <-  read_tsv(flow_url, col_names = FALSE, col_types = "ciDdc", skip = max(hdr_flow)+2)
colnames(dat_flow) <- unlist(strsplit(tolower(flow_raw[max(hdr_flow)+1]), split = "\\s+"))
head(dat_flow)

## ----trim_dat_flow-------------------------------------------------------
## keep only relevant columns
dat_flow <- dat_flow[c("datetime", grep("[0-9]$", colnames(dat_flow), value = TRUE))]
## nicer column names
colnames(dat_flow) <- c("date","flow")
## convert cubic feet to cubic meters
dat_flow$flow <- dat_flow$flow / 35.3147
## flow by year & month
dat_flow$year <- as.integer(format(dat_flow$date,"%Y"))
dat_flow$month <- as.integer(format(dat_flow$date,"%m"))
dat_flow <- dat_flow[,c("year","month","flow")]

## ----wtr_flow------------------------------------------------------------
## autumn flows in year t
flow_aut <- subset(dat_flow, (month>=10 & month<=12)
                   & year >= yr_frst & year <= yr_last-age_min+n_fore)
## spring flows in year t+1
flow_spr <- subset(dat_flow, (month>=1 & month<=3)
                   & year >= yr_frst+flow_lag & year <= yr_last-age_min+n_fore+flow_lag)
## change spr year index to match aut
flow_spr[,"year"] <- flow_spr[,"year"] - flow_lag
## combined flows indexed to brood year and calculate max flow over time period
dat_flow_wtr <- aggregate(flow ~ year, data = rbind(flow_aut,flow_spr), max)
dat_flow_wtr[,"flow"] <- round(dat_flow_wtr[,"flow"], 1) 
## change year index to brood year
dat_flow_wtr[,"year"] <- dat_flow_wtr[,"year"] 
## for plotting purpose later
colnames(dat_flow_wtr)[2] <- "Winter"

## ----sum_flow------------------------------------------------------------
## summer flows in year t
flow_sum <- subset(dat_flow, (month>=6 & month<=9)
                   & year >= yr_frst+flow_lag & year <= yr_last-age_min+n_fore+flow_lag)
## change year index to brood year
flow_sum[,"year"] <- flow_sum[,"year"] - flow_lag
## combined flows indexed to brood year and calculate max flow over time period
dat_flow_sum <- aggregate(flow ~ year, data = flow_sum, min)
dat_flow_sum <- round(dat_flow_sum, 2)
## for plotting purpose later
colnames(dat_flow_sum)[2] <- "Summer"

## ----get_NPGO_metadata---------------------------------------------------
## URL for NPGO data
url_NPGO <- "http://www.o3d.org/npgo/npgo.php"
## raw NPGO data 
NPGO_raw <- read_lines(url_NPGO)
## line with data headers
hdr_NPGO <- which(lapply(NPGO_raw,grep,pattern="YEAR")==1, arr.ind = TRUE)
## print PDO metadata
print(NPGO_raw[seq(hdr_NPGO)],quote = FALSE)

## ----get_NPGO------------------------------------------------------------
## NPGO data for years of interest
dat_NPGO <- read_table(url_NPGO, col_names = FALSE,
                       skip = hdr_NPGO + (yr_frst-1950)*12,
                       n_max = (n_yrs-1)*12)
colnames(dat_NPGO) <- c("year","month","NPGO")
## select only years of interest indexed by brood year 
dat_NPGO <- dat_NPGO[dat_NPGO$year >= yr_frst+marine_lag &
                     dat_NPGO$year <= yr_last-age_min+n_fore+marine_lag,]
dat_NPGO <- aggregate(dat_NPGO$NPGO, by = list(year = dat_NPGO$year), FUN = mean)
dat_NPGO <- data.frame(year = seq(yr_frst,yr_last-age_min+n_fore), NPGO = dat_NPGO[,2])
dat_NPGO[,"NPGO"] <- round(dat_NPGO[,"NPGO"], 2)
dat_NPGO

## ----get_CUI_metadata----------------------------------------------------
## URL for CUI data
url_CUI <- "https://www.pfeg.noaa.gov/products/PFELData/upwell/daily/p06dayac.all"
## raw CUI data from PFEL
CUI_raw <- read_lines(url_CUI)
## line with data headers
hdr_CUI <- which(lapply(CUI_raw,grep,pattern="YYYYMMDD")==1, arr.ind = TRUE)
## print CUI metadata
print(CUI_raw[seq(hdr_CUI-1)],quote = FALSE)

## ----get_CUI-------------------------------------------------------------
## get daily CUI data
dat_CUI <- read_table(url_CUI, col_names = TRUE, skip = hdr_CUI-1)
## extract year from date
dat_CUI$yr <- gsub("[0-9]{4}$","",dat_CUI$YYYYMMDD)
## select only years of interest
cui <- dat_CUI[dat_CUI$yr >= yr_frst+marine_lag & dat_CUI$yr <= yr_last-age_min+n_fore+marine_lag,]
## calculate cumulative upwelling by year
cum_CUI <- tapply(cui$Index, cui$yr, cumsum)
## calc STI for each year
dat_STI <- data.frame(year = seq(yr_frst,yr_last-age_min+n_fore),STI = sapply(cum_CUI,get_STI))

## ----get_hatchery_releases-----------------------------------------------
dat_hrel <- read_csv(paste0(ex_url,fn_hrel)) 
dat_hrel <- subset(dat_hrel, dat_hrel$year <= max(dat_STI$year))
dat_hrel[,2] <- dat_hrel[,2]/1000
dat_hrel

## ----combine_covars------------------------------------------------------
## covariate(s)
dat_cvrs <- Reduce(function(...) merge(..., all = TRUE),
                   list(dat_flow_sum,dat_NPGO,dat_flow_wtr,dat_STI,dat_hrel))
## drop year col
dat_cvrs <- dat_cvrs[,-1] 
dat_cvrs <- matrix(rnorm(37*5), 37, 5) 
## transform the covariates to z-scores
scl_cvrs <- as.matrix(scale(dat_cvrs)) 
## total number of covariates
n_cov <- dim(scl_cvrs)[2] 

## ----JAGS_RK_AR----------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(0.1,20);
  mu_Rkr_a <- log(alpha);
  E_Rkr_a <- mu_Rkr_a + sigma_r/(2 - 2*phi^2);
  
  ## strength of dens depend
  beta ~ dunif(0,0.01);
  
  ## AR(1) coef for proc errors
  phi ~ dunif(-0.999,0.999);
  
  ## process variance for recruits model
  sd_r ~ dunif(0.001,10);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,10);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(1,100);
  for(t in 1:(n_yrs-age_min+n_fore)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
#  ttl_run_mu ~ dunif(1,5);
#  ttl_run_tau ~ dunif(1,20);
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
#    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
#    tot_Run[i] <- exp(ln_tot_Run[i]);
    tot_Run_early[i] ~ dnorm(0,1e-5) I(0, 5e4);
    tot_Run[i] <- tot_Run_early[i];
  }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## 1st brood yr requires different innovation
  ## predicted recruits in BY t
  ln_Rkr_a[1] <- mu_Rkr_a;
  E_ln_Rec[1] <- ln_Rkr_a[1] + ln_Sp[1] - beta*Sp[1] + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
#    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
    Rec[1,a] <- tot_Rec[1] * pi_vec[1,a];
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min+n_fore)) {
    ## predicted recruits in BY t
    ln_Rkr_a[t] <- mu_Rkr_a; 
    E_ln_Rec[t] <- ln_Rkr_a[t] + ln_Sp[t] - beta*Sp[t] + phi*res_ln_Rec[t-1];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
#      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
      Rec[t,a] <- tot_Rec[t] * pi_vec[t,a];
    }
  } ## end t loop over year
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip+n_fore)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs+n_fore)) {
    ## obs model for spawners
    Sp[t] <- max(10,tot_Run[t] - dat_harv[t]);
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_RK_AR.txt")

## ----JAGS_RK_cov_AR------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(0.1,15);
  mu_Rkr_a <- log(alpha);
  E_Rkr_a <- mu_Rkr_a + sigma_r/(2 - 2*phi^2);
  
  ## strength of dens depend
  beta ~ dunif(0,0.01);
  
  ## covariate effects
  # for(i in 1:n_cov) { gamma[i] ~ dnorm(0,0.01) }
  gamma ~ dnorm(0,0.01);

  ## AR(1) coef for proc errors
  phi ~ dunif(-0.999,0.999);
  
  ## process variance for recruits model
  sd_r ~ dunif(0.001,10);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,10);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(1,100);
  for(t in 1:(n_yrs-age_min+n_fore)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
#  ttl_run_mu ~ dunif(1,5);
#  ttl_run_tau ~ dunif(1,20);
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
#    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
#    tot_Run[i] <- exp(ln_tot_Run[i]);
    tot_Run_early[i] ~ dnorm(0,1e-5) I(0, 5e4);
    tot_Run[i] <- tot_Run_early[i];
  }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## 1st brood yr requires different innovation
  ## predicted recruits in BY t
  # covar[1] <- inprod(gamma,scl_cvrs[1,]);
  covar[1] <- gamma*mod_cvrs[1];
  ln_Rkr_a[1] <- mu_Rkr_a + covar[1]; 
  E_ln_Rec[1] <- ln_Rkr_a[1] + ln_Sp[1] - beta*Sp[1] + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
#    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
    Rec[1,a] <- tot_Rec[1] * pi_vec[1,a];
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min+n_fore)) {
    ## predicted recruits in BY t
    # covar[t] <- inprod(gamma, scl_cvrs[t,]);
    covar[t] <- gamma*mod_cvrs[t];
    ln_Rkr_a[t] <- mu_Rkr_a + covar[t]; 
    E_ln_Rec[t] <- ln_Rkr_a[t] + ln_Sp[t] - beta*Sp[t] + phi*res_ln_Rec[t-1];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
#      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
      Rec[t,a] <- tot_Rec[t] * pi_vec[t,a];
    }
  } ## end t loop over year
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip+n_fore)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs+n_fore)) {
    ## obs model for spawners
    Sp[t] <- max(1,tot_Run[t] - dat_harv[t]);
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_RK_cov_AR.txt")

## ----JAGS_BH_AR----------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(1,100);
  mu_BH_a <- log(alpha);
  E_BH_a <- mu_BH_a + sigma_r/(2 - 2*phi^2);
  
  ## strength of dens depend
  beta ~ dunif(0,0.01);
  
  ## AR(1) coef for proc errors
  phi ~ dunif(-0.999,0.999);
  
  ## process variance for recruits model
  sd_r ~ dunif(0.001,10);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,10);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
#  ttl_run_mu ~ dunif(1,5);
#  ttl_run_tau ~ dunif(1,20);
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
#    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
#    tot_Run[i] <- exp(ln_tot_Run[i]);
    tot_Run_early[i] ~ dnorm(0,1e-5) I(0, 5e4);
    tot_Run[i] <- tot_Run_early[i];
  }
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(1,100);
  for(t in 1:(n_yrs-age_min+n_fore)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## 1st brood yr requires different innovation
  ## predicted recruits in BY t
  ln_BH_a[1] <- mu_BH_a; 
  E_ln_Rec[1] <- ln_BH_a[1] + ln_Sp[1] - log(1 + beta*Sp[1]) + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
#    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
    Rec[1,a] <- tot_Rec[1] * pi_vec[1,a];
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min+n_fore)) {
    ## predicted recruits in BY t
    ln_BH_a[t] <- mu_BH_a; 
    E_ln_Rec[t] <- ln_BH_a[t] + ln_Sp[t] - log(1 + beta*Sp[t]) + phi*res_ln_Rec[t-1];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
#      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
      Rec[t,a] <- tot_Rec[t] * pi_vec[t,a];
    }
  } ## end t loop over year
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip+n_fore)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs+n_fore)) {
    ## obs model for spawners
    Sp[t] <- max(10,tot_Run[t] - dat_harv[t]);
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_BH_AR.txt")

## ----JAGS_BH_cov_AR------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(1,100);
  mu_BH_a <- log(alpha);
  E_BH_a <- mu_BH_a + sigma_r/2;
  
  ## strength of dens depend
  beta ~ dunif(0,0.01);
  
  ## covariate effects
  # for(i in 1:n_cov) { gamma[i] ~ dnorm(0,0.01) }
  gamma ~ dnorm(0,0.01);

  ## AR(1) coef for proc errors
  phi ~ dunif(-0.999,0.999);
  
  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## process variance for recruits model
  sd_r ~ dunif(0.001,10);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,10);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
#  ttl_run_mu ~ dunif(1,5);
#  ttl_run_tau ~ dunif(1,20);
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
#    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
#    tot_Run[i] <- exp(ln_tot_Run[i]);
    tot_Run_early[i] ~ dnorm(0,1e-5) I(0, 5e4);
    tot_Run[i] <- tot_Run_early[i];
  }
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(1,100);
  for(t in 1:(n_yrs-age_min+n_fore)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## predicted recruits in BY t
  # covar[1] <- inprod(gamma,scl_cvrs[1,]);
  covar[1] <- gamma*mod_cvrs[1];
  ln_BH_a[1] <- mu_BH_a + covar[1]; 
  E_ln_Rec[1] <- ln_BH_a[1] + ln_Sp[1] - log(1 + beta*Sp[1]) + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
#    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
    Rec[1,a] <- tot_Rec[1] * pi_vec[1,a];
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min+n_fore)) {
    ## predicted recruits in BY t
    # covar[t] <- inprod(gamma, scl_cvrs[t,]);
    covar[t] <- gamma*mod_cvrs[t];
    ln_BH_a[t] <- mu_BH_a + covar[t]; 
    E_ln_Rec[t] <- ln_BH_a[t] + ln_Sp[t] - log(1 + beta*Sp[t]) + phi*res_ln_Rec[t-1];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
#      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
      Rec[t,a] <- tot_Rec[t] * pi_vec[t,a];
    }
  } ## end t loop over year
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip+n_fore)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs+n_fore)) {
    ## obs model for spawners
    Sp[t] <- max(10,tot_Run[t] - dat_harv[t]);
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t], ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_BH_cov_AR.txt")

## ----jags_setup----------------------------------------------------------
## 1. data to pass to JAGS
dat_jags <- c("ln_dat_esc", "dat_age", "dat_harv",
              "A", "age_min", "age_max", "age_skip",
              "n_yrs", "n_fore") 

## 2. model params/states for JAGS to return
par_jags <- c("lp_age","lp_esc")

## 3. MCMC control params
## MCMC parameters
mcmc_chains <- 4
mcmc_length <- 3e4
mcmc_burn <- 5e3
mcmc_thin <- 20
## total number of MCMC samples
mcmc_samp <- (mcmc_length-mcmc_burn)*mcmc_chains/mcmc_thin

## ----JAGS_IO_1, message = FALSE, warning = FALSE, cache = TRUE-----------
## empty list for fits
mod_fits <- vector("list", 2*(1+n_cov))

library(rjags)

init_vals_cov_AR <- function() {
	list(alpha = 10,
	     beta = 1/exp(mean(ln_dat_esc, na.rm = TRUE)),
	     gamma = 0,
	     pi_tau = 1,
	     pi_eta = rep(1,A),
	     pi_vec = matrix(c(0.015,0.35,0.465,0.15,0.015,0.005), n_yrs-age_min+n_fore, A,
	                     byrow = TRUE),
	     Rec_mu = log(1000),
	     Rec_sig = 0.1,
	     tot_ln_Rec = rep(log(1000), n_yrs - age_min + n_fore),
	     tot_Run_early = rep(8000, age_min+age_skip),
#	     ttl_run_mu = 1,
#	     ttl_run_tau  = 1,
	     innov_1 = 0,
	     phi = 0.5)
}

dat_jags <- list(dat_age=dat_age,
                 ln_dat_esc=ln_dat_esc,
                 dat_harv=dat_harv,
                 mod_cvrs=scl_cvrs[,1],
                 A=A,
                 age_min=age_min,
                 age_max=age_max,
                 age_skip=age_skip,
                 n_yrs=n_yrs,
                 n_fore=n_fore) 

par_jags <- c("alpha","E_Rkr_a","mu_Rkr_a",
              "beta",
              "gamma",
              "lp_age","lp_esc",
              "Sp","Rec","tot_ln_Rec","ln_RS",
              "pi_vec",
              "sigma_r","sigma_s","res_ln_Rec")

jm <- jags.model("IPM_RK_cov_AR.txt", data=dat_jags, 
                 n.chains=4, n.adapt=5000, inits = init_vals_cov_AR)

ss <- coda.samples(jm, par_jags, n.iter=25000, thin=20)

## ----db2-----------------------------------------------------------------
cn <- colnames(ss[[1]])
df <- data.frame(par=cn,
                 q0=rep(NA, length(cn)),
                 q25=rep(NA, length(cn)),
                 q50=rep(NA, length(cn)),
                 q75=rep(NA, length(cn)),
                 q100=rep(NA, length(cn)))
for(i in 1:length(cn)) {
  df[i, -1] <- apply(sapply(ss[,i], quantile), 1, mean)
}
df

