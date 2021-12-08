# Bang et al (2021) Neurocomputational mechanisms of confidence in self and other
# 
# Wrapper function for fitting computational models in Stan
# 
# Requires packages specified below and installation of Stan (https://mc-stan.org/users/interfaces/rstan.html)
# 
# Dan Bang danbang.db@gmail.com 2021

## FRESH MEMORY
rm(list=ls())

## SPECIFIFY LOCAL BASE DIRECTORY
base_dir <- "/Users/dbang/Dropbox/Ego/Matlab/ucl/self_other/Repository" # laptop
# base_dir <- "C:/Users/dbang/Dropbox/Ego/Matlab/ucl/self_other/Repository" # desktop

## LOAD PACKAGES
require(R.matlab) # for reading (readMAT) and writing (writeMAT) MAT-files (https://cran.r-project.org/web/packages/R.matlab/R.matlab.pdf)
require(matlab) # emulate matlab commands (https://cran.r-project.org/web/packages/matlab/matlab.pdf)
require(shinystan) # GUI for MCMC diagnostics and plots (https://cran.r-project.org/web/packages/shinystan/shinystan.pdf)
require(rstan) # [R] interface for STAN (https://cran.r-project.org/web/packages/rstan/rstan.pdf)
require(parallel) # for parallelisation (https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)

## ADD CUSTOM FUNCTIONS (https://github.com/dollbb/estRLParam/tree/master/Rcode)
source(paste(base_dir,filesep,'Modelling',filesep,'Utils',filesep,'diagnoseMCMC.R',sep=""))
source(paste(base_dir,filesep,'Modelling',filesep,'Utils',filesep,'waic1.R',sep=""))
source(paste(base_dir,filesep,'Modelling',filesep,'Utils',filesep,'waic2.R',sep=""))
    
## DIRECTORIES
dir_fits <- paste(base_dir,filesep,'Modelling',filesep,'Stan_fits',sep="") # Stan fits go here
dir_data <- paste(base_dir,filesep,'Data',filesep,'Behaviour',filesep,'Scan',filesep,'Task',sep="") # behavioural data
dir_stan <- paste(base_dir,filesep,'Modelling',sep="") # working directory
dir_modz <- paste(base_dir,filesep,'Modelling',filesep,'Stan_models',sep="") # Stan models
dir_file <- paste(base_dir,filesep,'Modelling',filesep,'Files',sep="") # relevant files

## PREPARE RSTAN 
rstan_options(auto_write = TRUE) # save bare version of a compiled STAN program to harddrive to avoid unneccessary recompilations
options(mc.cores = parallel::detectCores()) # execute in parallel if applicable

## SPECIFY MODELS FOR FITTING (ABBREVIATIONS SAME AS IN PAPER)
models <- c('S_B1T1P0','S_B1T2P0','S_B2T1P0','S_B2T2P0',
            'S_B1T1P1','S_B1T2P1','S_B2T1P1','S_B2T2P1',
            'Q_B1T1P0','Q_B1T2P0','Q_B2T1P0','Q_B2T2P0',
            'Q_B1T1P1','Q_B1T2P1','Q_B2T1P1','Q_B2T2P1',
            'T_B1T1P0','T_B1T2P0','T_B2T1P0','T_B2T2P0',
            'T_B1T1P1','T_B1T2P1','T_B2T1P1','T_B2T2P1')

## SUBJECTS
subjects <- c('22','23','24','25','26','27','29','30','31','32','34','35','36','37','38','39','40','41','42','44','45')

## SAMPLING OPTIONS
Nchains <- 1 # number of chains
Niter <- 5000 # number of samples
Nreruns <- 4 # number of re-runs
warmup <- Niter/10 # burn in
outsamp <- 500 # samples from posterior
seedz <- c('1','2','3','4') # seed names

## PARAMETER PRIORS FOR STAN (NOT ALL PARAMETERS USED FOR ALL MODELS)
# Note different notation from paper (y: sensory evidence [paper, x]; d: direction [paper, k]; k: other's noise [paper, sigma_j])
# bounds on sensory evidence to avoid STAN blowing up
y_lb <- -1
y_ub <- 1
# group priors and bounds on sensory noise
mu_sigma_prior_mu <- .2
mu_sigma_prior_sd <- .2
mu_sigma_lb <- .05
mu_sigma_ub <- 2
sd_sigma_lb <- 0
sd_sigma_ub <- .5
# group priors and bounds on softmax bias
mu_beta0_prior_mu <- 0
mu_beta0_prior_sd <- 10
mu_beta0_lb <- -10
mu_beta0_ub <- 10
sd_beta0_lb <- .01
sd_beta0_ub <- 10
# group priors and bounds on softmax temperature
mu_beta1_prior_mu <- .2
mu_beta1_prior_sd <- 2
mu_beta1_lb <- .01
mu_beta1_ub <- 2
sd_beta1_lb <- .01
sd_beta1_ub <- 2
# group priors and bounds on learning rate
mu_alpha_prior_mu <- .1
mu_alpha_prior_sd <- .1
mu_alpha_lb <- .0001
mu_alpha_ub <- 1
sd_alpha_lb <- 0
sd_alpha_ub <- 1
# group priors and bounds on state spacing
mu_factor_prior_mu <- 2
mu_factor_prior_sd <- 4
mu_factor_lb <- .01
mu_factor_ub <- 10
sd_factor_lb <- .01
sd_factor_ub <- 4
# when softmax bias is fitted separately for self- and other-trials
mu_beta0_s_prior_mu <- mu_beta0_prior_mu
mu_beta0_s_prior_sd <- mu_beta0_prior_sd
mu_beta0_s_lb <- mu_beta0_lb
mu_beta0_s_ub <- mu_beta0_ub
sd_beta0_s_lb <- sd_beta0_lb
sd_beta0_s_ub <- sd_beta0_ub
mu_beta0_o_prior_mu <- mu_beta0_prior_mu
mu_beta0_o_prior_sd <- mu_beta0_prior_sd
mu_beta0_o_lb <- mu_beta0_lb
mu_beta0_o_ub <- mu_beta0_ub
sd_beta0_o_lb <- sd_beta0_lb
sd_beta0_o_ub <- sd_beta0_ub
mu_beta1_s_prior_mu <- mu_beta1_prior_mu
mu_beta1_s_prior_sd <- mu_beta1_prior_sd
mu_beta1_s_lb <- mu_beta1_lb
mu_beta1_s_ub <- mu_beta1_ub
sd_beta1_s_lb <- sd_beta1_lb
sd_beta1_s_ub <- sd_beta1_ub
mu_beta1_o_prior_mu <- mu_beta1_prior_mu
mu_beta1_o_prior_sd <- mu_beta1_prior_sd
mu_beta1_o_lb <- mu_beta1_lb
mu_beta1_o_ub <- mu_beta1_ub
sd_beta1_o_lb <- sd_beta1_lb
sd_beta1_o_ub <- sd_beta1_ub
# when temperature is fitted separately for self-trials and other-trials
mu_beta1_s_prior_mu <- mu_beta1_prior_mu
mu_beta1_s_prior_sd <- mu_beta1_prior_sd
mu_beta1_s_lb <- mu_beta1_lb
mu_beta1_s_ub <- mu_beta1_ub
sd_beta1_s_lb <- sd_beta1_lb
sd_beta1_s_ub <- sd_beta1_ub
mu_beta1_o_prior_mu <- mu_beta1_prior_mu
mu_beta1_o_prior_sd <- mu_beta1_prior_sd
mu_beta1_o_lb <- mu_beta1_lb
mu_beta1_o_ub <- mu_beta1_ub
sd_beta1_o_lb <- sd_beta1_lb
sd_beta1_o_ub <- sd_beta1_ub
mu_beta1_s_prior_mu <- mu_beta1_prior_mu
mu_beta1_s_prior_sd <- mu_beta1_prior_sd
mu_beta1_s_lb <- mu_beta1_lb
mu_beta1_s_ub <- mu_beta1_ub
sd_beta1_s_lb <- sd_beta1_lb
sd_beta1_s_ub <- sd_beta1_ub
mu_beta1_o_prior_mu <- mu_beta1_prior_mu
mu_beta1_o_prior_sd <- mu_beta1_prior_sd
mu_beta1_o_lb <- mu_beta1_lb
mu_beta1_o_ub <- mu_beta1_ub
sd_beta1_o_lb <- sd_beta1_lb
sd_beta1_o_ub <- sd_beta1_ub

## GENERAL DETAILS
Ntrials <- 40 # number of trials
Nblocks <- 3 # number of blocks
Nsubjects <- length(subjects) # number of subjects

# INITIALISE BIG MATRICES FOR HIERARHCIAL BEHAVIOURAL DATA
g_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # gamble: fill 3 x 40 x subjects matrix with zeros
c_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # choice: fill 3 x 40 x subjects matrix with zeros
d_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # direction: fill 3 x 40 x subjects matrix with zeros
rg_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # reward risky: fill 3 x 40 x subjects matrix with zeros
rs_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # reward safe: fill 3 x 40 x subjects matrix with zeros
self_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # self-trial: fill 3 x 40 x subjects matrix with zeros
acc_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # accuracy: fill 3 x 40 x subjects matrix with zeros
theta_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # theta: fill 3 x 40 x subjects matrix with zeros
dtheta_matrix <- array(0, c(Nblocks, Ntrials, length(subjects))) # theta: fill 3 x 40 x subjects matrix with zeros
subNtrials <- array(0, c(Nblocks, length(subjects))) # number of error trials: fill 3 x subjects matrix with zeros
sigma_mu_vector <- array(0, length(subjects)) # initialise vector with subject priors on sensory noise
sigma_sd_vector <- array(0, length(subjects)) # initialise vector with subject priors on sensory noise

## COLLATE BEHAVIOURAL DATA ACROSS SUBJECTS AND RUNS
for (subi in 1:Nsubjects){
  
  # load priors on sensory noise for current subject (from calibration session)
  setwd(dir_file)
  TMP <- readMat('noise.mat')
  sigma_mu_vector[subi] <- TMP$sigma.mu[subi]
  sigma_sd_vector[subi] <- TMP$sigma.sd[subi]
  
  # load data file for current subject
  setwd(dir_data)
  DATA = readMat(paste('mriData_sub_',subjects[subi],'.mat',sep=""))
  dat = DATA$DATA
  
  # initalise data arrays for current subject
  k = array(0, c(Nblocks, Ntrials))
  theta = array(0, c(Nblocks, Ntrials))
  dtheta = array(0, c(Nblocks, Ntrials))
  self = array(0, c(Nblocks, Ntrials))
  rs = array(0, c(Nblocks, Ntrials))
  rg = array(0, c(Nblocks, Ntrials))
  acc = array(0, c(Nblocks, Ntrials))
  resp = array(0, c(Nblocks, Ntrials))
  RT = array(0, c(Nblocks, Ntrials))
  gamble = array(0, c(Nblocks, Ntrials))
  gamble_RT = array(0, c(Nblocks, Ntrials))
  err = array(0, c(Nblocks, Ntrials))
  
  # loop through blocks
  for (b in 1:Nblocks) {
    
    # get names of data structure fields
    trialInd = dimnames(dat[,,b]$trials)[[1]]
    typeIind = dimnames(dat[,,b]$typeI)[[1]]
    typeIIind = dimnames(dat[,,b]$typeII)[[1]]
    
    # retrieve block data using data structure field names as indices
    ind = which(match(trialInd, "d")==1)
    k[b,] = dat[,,b]$trials[[ind]]
    ind = which(match(trialInd, "theta")==1)
    theta[b,] = dat[,,b]$trials[[ind]]
    ind = which(match(trialInd, "self")==1)
    self[b,] = dat[,,b]$trials[[ind]]
    ind = which(match(trialInd, "rs")==1)
    rs[b,] = dat[,,b]$trials[[ind]]
    ind = which(match(trialInd, "rg")==1)
    rg[b,] = dat[,,b]$trials[[ind]]
    ind = which(match(typeIind, "response")==1)
    resp[b,] = dat[,,b]$typeI[[ind]]
    ind = which(match(typeIind, "correct")==1)
    acc[b,] = dat[,,b]$typeI[[ind]]
    ind = which(match(typeIind, "response.time")==1)
    RT[b,] = dat[,,b]$typeI[[ind]]
    ind = which(match(typeIIind, "response")==1)
    gamble[b,] = dat[,,b]$typeII[[ind]]
    ind = which(match(typeIIind, "response.time")==1)
    gamble_RT[b,] = dat[,,b]$typeII[[ind]]
    dtheta[b,] = k[b,]*theta[b,]
    
    # how many "useful" trials do we have on this block?
    err[b,] = is.nan(resp[b,]) | is.nan(RT[b,]) | is.nan(gamble_RT[b,]) | is.infinite(gamble_RT[b,])
    subNtrials[b,subi] <- sum(!err[b,])
    
    # add subject and block data to group matrix
    d_matrix[b,1:subNtrials[b,subi],subi] <- k[b,!err[b,]]
    c_matrix[b,1:subNtrials[b,subi],subi] <- resp[b,!err[b,]]-1
    g_matrix[b,1:subNtrials[b,subi],subi] <- gamble[b,!err[b,]]
    rg_matrix[b,1:subNtrials[b,subi],subi] <- rg[b,!err[b,]]
    rs_matrix[b,1:subNtrials[b,subi],subi] <- rs[b,!err[b,]]
    self_matrix[b,1:subNtrials[b,subi],subi] <- self[b,!err[b,]]
    acc_matrix[b,1:subNtrials[b,subi],subi] <- acc[b,!err[b,]]
    theta_matrix[b,1:subNtrials[b,subi],subi] <- theta[b,!err[b,]]
    dtheta_matrix[b,1:subNtrials[b,subi],subi] <- dtheta[b,!err[b,]]
    
    # randomly assign decisions on other-trials (STAN requires a value, doesn't affect fit)
    for (i in 1:40)
      if (self_matrix[b,i,subi]==0) {
        c_matrix[b,i,subi]= sample(0:1,1)
      }

  }
}

## LOOP THROUGH MODELS
for (m in 1:length(models)) {
  
  # stimulus space
  setwd(dir_file)
  # stimulus space
  thetaz <- readMat('thetaz.mat')
  # linear spacing
  v_thetaz <- thetaz$v.thetaz.comp.linear[1,] # smaller state space for approximate inference
  v_dthetaz <- thetaz$v.dthetaz.comp.linear[1,] # smaller state space for approximate inference
  n_thetaz <- length(v_thetaz) # number of states
  n_dthetaz <- length(v_dthetaz) # number of states
  n_thetaz_real <- length(v_thetaz) # number of states
  n_dthetaz_real <- length(v_dthetaz) # number of states
  # biased spacing
  spaceSizeFull <- n_thetaz;
  spaceSizeFullN <- n_thetaz;
  spaceSizeFullSignedN <- spaceSizeFullN*2;
  spaceMin <- .001;
  spaceMax <- .500;
  
  # data structure
  setwd(dir_stan)
  data <- list(Nblocks, maxTrials=Ntrials, Ns=Nsubjects, subNtrials=subNtrials,
               g=g_matrix, c=c_matrix, theta=theta_matrix, dtheta=dtheta_matrix, acc=acc_matrix, rg=rg_matrix, rs=rs_matrix, self=self_matrix, 
               noise=sigma_mu_vector, noise_sd=sigma_sd_vector, 
               n_thetaz, n_dthetaz, n_thetaz_real, n_dthetaz_real, v_thetaz, v_dthetaz,
               spaceSizeFull, spaceSizeFullN, spaceSizeFullSignedN, spaceMin, spaceMax,
               y_lb, y_ub, 
               mu_sigma_lb, mu_sigma_ub, sd_sigma_lb, sd_sigma_ub, mu_sigma_prior_mu, mu_sigma_prior_sd, 
               mu_beta0_prior_mu, mu_beta0_prior_sd, mu_beta0_lb, mu_beta0_ub, sd_beta0_lb, sd_beta0_ub, 
               mu_beta0_s_prior_mu, mu_beta0_s_prior_sd, mu_beta0_s_lb, mu_beta0_s_ub, sd_beta0_s_lb, sd_beta0_s_ub,
               mu_beta0_o_prior_mu, mu_beta0_o_prior_sd, mu_beta0_o_lb, mu_beta0_o_ub, sd_beta0_o_lb, sd_beta0_o_ub,
               mu_beta1_prior_mu, mu_beta1_prior_sd, mu_beta1_lb, mu_beta1_ub, sd_beta1_lb, sd_beta1_ub,
               mu_beta1_s_prior_mu, mu_beta1_s_prior_sd, mu_beta1_s_lb, mu_beta1_s_ub, sd_beta1_s_lb, sd_beta1_s_ub,
               mu_beta1_o_prior_mu, mu_beta1_o_prior_sd, mu_beta1_o_lb, mu_beta1_o_ub, sd_beta1_o_lb, sd_beta1_o_ub,
               mu_alpha_prior_mu, mu_alpha_prior_sd, mu_alpha_lb, mu_alpha_ub, sd_alpha_lb, sd_alpha_ub,
               mu_factor_prior_mu, mu_factor_prior_sd, mu_factor_lb, mu_factor_ub, sd_factor_lb, sd_factor_ub)
  
  ## FIT N TIMES (Nreruns) USING DIFFERENT STARTING POINTS
  for (iter in 1:Nreruns) {
      
      ## INITIALIZE CURRENT MODEL: PARAMETERS FOR FITTING AND MODEL SAMPLES (PREFIX: pp)
      # seed
      seed <- sample(1:10000,1) # random seed
      # model
      if (models[m] == 'S_B1T1P1') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1','beta1','mu_factor','factor','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1','mu_factor')
        subject_param <- c('sigma','beta0','beta1','factor')
      } else if (models[m] == 'S_B1T2P1') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_factor','factor','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta0','beta0','mu_beta1_o','mu_factor')
        subject_param <- c('sigma','beta0','beta1_s','beta1_o','factor')
      } else if (models[m] == 'S_B2T1P1') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1','beta1','mu_factor','factor','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1','mu_factor')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1','factor')
      } else if (models[m] == 'S_B2T2P1') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_factor','factor','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1_s','mu_beta1_o','mu_factor')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1_s','beta1_o','factor')
      } else if (models[m] == 'Q_B1T1P1') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1','beta1','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0','beta1','alpha','factor')
      } else if (models[m] == 'Q_B1T2P1') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1_s','mu_beta1_o','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0','beta1_s','beta1_o','alpha','factor')
      } else if (models[m] == 'Q_B2T1P1') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1','beta1','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1','alpha','factor')
      } else if (models[m] == 'Q_B2T2P1') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1_s','mu_beta1_o','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1_s','beta1_o','alpha','factor')
      } else if (models[m] == 'T_B1T1P1') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1','beta1','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0','beta1','alpha','factor')
      } else if (models[m] == 'T_B1T2P1') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1_s','mu_beta1_o','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0','beta1_s','beta1_o','alpha','factor')
      } else if (models[m] == 'T_B2T1P1') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1','beta1','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1','alpha','factor')
      } else if (models[m] == 'T_B2T2P1') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','mu_factor','factor','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1_s','mu_beta1_o','mu_alpha','mu_factor')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1_s','beta1_o','alpha','factor')
      } else if (models[m] == 'S_B1T1P0') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1','beta1','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1')
        subject_param <- c('sigma','beta0','beta1')
      } else if (models[m] == 'S_B1T2P0') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta0','beta0','mu_beta1_o')
        subject_param <- c('sigma','beta0','beta1_s','beta1_o')
      } else if (models[m] == 'S_B2T1P0') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1','beta1','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1')
      } else if (models[m] == 'S_B2T2P0') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','ppx','ppCs','ppCo','ppEVr','ppPgamble','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1_s','mu_beta1_o')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1_s','beta1_o')
      } else if (models[m] == 'Q_B1T1P0') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1','beta1','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1','mu_alpha')
        subject_param <- c('sigma','beta0','beta1','alpha')
      } else if (models[m] == 'Q_B1T2P0') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1_s','mu_beta1_o','mu_alpha')
        subject_param <- c('sigma','beta0','beta1_s','beta1_o','alpha')
      } else if (models[m] == 'Q_B2T1P0') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1','beta1','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1','mu_alpha')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1','alpha')
      } else if (models[m] == 'Q_B2T2P0') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1_s','mu_beta1_o','mu_alpha')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1_s','beta1_o','alpha')
      } else if (models[m] == 'T_B1T1P0') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1','beta1','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1','mu_alpha')
        subject_param <- c('sigma','beta0','beta1','alpha')
      } else if (models[m] == 'T_B1T2P0') {
        parameters <- c('sigma','mu_beta0','beta0','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0','mu_beta1_s','mu_beta1_o','mu_alpha')
        subject_param <- c('sigma','beta0','beta1_s','beta1_o','alpha')
      } else if (models[m] == 'T_B2T1P0') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1','beta1','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1','mu_alpha')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1','alpha')
      } else if (models[m] == 'T_B2T2P0') {
        parameters <- c('sigma','mu_beta0_s','mu_beta0_o','beta0_s','beta0_o','mu_beta1_s','mu_beta1_o','beta1_s','beta1_o','mu_alpha','alpha','ppk','ppx','ppCs','ppCo','ppEVr','ppPgamble','ppd','ppu','lik','log_lik')
        group_param <- c('mu_beta0_s','mu_beta0_o','mu_beta1_s','mu_beta1_o','mu_alpha')
        subject_param <- c('sigma','beta0_s','beta0_o','beta1_s','beta1_o','alpha')
      }
      
      ## FIT MODEL
      filz <- stan_model(file = paste(dir_modz,filesep,'stan_',models[m],'.stan',sep=""))
      fitz <- vb(filz, data = data, pars = parameters, include = TRUE, 
                 seed = seed, iter=Niter, elbo_samples=200, eval_elbo=100, output_samples=outsamp, tol_rel_obj=.0001)
      
      ## GET DIAGNOSTICS FOR STAN FIT OBJECT
      samp = extract(fitz)
      
      ## EVALUATE FIT
      evlz = waic1(fitz)
      groupWAIC <- evlz$waic
      groupLOO <- -2*evlz$elpdLoo
      groupWAICStatsTotal <- evlz$total
      groupWAICStatsSe <- evlz$se
      subjectWAIC <- NULL
      subjectLOO <- NULL
      for (i in 1:length(subjects)) {
        lik <- (samp$lik[,i,])
        lik <- t(lik)
        keep <- lik[,1] < 1
        lik <- lik[keep,]
        log_lik <- log(lik)
        log_lik <- t(log_lik)
        lala <- 1;
        evlz <- waic2(log_lik)
        subjectWAIC[i] <- evlz$waic
        subjectLOO[i] <- -2*evlz$elpdLoo
      }
      
      ## EXTRACT GROUP-LEVEL PARAMETERS
      group_param_MAP <- NULL
      group_param_upper_CI <- NULL
      group_param_lower_CI <- NULL
      q95 <- c(0.025,0.975)
      for (i in 1:length(group_param)){
        sample_vector = as.numeric(samp[group_param[i]][[1]])
        group_param_MAP[group_param[i]] <- density(sample_vector)$x[which(density(sample_vector)$y ==
                                                                            max(density(sample_vector)$y))]
        group_param_upper_CI[group_param[i]] <- quantile(sample_vector,probs=q95[2])
        group_param_lower_CI[group_param[i]] <- quantile(sample_vector,probs=q95[1])
      }
      
      ## EXTRACT SUBJECT-LEVEL PARAMETERS
      subject_param <- NULL
      subject_predict <- NULL
      if (models[m] == 'S_B1T1P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean) # PP: sensory evidence - average across N samples (outsamp)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean) # PP: confidence self - average across N samples (outsamp)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean) # PP: confidence other - average across N samples (outsamp)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean) # PP: expected value - average across N samples (outsamp)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean) # PP: gamble - average across N samples (outsamp)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean) # likelihood - average across N samples (outsamp)
        ppPgamble <- (samp$ppPgamble) # PP: gamble - all N samples (outsamp)
      } else if (models[m] == 'S_B1T2P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'S_B2T1P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'S_B2T2P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B1T1P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B1T2P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B2T1P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B2T2P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B1T1P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B1T2P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B2T1P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B2T2P1') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_param$factor <- colMeans(samp$factor)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'S_B1T1P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'S_B1T2P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'S_B2T1P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'S_B2T2P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B1T1P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B1T2P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B2T1P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'Q_B2T2P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B1T1P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B1T2P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0 <- colMeans(samp$beta0)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B2T1P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1 <- colMeans(samp$beta1)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      } else if (models[m] == 'T_B2T2P0') { 
        subject_param$sigma <- colMeans(samp$sigma)
        subject_param$beta0_s <- colMeans(samp$beta0_s)
        subject_param$beta0_o <- colMeans(samp$beta0_o)
        subject_param$beta1_s <- colMeans(samp$beta1_s)
        subject_param$beta1_o <- colMeans(samp$beta1_o)
        subject_param$alpha <- colMeans(samp$alpha)
        subject_predict$ppk <- apply(samp$ppk,c(2,3),mean)
        subject_predict$ppx <- apply(samp$ppx,c(2,3),mean)
        subject_predict$ppCs <- apply(samp$ppCs,c(2,3),mean)
        subject_predict$ppCo <- apply(samp$ppCo,c(2,3),mean)
        subject_predict$ppEVr <- apply(samp$ppEVr,c(2,3),mean)
        subject_predict$ppPgamble <- apply(samp$ppPgamble,c(2,3),mean)
        subject_predict$ppd <- apply(samp$ppd,c(2,3),mean)
        subject_predict$ppu <- apply(samp$ppu,c(2,3),mean)
        subject_predict$lik <- apply(samp$lik,c(2,3),mean)
        ppPgamble <- (samp$ppPgamble)
      }
      
      lala <- 1;
      
      ## SAVE OUTPUT
      writeMat(paste('Stan_fits',filesep,'model_',models[m],'_seed_',seedz[iter],'_vb.mat',sep=""),
               group_param_names=group_param, group_param_MAP=group_param_MAP, group_param_upper_CI=group_param_upper_CI, group_param_lower_CI=group_param_lower_CI,
               subject_param=subject_param,subject_predict=subject_predict, ppPgamble = ppPgamble,
               groupWAIC=groupWAIC, subjectWAIC=subjectWAIC, groupLOO=groupLOO, subjectLOO=subjectLOO, seed=seed)
      
    }
  
}