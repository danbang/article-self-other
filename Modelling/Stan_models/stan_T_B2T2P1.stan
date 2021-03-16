data {
  
  // data
  int<lower=1> Ns;
  int<lower=1> Nblocks;
  int<lower=1> maxTrials;
  int<lower=1> n_thetaz;
  int<lower=1> n_dthetaz;
  real spaceSizeFull;
  int<lower=1> spaceSizeFullN;
  int<lower=1> spaceSizeFullSignedN;
  real spaceMin;
  real spaceMax;
  int<lower=1> subNtrials[Nblocks,Ns];
  int<lower=0,upper=1> g[Nblocks,maxTrials,Ns];
  int<lower=0,upper=1> c[Nblocks,maxTrials,Ns];
  int rg[Nblocks,maxTrials,Ns];
  int rs[Nblocks,maxTrials,Ns];
  int<lower=0,upper=1> self[Nblocks,maxTrials,Ns];
  real<lower=0,upper=1> theta[Nblocks,maxTrials,Ns];
  real<lower=-1,upper=1> dtheta[Nblocks,maxTrials,Ns];
  int<lower=0,upper=1> acc[Nblocks,maxTrials,Ns];
  real noise[Ns];
  real noise_sd[Ns];

  // parameters
  real mu_sigma_lb;
  real mu_sigma_ub;
  real sd_sigma_lb;
  real sd_sigma_ub;
  real mu_sigma_prior_mu;
  real mu_sigma_prior_sd;
  real y_lb;
  real y_ub;
  real mu_beta0_lb;
  real mu_beta0_ub;
  real sd_beta0_lb;
  real sd_beta0_ub;
  real mu_beta0_prior_mu;
  real mu_beta0_prior_sd;
  real mu_beta0_s_lb;
  real mu_beta0_s_ub;
  real mu_beta0_o_lb;
  real mu_beta0_o_ub;
  real mu_beta0_s_prior_mu;
  real mu_beta0_s_prior_sd;
  real mu_beta0_o_prior_mu;
  real mu_beta0_o_prior_sd;
  real mu_beta1_lb;
  real mu_beta1_ub;
  real sd_beta1_lb;
  real sd_beta1_ub;
  real mu_beta1_prior_mu;
  real mu_beta1_prior_sd;
  real mu_beta1_s_lb;
  real mu_beta1_s_ub;
  real mu_beta1_o_lb;
  real mu_beta1_o_ub;
  real mu_beta1_s_prior_mu;
  real mu_beta1_s_prior_sd;
  real mu_beta1_o_prior_mu;
  real mu_beta1_o_prior_sd;
  real mu_alpha_lb;
  real mu_alpha_ub;
  real sd_alpha_lb;
  real sd_alpha_ub;
  real mu_alpha_prior_mu;
  real mu_alpha_prior_sd;
  real mu_factor_lb;
  real mu_factor_ub;
  real sd_factor_lb;
  real sd_factor_ub;
  real mu_factor_prior_mu;
  real mu_factor_prior_sd;
  
}

parameters {
  
    // group level
    real<lower=mu_alpha_lb,upper=mu_alpha_ub> mu_alpha;
    real<lower=sd_alpha_lb,upper=sd_alpha_ub> sd_alpha;
    real<lower=sd_beta0_lb,upper=sd_beta0_ub> sd_beta0;
    real<lower=mu_beta0_s_lb,upper=mu_beta0_s_ub> mu_beta0_s;
    real<lower=mu_beta0_o_lb,upper=mu_beta0_o_ub> mu_beta0_o;
    real<lower=mu_beta1_s_lb,upper=mu_beta1_s_ub> mu_beta1_s;
    real<lower=mu_beta1_o_lb,upper=mu_beta1_o_ub> mu_beta1_o;
    real<lower=sd_beta1_lb,upper=sd_beta1_ub> sd_beta1;
    real<lower=mu_factor_lb,upper=mu_factor_ub> mu_factor;
    real<lower=sd_factor_lb,upper=sd_factor_ub> sd_factor;

    // subject level
    row_vector<lower=mu_sigma_lb,upper=mu_sigma_ub>[Ns] sigma;
    row_vector<lower=mu_alpha_lb,upper=mu_alpha_ub>[Ns] alpha;
    row_vector<lower=mu_beta0_s_lb,upper=mu_beta0_s_ub>[Ns] beta0_s;
    row_vector<lower=mu_beta0_o_lb,upper=mu_beta0_o_ub>[Ns] beta0_o;
    row_vector<lower=mu_beta1_s_lb,upper=mu_beta1_s_ub>[Ns] beta1_s;
    row_vector<lower=mu_beta1_o_lb,upper=mu_beta1_o_ub>[Ns] beta1_o;
    row_vector<lower=mu_factor_lb,upper=mu_factor_ub>[Ns] factor;
    
    // sensation
    real<lower=y_lb,upper=y_ub> y[Nblocks,maxTrials,Ns];

}

transformed parameters {

  real aa;
	real ab;

	aa = mu_alpha * pow(sd_alpha,-2);
	ab = pow(sd_alpha,-2) - aa;

}
model {
  
  // group-level priors
  mu_alpha ~ normal(mu_alpha_prior_mu,mu_alpha_prior_sd);
  mu_beta0_s ~ normal(mu_beta0_prior_mu,mu_beta0_prior_sd);
  mu_beta0_o ~ normal(mu_beta0_prior_mu,mu_beta0_prior_sd);
  mu_beta1_s ~ normal(mu_beta1_prior_mu,mu_beta1_prior_sd);
  mu_beta1_o ~ normal(mu_beta1_prior_mu,mu_beta1_prior_sd);
  mu_factor ~ normal(mu_factor_prior_mu,mu_factor_prior_sd);

  // subject-level parameters
  for (s in 1:Ns) {
    
        // subject-level priors
        sigma[s] ~ normal(noise[s],noise_sd[s]);
        alpha[s] ~ beta(aa,ab);
        beta0_s[s] ~ normal(mu_beta0_s, sd_beta0);
        beta0_o[s] ~ normal(mu_beta0_o, sd_beta0);
        beta1_s[s] ~ normal(mu_beta1_s, sd_beta1);
        beta1_o[s] ~ normal(mu_beta1_o, sd_beta1);
        factor[s] ~ normal(mu_factor, sd_factor);
        
        for (b in 1:Nblocks) {
          
        // initialise k
        real k[maxTrials+1];
        
        // initialise variables for constructing state space
        real stateBase[spaceSizeFullN];
        real R[2];
        real dR;
        real A[spaceSizeFullN];
        real spaceFull[spaceSizeFullN];
        real spaceFullSigned[spaceSizeFullSignedN];
        
        // create state space
        for (jj in 1:spaceSizeFullN) {
          stateBase[jj]= (jj/spaceSizeFull)^factor[s];
        }
        R[1]= spaceMin;
        R[2]= spaceMax;
        dR= R[2]-R[1];
        for (jj in 1:spaceSizeFullN) {
          A[jj]= stateBase[jj];
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]-min(A);
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]/max(A);
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]*dR;
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]+R[1];
        }
        for (jj in 1:spaceSizeFullN) {
          spaceFull[jj]= A[jj];
        }
        // signed
        for (jj in 1:spaceSizeFullN) {
          spaceFullSigned[jj]= -1*spaceFull[spaceSizeFullN+1-jj];
        }
        for (jj in 1:spaceSizeFullN) {
          spaceFullSigned[spaceSizeFullN+jj]= spaceFull[jj];
        }
            
        // initialise k
        k[1] = sigma[s];
          
            for (i in 1:subNtrials[b,s]) {
              
                // initialise trial variables
                real EVr;
                real PmCONDx_numerator[n_dthetaz];
                real PmCONDx_denomintor[n_dthetaz];
                real PMCONDx[n_dthetaz];
                real Pstate1CONDx_s;
                real Pstate2CONDx_s;
                real choice;
                real PcorCONDx_s;
                real PcorCONDm_o[n_dthetaz];
                real Ptemporary[n_dthetaz];
                real PcorCONDx_o;
                real pgamble;
                real p;
                real r;
                real d;
                real u;
                real term1;
                real term2;
                real term3;
                real gradient_v[n_dthetaz];
                real gradient;
                
                // draw sensory sample
                y[b,i,s] ~ normal(dtheta[b,i,s],sigma[s]);

                // calculate posterior over states
                for (kk in 1:n_dthetaz) {
                    PmCONDx_numerator[kk]= exp(normal_lpdf(y[b,i,s]|spaceFullSigned[kk],sigma[s]));
                }
                for (kk in 1:n_dthetaz) {
                    PmCONDx_denomintor[kk]= sum(PmCONDx_numerator);
                }
                for (kk in 1:n_dthetaz) {
                    PMCONDx[kk]= PmCONDx_numerator[kk]/PmCONDx_denomintor[kk];
                }
                
                // posterior probability of LEFT and RIGHT
                Pstate1CONDx_s= sum(PMCONDx[1:n_thetaz])/sum(PMCONDx);
                Pstate2CONDx_s= sum(PMCONDx[n_thetaz+1:n_thetaz*2])/sum(PMCONDx);
            
                // generate choice
                if (Pstate1CONDx_s > Pstate2CONDx_s) {
                  choice = 0;
                } else if (Pstate1CONDx_s < Pstate2CONDx_s) {
                  choice = 1;
                }
                
                // generate confidence in self
                if (choice == 0) {
                  PcorCONDx_s = Pstate1CONDx_s;
                } else if (choice == 1) {
                  PcorCONDx_s = Pstate2CONDx_s;
                }
                
                // generate confidence in other
                for (kk in 1:n_dthetaz) {
                  if (spaceFullSigned[kk]<0) {
                    PcorCONDm_o[kk]= normal_cdf(0,spaceFullSigned[kk],k[i]);
                  } else if (spaceFullSigned[kk]>0) {
                    PcorCONDm_o[kk]= 1-normal_cdf(0,spaceFullSigned[kk],k[i]);
                  }
                }
                for (kk in 1:n_dthetaz) {
                  Ptemporary[kk] = PMCONDx[kk] * PcorCONDm_o[kk];
                }
                PcorCONDx_o = sum(Ptemporary);
                
                // remove 0s and 1s to avoid infinite values
                if (PcorCONDx_s<.01) {
                    PcorCONDx_s=.01;
                } else if (PcorCONDx_s>.99) {
                    PcorCONDx_s=.99;
                } else if (PcorCONDx_o<.01) {
                    PcorCONDx_o=.01;
                } else if (PcorCONDx_o>.99) {
                    PcorCONDx_o=.99;
                }

                // expected value of risky option
                if (self[b,i,s] == 1) {
                  EVr = rg[b,i,s]*PcorCONDx_s - rg[b,i,s]*(1-PcorCONDx_s);
                } else if (self[b,i,s] == 0) {
                  EVr = rg[b,i,s]*PcorCONDx_o - rg[b,i,s]*(1-PcorCONDx_o);
                }  

                // sigmoid
                if (self[b,i,s] == 1) {
                  pgamble = inv_logit(beta0_s[s] + beta1_s[s]*(EVr - rs[b,i,s]));
                } else if (self[b,i,s] == 0) {
                  pgamble = inv_logit(beta0_o[s] + beta1_o[s]*(EVr - rs[b,i,s]));
                }
                
                // remove 0s and 1s to avoid infinite values
                if (pgamble<.01) {
                    pgamble=.01;
                } else if (pgamble>.99) {
                    pgamble=.99;
                }
                
                // generate choice
                if (self[b,i,s] == 1) {
                    c[b,i,s] ~  bernoulli_logit(100*(y[b,i,s]));
                }

                // generate gamble
                g[b,i,s] ~  bernoulli(pgamble);
                
                // update k
                if (self[b,i,s] == 1) {
                  k[i+1] = k[i];
                } else if (self[b,i,s] == 0) {
                // compute gradient
                for (kk in 1:n_dthetaz) {
                  if (spaceFullSigned[kk]<0) {
                    term1= spaceFullSigned[kk]/k[i];
                    term2= PmCONDx_numerator[kk]/PmCONDx_denomintor[kk];
                    term3= exp(normal_lpdf(0|spaceFullSigned[kk],k[i]));
                    gradient_v[kk]= term1*term2*term3;
                  } else if (spaceFullSigned[kk]>0) {
                    term1= -spaceFullSigned[kk]/k[i];
                    term2= PmCONDx_numerator[kk]/PmCONDx_denomintor[kk];
                    term3= exp(normal_lpdf(0|spaceFullSigned[kk],k[i]));
                    gradient_v[kk]= term1*term2*term3;
                  }
                }
                gradient= sum(gradient_v);
                // use gradient for updating
                p= PcorCONDx_o;
                if (acc[b,i,s]==0) {
                  r= 0;
                } else if (acc[b,i,s]==1) {
                  r= 1;
                }
                d= r-p;
                u= (alpha[s]*d)/gradient;
                k[i+1] = k[i]+u;
                
                // print("subject=",s)
                // print("block=",b)
                // print("trial=",i)
                // print("ss=",sigma[s])
                // print("sj=",k[i])
                // print("k=",dtheta[b,i,s])
                // print("x=",y[b,i,s])
                // print("d=",d)
                // print("u=",u)
                // print("gradient=",gradient)

                  // quality control on k
                  // remove low values to avoid blow-up
                  if (k[i+1]<.05) {
                      k[i+1]=.05;
                  }
                  // remove high values to avoid blow-up
                  if (k[i+1]>1) {
                      k[i+1]=1;
                  }
                  // remove nans
                  if (is_nan(k[i+1])==1) {
                      k[i+1]= k[i];
                  }
                
                }
              
            }
        }
  }
}

generated quantities {
  
  matrix[Ns,Nblocks*maxTrials] ppk;
  matrix[Ns,Nblocks*maxTrials] ppx;
  matrix[Ns,Nblocks*maxTrials] ppCs;
  matrix[Ns,Nblocks*maxTrials] ppCo;
  matrix[Ns,Nblocks*maxTrials] ppEVr;
  matrix[Ns,Nblocks*maxTrials] ppPgamble;
  matrix[Ns,Nblocks*maxTrials] ppd;
  matrix[Ns,Nblocks*maxTrials] ppu;
  matrix[Ns,Nblocks*maxTrials] ppgradient;
  matrix[Ns,Nblocks*maxTrials] lik;
  real log_lik[Ns];
  
  for (s in 1:Ns) {
    
    // trial counter
    int j;
    j = 1;
    
    // initialise
    log_lik[s] = 0;
    
    for (b in 1:Nblocks) {
      
        // initialise k
        real k[maxTrials+1];
        
        // initialise variables for constructing state space
        real stateBase[spaceSizeFullN];
        real R[2];
        real dR;
        real A[spaceSizeFullN];
        real spaceFull[spaceSizeFullN];
        real spaceFullSigned[spaceSizeFullSignedN];
        
        // create state space
        for (jj in 1:spaceSizeFullN) {
          stateBase[jj]= (jj/spaceSizeFull)^factor[s];
        }
        R[1]= spaceMin;
        R[2]= spaceMax;
        dR= R[2]-R[1];
        for (jj in 1:spaceSizeFullN) {
          A[jj]= stateBase[jj];
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]-min(A);
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]/max(A);
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]*dR;
        }
        for (jj in 1:spaceSizeFullN) {
          A[jj]= A[jj]+R[1];
        }
        for (jj in 1:spaceSizeFullN) {
          spaceFull[jj]= A[jj];
        }
        // signed
        for (jj in 1:spaceSizeFullN) {
          spaceFullSigned[jj]= -1*spaceFull[spaceSizeFullN+1-jj];
        }
        for (jj in 1:spaceSizeFullN) {
          spaceFullSigned[spaceSizeFullN+jj]= spaceFull[jj];
        }
            
        // initialise k
        k[1] = sigma[s];
    
       for (i in 1:maxTrials) {
          
          if (i <= subNtrials[b,s]) {
                
                // initialise trial variables
                real EVr;
                real PmCONDx_numerator[n_dthetaz];
                real PmCONDx_denomintor[n_dthetaz];
                real PMCONDx[n_dthetaz];
                real Pstate1CONDx_s;
                real Pstate2CONDx_s;
                real choice;
                real PcorCONDx_s;
                real PcorCONDm_o[n_dthetaz];
                real Ptemporary[n_dthetaz];
                real PcorCONDx_o;
                real pgamble;
                real p;
                real r;
                real d;
                real u;
                real x;
                real term1;
                real term2;
                real term3;
                real gradient_v[n_dthetaz];
                real gradient;
                
                
                // draw sensory sample
                x = y[b,i,s];
                ppx[s,j] = x;
                ppk[s,j] = k[i];
                
                // calculate posterior over states
                for (kk in 1:n_dthetaz) {
                    PmCONDx_numerator[kk]= exp(normal_lpdf(x|spaceFullSigned[kk],sigma[s]));
                }
                for (kk in 1:n_dthetaz) {
                    PmCONDx_denomintor[kk]= sum(PmCONDx_numerator);
                }
                for (kk in 1:n_dthetaz) {
                    PMCONDx[kk]= PmCONDx_numerator[kk]/PmCONDx_denomintor[kk];
                }
                
                // posterior probability of LEFT and RIGHT
                Pstate1CONDx_s= sum(PMCONDx[1:n_thetaz])/sum(PMCONDx);
                Pstate2CONDx_s= sum(PMCONDx[n_thetaz+1:n_thetaz*2])/sum(PMCONDx);
            
                // generate choice
                if (Pstate1CONDx_s > Pstate2CONDx_s) {
                  choice = 0;
                } else if (Pstate1CONDx_s < Pstate2CONDx_s) {
                  choice = 1;
                }
                
                // generate confidence in self
                if (choice == 0) {
                  PcorCONDx_s = Pstate1CONDx_s;
                } else if (choice == 1) {
                  PcorCONDx_s = Pstate2CONDx_s;
                }
                ppCs[s,j] = PcorCONDx_s;
                
                // generate confidence in other
                for (kk in 1:n_dthetaz) {
                  if (spaceFullSigned[kk]<0) {
                    PcorCONDm_o[kk]= normal_cdf(0,spaceFullSigned[kk],k[i]);
                  } else if (spaceFullSigned[kk]>0) {
                    PcorCONDm_o[kk]= 1-normal_cdf(0,spaceFullSigned[kk],k[i]);
                  }
                }
                for (kk in 1:n_dthetaz) {
                  Ptemporary[kk] = PMCONDx[kk] * PcorCONDm_o[kk];
                }
                PcorCONDx_o = sum(Ptemporary);
                ppCo[s,j] = PcorCONDx_o;
                
                // remove 0s and 1s to avoid infinite values
                if (PcorCONDx_s<.01) {
                    PcorCONDx_s=.01;
                } else if (PcorCONDx_s>.99) {
                    PcorCONDx_s=.99;
                } else if (PcorCONDx_o<.01) {
                    PcorCONDx_o=.01;
                } else if (PcorCONDx_o>.99) {
                    PcorCONDx_o=.99;
                }

                // expected value of risky option
                if (self[b,i,s] == 1) {
                  EVr = rg[b,i,s]*PcorCONDx_s - rg[b,i,s]*(1-PcorCONDx_s);
                } else if (self[b,i,s] == 0) {
                  EVr = rg[b,i,s]*PcorCONDx_o - rg[b,i,s]*(1-PcorCONDx_o);
                }
                ppEVr[s,j] = EVr;

                // sigmoid
                if (self[b,i,s] == 1) {
                  pgamble = inv_logit(beta0_s[s] + beta1_s[s]*(EVr - rs[b,i,s]));
                } else if (self[b,i,s] == 0) {
                  pgamble = inv_logit(beta0_o[s] + beta1_o[s]*(EVr - rs[b,i,s]));
                }
                
                // remove 0s and 1s to avoid infinite values
                if (pgamble<.01) {
                    pgamble=.01;
                } else if (pgamble>.99) {
                    pgamble=.99;
                }
                ppPgamble[s,j] = pgamble;

                // gamble likelihood
                lik[s,j] = exp(bernoulli_lpmf(g[b,i,s]| pgamble));
                log_lik[s] = log_lik[s] + (bernoulli_lpmf(g[b,i,s]| pgamble));

                // update k
                if (self[b,i,s] == 1) {
                    k[i+1] = k[i];
                    ppd[s,j] = 0;
                    ppu[s,j] = 0;
                } else if (self[b,i,s] == 0) {
                // compute gradient
                for (kk in 1:n_dthetaz) {
                  if (spaceFullSigned[kk]<0) {
                    term1= spaceFullSigned[kk]/k[i];
                    term2= PmCONDx_numerator[kk]/PmCONDx_denomintor[kk];
                    term3= exp(normal_lpdf(0|spaceFullSigned[kk],k[i]));
                    gradient_v[kk]= term1*term2*term3;
                  } else if (spaceFullSigned[kk]>0) {
                    term1= -spaceFullSigned[kk]/k[i];
                    term2= PmCONDx_numerator[kk]/PmCONDx_denomintor[kk];
                    term3= exp(normal_lpdf(0|spaceFullSigned[kk],k[i]));
                    gradient_v[kk]= term1*term2*term3;
                  }
                }
                gradient= sum(gradient_v);
                // use gradient for updating
                p= PcorCONDx_o;
                if (acc[b,i,s]==0) {
                  r= 0;
                } else if (acc[b,i,s]==1) {
                  r= 1;
                }
                d= r-p;
                u= (alpha[s]*d)/gradient;
                k[i+1] = k[i]+u;
                ppd[s,j]= d;
                ppu[s,j]= u;
                ppgradient[s,j]= gradient;

                  // quality control on k
                  // remove low values to avoid blow-up
                  if (k[i+1]<.05) {
                      k[i+1]=.05;
                  }
                  // remove high values to avoid blow-up
                  if (k[i+1]>1) {
                      k[i+1]=1;
                  }
                  // remove nans
                  if (is_nan(k[i+1])==1) {
                      k[i+1]= k[i];
                  }

                }

                j = j+1;

          } else {
            
                ppk[s,j] = 666;
                ppx[s,j] = 666;
                ppCs[s,j] = 666;
                ppCo[s,j] = 666;
                ppEVr[s,j] = 666;
                ppPgamble[s,j] = 666;
                ppd[s,j] = 666;
                ppu[s,j] = 666;
                ppgradient[s,j] = 666;
                lik[s,j] = 666;
                log_lik[s] = log_lik[s] + log(.5);
                j = j+1;
                
          }
                
      }
    }
  }
  
}
