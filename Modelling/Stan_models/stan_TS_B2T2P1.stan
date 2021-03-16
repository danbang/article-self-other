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
  real n_thetaz_real;
  real n_dthetaz_real;
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
  real mu_rho_lb;
  real mu_rho_ub;
  real sd_rho_lb;
  real sd_rho_ub;
  real mu_rho_prior_mu;
  real mu_rho_prior_sd;
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
    real<lower=mu_rho_lb,upper=mu_rho_ub> mu_rho;
    real<lower=sd_rho_lb,upper=sd_rho_ub> sd_rho;
    real<lower=mu_beta0_s_lb,upper=mu_beta0_s_ub> mu_beta0_s;
    real<lower=mu_beta0_o_lb,upper=mu_beta0_o_ub> mu_beta0_o;
    real<lower=sd_beta0_lb,upper=sd_beta0_ub> sd_beta0;
    real<lower=mu_beta1_s_lb,upper=mu_beta1_s_ub> mu_beta1_s;
    real<lower=mu_beta1_o_lb,upper=mu_beta1_o_ub> mu_beta1_o;
    real<lower=sd_beta1_lb,upper=sd_beta1_ub> sd_beta1;
    real<lower=mu_factor_lb,upper=mu_factor_ub> mu_factor;
    real<lower=sd_factor_lb,upper=sd_factor_ub> sd_factor;

    // subject level
    row_vector<lower=mu_sigma_lb,upper=mu_sigma_ub>[Ns] sigma;
    row_vector<lower=mu_alpha_lb,upper=mu_alpha_ub>[Ns] alpha;
    row_vector<lower=mu_rho_lb,upper=mu_rho_ub>[Ns] rho;
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
  mu_rho ~ normal(mu_rho_prior_mu,mu_rho_prior_sd);
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
        rho[s] ~ normal(mu_rho,sd_rho);
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
                real Pd;
                real Pk;
                real Px_sCONDstate_d1[n_thetaz];
                real Px_sCONDstate_d2[n_thetaz];
                real Px_sCONDstate_d1_a1[n_thetaz];
                real Px_sCONDstate_d1_a2[n_thetaz];
                real Px_sCONDstate_d2_a1[n_thetaz];
                real Px_sCONDstate_d2_a2[n_thetaz];
                real curr_mu_d1_a1;
                real curr_mu_d1_a2;
                real curr_mu_d2_a1;
                real curr_mu_d2_a2;
                real curr_sigma_d1_a1;
                real curr_sigma_d1_a2;
                real curr_sigma_d2_a1;
                real curr_sigma_d2_a2;
                real Pa_oCONDstate_x_s_d1_a1[n_thetaz];
                real Pa_oCONDstate_x_s_d1_a2[n_thetaz];
                real Pa_oCONDstate_x_s_d2_a1[n_thetaz];
                real Pa_oCONDstate_x_s_d2_a2[n_thetaz];
                real tmp_d1_a1[n_thetaz];
                real tmp_d1_a2[n_thetaz];
                real tmp_d2_a1[n_thetaz];
                real tmp_d2_a2[n_thetaz];
                real P_d1_a1;
                real P_d1_a2;
                real P_d2_a1;
                real P_d2_a2;
                real normalizer;
                real P_d1_a1_x_s;
                real P_d1_a2_x_s;
                real P_d2_a1_x_s;
                real P_d2_a2_x_s;
                real Px[n_dthetaz];
                real term1;
                real term2;
                real term3;
                real term4;
                real nominator;
                real denominator;
                real combine[n_thetaz];
                real gradientNeg;
                real gradientPos;
                real gradient;
                
                // draw sensory sample //
                y[b,i,s] ~ normal(dtheta[b,i,s],sigma[s]);
                
                // subjective priors
                Pd = .5;
                Pk = 1/n_thetaz_real;
                
                // likelihood of y_s given state
                for (kk in 1:n_thetaz) {
                    Px_sCONDstate_d1[kk] = exp(normal_lpdf(y[b,i,s]|-1*spaceFull[kk],sigma[s]));
                    Px_sCONDstate_d2[kk] = exp(normal_lpdf(y[b,i,s]|+1*spaceFull[kk],sigma[s]));
                }
                
                // augment with actions for later calculation
                Px_sCONDstate_d1_a1 = Px_sCONDstate_d1;
                Px_sCONDstate_d1_a2 = Px_sCONDstate_d1;
                Px_sCONDstate_d2_a1 = Px_sCONDstate_d2;
                Px_sCONDstate_d2_a2 = Px_sCONDstate_d2;
                
                // posterior probability of LEFT and RIGHT
                Pstate1CONDx_s= sum(Px_sCONDstate_d1)/(sum(Px_sCONDstate_d1)+sum(Px_sCONDstate_d2));
                Pstate2CONDx_s= sum(Px_sCONDstate_d2)/(sum(Px_sCONDstate_d1)+sum(Px_sCONDstate_d2));
            
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
                
                // probability of other action given sensation and state and coherence
                for (kk in 1:n_thetaz) {
                      
                      curr_mu_d1_a1 = -1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (y[b,i,s] - -1*spaceFull[kk]);
                      curr_sigma_d1_a1 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d1_a1[kk] = 0 + 1*normal_cdf(0,curr_mu_d1_a1,curr_sigma_d1_a1);
                      
                      curr_mu_d1_a2 = -1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (y[b,i,s] - -1*spaceFull[kk]);
                      curr_sigma_d1_a2 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d1_a2[kk] = 1 + -1*normal_cdf(0,curr_mu_d1_a2,curr_sigma_d1_a2);
                      
                      curr_mu_d2_a1 = +1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (y[b,i,s] - +1*spaceFull[kk]);
                      curr_sigma_d2_a1 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d2_a1[kk] = 0 + 1*normal_cdf(0,curr_mu_d2_a1,curr_sigma_d2_a1);
                      
                      curr_mu_d2_a2 = +1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (y[b,i,s] - +1*spaceFull[kk]);
                      curr_sigma_d2_a2 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d2_a2[kk] = 1 + -1*normal_cdf(0,curr_mu_d2_a2,curr_sigma_d2_a2);

                }
                
                // probability of each other action and state
                for (kk in 1:n_thetaz) {
                  tmp_d1_a1[kk] = Pd * Pk * (Px_sCONDstate_d1_a1[kk] * Pa_oCONDstate_x_s_d1_a1[kk]);
                  tmp_d1_a2[kk] = Pd * Pk * (Px_sCONDstate_d1_a2[kk] * Pa_oCONDstate_x_s_d1_a2[kk]);
                  tmp_d2_a1[kk] = Pd * Pk * (Px_sCONDstate_d2_a1[kk] * Pa_oCONDstate_x_s_d2_a1[kk]);
                  tmp_d2_a2[kk] = Pd * Pk * (Px_sCONDstate_d2_a2[kk] * Pa_oCONDstate_x_s_d2_a2[kk]);
                }
                P_d1_a1 = sum(tmp_d1_a1);
                P_d1_a2 = sum(tmp_d1_a2);
                P_d2_a1 = sum(tmp_d2_a1);
                P_d2_a2 = sum(tmp_d2_a2);
                
                // normalizing factor
                normalizer = Pd * Pk * (sum(Px_sCONDstate_d1)+sum(Px_sCONDstate_d2));
                
                //  probability of each other action and state given sensation
                P_d1_a1_x_s = P_d1_a1 / normalizer;
                P_d1_a2_x_s = P_d1_a2 / normalizer;
                P_d2_a1_x_s = P_d2_a1 / normalizer;
                P_d2_a2_x_s = P_d2_a2 / normalizer;
                
                // finally probability other is correct
                PcorCONDx_o = P_d1_a1_x_s + P_d2_a2_x_s;
                
                // remove 0s and 1s to avoid infinite values
                if (PcorCONDx_s<.01) {
                    PcorCONDx_s=.01;
                } else if (PcorCONDx_s>.99) {
                    PcorCONDx_s=.99;
                } 
                if (PcorCONDx_o<.01) {
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
                
                // Gradient
                // Calculate separately for d=-1 and d=1
                
                // d=-1
                // term 1
                for (kk in 1:n_dthetaz) {
                  Px[kk]= exp(normal_lpdf(y[b,i,s]|spaceFullSigned[kk],sigma[s])) * Pk * Pd;
                  // Px[kk]= Pd * Pk * exp(normal_lpdf(y[b,i,s]|-1*spaceFullSigned[kk],sigma[s]));
                } 
                term1= Pd/sum(Px);
                // terms 2-4 are found by summing over theta 
                for (kk in 1:n_thetaz) {
                  // term 2
                  term2= Pk;
                  // term 3
                  term3= exp(normal_lpdf(y[b,i,s]|-spaceFull[kk],sigma[s]));
                  // term 4
                  // calculate derivative for case where k=r=-1, dP(r=-1|k=-1, theta, x_s) / dsigma
                  nominator= -spaceFull[kk] * exp( (-(spaceFull[kk]- (k[i]/sigma[s])* rho[s]* (y[b,i,s]+spaceFull[kk]) )^2) / (2* (1-rho[s]^2)* k[i]^2));
                  denominator= sqrt(2* pi()* (1-rho[s]^2))* k[i]^2;
                  term4= nominator/denominator;
                  combine[kk]= term2* term3* term4;
                }
                gradientNeg= term1* sum(combine);
                
                // d=+1
                // term 1
                for (kk in 1:n_dthetaz) {
                  Px[kk]= exp(normal_lpdf(y[b,i,s]|spaceFullSigned[kk],sigma[s])) * Pk * Pd;
                } 
                term1= Pd/sum(Px);
                // terms 2-4 are found by summing over theta 
                for (kk in 1:n_thetaz) {
                  // term 2
                  term2= Pk;
                  // term 3
                  term3= exp(normal_lpdf(y[b,i,s]|spaceFull[kk],sigma[s]));
                  // term 4
                  // calculate derivative for case where k=r=-1, dP(r=-1|k=-1, theta, x_s) / dsigma
                  nominator= -spaceFull[kk] * exp( (-(spaceFull[kk]+ (k[i]/sigma[s])* rho[s]* (y[b,i,s]-spaceFull[kk]) )^2) / (2* (1-rho[s]^2)* k[i]^2));
                  denominator= sqrt(2* pi()* (1-rho[s]^2))* k[i]^2;
                  term4= nominator/denominator;
                  combine[kk]= term2* term3* term4;
                }
                gradientPos= term1* sum(combine);
              
                // sum gradients for d=-1 and d=1
                gradient= gradientNeg+gradientPos;
                
                // use gradient for updating
                p= PcorCONDx_o;
                r= acc[b,i,s];
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
                // print("rho=",rho[s])
                // print("alpha=",alpha[s])
                // print("p=",p)
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
                real Pd;
                real Pk;
                real Px_sCONDstate_d1[n_thetaz];
                real Px_sCONDstate_d2[n_thetaz];
                real Px_sCONDstate_d1_a1[n_thetaz];
                real Px_sCONDstate_d1_a2[n_thetaz];
                real Px_sCONDstate_d2_a1[n_thetaz];
                real Px_sCONDstate_d2_a2[n_thetaz];
                real curr_mu_d1_a1;
                real curr_mu_d1_a2;
                real curr_mu_d2_a1;
                real curr_mu_d2_a2;
                real curr_sigma_d1_a1;
                real curr_sigma_d1_a2;
                real curr_sigma_d2_a1;
                real curr_sigma_d2_a2;
                real Pa_oCONDstate_x_s_d1_a1[n_thetaz];
                real Pa_oCONDstate_x_s_d1_a2[n_thetaz];
                real Pa_oCONDstate_x_s_d2_a1[n_thetaz];
                real Pa_oCONDstate_x_s_d2_a2[n_thetaz];
                real tmp_d1_a1[n_thetaz];
                real tmp_d1_a2[n_thetaz];
                real tmp_d2_a1[n_thetaz];
                real tmp_d2_a2[n_thetaz];
                real P_d1_a1;
                real P_d1_a2;
                real P_d2_a1;
                real P_d2_a2;
                real normalizer;
                real P_d1_a1_x_s;
                real P_d1_a2_x_s;
                real P_d2_a1_x_s;
                real P_d2_a2_x_s;
                real Px[n_dthetaz];
                real term1;
                real term2;
                real term3;
                real term4;
                real nominator;
                real denominator;
                real combine[n_thetaz];
                real gradientNeg;
                real gradientPos;
                real gradient;
                
                // draw sensory sample
                x = y[b,i,s];
                ppx[s,j] = x;
                ppk[s,j] = k[i];
                
                // subjective priors
                Pd = .5;
                Pk = 1/n_thetaz_real;
                
                // likelihood of y_s given state
                for (kk in 1:n_thetaz) {
                    Px_sCONDstate_d1[kk] = exp(normal_lpdf(x|-1*spaceFull[kk],sigma[s]));
                    Px_sCONDstate_d2[kk] = exp(normal_lpdf(x|+1*spaceFull[kk],sigma[s]));
                }
                
                // augment with actions for later calculation
                Px_sCONDstate_d1_a1 = Px_sCONDstate_d1;
                Px_sCONDstate_d1_a2 = Px_sCONDstate_d1;
                Px_sCONDstate_d2_a1 = Px_sCONDstate_d2;
                Px_sCONDstate_d2_a2 = Px_sCONDstate_d2;
                
                // posterior probability of LEFT and RIGHT
                Pstate1CONDx_s= sum(Px_sCONDstate_d1)/(sum(Px_sCONDstate_d1)+sum(Px_sCONDstate_d2));
                Pstate2CONDx_s= sum(Px_sCONDstate_d2)/(sum(Px_sCONDstate_d1)+sum(Px_sCONDstate_d2));
            
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
                
                // probability of other action given sensation and state and coherence
                for (kk in 1:n_thetaz) {
                      
                      curr_mu_d1_a1 = -1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (x - -1*spaceFull[kk]);
                      curr_sigma_d1_a1 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d1_a1[kk] = 0 + 1*normal_cdf(0,curr_mu_d1_a1,curr_sigma_d1_a1);
                      
                      curr_mu_d1_a2 = -1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (x - -1*spaceFull[kk]);
                      curr_sigma_d1_a2 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d1_a2[kk] = 1 + -1*normal_cdf(0,curr_mu_d1_a2,curr_sigma_d1_a2);
                      
                      curr_mu_d2_a1 = +1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (x - +1*spaceFull[kk]);
                      curr_sigma_d2_a1 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d2_a1[kk] = 0 + 1*normal_cdf(0,curr_mu_d2_a1,curr_sigma_d2_a1);
                      
                      curr_mu_d2_a2 = +1*spaceFull[kk] + (k[i]/sigma[s]) * rho[s] * (x - +1*spaceFull[kk]);
                      curr_sigma_d2_a2 = sqrt((1-rho[s]^2) * (k[i]^2));
                      Pa_oCONDstate_x_s_d2_a2[kk] = 1 + -1*normal_cdf(0,curr_mu_d2_a2,curr_sigma_d2_a2);

                }
                
                // probability of each other action and state
                for (kk in 1:n_thetaz) {
                  tmp_d1_a1[kk] = Pd * Pk * (Px_sCONDstate_d1_a1[kk] * Pa_oCONDstate_x_s_d1_a1[kk]);
                  tmp_d1_a2[kk] = Pd * Pk * (Px_sCONDstate_d1_a2[kk] * Pa_oCONDstate_x_s_d1_a2[kk]);
                  tmp_d2_a1[kk] = Pd * Pk * (Px_sCONDstate_d2_a1[kk] * Pa_oCONDstate_x_s_d2_a1[kk]);
                  tmp_d2_a2[kk] = Pd * Pk * (Px_sCONDstate_d2_a2[kk] * Pa_oCONDstate_x_s_d2_a2[kk]);
                }
                P_d1_a1 = sum(tmp_d1_a1);
                P_d1_a2 = sum(tmp_d1_a2);
                P_d2_a1 = sum(tmp_d2_a1);
                P_d2_a2 = sum(tmp_d2_a2);
                
                // normalizing factor
                normalizer = Pd * Pk * (sum(Px_sCONDstate_d1)+sum(Px_sCONDstate_d2));
                
                //  probability of each other action and state given sensation
                P_d1_a1_x_s = P_d1_a1 / normalizer;
                P_d1_a2_x_s = P_d1_a2 / normalizer;
                P_d2_a1_x_s = P_d2_a1 / normalizer;
                P_d2_a2_x_s = P_d2_a2 / normalizer;
                
                // finally probability other is correct
                PcorCONDx_o = P_d1_a1_x_s + P_d2_a2_x_s;
                ppCo[s,j] = PcorCONDx_o;
                
                // remove 0s and 1s to avoid infinite values
                if (PcorCONDx_s<.01) {
                    PcorCONDx_s=.01;
                } else if (PcorCONDx_s>.99) {
                    PcorCONDx_s=.99;
                } 
                if (PcorCONDx_o<.01) {
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

                // Update k-value
                if (self[b,i,s] == 1) {
                  k[i+1] = k[i];
                  ppd[s,j] = 0;
                  ppu[s,j] = 0;
                } else if (self[b,i,s] == 0) {

                  // Gradient
                  // Calculate separately for d=-1 and d=1
                  
                  // d=-1
                  // term 1
                  for (kk in 1:n_dthetaz) {
                    Px[kk]= exp(normal_lpdf(x|spaceFullSigned[kk],sigma[s])) * Pk * Pd; 
                  } 
                  term1= Pd/sum(Px);
                  // terms 2-4 are found by summing over theta 
                  for (kk in 1:n_thetaz) {
                    // term 2
                    term2= Pk;
                    // term 3
                    term3= exp(normal_lpdf(x|-spaceFull[kk],sigma[s]));
                    // term 4
                    // calculate derivative for case where k=r=-1, dP(r=-1|k=-1, theta, x_s) / dsigma
                    nominator= -spaceFull[kk] * exp( (-(spaceFull[kk]- (k[i]/sigma[s])* rho[s]* (x+spaceFull[kk]) )^2) / (2* (1-rho[s]^2)* k[i]^2));
                    denominator= sqrt(2* pi()* (1-rho[s]^2))* k[i]^2;
                    term4= nominator/denominator;
                    combine[kk]= term2* term3* term4;
                  }
                  gradientNeg= term1* sum(combine);
                  
                  // d=+1
                  // term 1
                  for (kk in 1:n_dthetaz) {
                    Px[kk]= exp(normal_lpdf(x|spaceFullSigned[kk],sigma[s])) * Pk * Pd;
                  } 
                  term1= Pd/sum(Px);
                  // terms 2-4 are found by summing over theta 
                  for (kk in 1:n_thetaz) {
                    // term 2
                    term2= Pk;
                    // term 3
                    term3= exp(normal_lpdf(x|spaceFull[kk],sigma[s]));
                    // term 4
                    // calculate derivative for case where k=r=-1, dP(r=-1|k=-1, theta, x_s) / dsigma
                    nominator= -spaceFull[kk] * exp( (-(spaceFull[kk]+ (k[i]/sigma[s])* rho[s]* ((x-spaceFull[kk])) )^2) / (2* (1-rho[s]^2)* k[i]^2));
                    denominator= sqrt(2* pi()* (1-rho[s]^2))* k[i]^2;
                    term4= nominator/denominator;
                    combine[kk]= term2* term3* term4;
                  }
                  gradientPos= term1* sum(combine);
                
                  // sum gradients for d=-1 and d=1
                  gradient= gradientNeg+gradientPos;
                  
                  // use gradient for updating
                  p= PcorCONDx_o;
                  r= acc[b,i,s];
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
