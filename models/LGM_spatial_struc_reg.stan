# LGM data gen
# 19 June 2022
# Zachary Roman <zjr159811@gmail.com>
# # # # # # # # # # # # # # # #

data {
 int<lower =0 > N; // N
 int<lower =0 > Kx; // N measurement occ.
 matrix [N, Kx] x; // X matrix of measurements
 vector[N] Y; // Endogenous variable predicted by growth
 matrix<lower = 0, upper = 1>[N, N] W; // Contiguity matrix
}

parameters {
 vector<lower=0>[Kx] sigmax;
 vector<lower=0>[2] sigmaxi;
 real<lower=0> sigma;  // Structural Error
 cholesky_factor_corr[2] L1;
 matrix [N ,2] zi;
 vector[2] muxi1; // Means of Int and Slope
 real beta1;
 real beta2;
 real rho;
}

model {
 matrix [N, Kx] mux ; // E [ x | xi ]
 matrix[N,2] muxivector; // vector with copypaste means
 matrix [N ,2] xi; // Xi[,1] = intercept, xi[,2] = slope
 
 
 for(n in 1:N){for(j in 1:2){muxivector[n,j]=muxi1[j];}}

 xi = muxivector + zi * diag_pre_multiply(sigmaxi,L1)'; 

 mux [,1] = 1*xi[,1] + 0*xi[,2]; // Linear growth model
 mux [,2] = 1*xi[,1] + 1*xi[,2]; // Ugly coding, but human readable
 mux [,3] = 1*xi[,1] + 2*xi[,2];
 mux [,4] = 1*xi[,1] + 3*xi[,2];
 mux [,5] = 1*xi[,1] + 4*xi[,2];

for(z in 1:Kx){x[,z] ~ normal(mux[,z], sigmax[z]);}
  to_vector(zi) ~ normal(0,1);     
  sigmax ~ cauchy(0,2.5); // Res. Var. observed measurement occ.   
  sigma ~ cauchy(0,2.5); // Res. Var. structural equation
  sigmaxi ~ cauchy(0,2.5); // Res. Var. exogenous factors
  L1 ~ lkj_corr_cholesky(2); // Cholesky decomposition for efficiency
  muxi1 ~ normal(0,1); // Mean of factors
  rho ~ uniform(-1,1); // Spatial auto-regressive effect
  beta1 ~ normal(0,1); // Structural slopes
  beta2 ~ normal(0,1);
  
  // Structural regression
  // rho*(W*Y) is the saptial lag of the endogenous variable.
  Y ~ normal(beta1*xi[,1] + beta2*xi[,2] + rho*(W*Y), sigma);
}

generated quantities{
  matrix[2,2] phi;
  phi = diag_pre_multiply(sigmaxi,L1)*diag_pre_multiply(sigmaxi,L1)';

}







