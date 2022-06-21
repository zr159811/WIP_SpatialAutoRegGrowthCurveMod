data {
 int<lower =0 > N;
 int<lower =0 > Kx;
 matrix [N, Kx] x;
}

parameters {
 vector<lower=0>[Kx] sigmax;
 vector<lower=0>[2] sigmaxi;
 cholesky_factor_corr[2] L1;
 matrix [N ,2] zi;
 vector[2] muxi1; // Means of Int and Slope
 matrix[N,2] Xi0; // Factor Scores
}

model {
 matrix [N, Kx] mux ; // E [ x | xi ]
 matrix[N,2] muxivector; // vector with copypaste means
 matrix [N ,2] xi; // Xi[,1] = intercept, xi[,2] = slope
 
 
 for(n in 1:N){for(j in 1:2){muxivector[n,j]=muxi1[j];}}

 xi = muxivector + zi * diag_pre_multiply(sigmaxi,L1)'; 

 mux [,1] = 1*xi[,1] + 0*xi[,2];
 mux [,2] = 1*xi[,1] + 1*xi[,2];
 mux [,3] = 1*xi[,1] + 2*xi[,2];
 mux [,4] = 1*xi[,1] + 3*xi[,2];
 mux [,5] = 1*xi[,1] + 4*xi[,2];
 
 to_vector(Xi0) ~ normal(0,1);
 

for(z in 1:Kx){x[,z] ~ normal(mux[,z], sigmax[z]);}
  to_vector(zi) ~ normal(0,1);     
  sigmax ~ cauchy(0,2.5);          
  sigmaxi ~ cauchy(0,2.5);
  L1 ~ lkj_corr_cholesky(2);
  muxi1 ~ normal(0,1); //Mean of factors
}

generated quantities{
  matrix[2,2] phi;
  phi = diag_pre_multiply(sigmaxi,L1)*diag_pre_multiply(sigmaxi,L1)';

}







