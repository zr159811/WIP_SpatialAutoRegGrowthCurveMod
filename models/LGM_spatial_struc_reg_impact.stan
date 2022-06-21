data {
 int<lower =0 > N;
 int<lower =0 > Kx;
 matrix [N, Kx] x;
 vector[N] Y;
 matrix<lower = 0, upper = 1>[N, N] W;
 //matrix[N, N] I; 
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


for(i in 1:Kx){
 mux [,i] = 1*xi[,1] + (i-1)*xi[,2];// Linear
}


for(z in 1:Kx){x[,z] ~ normal(mux[,z], sigmax[z]);}
  to_vector(zi) ~ normal(0,1);     
  sigmax ~ cauchy(0,2.5);    
  sigma ~ cauchy(0,2.5); 
  sigmaxi ~ cauchy(0,2.5);
  L1 ~ lkj_corr_cholesky(2);
  muxi1[1] ~ normal(0,1); //Mean of factors
  muxi1[2] ~ normal(0,1); //Mean of factors 
  rho ~ normal(0,0.33) T[-1,1]; 
  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  
  // Structural regression
  Y ~ normal(beta1*xi[,1] + beta2*xi[,2] + rho*(W*Y), sigma);
}

generated quantities{
  matrix[2,2] phi;
  //matrix[N,N] CPDM;
  //real direct;
  //real indirect;
  //real total;
  //matrix[N,N] IrhoW;
  phi = diag_pre_multiply(sigmaxi,L1)*diag_pre_multiply(sigmaxi,L1)';

  // Spatial spillover (impact) calculations
  // b1[3] is effect of xi[2] and is the effect of interest
  // xi[1] is the moderator
  // b[4] is the interaction effect of xi[1] and xi[2]
  //IrhoW = inverse(I - rho*W);

  //CPDM =  IrhoW * I *(beta);
  //direct = mean(diagonal(CPDM));
  // Hacky way of computing mean of an off diagonal
  //indirect =  sum(add_diag(CPDM,0))/(N * (N - 1));  
  //total = (indirect + direct);


}







