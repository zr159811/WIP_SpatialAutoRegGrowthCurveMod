# LGM data gen
# 19 June 2022
# Zachary Roman <zjr159811@gmail.com>
# # # # # # # # # # # # # # # #
# Analysis
# # # # # # # #
library(rstan)

options(mc.cores = parallel::detectCores())
source("Spatial_LGM_data_gen.R")

# LGM Parameters
# Av_Int <-  1 # Average intercept
# Av_Slope <-  0.45 # Average slope
# Var_Int <- 1 # Variance of intercept
# Var_Slope <- 1  # Variance of slope
# cov_int_slope <- 0.3 # Covariance of slope and intercept
# Nx <- Nx # Number of measurement occasions
# rho <- rho # Spatial autocorr magnitude
# Beta <-  matrix(c(0.2,0.6)) # Relationship between slope predicting Y


datenstan <- RspatLGM(N1 = 20,
                      N2 = 20,
                      Nx = 5,
                      rho = -0.5)

# W check for sanity
if(sum(!rowSums(datenstan$W) == 1)==0){
  print("W matrix is properly Row Standardized")
}else{print("Issue with normalization of W")}


synt <- stanc("LGM_spatial_struc_reg.stan")
LGM <- stan_model( stanc_ret = synt, verbose = FALSE)

fit <- sampling(LGM,
                data = datenstan,
                chains = 2,
                iter = 2000,
                warmup = 1000,
                save_warmup = FALSE)

pars <- c("muxi1[1]",
          "muxi1[2]",
          "beta2",
          "beta1",
          "rho",
          "phi[1,1]",
          "phi[2,1]",
          "phi[2,2]")

summary(fit, par = pars)$summary


# Test effect of ignoring spatial dependence in Y

synt <- stanc("LGM_strucreg.stan")
LGM <- stan_model( stanc_ret = synt, verbose = FALSE)

fit2 <- sampling(LGM,
                data = datenstan,
                chains = 2,
                iter = 2000,
                warmup = 1000,
                save_warmup = FALSE)

pars <- c("muxi1[1]",
          "muxi1[2]",
          "beta2",
          "beta1",
          "phi[1,1]",
          "phi[2,1]",
          "phi[2,2]")

summary(fit2, par = pars)$summary


