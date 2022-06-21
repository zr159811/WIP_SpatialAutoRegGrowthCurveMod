
# LGM data gen
# 18 May 2022
# Zachary Roman <zjr159811@gmail.com>
# # # # # # # # # # # # # # # #
library(mvtnorm)
library(semPlot)
library(lavaan)
library(ggplot2)

Av_Int <-  10
Av_Slope <-  20

Var_Int <- 5
Var_Slope <- 10 

Nx <- 5

lx <- matrix(0,Nx,2)
lx[1:Nx,1] <- 1
lx[1:Nx,2] <- 0:(Nx-1) 

# td <- te <- .25 # rel of .80  
td <- te <- .25

N = 400
I = diag(N)
sigmaxi <- diag(2) # xi variance = 1
sigmaxi[1,1] <- Var_Int
sigmaxi[2,2] <- Var_Slope

sigmaxi[1,2] <- sigmaxi[2,1] <- 5 # covariance of slope and intercept

xi   <- rmvnorm(N,c(Av_Int,Av_Slope),sigmaxi)

dx    <- rmvnorm(N,rep(0,Nx),td*diag(Nx))

x <- data.frame(xi%*%t(lx) + dx)

mod <- "
INT =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5
SLOPE =~ 0*X1 + 1*X2 + 2*X3 + 3*X4 + 4*X5
"
summary(growth(mod,x))

# Quadratic
# # # # # # # #

# Number of measurement occasions
Nx <- 5

# Average: Int, Slope, Quad
MU <- c(2,2,4)

# Variance: Int, Slope, Quad
Sigma <- c(4,5,6)

# Covariances of LVs
Int_Slope <- 2
Int_Quad <- 3
Slope_Quad <- 4 

sigmaxi <- diag(3)

sigmaxi[1,2] <- sigmaxi[2,1] <- Int_Slope
sigmaxi[1,3] <- sigmaxi[3,1] <- Int_Quad
sigmaxi[2,3] <- sigmaxi[3,2] <- Slope_Quad


lx <- matrix(0,Nx,3)
lx[1:Nx,1] <- 1
lx[1:Nx,2] <- 0:(Nx-1) 
lx[1:Nx,3] <- (0:(Nx-1))^2

# td <- te <- .25 # rel of .80  
td <- te <- .25

N = 4000
I = diag(N)

diag(sigmaxi) <- Sigma 

xi   <- rmvnorm(N,MU,sigmaxi)

dx    <- rmvnorm(N,rep(0,Nx),td*diag(Nx))


x <- data.frame(xi%*%t(lx) + dx)

mod <- "
INT =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5
SLOPE =~ 0*X1 + 1*X2 + 2*X3 + 3*X4 + 4*X5
QUAD =~ 0*X1 + 1*X2 + 4*X3 + 9*X4 + 16*X5
"
summary(growth(mod,x))


# Linear with slope predicting an outcome variable
# # # # # # # # # 

# SD of structural regression
sd_zeta <- 2
Av_Int <-  10
Av_Slope <-  20
Var_Int <- 5
Var_Slope <- 10 
cov_int_slope <- 5
Nx <- 5

lx <- matrix(0,Nx,2)
lx[1:Nx,1] <- 1
lx[1:Nx,2] <- 0:(Nx-1) 

# td <- te <- .25 # rel of .80  
td <- te <- .25

N = 500
I = diag(N)
sigmaxi <- diag(2) # xi variance = 1
sigmaxi[1,1] <- Var_Int
sigmaxi[2,2] <- Var_Slope

sigmaxi[1,2] <- sigmaxi[2,1] <- cov_int_slope # covariance of slope and intercept

xi   <- rmvnorm(N,c(Av_Int,Av_Slope),sigmaxi)

dx    <- rmvnorm(N,rep(0,Nx),td*diag(Nx))

x <- data.frame(xi%*%t(lx) + dx)
Y<- 3*xi[,2] + rnorm(N,0,2)


mod <- "
INT =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5
SLOPE =~ 0*X1 + 1*X2 + 2*X3 + 3*X4 + 4*X5

Y ~ SLOPE
"
summary(mod <- growth(mod,data.frame(x,Y)))

semPaths(mod)


# Rstan Analysis
# # # # # # # # # #

library(rstan)
#options(mc.cores = parallel::detectCores())
#rstan_options(disable_march_warning = FALSE)

synt <- stanc("LGM_strucreg.stan")
LGM <- stan_model( stanc_ret = synt, verbose = FALSE)


datenstan <- list("N" = N,
                  "Kx" = Nx,
                  "x" = x,
                  "Y" = Y)

fit <- sampling(LGM,
                data = datenstan,
                chains = 2,
                iter = 1000,
                warmup = 500,
                save_warmup = FALSE)

# output from model


pars <-summary(fit, pars = c("phi",
                             "beta",
                             "muxi1",
                             "sigmax",
                             "sigmaxi"))$summary

pars


# Recovers nicely
# Moving to spatial adaption here
# # # # # # # # # # # # # # # # #


# # # # # # # #
# Spatial Data Generation

N1 = 20
N2 = 20
N = N1*N2

# Makes W matrix
x <- rep(1:N1)
y <- rep(1:N2)
dat <- as.matrix(expand.grid(x,y))
W_d1 <- as.matrix(dist(dat,method = "euclidean",diag = F), ncol = N1*N2)
W_d <- 1/W_d1
diag(W_d) <- rep(0, length(diag(W_d)))
diag(W_d) <- rep(0, length(diag(W_d)))
W <- ifelse(W_d >= 1, 1,0) # Cut of 1 = first order contiguity matrix

# LGM Parameters
Av_Int <-  10 # average intercept
Av_Slope <-  20 # average slope
Var_Int <- 5 # Variance of intercept
Var_Slope <- 10  # Variance of slope
cov_int_slope <- 5 # covariance of slope and intercept
Nx <- 5 # Number of measurement occasions
rho <- 0.6 # Spatial autocorr magnitude
beta <-  0.3 # Relationship between slope predicting Y

lx <- matrix(0,Nx,2)
lx[1:Nx,1] <- 1
lx[1:Nx,2] <- 0:(Nx-1) 

td <- te <- .25

I = diag(N)
sigmaxi <- diag(2) # xi variance = 1
sigmaxi[1,1] <- Var_Int
sigmaxi[2,2] <- Var_Slope
sigmaxi[1,2] <- sigmaxi[2,1] <- cov_int_slope # covariance of slope and intercept
xi   <- rmvnorm(N,c(Av_Int,Av_Slope),sigmaxi)
dx    <- rmvnorm(N,rep(0,Nx),td*diag(Nx))
x <- data.frame(xi%*%t(lx) + dx)
yz <- rnorm(N,0,1)
Y <- beta*xi[,2] + rho*(W%*%yz) + yz


str(Y)
# Analysis
# # # # # # # #

synt <- stanc("LGM_spatial_struc_reg.stan")
LGM <- stan_model( stanc_ret = synt, verbose = FALSE)


datenstan <- list("N" = N,
                  "Kx" = Nx,
                  "x" = x,
                  "Y" = as.vector(Y),
                  "Wy" = as.vector(W%*%Y))

fit <- sampling(LGM,
                data = datenstan,
                chains = 2,
                iter = 1000,
                warmup = 500,
                save_warmup = FALSE)





