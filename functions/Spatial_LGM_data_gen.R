# LGM data gen
# 19 May 2022
# Zachary Roman <zjr159811@gmail.com>
# # # # # # # # # # # # # # # #

RspatLGM <- function(N1 = 10, N2 = 10, Nx = 5, rho = 0.45, Av_Int = 1, 
                     Av_Slope = 1, Var_Int = 1, Var_Slope = 1,
                     cov_int_slope = 0.3){
require(mvtnorm)
require(semPlot)
require(lavaan)
require(ggplot2)
require(rstan)
N1 = N1
N2 = N2
N = N1*N2

# Makes W matrix
x <- rep(1:N1)
y <- rep(1:N2)
dat <- as.matrix(expand.grid(x,y))
W_d1 <- as.matrix(dist(dat,method = "euclidean",diag = F), ncol = N1*N2)
W_d <- 1/W_d1
diag(W_d) <- rep(0, length(diag(W_d)))
diag(W_d) <- rep(0, length(diag(W_d)))
W_c <- ifelse(W_d >= 1, 1,0) # Cut of 1 = first order contiguity matrix

# Normalizes W
W <- W_c/rowSums(W_c)

# Sets up "factor loadings" (i.e., in this context growth rate)
lx <- matrix(0,Nx,2)
lx[1:Nx,1] <- 1
lx[1:Nx,2] <- 0:(Nx-1) 

td <- te <- .25 

I = diag(N)
sigmaxi <- diag(2) 
sigmaxi[1,1] <- Var_Int
sigmaxi[2,2] <- Var_Slope
sigmaxi[1,2] <- sigmaxi[2,1] <- cov_int_slope # covariance of slope and intercept
xi   <- rmvnorm(N,c(Av_Int,Av_Slope),sigmaxi)
dx    <- rmvnorm(N,rep(0,Nx),td*diag(Nx))
x <- data.frame(xi%*%t(lx) + dx)
zeta <- rnorm(N,0,0.25)

# data gen model, see Lesage and Pace intro pg. 16 for more information
# on the data generating model in observed spatial-autoreg models.
Y <- (solve(I - rho*W))%*%xi%*%Beta + (solve(I-rho*W))%*%zeta

datenstan <- list("N" = N,
                  "Kx" = Nx,
                  "x" = x,
                  "Y" = as.vector(Y),
                  "W" = W,
                  "I" = I)

return(datenstan)

}
