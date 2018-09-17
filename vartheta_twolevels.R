# Investigating the variance of the treatment effect estimator
# under different models and designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ltsa)
library(plyr)

vartheta_mean <- function(Vi, Xmat){
  # Calculates the variance of the treatment effect, theta, for a model at the
  # cluster mean level with a particular treatment schedule
  #
  # Inputs:
  # Vi - a T x T covariance matrix for one cluster
  # Xmat - a K x T matrix of the treatment schedule (note: all elements either 0 or 1)
  
  K <- nrow(Xmat)
  Tp <- ncol(Xmat)
  Xvec <- as.vector(t(Xmat))
  Vi_inv <- solve(Vi)
  var <- 1/(t(Xvec) %*% (diag(1,K) %x% Vi_inv) %*% Xvec - 
            colSums(Xmat) %*% Vi_inv %*% (matrix(colSums(Xmat),nrow=Tp, ncol=1))/K)
  return(var)
}

vartheta_ind_vec <- function(Vi, Xmat, Toeplitz=TRUE){
  # Calculates the variance of the treatment effect, theta, for a model at the
  # individual level with a particular treatment schedule
  #
  # Inputs:
  # Vi - a Tm x Tm variance matrix for one cluster
  # Xmat - a vector of K x T matrices of the treatment schedule (note: all elements either 0 or 1)
  
  # If continuous time matrix, use Toeplitz inversion algorithm
  if(Toeplitz){
    Vi_inv <- TrenchInverse(Vi)
  } else{
    Vi_inv <- solve(Vi)
  }
  
  vars <- laply(Xmat, vartheta, Vi_inv)
  return(vars)
}

vartheta <- function(Xmat, Vi_inv) {
  # Returns variance of treatment effect estimator for an inverse covariance
  # matrix and a design matrix
  
  K <- nrow(Xmat)
  Tp <- ncol(Xmat)
  m <- nrow(Vi_inv)/Tp
  
  Q <- Xmat %x% t(rep(1,m))
  B <- colSums(Xmat) %x% rep(1,m)
  C <- diag(Tp) %x% rep(1,m)
  term1 <- sum(diag(Q %*% Vi_inv %*% t(Q)))
  term2 <- t(B) %*% Vi_inv %*% C
  term3 <- solve(t(C) %*% Vi_inv %*% C)
  term4 <- t(C) %*% Vi_inv %*% B
  var <- 1/(term1 - (1/K)*term2 %*% term3 %*% term4)
  return(var)
}

# Variance matrices under different models
# Uniform correlation (Hussey & Hughes), discrete time decay, continuous time decay

HHVi <- function(Tp, m, rho0, meanlvl=TRUE){
  # Constructs the covariance matrix for a single cluster, Vi, under the
  # Hussey & Hughes model (2007), at either the cluster mean level or
  # at the individual level
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of subjects per cluster-period
  # rho0 - base correlation between a pair of subjects' outcomes
  # meanlvl - boolean for whether to construct the covariance matrix for the
  #           cluster mean level (default) or the individual level

  totalvar <- 1
  tau2 <- rho0*totalvar
  sig2E <- totalvar - tau2
  sig2 <- sig2E/m # subject-specific variance at mean level
  if(meanlvl==TRUE){
    Vi <- diag(sig2,Tp) + tau2*matrix(1, nrow=Tp, ncol=Tp)
  }
  else{
    Vi <- diag(sig2E,Tp*m) + tau2*matrix(1, nrow=Tp*m, ncol=Tp*m)
  }
  return(Vi)
}

expdecayVi <- function(r, Tp, m, rho0, meanlvl=TRUE){
  # Constructs the covariance matrix for a single cluster, Vi, under the
  # discrete-time decay model (Kasza et al 2017), at either the cluster-
  # period mean level or at the individual level
  #
  # Inputs:
  # r - proportionate reduction in correlation
  # Tp - number of time periods
  # m - number of subjects per cluster-period
  # rho0 - base correlation between a pair of subjects' outcomes
  # meanlvl - boolean for whether to construct the covariance matrix for the
  #           cluster-period mean level (default) or the individual level

  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  if(meanlvl==TRUE){
    Vi <- diag(sig2,Tp) +
            sig2CP*(r^abs(matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=FALSE) -
                          matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=TRUE)))    
  }
  else{
    Vi <- diag(sig2E,Tp*m) +
            sig2CP*(r^abs(matrix(rep(1:Tp, each=m), nrow=Tp*m, ncol=Tp*m, byrow=FALSE) -
                          matrix(rep(1:Tp, each=m), nrow=Tp*m, ncol=Tp*m, byrow=TRUE)))
  }
  return(Vi)
}

expdecayVicont <- function(r, Tp, m, rho0, meanlvl=FALSE){
  # Constructs the covariance matrix for a single cluster, Vi, under the
  # continuous-time correlation decay model, at either the cluster-
  # period mean level or at the individual level
  #
  # Inputs:
  # r - proportionate reduction in correlation
  # Tp - number of time periods
  # m - number of subjects per cluster-period
  # rho0 - base correlation between a pair of subjects' outcomes
  # meanlvl - boolean for whether to construct the covariance matrix for the
  #           cluster-period mean level or the individual level (default)
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  if(meanlvl==TRUE){
    if(r==1){
      Vi <- expdecayVi(r=r, Tp=Tp, m=m, rho0=rho0, meanlvl=meanlvl) # Reverts to Hussey and Hughes
    }
    else{
      vars <- (1/m^2)*(m + 2*((r^(1/m))*(-m*(r^(1/m)) + m + r - 1))/((r^(1/m)) - 1)^2) # to be multiplied by sig2CP and added to sig2
      covars <- (1/m^2)*(((r^((1/m)-1))*(r-1)^2)/((r^(1/m)) - 1)^2) # to be multiplied by sig2CP and r^{|j-j'|}
      mat <- matrix(covars, nrow=Tp, ncol=Tp) - diag(covars,Tp) + diag(vars,Tp)
      Vi <- diag(sig2,Tp) +
            (sig2CP*(r^abs(matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=FALSE) -
                           matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=TRUE)))*mat)  
    }
  }
  else{
      Vi <- diag(sig2E,Tp*m) +
            sig2CP*(r^(abs(matrix(rep(1:(Tp*m)), nrow=Tp*m, ncol=Tp*m, byrow=FALSE) -
                           matrix(rep(1:(Tp*m)), nrow=Tp*m, ncol=Tp*m, byrow=TRUE))/m))
  }
  return(Vi)
}

# Design matrices for different trial designs
# Stepped wedge (SW), cluster randomised crossover (CRXO), parallel

SWdesmat <- function(Tp, N=Tp-1){
  Xsw <- matrix(data=0, ncol=Tp, nrow=(Tp-1))
  for(i in 1:(Tp-1)){
    Xsw[i,(i+1):Tp] <- 1
  }
  if(N%%(Tp-1) != 0) stop('N must be a multiple of Tp-1')
  Xsw %x% rep(1, N/(Tp-1))
}

crxodesmat <- function(Tp, N=Tp){
  if(N%%2 != 0) stop('N must be even')
  Xcrxo <- matrix(data=0, ncol=Tp, nrow=2)
  Xcrxo[1, seq(1,Tp,2)] <- 1
  Xcrxo[2, seq(2,Tp,2)] <- 1
  Xcrxo %x% rep(1, N/2)
}

plleldesmat <- function(Tp, N=Tp){
  if(N%%2 != 0) stop('N must be even')
  Xpllel <- matrix(data=0, ncol=Tp, nrow=2)
  Xpllel[1,] <- 1
  Xpllel %x% rep(1, N/2)
}


var_ct_mean_results <- function(Tp, N, m, rho0){
  # Get variances using cluster-mean-level covariance matrices
  rs <- seq(0.5, 1, 0.01)
  # Specify the covariance matrices under the continuous-time model
  ctmeanvarmat <- llply(rs, expdecayVicont, Tp, m, rho0, meanlvl=TRUE)
  
  # Get the variance of the treatment effect estimator under the
  # different models and designs
  SWXmat <- SWdesmat(Tp, N)
  crxoXmat <- crxodesmat(Tp, N)
  pllelXmat <- plleldesmat(Tp, N)
  Xmats <- list(SWXmat, crxoXmat, pllelXmat)
  ctres <- laply(ctmeanvarmat, vartheta_ind_vec, Xmat=Xmats)
  varvals <- data.frame(decay = 1-rs,
                        ctmeanSW = ctres[,1],
                        ctmeancrxo = ctres[,2],
                        ctmeanpllel = ctres[,3])
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2] # get numbers after decimal point
  save(varvals, file=paste0("plots/vars_ct_mean_T", Tp, "_N", N, "_m", m, "_rho", rho0char, ".Rda"))
}

# Generate main variance results
generate_var_results  <- function(Tp, N, m, rho0_CD, rho0_UC) {
  # Calculates the variance of the treatment effect estimator under the models:
  #    continuous time (ct), discrete time (dt), Hussey & Hughes (HH)
  # with trial designs:
  #    stepped wedge (SW),
  #    cluster randomised crossover (CRXO),
  #    parallel (pllel)
  #
  # Inputs:
  # Tp - number of time periods in the trial
  # N - number of clusters (assumed common across designs)
  # m - number of subjects measured in each cluster-period
  # rho0_CD - base correlation between a pair of subjects' outcomes, under
  #           decaying correlation models
  # rho0_UC - correlation between any pair of subjects' outcomes, under
  #           uniform correlation model
  # Example usage: vals <- generate_var_results(Tp=4, N=6, m=50, rho0_CD=0.023, rho0_UC=0.019)
  
  # Set vector of r values (Decay = 1-r)
  rs <- seq(0.5, 1, 0.01)
  # Specify the covariance matrices under the different models
  ctvarmat <- llply(rs, expdecayVicont, Tp, m, rho0_CD, meanlvl=FALSE)
  dtvarmat <- llply(rs, expdecayVi, Tp, m, rho0_CD, meanlvl=TRUE)
  HHvarmat <- HHVi(Tp, m, rho0_UC, meanlvl=TRUE)

  # Get the variances of the treatment effect estimator under the
  # different models and designs
  SWXmat <- SWdesmat(Tp, N)
  crxoXmat <- crxodesmat(Tp, N)
  pllelXmat <- plleldesmat(Tp, N)
  Xmats <- list(SWXmat, crxoXmat, pllelXmat)
  ctres <- laply(ctvarmat, vartheta_ind_vec, Xmat=Xmats)
  varvals <- data.frame(decay=1-rs,
                        ctSW = ctres[,1],
                        ctcrxo = ctres[,2],
                        ctpllel = ctres[,3],
                        dtSW = laply(dtvarmat, vartheta_mean, Xmat=SWXmat),
                        dtcrxo = laply(dtvarmat, vartheta_mean, Xmat=crxoXmat),
                        dtpllel = laply(dtvarmat, vartheta_mean, Xmat=pllelXmat),
                        HHSW = vartheta_mean(Vi=HHvarmat, Xmat=SWXmat),
                        HHcrxo = vartheta_mean(Vi=HHvarmat, Xmat=crxoXmat),
                        HHpllel = vartheta_mean(Vi=HHvarmat, Xmat=pllelXmat))

  # Save results to R data file
  rho0CDchar <- strsplit(as.character(rho0_CD),"\\.")[[1]][2] # get numbers after decimal point
  rho0UCchar <- strsplit(as.character(rho0_UC),"\\.")[[1]][2]
  save(varvals, file=paste0("plots/vars_T", Tp, "_N", N, "_m", m, "_rhoCD", rho0CDchar, "_rhoUC", rho0UCchar, ".Rda"))
  return(varvals)
}
