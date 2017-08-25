# Simulates arrival times and calculates variance of treatment effect
# Note: Does not assume evenly spaced arrival times
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

Tp <- 4
m <- 2
rho0 <- 0.035
r <- 0.6
Xmat <- SWdesmat(Tp)

Vicontsim <- function(r, Tp, m, rho0){ # Include meanlvl option? # Include option to generate matrix of times?
  # Constructs the variance matrix for a single cluster, Vi, under the
  # continuous time model at the individual level with arrival times
  # randomly drawn from a uniform distribution, i.e. WITHOUT the assumption
  # of evenly placed arrival times.
  #
  # Inputs:
  # r - base correlation term; rate of decay is given by (1-r)
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster-period random effects

  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP

  # Generate arrival times
  x <- runif(Tp*m, 0, 1) # Generate vector of fractional times
  # Note: Does not sample from extreme values (should allow this)
  j <- rep(1:Tp, each=m) # Create vector of time period indices and add to fractional times
  times <- sort(j + x) # Combine to create arrival times
  
  A <- matrix(raw(), Tp*m, Tp*m)
  Vi <- diag(sig2E, Tp*m) +
        sig2CP*(r^(abs(matrix(times[col(A)], Tp*m, Tp*m) -
                       matrix(times[row(A)], Tp*m, Tp*m))))
  return(Vi)
}

V <- Vicontsim(r, Tp, m, rho0)
var <- vartheta_ind(Vi = V, Xmat = Xmat, Toeplitz = FALSE)
