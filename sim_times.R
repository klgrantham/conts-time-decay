# Simulates arrival times and calculates variance of treatment effect
# Note: Does not assume evenly spaced arrival times
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)

source('vartheta_twolevels.R')

Vicontsim <- function(Tp, m, rho0, type="uniform"){
  # Returns a function for the variance matrix for a single cluster, Vi,
  # in terms of the base correlation, r (where rate of decay is 1-r),
  # under the continuous time model at the individual level with arrival times
  # randomly drawn from a uniform distribution, i.e. WITHOUT the assumption
  # of evenly spaced arrival times.
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster-period random effects

  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP

  # Generate arrival times
  if (type=="uniform"){
    x <- runif(Tp*m, 0, 1) # Generate vector of fractional times
    # Note: Does not sample from extreme values (should allow this)
    j <- rep(1:Tp, each=m) # Create vector of time period indices and add to fractional times
    times <- sort(j + x) # Combine to create arrival times
  } else if (type=="exponential"){
    x <- matrix(rexp(Tp*m), m, Tp)
    colmaxs <- apply(x, 2, max)
    y <- x %*% diag(1/colmaxs)
    colmins <- apply(y, 2, min)
    z <- t(t(y) - colmins)
    j <- matrix(1:Tp, m, Tp, byrow=TRUE)
    times <- sort(as.vector(z + j))
  }
  
  Vi <- function(r){
    A <- matrix(raw(), Tp*m, Tp*m)
    diag(sig2E, Tp*m) +
        sig2CP*(r^(abs(matrix(times[col(A)], Tp*m, Tp*m) -
                       matrix(times[row(A)], Tp*m, Tp*m))))
  }
  return(Vi)
}

variances_sim <- function(rs, Tp, m, rho0, type="uniform"){
  # Calculates the variance of the treatment effect under the
  # continuous time model at the individual level with simulated
  # arrival times from a uniform distribution, with trial designs:
  #    - stepped wedge (SW)
  #    - cluster randomised crossover (CRXO)
  #    - parallel (pllel)
  #    - parallel with baseline (pllelbase)
  
  Vs <- Vicontsim(Tp, m, rho0, type)
  varmats <- llply(rs, Vs) # Creates a list of covariance matrices for the values in rs
  scalefactor <- Tp/(Tp-1)
  vals <- data.frame(decay = 1-rs,
                     SW = laply(varmats, vartheta_ind, Xmat=SWdesmat(Tp), Toeplitz=FALSE),
                     crxo = scalefactor*laply(varmats, vartheta_ind, Xmat=crxodesmat(Tp), Toeplitz=FALSE),
                     pllel = scalefactor*laply(varmats, vartheta_ind, Xmat=plleldesmat(Tp), Toeplitz=FALSE),
                     pllelbase = scalefactor*laply(varmats, vartheta_ind, Xmat=pllelbasedesmat(Tp), Toeplitz=FALSE))
  return(vals)
}

sim_results <- function(nsims, rs, Tp, m, rho0, type="uniform"){
  # Calculates the variance of the treatment effect under the model:
  #    - continuous time (ct)
  # using simulated arrival times, with trial designs:
  #    - stepped wedge (SW)
  #    - cluster randomised crossover (CRXO)
  #    - parallel (pllel)
  #    - parallel with baseline (pllelbase)
  # and returns the average, minimum and maximum observed variances
  # for each design

  simvalslist <- replicate(nsims, variances_sim(rs, Tp, m, rho0, type), simplify=FALSE)
  
  nrows <- length(rs)
  SWvals <- matrix(numeric(0), nrows, nsims)
  crxovals <- matrix(numeric(0), nrows, nsims)
  pllelvals <- matrix(numeric(0), nrows, nsims)
  pllelbasevals <- matrix(numeric(0), nrows, nsims)
  for (i in 1:nsims){
    SWvals[,i] <- simvalslist[[i]]$SW
    crxovals[,i] <- simvalslist[[i]]$crxo
    pllelvals[,i] <- simvalslist[[i]]$pllel
    pllelbasevals[,i] <- simvalslist[[i]]$pllelbase
  }
  simvals <- data.frame(decay=1-rs,
                        SWmin = apply(SWvals, 1, min),
                        SWavg = apply(SWvals, 1, mean),
                        SWmax = apply(SWvals, 1, max),
                        crxomin = apply(crxovals, 1, min),
                        crxoavg = apply(crxovals, 1, mean),
                        crxomax = apply(crxovals, 1, max),
                        pllelmin = apply(pllelvals, 1, min),
                        pllelavg = apply(pllelvals, 1, mean),
                        pllelmax = apply(pllelvals, 1, max),
                        pllelbasemin = apply(pllelbasevals, 1, min),
                        pllelbaseavg = apply(pllelbasevals, 1, mean),
                        pllelbasemax = apply(pllelbasevals, 1, max))
  return(simvals)
}


sim_results_plots <- function(nsims=100, rs=seq(0.5, 1, 0.01), Tp, m, rho0, type="uniform"){

  simvals <- sim_results(nsims=nsims, rs=rs, Tp=Tp, m=m, rho0=rho0, type=type)
  save(simvals, file=paste0("plots/vars_sim_", type, "_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))

  # Convert for plotting
  simvals_long <- simvals %>%
                    gather(col, variance, -decay) %>%
                    separate(col, c("design", "measure"), sep=-4) %>%
                    spread(measure, variance)

  # Obtain results under evenly-spaced assumption
  load(paste0("plots/vars_T", Tp, "_m", m, ".Rda")); vars_assump <- varvals
  ctvarvals_long <- vars_assump %>%
                      select(decay, starts_with('ct')) %>%
                      gather(design, variances, ctSW:ctpllelbase, convert=TRUE) %>%
                      separate(design, c("model", "design"), sep=2) %>%
                      select(-model)
  varvals <- merge(simvals_long, ctvarvals_long)
  
  # Plot results and compare to results under evenly-spaced assumption
  pctall <- ggplot(data=varvals, aes(x=decay, group=design, colour=design)) +
    geom_line(aes(y=avg), size=1.0) +
    geom_line(aes(y=variances), size=1.0, colour='black', linetype="longdash") +
    geom_ribbon(aes(ymin=min, ymax=max, fill=design), alpha=0.3, show.legend=FALSE) +
    expand_limits(y=0) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    labs(title="Variances of treatment effect, continuous time",
         subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
    scale_color_hue(labels = c("CRXO", "Parallel", "Parallel w/ baseline", "SW")) +
    guides(colour=guide_legend("Design")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.title=element_text(size=16), legend.text=element_text(size=14),
          legend.position="bottom")
  ggsave(paste0("plots/conts_sim_", type, "_compare_T", Tp, "_m", m, ".pdf"), pctall, width=297, height=210, units="mm")
}

sim_results_plots(nsims=100, rs=seq(0.5, 1, 0.01), Tp=4, m=50, rho0=0.035, type="uniform")
sim_results_plots(nsims=100, rs=seq(0.5, 1, 0.01), Tp=4, m=50, rho0=0.035, type="exponential")
sim_results_plots(nsims=100, rs=seq(0.5, 1, 0.01), Tp=4, m=500, rho0=0.035, type="uniform")
sim_results_plots(nsims=100, rs=seq(0.5, 1, 0.01), Tp=4, m=500, rho0=0.035, type="exponential")
