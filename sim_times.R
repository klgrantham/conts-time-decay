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
  Xmats <- list(SWdesmat(Tp), crxodesmat(Tp), plleldesmat(Tp))
  ctres <- laply(varmats, vartheta_ind_vec, Xmat=Xmats, Toeplitz=FALSE)
  vals <- data.frame(decay = 1-rs,
                     SW = ctres[,1],
                     crxo = scalefactor*ctres[,2],
                     pllel = scalefactor*ctres[,3])
  return(vals)
}

sim_results <- function(nsims, rs, Tp, m, rho0, type="uniform"){
  # Calculates the variance of the treatment effect under the model:
  #    - continuous time (ct)
  # using simulated arrival times, with trial designs:
  #    - stepped wedge (SW)
  #    - cluster randomised crossover (CRXO)
  #    - parallel (pllel)
  # and returns the average, minimum and maximum observed variances
  # for each design
  
  simvalslist <- replicate(nsims, variances_sim(rs, Tp, m, rho0, type), simplify=FALSE)
  # Save sim results
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2]
  save(simvalslist, file=paste0("plots/vars_nsims_", nsims, "_", type, "_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))
  
  nrows <- length(rs)
  SWvals <- matrix(numeric(0), nrows, nsims)
  crxovals <- matrix(numeric(0), nrows, nsims)
  pllelvals <- matrix(numeric(0), nrows, nsims)
  for (i in 1:nsims){
    SWvals[,i] <- simvalslist[[i]]$SW
    crxovals[,i] <- simvalslist[[i]]$crxo
    pllelvals[,i] <- simvalslist[[i]]$pllel
  }
  SWvals_sorted <- t(apply(SWvals, 1, sort))
  crxovals_sorted <- t(apply(crxovals, 1, sort))
  pllelvals_sorted <- t(apply(pllelvals, 1, sort))
  ub <- 0.975; ub_index <- round(ub*nsims)
  lb <- 0.025; lb_index <- round(lb*nsims)
  simvals <- data.frame(decay=1-rs,
                        SW.lower = SWvals_sorted[,lb_index],
                        SW.avg = apply(SWvals, 1, mean),
                        SW.upper = SWvals_sorted[,ub_index],
                        crxo.lower = crxovals_sorted[,lb_index],
                        crxo.avg = apply(crxovals, 1, mean),
                        crxo.upper = crxovals_sorted[,ub_index],
                        pllel.lower = pllelvals_sorted[,lb_index],
                        pllel.avg = apply(pllelvals, 1, mean),
                        pllel.upper = pllelvals_sorted[,ub_index])
  return(simvals)
}

##
Tp <- 4
m <- 50
rho0 <- 0.04
rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2]

simvals <- sim_results(nsims=1000, rs=seq(0.5, 1, 0.01), Tp=Tp, m=m, rho0=rho0, type="uniform")
save(simvals, file=paste0("plots/vars_nsims_1000_summary_uniform_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))

# TODO: Run these results
simvalsexp <- sim_results(nsims=1000, rs=seq(0.5, 1, 0.01), Tp=Tp, m=m, rho0=rho0, type="exponential")
save(simvalsexp, file=paste0("plots/vars_nsims_1000_summary_exponential_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))

sim_results_plots <- function(simvals, Tp, m, rho0, type="uniform"){
  
  # Convert for plotting
  simvals_long <- simvals %>%
    gather(Design, Variance, -decay) %>%
    separate(Design, c("Design", "measure")) %>%
    spread(measure, Variance)
  
  # Obtain results under evenly-spaced assumption
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2]
  load(paste0("plots/vars_T", Tp, "_m", m, "_rho", rho0char, ".Rda")); vars_assump <- varvals
  ctvarvals_long <- vars_assump %>%
    select(decay, starts_with('ct'), -ends_with('base')) %>%
    gather(Design, Variance, -decay, convert=TRUE) %>%
    separate(Design, c("model", "Design"), sep=2) %>%
    select(-model)
  varvals <- merge(simvals_long, ctvarvals_long)
  
  # TODO: Add comparison to ctvarvals_long back in
  # Plot results and compare to results under evenly-spaced assumption
  p <- ggplot(data=simvals_long, aes(x=decay, colour=Design, linetype=Design)) +
    geom_line(aes(y=avg), size=1.2) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Design), alpha=0.3, show.legend=FALSE) +
    expand_limits(y=0) +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                       labels = c("CRXO", "Parallel", "SW")) +
    scale_linetype_manual(values = c("twodash", "dashed", "solid"),
                          labels = c("CRXO", "Parallel", "SW")) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    labs(title="Variance of treatment effect estimator, continuous-time correlation decay",
         subtitle=paste0("Simulated measurement times,", type, " distribution")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=14),
          legend.position="bottom")
  ggsave(paste0("plots/conts_sim_", type, "_compare_T", Tp, "_m", m, ".pdf"), p, width=297, height=210, units="mm")
}
