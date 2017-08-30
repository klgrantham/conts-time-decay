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

Tp <- 4
m <- 50
rho0 <- 0.035

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

get_variance_sim <- function(r, Tp, m, rho0, nsims){
  # Calculates the variance of the treatment effect under the model:
  #    - continuous time (ct)
  # using simulated arrival times, with trial designs:
  #    - stepped wedge (SW)
  #    - cluster randomised crossover (CRXO)
  #    - parallel (pllel)
  #    - parallel with baseline (pllelbase)
  #
  # Inputs:
  # r - base correlation term
  # Tp - number of time periods in the trial
  # m - number of subjects measured in each time period
  # rho0 - base correlation between a pair of subjects
  #
  # Example usage: val <- get_variance_sim(r=0.6, Tp=4, m=50, rho0=0.035)

  # Generate covariance matrix using randomly generated arrival times
  ctvarmat <- replicate(nsims, Vicontsim(r, Tp, m, rho0), simplify=FALSE)

  # Get the variance of the treatment effect under the different designs
  # Scale the non-SW variances by (Tp/(Tp-1)) to account for uneven clusters across designs
  scalefactor <- Tp/(Tp-1)
  vals <- data.frame(SW = laply(ctvarmat, vartheta_ind, Xmat=SWdesmat(Tp), Toeplitz=FALSE),
                     crxo = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=crxodesmat(Tp), Toeplitz=FALSE),
                     pllel = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=plleldesmat(Tp), Toeplitz=FALSE),
                     pllelbase = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=pllelbasedesmat(Tp), Toeplitz=FALSE))
  sumvals <- data.frame(decay=1-r,
                        SWmin = min(vals$SW), SWavg = mean(vals$SW), SWmax = max(vals$SW),
                        crxomin = min(vals$crxo), crxoavg = mean(vals$crxo), crxomax = max(vals$crxo),
                        pllelmin = min(vals$pllel), pllelavg = mean(vals$pllel), pllelmax = max(vals$pllel),
                        pllelbasemin = min(vals$pllelbase), pllelbaseavg = mean(vals$pllelbase),
                        pllelbasemax = max(vals$pllelbase))
  return(sumvals)
}

rs <- seq(0.5, 1, 0.01)
simvalslist <- llply(rs, get_variance_sim, Tp, m, rho0, nsims=100)
simvals <- do.call("rbind", simvalslist)

# Convert for plotting
simvals_long <- simvals %>%
                  gather(col, variance, -decay) %>%
                  separate(col, c("design", "measure"), sep=-4) %>%
                  spread(measure, variance)

pctsim <- ggplot(data=simvals_long, aes(x=decay, group=design, colour=design)) +
  geom_line(aes(y=avg), size=1.0) +
  geom_ribbon(aes(ymin=min, ymax=max, fill=design), alpha=0.3, show.legend=FALSE) +
  expand_limits(y=0) +
  xlab("Decay (1-r)") +
  ylab("Variance") +
  labs(title="Variances of treatment effect, continuous time, simulation",
       subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
  scale_color_hue(labels = c("CRXO", "Parallel", "ParallelBase", "SW")) +
  guides(colour=guide_legend("Design")) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=20),
        plot.subtitle=element_text(hjust=0.5, size=18),
        axis.title=element_text(size=16), axis.text=element_text(size=16),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        legend.position="bottom")
ggsave(paste0("plots/conts_sim_T", Tp, "_m", m, ".pdf"), pctsim, width=297, height=210, units="mm")

# Compare to results with uniform arrival time assumption
load("plots/vars_T4_m50.Rda"); vars_T4_m50 <- varvals
ctvarvals_long <- vars_T4_m50 %>%
              select(decay, starts_with('ct')) %>%
              gather(design, variances, ctSW:ctpllelbase, convert=TRUE) %>%
              separate(design, c("model", "design"), sep=2) %>%
              select(-model)
varvals <- merge(simvals_long, ctvarvals_long)

pctall <- ggplot(data=varvals, aes(x=decay, group=design, colour=design)) +
  geom_line(aes(y=avg), size=1.0) +
  geom_line(aes(y=variances), size=1.0, colour='black') +
  geom_ribbon(aes(ymin=min, ymax=max, fill=design), alpha=0.3, show.legend=FALSE) +
  expand_limits(y=0) +
  xlab("Decay (1-r)") +
  ylab("Variance") +
  labs(title="Variances of treatment effect, continuous time",
       subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
  scale_color_hue(labels = c("CRXO", "Parallel", "ParallelBase", "SW")) +
  guides(colour=guide_legend("Design")) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=20),
        plot.subtitle=element_text(hjust=0.5, size=18),
        axis.title=element_text(size=16), axis.text=element_text(size=16),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        legend.position="bottom")
ggsave(paste0("plots/conts_sim_compare_T", Tp, "_m", m, ".pdf"), pctall, width=297, height=210, units="mm")

