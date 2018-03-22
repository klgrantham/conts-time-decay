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
  # Calculates the variance of the treatment effect estimator under the model:
  #    - continuous time (ct)
  # using simulated arrival times, with trial designs:
  #    - stepped wedge (SW)
  #    - cluster randomised crossover (CRXO)
  #    - parallel (pllel)
  # and returns the median and lower and upper bound of the 95% CI of the
  # observed variances for each design
  
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
  simvals <- data.frame(decay=1-rs,
                        SW.lower = apply(SWvals, 1, quantile, probs=0.025),
                        SW.median = apply(SWvals, 1, median),
                        SW.upper = apply(SWvals, 1, quantile, probs=0.975),
                        crxo.lower = apply(crxovals, 1, quantile, probs=0.025),
                        crxo.median = apply(crxovals, 1, median),
                        crxo.upper = apply(crxovals, 1, quantile, probs=0.975),
                        pllel.lower = apply(pllelvals, 1, quantile, probs=0.025),
                        pllel.median = apply(pllelvals, 1, median),
                        pllel.upper = apply(pllelvals, 1, quantile, probs=0.975))
  return(simvals)
}

##
Tp <- 4
m <- 50
rho0 <- 0.04
rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2]

simvals <- sim_results(nsims=1000, rs=seq(0.5, 1, 0.01), Tp=Tp, m=m, rho0=rho0, type="uniform")
save(simvals, file=paste0("plots/vars_nsims_1000_summary_uniform_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))

simvalsexp <- sim_results(nsims=1000, rs=seq(0.5, 1, 0.01), Tp=Tp, m=m, rho0=rho0, type="exponential")
save(simvalsexp, file=paste0("plots/vars_nsims_1000_summary_exponential_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))

# Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

make_1x2_multiplot <- function(p1, p2, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

sim_results_plots <- function(simvals, Tp, m, rho0, title){
  
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
  
  # Plot results and compare to results under evenly-spaced assumption
  p <- ggplot(data=varvals, aes(x=decay, colour=Design, linetype=Design)) +
    geom_line(aes(y=median), size=1.0) +
    geom_line(aes(y=Variance), size=0.5, colour="black", show.legend=FALSE) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Design), alpha=0.3, show.legend=FALSE, colour=NA) +
    expand_limits(y=0) +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                       labels = c("CRXO", "Parallel", "SW")) +
    scale_linetype_manual(values = c("twodash", "dashed", "solid"),
                          labels = c("CRXO", "Parallel", "SW")) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    labs(title=title) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=11),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
}

load('plots/vars_nsims_1000_summary_uniform_T4_m50_rho04.Rda'); simvals_unif <- simvals
load('plots/vars_nsims_1000_summary_exponential_T4_m50_rho04.Rda'); simvals_exp <- simvalsexp
p1 <- sim_results_plots(simvals_unif, Tp=4, m=50, rho0=0.04, title="Uniformly-distributed simulated measurement times")
p2 <- sim_results_plots(simvals_exp, Tp=4, m=50, rho0=0.04, title="Exponentially-distributed simulated measurement times")
mylegend <- g_legend(p1)
p1to2 <- make_1x2_multiplot(p1, p2, mylegend,
                            title="Variance of treatment effect estimator, continuous-time correlation decay")
ggsave(paste0("plots/conts_sim_unif_exp_compare_T4_m50_rho04.jpg"), p1to2, width=9, height=4, units="in", dpi=600)
