# Comparisons of properties of continuous time model with
# non-uniform correlation
# 
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)

source('vartheta_twolevels.R')

comparison_plots <- function(Tp, m, rho0){
  # Calculates the variance of the treatment effect under the models:
  #    - continuous time (ct)
  #    - discrete time (dt)
  #    - Hussey & Hughes (HH)
  # with trial designs:
  #    - stepped wedge (SW)
  #    - cluster randomised crossover (CRXO)
  #    - parallel (pllel)
  #    - parallel with baseline (pllelbase)
  # and saves a set of comparison plots
  #
  # Inputs:
  # Tp - number of time periods in the trial
  # m - number of subjects measured in each time period
  # rho0 - base correlation between a pair of subjects
  #
  # Example usage: vals <- comparison_plots(Tp=4, m=50, rho0=0.035)
  
  rs <- seq(0.5, 1, 0.01) # Note: Could only calculate from 0.5 to 1
  # Specify the variance matrices under the different models
  ctvarmat <- llply(rs, expdecayVicont, Tp, m, rho0, meanlvl=FALSE)
  dtvarmat <- llply(rs, expdecayVi, Tp, m, rho0, meanlvl=TRUE)
  HHvarmat <- HHVi(Tp, m, rho0, meanlvl=TRUE)
  
  # Get the variances of the treatment effect under the
  # different models and designs
  # Note: Still need to expand the Xmats according to nclust
  # Scale the non-SW variances by (Tp/(Tp-1)) to account for uneven clusters across designs
  scalefactor <- Tp/(Tp-1)
  varvals <- data.frame(decay = 1-rs,
              ctSW = laply(ctvarmat, vartheta_ind, Xmat=SWdesmat(Tp)),
              ctcrxo = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=crxodesmat(Tp)),
              ctpllel = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=plleldesmat(Tp)),
              ctpllelbase = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=pllelbasedesmat(Tp)),
              dtSW = laply(dtvarmat, vartheta_mean, Xmat=SWdesmat(Tp)),
              dtcrxo = scalefactor*laply(dtvarmat, vartheta_mean, Xmat=crxodesmat(Tp)),
              dtpllel = scalefactor*laply(dtvarmat, vartheta_mean, Xmat=plleldesmat(Tp)),
              dtpllelbase = scalefactor*laply(dtvarmat, vartheta_mean, Xmat=pllelbasedesmat(Tp)),
              HHSW = vartheta_mean(Vi=HHvarmat, Xmat=SWdesmat(Tp)),
              HHcrxo = scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=crxodesmat(Tp)),
              HHpllel = scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=plleldesmat(Tp)),
              HHpllelbase = scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=pllelbasedesmat(Tp)))
  
  # Variance, continuous time, all designs
  ctvarvals <- varvals %>%
    select(decay, starts_with('ct')) %>%
    filter(decay<=0.5)
  ctvarvals_long <- gather(data=ctvarvals, key=design, value=variances,
                           ctSW:ctpllelbase, convert=TRUE)
  
  pctvar <- ggplot(data=ctvarvals_long, aes(x=decay, y=variances, group=design, colour=design)) +
    geom_line(size=1.0) +
    expand_limits(y=0) +
    xlab("Decay (1-r)") +
    ylab("Variance") +
    labs(title="Variances of treatment effect, continuous time",
         subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.title=element_text(size=16))
  
  # Relative variance, continuous vs discrete time, all designs
  ctvdtvarvals <- varvals %>%
    mutate(ratioSW=ctSW/dtSW, ratiocrxo=ctcrxo/dtcrxo,
           ratiopllel=ctpllel/dtpllel, ratiopllelbase=ctpllelbase/dtpllelbase) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvdtvarvals_long <- gather(data=ctvdtvarvals, key=design, value=relative_variance,
                              ratioSW:ratiopllelbase, convert=TRUE)
  
  pctvdtvar <- ggplot(data=ctvdtvarvals_long, aes(x=decay, y=relative_variance, group=design, colour=design)) +
    geom_line(size=1.0) +
    xlab("Decay (1-r)") +
    ylab("Relative variance") +
    labs(title="Relative variances of treatment effect, conts vs discrete",
         subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.title=element_text(size=16))
  
  # Relative variance, continuous time vs HH, all designs
  ctvHHvarvals <- varvals %>%
    mutate(ratioSW=ctSW/HHSW, ratiocrxo=ctcrxo/HHcrxo,
           ratiopllel=ctpllel/HHpllel, ratiopllelbase=ctpllelbase/HHpllelbase) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=design, value=relative_variance,
                              ratioSW:ratiopllelbase, convert=TRUE)
  
  pctvHHvar <- ggplot(data=ctvHHvarvals_long, aes(x=decay, y=relative_variance, group=design, colour=design)) +
    geom_line(size=1.0) +
#    expand_limits(y=1.5) +
    xlab("Decay (1-r)") +
    ylab("Relative variance") +
    labs(title="Relative variances of treatment effect, conts vs HH",
         subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.title=element_text(size=16))
  
  # Variance, continuous and discrete time, each design separately
  ctdtvarvals <- varvals %>%
    select(decay, starts_with('ct'), starts_with('dt')) %>%
    gather(col, variance, -decay) %>%
    separate(col, c("model", "design"), sep=2) %>%
    filter(decay<=0.5)
  
  pgen <- function(df, designtype, title){
    p <- ggplot(subset(df, design==designtype), aes(x=decay, y=variance, colour=model)) +
                geom_line(size=1.0) +
                xlab("Decay (1-r)") +
                ylab("Variance") +
                ggtitle(title) +
                theme_bw() +
                theme(plot.title=element_text(hjust=0.5))
    return(p)
  }
  p1 <- pgen(ctdtvarvals, "SW", "SW design")
  p2 <- pgen(ctdtvarvals, "crxo", "CRXO design")
  p3 <- pgen(ctdtvarvals, "pllel", "Parallel design")
  p4 <- pgen(ctdtvarvals, "pllelbase", "Parallel w/ baseline design")
  pctdtvarmult <- grid.arrange(p1, p2, p3, p4, ncol=2)
  
  # Save plots to pdf
  ggsave(paste0("plots/conts_T", Tp, "_m", m, ".pdf"), pctvar) # TODO: Include rho0 values in filenames
  ggsave(paste0("plots/conts_vs_disc_T", Tp, "_m", m, ".pdf"), pctvdtvar)
  ggsave(paste0("plots/conts_vs_HH_T", Tp, "_m", m, ".pdf"), pctvHHvar)
  ggsave(paste0("plots/conts_disc_multi_T", Tp, "_m", m, ".pdf"), pctdtvarmult)
  
  save(varvals, file=paste0("plots/vars_T", Tp, "_m", m, ".Rda"))
  
  return(varvals)
}
