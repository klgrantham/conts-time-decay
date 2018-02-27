# Generate plots of variance ratios for WIP presentation
# 
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)

source('vartheta_twolevels.R')

load("plots/vars_T4_m500_rho035.Rda"); vars_T4_m500_rho035 <- varvals
load("plots/vars_T8_m500_rho035.Rda"); vars_T8_m500_rho035 <- varvals
load("plots/vars_T4_m500_rho04.Rda"); vars_T4_m500_rho040 <- varvals
load("plots/vars_T4_m500_rho005.Rda"); vars_T4_m500_rho005 <- varvals
load("plots/vars_T4_m50_rho035.Rda"); vars_T4_m50_rho035 <- varvals


# Convert continuous vs discrete time results to long format
long_rel_dt_SW <- function(df){
  ctvdtvarvals <- df %>%
    mutate(ratioSW=ctSW/dtSW) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvdtvarvals_long <- gather(data=ctvdtvarvals, key=design, value=relative_variance,
                              ratioSW, convert=TRUE)
  return(ctvdtvarvals_long)
}

# Convert continuous time vs HH results to long format
long_rel_HH_SW <- function(df){
  ctvHHvarvals <- df %>%
    mutate(ratioSW=ctSW/HHSW) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=design, value=relative_variance,
                              ratioSW, convert=TRUE)
  return(ctvHHvarvals_long)
}

# General plotting function using long format data frame
plot_variances <- function(df.long, ylabel, ylimits, color="#00BFC4"){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value)) +
    geom_line(size=2.0, color=color) +
    geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash") +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab(ylabel) +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20))
  return(p)
}

plot_var_ratio <- function(df.long, ylimits, color="#00BFC4"){
  p <- ggplot(data=df.long, aes(x=decay)) +
    geom_line(aes(y=relative_variance), size=2.0, color=color) +
    geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash") +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab("Variance ratio") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20))
  return(p)
}

# Functions for comparing all designs

plot_var_ratios <- function(df.long, ylimits, color="#00BFC4"){
  p <- ggplot(data=df.long, aes(x=decay, group=design, colour=design)) +
    geom_line(aes(y=relative_variance), size=2.0) +
    geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash") +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab("Variance ratio") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
          legend.position="none")
  return(p)
}

# Convert continuous time results to long format
long_ct <- function(df){
  ctvarvals <- df %>%
    select(decay, starts_with('ct')) %>%
    filter(decay<=0.5)
  ctvarvals_long <- gather(data=ctvarvals, key=design, value=variances,
                           starts_with('ct'), convert=TRUE)
  return(ctvarvals_long)
}

# Convert discrete time results to long format
long_dt <- function(df){
  dtvarvals <- df %>%
    select(decay, starts_with('dt')) %>%
    filter(decay<=0.5)
  dtvarvals_long <- gather(data=dtvarvals, key=design, value=variances,
                           dtSW:dtpllelbase, convert=TRUE)
  return(dtvarvals_long)
}

# Convert continuous vs discrete time results to long format
long_rel_dt <- function(df){
  ctvdtvarvals <- df %>%
    mutate(ratioSW=ctSW/dtSW, ratiocrxo=ctcrxo/dtcrxo,
           ratiopllel=ctpllel/dtpllel, ratiopllelbase=ctpllelbase/dtpllelbase) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvdtvarvals_long <- gather(data=ctvdtvarvals, key=design, value=relative_variance,
                              ratioSW:ratiopllelbase, convert=TRUE)
  return(ctvdtvarvals_long)
}

# Convert continuous time vs HH results to long format
long_rel_HH <- function(df){
  ctvHHvarvals <- df %>%
    mutate(ratioSW=ctSW/HHSW, ratiocrxo=ctcrxo/HHcrxo,
           ratiopllel=ctpllel/HHpllel, ratiopllelbase=ctpllelbase/HHpllelbase) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=design, value=relative_variance,
                              ratioSW:ratiopllelbase, convert=TRUE)
  return(ctvHHvarvals_long)
}

sim_results_plots <- function(simvals, Tp, m, rho0, type="uniform"){
  
  # Convert for plotting
  simvals_long <- simvals %>%
    gather(col, variance, -decay) %>%
    separate(col, c("design", "measure")) %>%
    spread(measure, variance)
  
  # Obtain results under evenly-spaced assumption
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2]
  load(paste0("plots/vars_T", Tp, "_m", m, "_rho", rho0char, ".Rda")); vars_assump <- varvals
  ctvarvals_long <- vars_assump %>%
    select(decay, starts_with('ct')) %>%
    gather(design, variances, ctSW:ctpllelbase, convert=TRUE) %>%
    separate(design, c("model", "design"), sep=2) %>%
    select(-model)
  varvals <- merge(simvals_long, ctvarvals_long)
  
  # Plot results and compare to results under evenly-spaced assumption
  pctall <- ggplot(data=varvals, aes(x=decay, group=design, colour=design)) +
    geom_line(aes(y=avg), size=1.0) +
    geom_line(aes(y=variances), size=1.0, colour='black', linetype="dotted") +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=design), alpha=0.3, show.legend=FALSE) +
    expand_limits(y=0) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
          legend.position="none")
  ggsave(paste0("plots/conts_sim_", type, "_compare_T", Tp, "_m", m, ".png"),
         pctall, width=7.99, height=5.67, units="in", dpi=900)
}

## SW design only

ctvdtT4 <- long_rel_dt_SW(vars_T4_m500_rho035)
ctvHHT4 <- long_rel_HH_SW(vars_T4_m500_rho035)
ctvHHT8 <- long_rel_HH_SW(vars_T8_m500_rho035)

# Compare ct to dt, T=4, m=500, rho0=0.035
plot_var_ratio(ctvdtT4, c(0.0, 1.0), "#C77CFF")
ggsave(paste0("plots/ctvdt_T4_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4, m=500, rho0=0.035
plot_var_ratio(ctvHHT4, c(0.0, 6.2), "#C77CFF")
ggsave(paste0("plots/ctvHH_T4_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4 vs T=8, m=500, rho0=0.035
ctvHHT4T8 <- data.frame(ctvHHT4, relative_variance_T8 = ctvHHT8$relative_variance)
p <- plot_var_ratio(ctvHHT4T8, c(0.0, 6.2), "#C77CFF")
p + geom_line(aes(y=relative_variance_T8), size=2.0, color="#7f00de")
ggsave(paste0("plots/ctvHH_T4_8_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct (rho0=0.04) to HH (rho0=0.035)
vars_T4_m500_rho035_040 <- data.frame(select(vars_T4_m500_rho035, decay, HHSW),
                                      ctSW = vars_T4_m500_rho040$ctSW)
ctvHHT4_rho035_040 <- long_rel_HH_SW(vars_T4_m500_rho035_040)
plot_var_ratio(ctvHHT4_rho035_040, c(0.0, 6.2), "#F8766D")
ggsave(paste0("plots/ctvHH_T4_m500_rho035_040.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4, m=500, rho0=0.035 vs 0.005
ctvHHT4rho005 <- long_rel_HH_SW(vars_T4_m500_rho005)
ctvHHT4rho035_005 <- data.frame(ctvHHT4, relative_variance_rho005 = ctvHHT4rho005$relative_variance)
p <- plot_var_ratio(ctvHHT4rho035_005, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_rho005), size=2.0, color="#00BFC4")
ggsave(paste0("plots/ctvHH_T4_m500_rho035_005.png"), width=7.99, height=5.67, units="in", dpi=900)

# m = 500 vs m = 50
ctvHHT4m50 <- long_rel_HH_SW(vars_T4_m50_rho035)
ctvHHT4m500_50 <- data.frame(ctvHHT4, relative_variance_m50 = ctvHHT4m50$relative_variance)
p <- plot_var_ratio(ctvHHT4m500_50, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_m50), size=2.0, color="#00BFC4")
ggsave(paste0("plots/ctvHH_T4_m500_50_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

## Compare all designs

# Plot relative variance, continuous vs HH, all designs
ylims <- c(0.0,5.0)
p1 <- plot_var_ratios(df.long=long_rel_HH(vars_T4_m500_rho035), ylimits=ylims)
ggsave(paste0("plots/ctvHH_T4_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot relative variance, continuous vs HH, all designs
ylims <- c(0.0,6.3)
p1 <- plot_var_ratios(df.long=long_rel_HH(vars_T4_m500_rho035), ylimits=ylims)
ggsave(paste0("plots/ctvHH_T4_m500_rho035_samescale.png"), width=7.99, height=5.67, units="in", dpi=900)
p2 <- plot_var_ratios(df.long=long_rel_HH(vars_T8_m500_rho035), ylimits=ylims)
ggsave(paste0("plots/ctvHH_T8_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot relative variance, continuous vs HH, all designs
ylims <- c(0.0,1.0)
p3 <- plot_var_ratios(df.long=long_rel_dt(vars_T4_m500_rho035), ylimits=ylims)
ggsave(paste0("plots/ctvdt_T4_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)
p4 <- plot_var_ratios(df.long=long_rel_dt(vars_T8_m500_rho035), ylimits=ylims)
ggsave(paste0("plots/ctvdt_T8_m500_rho035.png"), width=7.99, height=5.67, units="in", dpi=900)

## Simulated arrival times

# Uniform distribution
load("plots/vars_nsims_1000_summary_uniform_T4_m50_rho035.Rda"); simvalsunif <- simvals
sim_results_plots(simvalsunif, Tp=4, m=50, rho0=0.035, type="uniform")

# Exponential distribution
load("plots/vars_nsims_1000_summary_exponential_T4_m50_rho035.Rda"); simvalsexp <- simvals
sim_results_plots(simvalsexp, Tp=4, m=50, rho0=0.035, type="exponential")


## Not included

# Visualize simulated times
# Uniform distribution
Tp <- 4; m <- 50
x <- runif(Tp*m, 0, 1) # Generate vector of fractional times
j <- rep(1:Tp, each=m) # Create vector of time period indices and add to fractional times
uniftimes <- sort(j + x) # Combine to create arrival times

# Exponential distribution
x <- matrix(rexp(Tp*m), m, Tp)
colmaxs <- apply(x, 2, max)
y <- x %*% diag(1/colmaxs)
colmins <- apply(y, 2, min)
z <- t(t(y) - colmins)
j <- matrix(1:Tp, m, Tp, byrow=TRUE)
exptimes <- sort(as.vector(z + j))

stripchart(uniftimes, method="jitter", pch=1)
stripchart(exptimes, method="jitter", pch=1)

hist(uniftimes,main='Simulated measurement times',xlab='Trial period', ylim=c(0,30))
stripchart(uniftimes,add=TRUE,at=-0.5)

hist(uniftimes[1:50],main='Simulated measurement times',xlab='Trial period', ylim=c(0,30))

unifsimtimes <- data.frame(position=1:200, const=1, period=as.factor(rep(1:Tp, each=m)), time=uniftimes)
ggplot(unifsimtimes, aes(x=time, y=period)) +
  geom_jitter(height=0.1)

expsimtimes <- data.frame(const=1, period=as.factor(rep(1:Tp, each=m)), time=exptimes)
ggplot(expsimtimes, aes(x=time, y=const, color=period)) +
  expand_limits(y=c(0,2)) +
  geom_jitter(height=0.05, alpha=0.5) +
  theme_bw()
# + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

ggplot(unifsimtimes, aes(x=time)) + geom_density(aes(fill=period), alpha = 0.4, trim=TRUE)
ggplot(expsimtimes, aes(x=time)) + geom_density(aes(fill=period), alpha = 0.4, trim=TRUE)

ggplot(unifsimtimes, aes(x=time)) + geom_histogram(aes(fill=period), alpha = 0.4, binwidth=0.1)


# Several strip plots of times
simtimesunif <- function(){
  x <- runif(Tp*m, 0, 1) # Generate vector of fractional times
  j <- rep(1:Tp, each=m) # Create vector of time period indices and add to fractional times
  times <- sort(j + x) # Combine to create arrival times
  uniftimes <- data.frame(const=1, period=as.factor(rep(1:Tp, each=m)), time=times)
  return(uniftimes)
}

simtimesexp <- function(){
  x <- matrix(rexp(Tp*m), m, Tp)
  colmaxs <- apply(x, 2, max)
  y <- x %*% diag(1/colmaxs)
  colmins <- apply(y, 2, min)
  z <- t(t(y) - colmins)
  j <- matrix(1:Tp, m, Tp, byrow=TRUE)
  times <- sort(as.vector(z + j))
  exptimes <- data.frame(const=1, period=as.factor(rep(1:Tp, each=m)), time=times)
  return(exptimes)
}

ggplot(simtimesunif(), aes(x=time, y=const)) +
  expand_limits(y=c(0,2)) +
  geom_jitter(size=2.0, color="darkblue", height=0.05, alpha=0.5) +
  theme_void() +
  theme(legend.position = "none")
ggsave(paste0("plots/arrivaltimes_unif.png"), width=7.99, height=5.67, units="in", dpi=900)


ggplot(simtimesexp(), aes(x=time, y=const)) +
  expand_limits(y=c(0,2)) +
  geom_jitter(size=2.0, color="darkblue", height=0.05, alpha=0.5) +
  theme_void() +
  theme(legend.position = "none")
ggsave(paste0("plots/arrivaltimes_exp.png"), width=7.99, height=5.67, units="in", dpi=900)
#  scale_color_brewer(palette="Dark2")

# Evenly spaced times plot
eventimes <- data.frame(const=1, period=as.factor(rep(1:Tp, each=m)), time=seq(1, (Tp+1)-1/m, 1/m))
ggplot(eventimes, aes(x=time, y=const)) +
  expand_limits(y=c(0,2)) +
  geom_jitter(size=2.0, color="darkblue", height=0.05, alpha=0.5) +
  theme_void() +
  theme(legend.position = "none")
ggsave(paste0("plots/arrivaltimes_even.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct (rho0=0.04) to dt (rho0=0.035)
vars_ct_dt_T4_m500_rho035_040 <- data.frame(select(vars_T4_m500_rho035, decay, starts_with("dt")),
                                            select(vars_T4_m500_rho040, starts_with("ct")))
ctvdtT4_rho035_040 <- long_rel_dt(vars_ct_dt_T4_m500_rho035_040)
p <- plot_var_ratios(ctvdtT4_rho035_040, c(0.0, 1.0))
p + labs(title=bquote(paste("Relative variance, continuous (", rho[0]==0.04, ") vs discrete (", rho[0]==0.035, ")")),
         subtitle=bquote(paste(T==4, ", ", m==500))) +
    scale_color_hue(labels = c("CRXO", "Parallel", "Parallel w/ baseline", "SW")) +
    guides(colour=guide_legend("Design")) +
    theme(plot.title=element_text(hjust=0.5, size=20),
        plot.subtitle=element_text(hjust=0.5, size=18),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        legend.position = "bottom")
p + ggtitle("Relative variance, continuous vs discrete time")
ggsave("plots/ctvdt_T4_m500_rho035_040.png", width=7.99, height=5.67, units="in", dpi=900)