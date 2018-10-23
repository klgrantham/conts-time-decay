# Generate figures for continuous-time correlation decay paper
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)

source('vartheta_twolevels.R')

# Load results for different trial configurations
load("plots/vars_T4_N6_m50_rhoCD023_rhoUC.Rda"); vars_T4_m50 <- varvals
load("plots/vars_T8_N14_m50_rhoCD023_rhoUC.Rda"); vars_T8_m50 <- varvals
load("plots/vars_T4_N6_m150_rhoCD023_rhoUC.Rda"); vars_T4_m150 <- varvals
load("plots/vars_T8_N14_m150_rhoCD023_rhoUC.Rda"); vars_T8_m150 <- varvals

load("plots/vars_ct_mean_T4_N6_m50_rho023.Rda"); vars_ct_mean_T4_m50 <- varvals
load("plots/vars_ct_mean_T8_N14_m50_rho023.Rda"); vars_ct_mean_T8_m50 <- varvals
load("plots/vars_ct_mean_T4_N6_m150_rho023.Rda"); vars_ct_mean_T4_m150 <- varvals
load("plots/vars_ct_mean_T8_N14_m150_rho023.Rda"); vars_ct_mean_T8_m150 <- varvals

load("plots/vars_T4_N6_m50_rhoCD05_rhoUC.Rda"); vars_T4_m50_rho05 <- varvals
load("plots/vars_T4_N6_m10_rhoCD01_rhoUC.Rda"); vars_T4_m10_rho01 <- varvals

load("plots/vars_T4_N6_m50_rhoCD022_rhoUC.Rda"); vars_T4_m50_ct <- varvals
load("plots/vars_T8_N14_m50_rhoCD022_rhoUC.Rda"); vars_T8_m50_ct <- varvals
load("plots/vars_T4_N6_m150_rhoCD022_rhoUC.Rda"); vars_T4_m150_ct <- varvals
load("plots/vars_T8_N14_m150_rhoCD022_rhoUC.Rda"); vars_T8_m150_ct <- varvals


# Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Multiplot framework
make_2x2_multiplot <- function(p1, p2, p3, p4, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                p3 + theme(legend.position="none"),
                                p4 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

make_1x2_multiplot <- function(p1, p2, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

# General plotting function using long format data frame
compare_designs <- function(df.long, ylabel, ylimits, Tp, m){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value, colour=Design, linetype=Design)) +
    geom_line(size=1.2) +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                       labels = c("CRXO", "Parallel", "SW")) +
    scale_linetype_manual(values = c("twodash", "dashed", "solid"),
                          labels = c("CRXO", "Parallel", "SW")) +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab(ylabel) +
    labs(title=bquote(paste(.(Tp), " periods, ", .(m), " subjects/cluster-period"))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
  return(p)
}

# Convert continuous-time results to long format
long_ct <- function(df){
  ctvarvals <- df %>%
    select(decay, starts_with('ct'))
  ctvarvals_long <- gather(data=ctvarvals, key=Design, value=Variance,
                           -decay, convert=TRUE)
  return(ctvarvals_long)
}

# Convert uniform correlation results to long format
long_HH <- function(df){
  HHvarvals <- df %>%
    select(decay, starts_with('HH'))
  HHvarvals_long <- gather(data=HHvarvals, key=Design, value=Variance,
                           -decay, convert=TRUE)
  return(HHvarvals_long)
}

# Convert discrete-time results to long format
long_dt <- function(df){
  dtvarvals <- df %>%
    select(decay, starts_with('dt'))
  dtvarvals_long <- gather(data=dtvarvals, key=Design, value=Variance,
                           -decay, convert=TRUE)
  return(dtvarvals_long)
}

# Convert HH vs continuous-time results to long format
long_rel_HH_ct <- function(df){
  ctvHHvarvals <- df %>%
    mutate(ratioSW=HHSW/ctSW, ratiocrxo=HHcrxo/ctcrxo,
           ratiopllel=HHpllel/ctpllel) %>%
    select(decay, starts_with('ratio'))
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=Design, value=relative_variance,
                              ratioSW:ratiopllel, convert=TRUE)
  return(ctvHHvarvals_long)
}

# Convert discrete-time vs continuous-time results to long format
long_rel_dt_ct <- function(df){
  ctvdtvarvals <- df %>%
    mutate(ratioSW=dtSW/ctSW, ratiocrxo=dtcrxo/ctcrxo,
           ratiopllel=dtpllel/ctpllel) %>%
    select(decay, starts_with('ratio'))
  ctvdtvarvals_long <- gather(data=ctvdtvarvals, key=Design, value=relative_variance,
                              ratioSW:ratiopllel, convert=TRUE)
  return(ctvdtvarvals_long)
}

# Calculate power for a set of variances and a given effect size and sig level
powdf <- function(df, effsize, siglevel=0.05){
  powvals <- apply(df[-df$decay], MARGIN=2, pow, effsize, siglevel)
  powdf <- data.frame(decay=df$decay, powvals)
  return(powdf)
}


# Plot variances, continuous-time, all designs
ylims <- c(0.0,0.02)
p1 <- compare_designs(df.long=long_ct(vars_T4_m50),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=50)
p2 <- compare_designs(df.long=long_ct(vars_T4_m150),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=150)
p3 <- compare_designs(df.long=long_ct(vars_T8_m50),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=50)
p4 <- compare_designs(df.long=long_ct(vars_T8_m150),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=150)
mylegend <- g_legend(p1)
title <- expression(paste("Variance of treatment effect estimator, ", Var(hat(theta)[CCD])))
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
ggsave("plots/conts_T4_8_m50_150_rho023.jpg", p1to4, width=9, height=7, units="in", dpi=600)
ggsave("plots/conts_T4_8_m50_150_rho023.pdf", p1to4, width=9, height=7, units="in", dpi=600)

# Plot relative variance, HH vs continuous, all designs
ylims <- c(0.2,5.0)
p1 <- compare_designs(df.long=long_rel_HH_ct(vars_T4_m50),
                      ylabel="Relative variance", ylimits=ylims, Tp=4, m=50) +
      geom_hline(aes(yintercept=1)) + scale_y_log10(breaks=c(0.2,0.5,1.0,2.0,5.0))
p2 <- compare_designs(df.long=long_rel_HH_ct(vars_T4_m150),
                      ylabel="Relative variance", ylimits=ylims, Tp=4, m=150) +
      geom_hline(aes(yintercept=1)) + scale_y_log10(breaks=c(0.2,0.5,1.0,2.0,5.0))
p3 <- compare_designs(df.long=long_rel_HH_ct(vars_T8_m50),
                      ylabel="Relative variance", ylimits=ylims, Tp=8, m=50) +
      geom_hline(aes(yintercept=1)) + scale_y_log10(breaks=c(0.2,0.5,1.0,2.0,5.0))
p4 <- compare_designs(df.long=long_rel_HH_ct(vars_T8_m150),
                      ylabel="Relative variance", ylimits=ylims, Tp=8, m=150) +
      geom_hline(aes(yintercept=1)) + scale_y_log10(breaks=c(0.2,0.5,1.0,2.0,5.0))
mylegend <- g_legend(p1)
title <- expression(paste("Relative variance of treatment effect estimators, ",
                          Var(hat(theta)[UC])/Var(hat(theta)[CCD])))
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
ggsave("plots/HH_vs_conts_50_150_rhoCD023_rhoUC.jpg", p1to4, width=9, height=7, units="in", dpi=600)
ggsave("plots/HH_vs_conts_50_150_rhoCD023_rhoUC.pdf", p1to4, width=9, height=7, units="in", dpi=600)

# Plot relative variance, discrete vs continuous, all designs
ylims <- c(0.8,2.0)
p1 <- compare_designs(df.long=long_rel_dt_ct(vars_T4_m50),
                      ylabel="Relative variance", ylimits=ylims, Tp=4, m=50) +
      geom_hline(aes(yintercept=1)) + scale_y_continuous(breaks=c(0.8,1.2,1.6,2.0))
p2 <- compare_designs(df.long=long_rel_dt_ct(vars_T4_m150),
                      ylabel="Relative variance", ylimits=ylims, Tp=4, m=150) +
      geom_hline(aes(yintercept=1)) + scale_y_continuous(breaks=c(0.8,1.2,1.6,2.0))
p3 <- compare_designs(df.long=long_rel_dt_ct(vars_T8_m50),
                      ylabel="Relative variance", ylimits=ylims, Tp=8, m=50) +
      geom_hline(aes(yintercept=1)) + scale_y_continuous(breaks=c(0.8,1.2,1.6,2.0))
p4 <- compare_designs(df.long=long_rel_dt_ct(vars_T8_m150),
                      ylabel="Relative variance", ylimits=ylims, Tp=8, m=150) +
      geom_hline(aes(yintercept=1)) + scale_y_continuous(breaks=c(0.8,1.2,1.6,2.0))
mylegend <- g_legend(p1)
title <- expression(paste("Relative variance of treatment effect estimators, ",
                          Var(hat(theta)[DCD])/Var(hat(theta)[CCD])))
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
ggsave("plots/dt_vs_conts_50_150_rho023.jpg", p1to4, width=9, height=7, units="in", dpi=600)
ggsave("plots/dt_vs_conts_50_150_rho023.pdf", p1to4, width=9, height=7, units="in", dpi=600)

# Plot relative variance, continuous, mean vs individual level, all designs
vars_ind_mean_T4_m50 <- data.frame(long_ct(vars_T4_m50),
                            Variance_mean = long_ct(vars_ct_mean_T4_m50)$Variance)
vars_ind_mean_T4_m150 <- data.frame(long_ct(vars_T4_m150),
                             Variance_mean = long_ct(vars_ct_mean_T4_m150)$Variance)
vars_ind_mean_T8_m50 <- data.frame(long_ct(vars_T8_m50),
                             Variance_mean = long_ct(vars_ct_mean_T8_m50)$Variance)
vars_ind_mean_T8_m150 <- data.frame(long_ct(vars_T8_m150),
                             Variance_mean = long_ct(vars_ct_mean_T8_m150)$Variance)

vars_meanvsind_ratios_T4_m50 <- vars_ind_mean_T4_m50 %>%
  mutate(var_ratio = Variance_mean/Variance) %>%
  select(decay, Design, var_ratio)
vars_meanvsind_ratios_T4_m150 <- vars_ind_mean_T4_m150 %>%
  mutate(var_ratio = Variance_mean/Variance) %>%
  select(decay, Design, var_ratio)
vars_meanvsind_ratios_T8_m50 <- vars_ind_mean_T8_m50 %>%
  mutate(var_ratio = Variance_mean/Variance) %>%
  select(decay, Design, var_ratio)
vars_meanvsind_ratios_T8_m150 <- vars_ind_mean_T8_m150 %>%
  mutate(var_ratio = Variance_mean/Variance) %>%
  select(decay, Design, var_ratio)

p1 <- compare_designs(vars_meanvsind_ratios_T4_m50, ylabel="Relative variance",
                     ylimits=c(0.95,1.1), Tp=4, m=50) +
      geom_hline(aes(yintercept=1))
p2 <- compare_designs(vars_meanvsind_ratios_T4_m150, ylabel="Relative variance",
                      ylimits=c(0.95,1.1), Tp=4, m=150) +
      geom_hline(aes(yintercept=1))
p3 <- compare_designs(vars_meanvsind_ratios_T8_m50, ylabel="Relative variance",
                      ylimits=c(0.95,1.1), Tp=8, m=50) +
      geom_hline(aes(yintercept=1))
p4 <- compare_designs(vars_meanvsind_ratios_T8_m150, ylabel="Relative variance",
                      ylimits=c(0.95,1.1), Tp=8, m=150) +
      geom_hline(aes(yintercept=1))
mylegend <- g_legend(p1)
title <- expression(paste("Relative variance of treatment effect estimators, ",
                          Var(hat(theta)[CCD][mean])/Var(hat(theta)[CCD][ind])))
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
ggsave("plots/conts_meanvsind_ratio_50_150_rho023.jpg", p1to4, width=9, height=7, units="in", dpi=600)
ggsave("plots/conts_meanvsind_ratio_50_150_rho023.pdf", p1to4, width=9, height=7, units="in", dpi=600)

# Plot variance and relative variance, continuous, different rho0 values

compare_designs_1by2 <- function(df.long, ylabel, ylimits, title){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value, colour=Design, linetype=Design)) +
    geom_line(size=1.2) +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                       labels = c("CRXO", "Parallel", "SW")) +
    scale_linetype_manual(values = c("twodash", "dashed", "solid"),
                          labels = c("CRXO", "Parallel", "SW")) +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab(ylabel) +
    labs(title=title) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
  return(p)
}

# Tp=4, m=50, rho0=0.05
p1title <- expression(paste("Variance of treatment effect estimator, ", Var(hat(theta)[CCD])))
p1 <- compare_designs_1by2(df.long=long_ct(vars_T4_m50_rho05),
                           ylabel="Variance", ylimits=c(0.0,0.04), title=p1title)
p2title <- expression(paste("Relative variance, ", Var(hat(theta)[UC])/Var(hat(theta)[CCD])))
p2 <- compare_designs_1by2(df.long=long_rel_HH_ct(vars_T4_m50_rho05),
                           ylabel="Relative variance", ylimits=c(0.5,2.0), title=p2title) +
      geom_hline(aes(yintercept=1)) + scale_y_log10(breaks=c(0.5,1.0,2.0))
mylegend <- g_legend(p1)
title <- expression(paste("4 periods, 50 subjects/cluster-period, ", rho, "=0.05"))
p1to2 <- make_1x2_multiplot(p1, p2, mylegend, title=title)
ggsave("plots/vars_T4_m50_rho05.jpg", p1to2, width=9, height=4, units="in", dpi=600)
ggsave("plots/vars_T4_m50_rho05.pdf", p1to2, width=9, height=4, units="in", dpi=600)

# Tp=4, m=10, rho0=0.01
p1title <- expression(paste("Variance of treatment effect estimator, ", Var(hat(theta)[CCD])))
p1 <- compare_designs_1by2(df.long=long_ct(vars_T4_m10_rho01),
                           ylabel="Variance", ylimits=c(0.0,0.05), title=p1title)
p2title <- expression(paste("Relative variance, ", Var(hat(theta)[UC])/Var(hat(theta)[CCD])))
p2 <- compare_designs_1by2(df.long=long_rel_HH_ct(vars_T4_m10_rho01),
                           ylabel="Relative variance", ylimits=c(0.9,1.1), title=p2title) +
  geom_hline(aes(yintercept=1)) + scale_y_continuous(breaks=c(0.9,1.0,1.1))
mylegend <- g_legend(p1)
title <- expression(paste("4 periods, 10 subjects/cluster-period, ", rho, "=0.01"))
p1to2 <- make_1x2_multiplot(p1, p2, mylegend, title=title)
ggsave("plots/vars_T4_m10_rho01.jpg", p1to2, width=9, height=4, units="in", dpi=600)
ggsave("plots/vars_T4_m10_rho01.pdf", p1to2, width=9, height=4, units="in", dpi=600)


# Plot variances, continuous-time, continuous time parameterisation, all designs
ylims <- c(0.0,0.02)
p1 <- compare_designs(df.long=long_ct(vars_T4_m50_ct),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=50)
p2 <- compare_designs(df.long=long_ct(vars_T4_m150_ct),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=150)
p3 <- compare_designs(df.long=long_ct(vars_T8_m50_ct),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=50)
p4 <- compare_designs(df.long=long_ct(vars_T8_m150_ct),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=150)
mylegend <- g_legend(p1)
title <- expression(paste("Variance of treatment effect estimator, ", Var(hat(theta)[CCD])))
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
ggsave("plots/conts_T4_8_m50_150_rho022_ct.jpg", p1to4, width=9, height=7, units="in", dpi=600)
ggsave("plots/conts_T4_8_m50_150_rho022_ct.pdf", p1to4, width=9, height=7, units="in", dpi=600)

# Power plots
effsize <- 0.2
siglevel <- 0.05
pow_T4_m50 <- powdf(vars_T4_m50, effsize, siglevel)
pow_T4_m150 <- powdf(vars_T4_m150, effsize, siglevel)
pow_T8_m50 <- powdf(vars_T8_m50, effsize, siglevel)
pow_T8_m150 <- powdf(vars_T8_m150, effsize, siglevel)

# Power, continuous-time, all designs
ylims <- c(0,1)
p1 <- compare_designs(df.long=long_ct(pow_T4_m50),
                      ylabel="Power", ylimits=ylims, Tp=4, m=50) + geom_hline(aes(yintercept=1))
p2 <- compare_designs(df.long=long_ct(pow_T4_m150),
                      ylabel="Power", ylimits=ylims, Tp=4, m=150) + geom_hline(aes(yintercept=1))
p3 <- compare_designs(df.long=long_ct(pow_T8_m50),
                      ylabel="Power", ylimits=c(0.5,1), Tp=8, m=50) + geom_hline(aes(yintercept=1))
p4 <- compare_designs(df.long=long_ct(pow_T8_m150),
                      ylabel="Power", ylimits=c(0.5,1), Tp=8, m=150) + geom_hline(aes(yintercept=1))
mylegend <- g_legend(p1)
title <- paste0('Power to detect effect size of ', effsize, ', CCD')
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
ggsave("plots/power_conts_T4_8_m50_150_rho023.jpg", p1to4, width=9, height=7, units="in", dpi=600)
ggsave("plots/power_conts_T4_8_m50_150_rho023.pdf", p1to4, width=9, height=7, units="in", dpi=600)
