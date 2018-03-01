# Generate figures for continuous time decay paper
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
load("plots/vars_T4_m50_rho04.Rda"); vars_T4_m50 <- varvals
load("plots/vars_T8_m50_rho04.Rda"); vars_T8_m50 <- varvals
load("plots/vars_T4_m500_rho04.Rda"); vars_T4_m500 <- varvals
load("plots/vars_T8_m500_rho04.Rda"); vars_T8_m500 <- varvals

load("plots/vars_ct_mean_T4_m50_rho04.Rda"); vars_ct_mean_T4_m50 <- varvals
load("plots/vars_ct_mean_T8_m50_rho04.Rda"); vars_ct_mean_T8_m50 <- varvals
load("plots/vars_ct_mean_T4_m500_rho04.Rda"); vars_ct_mean_T4_m500 <- varvals
load("plots/vars_ct_mean_T8_m500_rho04.Rda"); vars_ct_mean_T8_m500 <- varvals


# Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Perhaps better to use facetwrap instead?
# Would have to merge all results and have a variable for the config
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

# Convert continuous time results to long format
long_ct <- function(df){
  ctvarvals <- df %>%
    select(decay, starts_with('ct'), -ends_with('base')) %>%
    filter(decay<=0.5)
  ctvarvals_long <- gather(data=ctvarvals, key=Design, value=Variance,
                           -decay, convert=TRUE)
  return(ctvarvals_long)
}

# Convert uniform results to long format
long_HH <- function(df){
  HHvarvals <- df %>%
    select(decay, starts_with('HH'), -ends_with('base')) %>%
    filter(decay<=0.5)
  HHvarvals_long <- gather(data=HHvarvals, key=Design, value=Variance,
                           -decay, convert=TRUE)
  return(HHvarvals_long)
}

# Convert discrete time results to long format
long_dt <- function(df){
  dtvarvals <- df %>%
    select(decay, starts_with('dt'), -ends_with('base')) %>%
    filter(decay<=0.5)
  dtvarvals_long <- gather(data=dtvarvals, key=Design, value=Variance,
                           -decay, convert=TRUE)
  return(dtvarvals_long)
}

# Convert HH vs continuous time results to long format
long_rel_HH_ct <- function(df){
  ctvHHvarvals <- df %>%
    mutate(ratioSW=HHSW/ctSW, ratiocrxo=HHcrxo/ctcrxo,
           ratiopllel=HHpllel/ctpllel) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=Design, value=relative_variance,
                              ratioSW:ratiopllel, convert=TRUE)
  return(ctvHHvarvals_long)
}

# Convert discrete time vs continuous time results to long format
long_rel_dt_ct <- function(df){
  ctvdtvarvals <- df %>%
    mutate(ratioSW=dtSW/ctSW, ratiocrxo=dtcrxo/ctcrxo,
           ratiopllel=dtpllel/ctpllel) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvdtvarvals_long <- gather(data=ctvdtvarvals, key=Design, value=relative_variance,
                              ratioSW:ratiopllel, convert=TRUE)
  return(ctvdtvarvals_long)
}


# Plot variances, continuous time, all designs
ylims <- c(0.0,0.06)
p1 <- compare_designs(df.long=long_ct(vars_T4_m50),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=50)
p2 <- compare_designs(df.long=long_ct(vars_T4_m500),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=500)
p3 <- compare_designs(df.long=long_ct(vars_T8_m50),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=50)
p4 <- compare_designs(df.long=long_ct(vars_T8_m500),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=500)
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Variance of treatment effect estimator, continuous-time decay")
ggsave(paste0("plots/conts_T4_8_m50_500_rho04.jpg"), p1to4, width=9, height=7, units="in", dpi=600)

# Plot relative variance, HH vs continuous, all designs
ylims <- c(0.0,4.0)
p1 <- compare_designs(df.long=long_rel_HH_ct(vars_T4_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=50) +
      geom_hline(aes(yintercept=1))
p2 <- compare_designs(df.long=long_rel_HH_ct(vars_T4_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=500) +
      geom_hline(aes(yintercept=1))
p3 <- compare_designs(df.long=long_rel_HH_ct(vars_T8_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=50) +
      geom_hline(aes(yintercept=1))
p4 <- compare_designs(df.long=long_rel_HH_ct(vars_T8_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=500) +
      geom_hline(aes(yintercept=1))
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Relative variance, uniform vs continuous-time decay")
ggsave(paste0("plots/HH_vs_conts_50_500_rho04.jpg"), p1to4, width=9, height=7, units="in", dpi=600)

# Plot relative variance, discrete vs continuous, all designs
ylims <- c(0.5,3.0)
p1 <- compare_designs(df.long=long_rel_dt_ct(vars_T4_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=50) +
      geom_hline(aes(yintercept=1))
p2 <- compare_designs(df.long=long_rel_dt_ct(vars_T4_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=500) +
      geom_hline(aes(yintercept=1))
p3 <- compare_designs(df.long=long_rel_dt_ct(vars_T8_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=50) +
      geom_hline(aes(yintercept=1))
p4 <- compare_designs(df.long=long_rel_dt_ct(vars_T8_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=500) +
      geom_hline(aes(yintercept=1))
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Relative variance, discrete- vs continuous-time decay")
ggsave(paste0("plots/dt_vs_conts_50_500_rho04.jpg"), p1to4, width=9, height=7, units="in", dpi=600)


# Plot relative variance, continuous, mean vs individual level, all designs
vars_ind_mean_T4_m50 <- data.frame(long_ct(vars_T4_m50),
                            Variance_mean = long_ct(vars_ct_mean_T4_m50)$Variance)
vars_ind_mean_T4_m500 <- data.frame(long_ct(vars_T4_m500),
                             Variance_mean = long_ct(vars_ct_mean_T4_m500)$Variance)
vars_ind_mean_T8_m50 <- data.frame(long_ct(vars_T8_m50),
                             Variance_mean = long_ct(vars_ct_mean_T8_m50)$Variance)
vars_ind_mean_T8_m500 <- data.frame(long_ct(vars_T8_m500),
                             Variance_mean = long_ct(vars_ct_mean_T8_m500)$Variance)

vars_indvsmean_ratios_T4_m50 <- vars_ind_mean_T4_m50 %>%
  mutate(var_ratio = Variance/Variance_mean) %>%
  select(decay, Design, var_ratio)
vars_indvsmean_ratios_T4_m500 <- vars_ind_mean_T4_m500 %>%
  mutate(var_ratio = Variance/Variance_mean) %>%
  select(decay, Design, var_ratio)
vars_indvsmean_ratios_T8_m50 <- vars_ind_mean_T8_m50 %>%
  mutate(var_ratio = Variance/Variance_mean) %>%
  select(decay, Design, var_ratio)
vars_indvsmean_ratios_T8_m500 <- vars_ind_mean_T8_m500 %>%
  mutate(var_ratio = Variance/Variance_mean) %>%
  select(decay, Design, var_ratio)

p1 <- compare_designs(vars_indvsmean_ratios_T4_m50, ylabel="Relative variance",
                     ylimits=c(0.5,1.0), Tp=4, m=50) +
      geom_hline(aes(yintercept=1))
p2 <- compare_designs(vars_indvsmean_ratios_T4_m500, ylabel="Relative variance",
                      ylimits=c(0.5,1.0), Tp=4, m=500) +
      geom_hline(aes(yintercept=1))
p3 <- compare_designs(vars_indvsmean_ratios_T8_m50, ylabel="Relative variance",
                      ylimits=c(0.5,1.0), Tp=8, m=50) +
      geom_hline(aes(yintercept=1))
p4 <- compare_designs(vars_indvsmean_ratios_T8_m500, ylabel="Relative variance",
                      ylimits=c(0.5,1.0), Tp=8, m=500) +
      geom_hline(aes(yintercept=1))
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Relative variance, continuous-time decay, mean vs individual level")
ggsave(paste0("plots/conts_meanvsind_ratio_50_500_rho04.jpg"), p1to4, width=9, height=7, units="in", dpi=600)
