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
load("plots/vars_T4_m50.Rda"); vars_T4_m50 <- varvals
load("plots/vars_T8_m50.Rda"); vars_T8_m50 <- varvals
load("plots/vars_T4_m100.Rda"); vars_T4_m100 <- varvals
load("plots/vars_T8_m100.Rda"); vars_T8_m100 <- varvals
load("plots/vars_T4_m200.Rda"); vars_T4_m200 <- varvals
load("plots/vars_T8_m200.Rda"); vars_T8_m200 <- varvals
load("plots/vars_T4_m400.Rda"); vars_T4_m400 <- varvals
load("plots/vars_T8_m400.Rda"); vars_T8_m400 <- varvals
load("plots/vars_T4_m500.Rda"); vars_T4_m500 <- varvals
load("plots/vars_T8_m500.Rda"); vars_T8_m500 <- varvals

# Extract legend
# Source: github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
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
                                 gp=gpar(fontsize=20)))
  return(p)
}

# General plotting function using long format data frame
compare_designs <- function(df.long, ylabel, ylimits, Tp, m, rho0){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value, group=design, colour=design)) +
    geom_line(size=1.0) +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab(ylabel) +
    labs(title=bquote(paste(T==.(Tp), ", ", m==.(m)))) + #, ", ", rho[0]==.(rho0)
    scale_color_hue(labels = c("CRXO", "Parallel", "Parallel w/ baseline", "SW")) +
    guides(colour=guide_legend("Design")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          legend.title=element_text(size=16), legend.text=element_text(size=14),
          legend.position="bottom")
  return(p)
}

# Convert continuous time results to long format
long_ct <- function(df){
  ctvarvals <- df %>%
    select(decay, starts_with('ct')) %>%
    filter(decay<=0.5)
  ctvarvals_long <- gather(data=ctvarvals, key=design, value=variances,
                           ctSW:ctpllelbase, convert=TRUE)
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

# Convert uniform results to long format
long_HH <- function(df){
  HHvarvals <- df %>%
    select(decay, starts_with('HH')) %>%
    filter(decay<=0.5)
  HHvarvals_long <- gather(data=HHvarvals, key=design, value=variances,
                           HHSW:HHpllelbase, convert=TRUE)
  return(HHvarvals_long)
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


rho0 <- 0.035

# Plot variances, continuous time, all designs
ylims <- c(0.0,0.06)
p1 <- compare_designs(df.long=long_ct(vars_T4_m50),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=50, rho0=rho0)
p2 <- compare_designs(df.long=long_ct(vars_T4_m500),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=500, rho0=rho0)
p3 <- compare_designs(df.long=long_ct(vars_T8_m50),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=50, rho0=rho0)
p4 <- compare_designs(df.long=long_ct(vars_T8_m500),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=500, rho0=rho0)
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Variance of treatment effect, continuous time")
ggsave(paste0("plots/conts_50_500.pdf"), p1to4, width=297, height=210, units="mm")

# Plot relative variance, continuous vs discrete time, all designs
ylims <- c(0.3,1.0)
p1 <- compare_designs(df.long=long_rel_dt(vars_T4_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=50, rho0=rho0)
p2 <- compare_designs(df.long=long_rel_dt(vars_T4_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=500, rho0=rho0)
p3 <- compare_designs(df.long=long_rel_dt(vars_T8_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=50, rho0=rho0)
p4 <- compare_designs(df.long=long_rel_dt(vars_T8_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=500, rho0=rho0)
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Relative variance, continuous vs discrete time")
ggsave(paste0("plots/conts_vs_disc_50_500.pdf"), p1to4, width=297, height=210, units="mm")

# Plot relative variance, continuous vs HH, all designs
ylims <- c(0.0,6.0)
p1 <- compare_designs(df.long=long_rel_HH(vars_T4_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=50, rho0=rho0)
p2 <- compare_designs(df.long=long_rel_HH(vars_T4_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=4, m=500, rho0=rho0)
p3 <- compare_designs(df.long=long_rel_HH(vars_T8_m50),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=50, rho0=rho0)
p4 <- compare_designs(df.long=long_rel_HH(vars_T8_m500),
                      ylabel="Variance ratio", ylimits=ylims, Tp=8, m=500, rho0=rho0)
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Relative variance, continuous time vs uniform")
ggsave(paste0("plots/conts_vs_HH_50_500.pdf"), p1to4, width=297, height=210, units="mm")

# Plot variances, discrete time, all designs
ylims <- c(0.0,0.06)
p1 <- compare_designs(df.long=long_dt(vars_T4_m50),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=50, rho0=rho0)
p2 <- compare_designs(df.long=long_dt(vars_T4_m500),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=500, rho0=rho0)
p3 <- compare_designs(df.long=long_dt(vars_T8_m50),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=50, rho0=rho0)
p4 <- compare_designs(df.long=long_dt(vars_T8_m500),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=500, rho0=rho0)
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Variance of treatment effect, discrete time")
ggsave(paste0("plots/disc_50_500.pdf"), p1to4, width=297, height=210, units="mm")

# Plot variances, uniform, all designs
ylims <- c(0.0,0.06)
p1 <- compare_designs(df.long=long_HH(vars_T4_m50),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=50, rho0=rho0)
p2 <- compare_designs(df.long=long_HH(vars_T4_m500),
                      ylabel="Variance", ylimits=ylims, Tp=4, m=500, rho0=rho0)
p3 <- compare_designs(df.long=long_HH(vars_T8_m50),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=50, rho0=rho0)
p4 <- compare_designs(df.long=long_HH(vars_T8_m500),
                      ylabel="Variance", ylimits=ylims, Tp=8, m=500, rho0=rho0)
mylegend <- g_legend(p1)
p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend,
                            title="Variance of treatment effect, uniform correlation")
ggsave(paste0("plots/HH_50_500.pdf"), p1to4, width=297, height=210, units="mm")