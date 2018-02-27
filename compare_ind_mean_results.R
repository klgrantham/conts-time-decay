# Compare individual level to cluster-mean level, continuous time model
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('conts_time_comparison_plots.R')

var_ct_mean_results <- function(Tp, m, rho0){
  # Get variances using cluster-mean-level covariance matrices
  rs <- seq(0.5, 1, 0.01)
  # Specify the variance matrices under the continuous time model
  ctmeanvarmat <- llply(rs, expdecayVicont, Tp, m, rho0, meanlvl=TRUE)
  
  # Get the variances of the treatment effect under the
  # different models and designs
  # Note: Still need to expand the Xmats according to nclust
  # Scale the non-SW variances by (Tp/(Tp-1)) to account for uneven clusters across designs
  scalefactor <- Tp/(Tp-1)
  varvals <- data.frame(decay = 1-rs,
                        ctmeanSW = laply(ctmeanvarmat, vartheta_ind, Xmat=SWdesmat(Tp)),
                        ctmeancrxo = scalefactor*laply(ctmeanvarmat, vartheta_ind, Xmat=crxodesmat(Tp)),
                        ctmeanpllel = scalefactor*laply(ctmeanvarmat, vartheta_ind, Xmat=plleldesmat(Tp)),
                        ctmeanpllelbase = scalefactor*laply(ctmeanvarmat, vartheta_ind, Xmat=pllelbasedesmat(Tp)))
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2] # get numbers after decimal point
  save(varvals, file=paste0("plots/vars_ct_mean_T", Tp, "_m", m, "_rho", rho0char, ".Rda"))
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

# Plot variances, individual and mean levels, continuous time model
plot_indvsmean <- function(df_long, Tp, m, rho0){
  pct_indvsmean_all <- ggplot(data=df_long, aes(x=decay, group=design, colour=design)) +
    geom_line(aes(y=variances), size=2.0) +
    geom_line(aes(y=variances_mean), size=2.0, linetype="dotted") +
    expand_limits(y=0) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    labs(title="Variance of treatment effect, continuous time",
         subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
    scale_color_hue(labels = c("CRXO", "Parallel", "Parallel w/ baseline", "SW")) +
    guides(colour=guide_legend("Design")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.title=element_text(size=16), legend.text=element_text(size=14),
          legend.position="bottom")
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2] # get numbers after decimal point
  ggsave(paste0("plots/conts_indvsmean_compare_T", Tp, "_m", m, "_rho", rho0char, ".pdf"),
         pct_indvsmean_all, width=297, height=210, units="mm")
  return(pct_indvsmean_all)
}

# General plotting function using long format data frame
compare_designs <- function(df.long, title, ylabel, ylimits, Tp, m, rho0){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value, group=design, colour=design)) +
    geom_line(size=2.0) +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab(ylabel) +
    labs(title=title,
         subtitle=bquote(paste(T==.(Tp), ", ", m==.(m), ", ", rho[0]==.(rho0)))) +
    scale_color_hue(labels = c("CRXO", "Parallel", "Parallel w/ baseline", "SW")) +
    guides(colour=guide_legend("Design")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20),
          plot.subtitle=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.title=element_text(size=16), legend.text=element_text(size=14),
          legend.position="bottom")
  return(p)
}

Tp <- 4
m <- 500
rho0 <- 0.035

generate_ind_mean_plots <- function(Tp, m, rho0){
  
  var_ct_mean_results(Tp=Tp, m=m, rho0=rho0)
  
  rho0char <- strsplit(as.character(rho0),"\\.")[[1]][2] # get numbers after decimal point
  load(paste0("plots/vars_ct_mean_T", Tp, "_m", m, "_rho", rho0char, ".Rda")); vars_ct_mean <- varvals
  load(paste0("plots/vars_T", Tp, "_m", m, "_rho", rho0char, ".Rda")); vars_ct_ind <- varvals
  
  vars_ct_mean_long <- long_ct(vars_ct_mean)
  vars_ct_ind_long <- long_ct(vars_ct_ind)
  
  vars_ind_mean <- data.frame(vars_ct_ind_long,
                              variances_mean = vars_ct_mean_long$variances)
  
  
  plot_indvsmean(vars_ind_mean, Tp=Tp, m=m, rho0=rho0)
  
  
  # Plot variance ratios, individual vs mean levels, continuous time model
  vars_indvsmean_ratios <- vars_ind_mean %>%
    mutate(var_ratio_indvsmean = variances/variances_mean) %>%
    select(decay, design, var_ratio_indvsmean)
  
  p <- compare_designs(vars_indvsmean_ratios, title="Relative variance, continuous time decay",
                       ylabel="Variance ratio", ylimits=c(0.0,1.0), Tp=Tp, m=m, rho0=rho0)
  p2 <- p + geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash")
  ggsave(paste0("plots/conts_indvsmean_ratio_T", Tp, "_m", m, "_rho", rho0char, ".pdf"),
         p2, width=297, height=210, units="mm")
  
  # Plot precision, individual vs mean levels, continuous time model
  prec_ind_mean <- vars_ind_mean %>%
    mutate(precision = 1/variances, precision_mean = 1/variances_mean,
           precision_ratio = precision/precision_mean) %>%
    select(decay, design, precision, precision_mean, precision_ratio)
  
  p <- compare_designs(prec_ind_mean, title="Relative precision, continuous time decay",
                       ylabel="Precision ratio", ylimits=c(1.0,1.5), Tp=Tp, m=m, rho0=rho0)
  p2 <- p + geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash")
  ggsave(paste0("plots/conts_indvsmean_precratio_T", Tp, "_m", m, "_rho", rho0char, ".pdf"),
         p2, width=297, height=210, units="mm")
  
  # Plot precision, mean vs individual levels, continuous time model
  prec_mean_ind <- vars_ind_mean %>%
    mutate(precision = 1/variances, precision_mean = 1/variances_mean,
           precision_ratio = precision_mean/precision) %>%
    select(decay, design, precision, precision_mean, precision_ratio)
  
  p <- compare_designs(prec_mean_ind, title="Relative precision, continuous time decay",
                       ylabel="Precision ratio", ylimits=c(0.5,1.0), Tp=Tp, m=m, rho0=rho0)
  p2 <- p + geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash")
  ggsave(paste0("plots/conts_meanvsind_precratio_T", Tp, "_m", m, "_rho", rho0char, ".pdf"),
         p2, width=297, height=210, units="mm")
}

generate_ind_mean_plots(Tp=4, m=500, rho0=0.035)
generate_ind_mean_plots(Tp=4, m=50, rho0=0.035)
generate_ind_mean_plots(Tp=8, m=500, rho0=0.035)
generate_ind_mean_plots(Tp=8, m=50, rho0=0.035)
