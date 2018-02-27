# Generate plots of variance ratios for YSC2017 presentation
# 
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)

source('vartheta_twolevels.R')

#vars_T4_m500_rho005 <- comparison_plots(Tp=4, m=500, rho0=0.005)
#vars_T4_m50_rho035 <- comparison_plots(Tp=4, m=50, rho0=0.035)

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

ctvdtT4 <- long_rel_dt_SW(vars_T4_m500_rho035)
ctvHHT4 <- long_rel_HH_SW(vars_T4_m500_rho035)
ctvHHT8 <- long_rel_HH_SW(vars_T8_m500_rho035)

# Compare ct to dt, T=4, m=500, rho0=0.035
plot_var_ratio(ctvdtT4, c(0.0, 1.0), "#C77CFF")
ggsave("plots/ctvdt_T4_m500_rho035.png", width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4, m=500, rho0=0.035
plot_var_ratio(ctvHHT4, c(0.0, 6.2), "#F8766D")
ggsave("plots/ctvHH_T4_m500_rho035.png", width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4 vs T=8, m=500, rho0=0.035
ctvHHT4T8 <- data.frame(ctvHHT4, relative_variance_T8 = ctvHHT8$relative_variance)
p <- plot_var_ratio(ctvHHT4T8, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_T8), size=2.0, color="#00BFC4")
ggsave("plots/ctvHH_T4_8_m500_rho035.png", width=7.99, height=5.67, units="in", dpi=900)

# Compare ct (rho0=0.04) to HH (rho0=0.035)
vars_T4_m500_rho035_040 <- data.frame(select(vars_T4_m500_rho035, decay, HHSW),
                                      ctSW = vars_T4_m500_rho040$ctSW)
ctvHHT4_rho035_040 <- long_rel_HH_SW(vars_T4_m500_rho035_040)
plot_var_ratio(ctvHHT4_rho035_040, c(0.0, 6.2), "#F8766D")
ggsave("plots/ctvHH_T4_m500_rho035_040.png", width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4, m=500, rho0=0.035 vs 0.005
ctvHHT4rho005 <- long_rel_HH_SW(vars_T4_m500_rho005)
ctvHHT4rho035_005 <- data.frame(ctvHHT4, relative_variance_rho005 = ctvHHT4rho005$relative_variance)
p <- plot_var_ratio(ctvHHT4rho035_005, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_rho005), size=2.0, color="#00BFC4")
ggsave("plots/ctvHH_T4_m500_rho035_005.png", width=7.99, height=5.67, units="in", dpi=900)

# compare ct to HH, T=4, m = 500 vs 50, rho0=0.035
ctvHHT4m50 <- long_rel_HH_SW(vars_T4_m50_rho035)
ctvHHT4m500_50 <- data.frame(ctvHHT4, relative_variance_m50 = ctvHHT4m50$relative_variance)
p <- plot_var_ratio(ctvHHT4m500_50, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_m50), size=2.0, color="#00BFC4")
ggsave("plots/ctvHH_T4_m500_50_rho035.png", width=7.99, height=5.67, units="in", dpi=900)
