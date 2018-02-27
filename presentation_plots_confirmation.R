# Generate plots of variance ratios for confirmation presentation
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
load("plots/vars_T4_m500_rho005.Rda"); vars_T4_m500_rho005 <- varvals


# Convert continuous time results to long format
long_ct <- function(df){
  ctvarvals <- df %>%
    select(decay, starts_with('ct'), -ends_with('base')) %>%
    filter(decay<=0.5)
  ctvarvals_long <- gather(data=ctvarvals, key=design, value=variance,
                           -decay, convert=TRUE)
  return(ctvarvals_long)
}

# Convert continuous time results to long format
long_HH <- function(df){
  HHvarvals <- df %>%
    select(decay, starts_with('HH'), -ends_with('base')) %>%
    filter(decay<=0.5)
  HHvarvals_long <- gather(data=HHvarvals, key=design, value=variance,
                           -decay, convert=TRUE)
  return(HHvarvals_long)
}

# Convert continuous time vs HH results to long format, SW design
long_rel_HH_SW <- function(df){
  ctvHHvarvals <- df %>%
    mutate(ratioSW=ctSW/HHSW) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=design, value=relative_variance,
                              ratioSW, convert=TRUE)
  return(ctvHHvarvals_long)
}

plot_var_ratio <- function(df.long, ylimits, color="#6666FF", alpha=1){
  p <- ggplot(data=df.long, aes(x=decay)) +
    geom_line(aes(y=relative_variance), size=2.0, color=color, alpha=alpha) +
    geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash") +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab("Variance ratio") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20))
  return(p)
}

plot_variance <- function(df.long, ylimits, color="#00BFC4"){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay)) +
    geom_line(aes(y=value), size=2.0, color=color) +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20))
  return(p)
}

# Functions for comparing all designs

# General plotting function using long format data frame
plot_variances <- function(df.long){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value, group=design, colour=design)) +
    geom_line(size=2.0) +
    scale_color_manual(values = c("#F8766D", "#00CC99", "#6666FF"), #c("#F8766D", "#7CAE00", "#C77CFF")
                       labels = c("CRXO", "Parallel", "SW")) +
    xlab("Decay (1 - r)") +
    ylab("Variance") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
          legend.position="none")
  return(p)
}

plot_var_ratios <- function(df.long, ylimits, color="#00BFC4"){
  p <- ggplot(data=df.long, aes(x=decay, group=design, colour=design)) +
    geom_line(aes(y=relative_variance), size=2.0) +
    geom_hline(yintercept=1.0, size=1.0, color="black", linetype="longdash") +
    scale_color_manual(values = c("#F8766D", "#00CC99", "#6666FF"),
                       labels = c("CRXO", "Parallel", "SW")) +
    expand_limits(y=ylimits) +
    xlab("Decay (1 - r)") +
    ylab("Variance ratio") +
    theme_bw() +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
          legend.position="none")
  return(p)
}

# Convert continuous time vs HH results to long format, all designs
long_rel_HH <- function(df){
  ctvHHvarvals <- df %>%
    mutate(ratioSW=ctSW/HHSW, ratiocrxo=ctcrxo/HHcrxo,
           ratiopllel=ctpllel/HHpllel) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=design, value=relative_variance,
                              -decay, convert=TRUE)
  return(ctvHHvarvals_long)
}


## SW design only
ctvHHT4 <- long_rel_HH_SW(vars_T4_m500)
ctvHHT8 <- long_rel_HH_SW(vars_T8_m500)

# Variances under ct, HH and ct&HH, T=4, m=500, rho0=0.04
ctSWvar <- subset(long_ct(vars_T4_m500), design=="ctSW")
plot_variance(ctSWvar, c(0,0.01), color="#FF9966")
ggsave(paste0("plots/ct_T4_m500_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

HHSWvar <- subset(long_HH(vars_T4_m500), design=="HHSW")
plot_variance(HHSWvar, c(0,0.01), color="#B33C00")
ggsave(paste0("plots/HH_T4_m500_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

ct_HHSWvar <- rbind(ctSWvar, HHSWvar)
p <- ggplot(data=ct_HHSWvar, aes(x=decay, y=variance, group=design, colour=design)) +
  geom_line(size=2.0) +
  scale_color_manual(values = c("#FF9966", "#B33C00")) + #c("#FF884D", "#00CC66")
  expand_limits(y=c(0,0.01)) +
  xlab("Decay (1 - r)") +
  ylab("Variance") +
  theme_bw() +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
        legend.position="none")
ggsave(paste0("plots/ct_HH_T4_m500_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

# Variance under ct & HH, SW design, T=4 (light) and T=8
varsT4 <- select(vars_T4_m500, c(decay, ctSWT4=ctSW, HHSWT4=HHSW))
varsT8 <- select(vars_T8_m500, c(ctSWT8=ctSW, HHSWT8=HHSW))
ct_HHSWvar_T4_8 <- cbind(varsT4, varsT8)
ct_HHSWvar_T4_8_long <- ct_HHSWvar_T4_8 %>% gather(key=design, value=variance, -decay, convert=TRUE)

p <- ggplot(data=ct_HHSWvar_T4_8_long, aes(x=decay, y=variance, group=design, colour=design, alpha=design)) +
  geom_line(size=2.0) +
  scale_color_manual(values = c("#FF9966", "#FF9966", "#B33C00", "#B33C00")) + #c("#FF884D", "#00CC66") #c("#99B3FF", "#3366FF")
  scale_alpha_manual(values = c(0.15, 1, 0.15, 1)) +
  expand_limits(y=c(0,0.01)) +
  xlab("Decay (1 - r)") +
  ylab("Variance") +
  theme_bw() +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
        legend.position="none")
ggsave(paste0("plots/ct_HH_T4_8_m500_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)


# Compare ct to HH, T=4, m=500, rho0=0.04
plot_var_ratio(ctvHHT4, c(0.0, 6.8), "#6666FF")
ggsave(paste0("plots/ctvHH_T4_m500_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4 vs T=8, m=500, rho0=0.04
ctvHHT4T8 <- data.frame(ctvHHT4, relative_variance_T8 = ctvHHT8$relative_variance)
p <- plot_var_ratio(ctvHHT4T8, c(0.0, 6.8), "#6666FF", alpha=0.15)
p + geom_line(aes(y=relative_variance_T8), size=2.0, color="#000099") #7f00de
ggsave(paste0("plots/ctvHH_T4_8_m500_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

# Compare ct to HH, T=4, m=500, rho0=0.04 vs 0.005
ctvHHT4rho005 <- long_rel_HH_SW(vars_T4_m500_rho005)
ctvHHT4rho04_005 <- data.frame(ctvHHT4, relative_variance_rho005 = ctvHHT4rho005$relative_variance)
p <- plot_var_ratio(ctvHHT4rho04_005, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_rho005), size=2.0, color="#00BFC4")
ggsave(paste0("plots/ctvHH_T4_m500_rho04_005_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

# m = 500 vs m = 50
ctvHHT4m50 <- long_rel_HH_SW(vars_T4_m50)
ctvHHT4m500_50 <- data.frame(ctvHHT4, relative_variance_m50 = ctvHHT4m50$relative_variance)
p <- plot_var_ratio(ctvHHT4m500_50, c(0.0, 6.2), "#F8766D")
p + geom_line(aes(y=relative_variance_m50), size=2.0, color="#00BFC4")
ggsave(paste0("plots/ctvHH_T4_m500_50_rho04_SW.png"), width=7.99, height=5.67, units="in", dpi=900)

## Compare all designs

# Plot variance under continuous, all designs
p <- plot_variances(long_ct(vars_T4_m500)) + expand_limits(y=c(0,0.06))
ggsave(paste0("plots/ct_T4_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot variance under HH, all designs
p <- plot_variances(long_HH(vars_T4_m500)) + expand_limits(y=c(0,0.06))
ggsave(paste0("plots/HH_T4_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot variance under ct & HH, all designs
ct_HHvar <- rbind(long_ct(vars_T4_m500), long_HH(vars_T4_m500))
p <- ggplot(data=ct_HHvar, aes(x=decay, y=variance, group=design, colour=design, alpha=design)) +
  geom_line(size=2.0) +
  scale_color_manual(values = c("#F8766D", "#00CC99", "#6666FF", "#F8766D", "#00CC99", "#6666FF"), #c("#F8766D", "#7CAE00", "#C77CFF")
                     labels = c("CRXO_CT", "Parallel_CT", "SW_CT", "CRXO_U", "Parallel_U", "SW_U")) +
  scale_alpha_manual(values = c(0.15, 0.15, 0.15, 1, 1, 1),
                     labels = c("CRXO_CT", "Parallel_CT", "SW_CT", "CRXO_U", "Parallel_U", "SW_U")) +
  expand_limits(y=c(0,0.06)) +
  xlab("Decay (1 - r)") +
  ylab("Variance") +
  theme_bw() +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
        legend.position="none")
ggsave(paste0("plots/ct_HH_T4_8_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)


# Plot relative variance, continuous vs HH, all designs
ylims <- c(0.0,4.5)
p1 <- plot_var_ratios(df.long=long_rel_HH(vars_T4_m500), ylimits=ylims)
ggsave(paste0("plots/ctvHH_T4_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot relative variance, continuous vs HH, all designs
ylims <- c(0.0,6.8)
p1 <- plot_var_ratios(df.long=long_rel_HH(vars_T4_m500), ylimits=ylims)
ggsave(paste0("plots/ctvHH_T4_m500_rho04_samescale.png"), width=7.99, height=5.67, units="in", dpi=900)
p2 <- plot_var_ratios(df.long=long_rel_HH(vars_T8_m500), ylimits=ylims)
ggsave(paste0("plots/ctvHH_T8_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot variance under continuous, all designs, T=4 (light), T=8
vars_T8_m500mod <- select(vars_T8_m500, c(decay, ctSWT8=ctSW, ctcrxoT8=ctcrxo, ctpllelT8=ctpllel))
ct_T4_8_var <- rbind(long_ct(vars_T4_m500), long_ct(vars_T8_m500mod))
p <- ggplot(data=ct_T4_8_var, aes(x=decay, y=variance, group=design, colour=design, alpha=design)) +
  geom_line(size=2.0) + # c("#F8766D", "#00CC99", "#6666FF")
  scale_color_manual(values = c("#F8766D", "#F8766D", "#00CC99", "#00CC99", "#6666FF", "#6666FF"), #c("#F8766D", "#7CAE00", "#C77CFF")
                     labels = c("CRXO_T4", "CRXO_T8", "Parallel_T4", "Parallel_T8", "SW_T4", "SW_T8")) +
  scale_alpha_manual(values = c(0.15, 1, 0.15, 1, 0.15, 1),
                     labels = c("CRXO_T4", "CRXO_T8", "Parallel_T4", "Parallel_T8", "SW_T4", "SW_T8")) +
  expand_limits(y=c(0,0.06)) +
  xlab("Decay (1 - r)") +
  ylab("Variance") +
  theme_bw() +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),legend.position="none")
ggsave(paste0("plots/ct_T8_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)

# Plot variance under HH, all designs, T=4 (light), T=8
vars_T8_m500mod <- select(vars_T8_m500, c(decay, HHSWT8=HHSW, HHcrxoT8=HHcrxo, HHpllelT8=HHpllel))
HH_T4_8_var <- rbind(long_HH(vars_T4_m500), long_HH(vars_T8_m500mod))
p <- ggplot(data=HH_T4_8_var, aes(x=decay, y=variance, group=design, colour=design, alpha=design)) +
  geom_line(size=2.0) + # c("#F8766D", "#00CC99", "#6666FF")
  scale_color_manual(values = c("#F8766D", "#F8766D", "#00CC99", "#00CC99", "#6666FF", "#6666FF"), #c("#F8766D", "#7CAE00", "#C77CFF")
                     labels = c("CRXO_T4", "CRXO_T8", "Parallel_T4", "Parallel_T8", "SW_T4", "SW_T8")) +
  scale_alpha_manual(values = c(0.15, 1, 0.15, 1, 0.15, 1),
                     labels = c("CRXO_T4", "CRXO_T8", "Parallel_T4", "Parallel_T8", "SW_T4", "SW_T8")) +
  expand_limits(y=c(0,0.06)) +
  xlab("Decay (1 - r)") +
  ylab("Variance") +
  theme_bw() +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),legend.position="none")
ggsave(paste0("plots/HH_T8_m500_rho04.png"), width=7.99, height=5.67, units="in", dpi=900)
