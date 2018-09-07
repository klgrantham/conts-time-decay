# Generate variance results for plotting
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('vartheta_twolevels.R')

rhoCD <- 0.023
rhoUC <- 0.019
generate_var_results(Tp=4, N=6, m=50, rho0_CD=rhoCD, rho0_UC=rhoUC)
generate_var_results(Tp=4, N=6, m=150, rho0_CD=rhoCD, rho0_UC=rhoUC)
generate_var_results(Tp=8, N=14, m=50, rho0_CD=rhoCD, rho0_UC=rhoUC)
generate_var_results(Tp=8, N=14, m=150, rho0_CD=rhoCD, rho0_UC=rhoUC)

# Change this function
var_ct_mean_results(Tp=4, N=6, m=50, rho0=rhoCD)
var_ct_mean_results(Tp=4, N=6, m=150, rho0=rhoCD)
var_ct_mean_results(Tp=8, N=14, m=50, rho0=rhoCD)
var_ct_mean_results(Tp=8, N=14, m=150, rho0=rhoCD)

generate_var_results(Tp=4, N=6, m=50, rho0_CD=0.05, rho0_UC=0.05)
generate_var_results(Tp=4, N=6, m=10, rho0_CD=0.01, rho0_UC=0.01)
