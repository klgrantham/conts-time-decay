# Generate variance results for plotting
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('vartheta_twolevels.R')

rho <- 0.023
generate_var_results(Tp=4, m=50, rho0=rho)
generate_var_results(Tp=4, m=150, rho0=rho)
generate_var_results(Tp=8, m=50, rho0=rho)
generate_var_results(Tp=8, m=150, rho0=rho)

var_ct_mean_results(Tp=4, m=50, rho0=rho)
var_ct_mean_results(Tp=4, m=150, rho0=rho)
var_ct_mean_results(Tp=8, m=50, rho0=rho)
var_ct_mean_results(Tp=8, m=150, rho0=rho)

generate_var_results(Tp=4, m=50, rho0=0.05)
generate_var_results(Tp=4, m=10, rho0=0.01)
