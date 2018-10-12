# Generate variance results for plotting
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('vartheta_twolevels.R')

# General variance results: CCD, DCD, UC
rhoCD <- 0.023
generate_var_results(Tp=4, N=6, m=50, rho0_CD=rhoCD)
generate_var_results(Tp=4, N=6, m=150, rho0_CD=rhoCD)
generate_var_results(Tp=8, N=14, m=50, rho0_CD=rhoCD)
generate_var_results(Tp=8, N=14, m=150, rho0_CD=rhoCD)

# CCD only, cluster-period mean level
var_ct_mean_results(Tp=4, N=6, m=50, rho0=rhoCD)
var_ct_mean_results(Tp=4, N=6, m=150, rho0=rhoCD)
var_ct_mean_results(Tp=8, N=14, m=50, rho0=rhoCD)
var_ct_mean_results(Tp=8, N=14, m=150, rho0=rhoCD)

# Under alternative rho values
generate_var_results(Tp=4, N=6, m=50, rho0_CD=0.05)
generate_var_results(Tp=4, N=6, m=10, rho0_CD=0.01)

# CCD with correlation parameters estimated from continuous-time
# time parameterisation
rhoCD <- 0.022
generate_var_results(Tp=4, N=6, m=50, rho0_CD=rhoCD)
generate_var_results(Tp=4, N=6, m=150, rho0_CD=rhoCD)
generate_var_results(Tp=8, N=14, m=50, rho0_CD=rhoCD)
generate_var_results(Tp=8, N=14, m=150, rho0_CD=rhoCD)
