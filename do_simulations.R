library(baydem)
library(magrittr)
library(Bchron)
library(testthat)

# Clear the workspace
rm(list=ls())

# Simulation results are created by calling, in turn, the following two scripts:
#
# (1) do_simulations.R [this script]
# (2) create_simulation_plots.R
#
# Aside from doing fits and Bayesian inference, do_simulations.R creates OxCal
# KDE models that are run using the OxCal web interface prior to creating the
# final publication results using create_simulation_plots.R.

# ------------------------------------------------------------------------------
# (0) Define functions used below
# ------------------------------------------------------------------------------

# Write files for doing KDE Plot fits with OxCal in the format expected by OxCal
write_KDE_Plot_model <- function(trc_m,sig_trc_m,file_name) {
  sink(file_name)
  cat('KDE_Plot()\n')
  cat('{\n')
  for(n in 1:length(trc_m)) {
    cat(paste0('R_Date(\"obs',
               as.character(n),
               '\",',
               as.character(trc_m[n]),
               ',',
               as.character(sig_trc_m[n]),');'))
    cat('\n')
  }
  cat('};\n')
  sink()
}

# Do a Gaussian mixture fit to samples drawn from the posterior densities of
# individually calibrated samples. The following function is based on:
# https://github.com/andrewcparnell/Bchron/blob/master/R/BchronDensityFast.R
#
# The only difference is that in the call to mclust::densityMclust a one
# dimensional variable/unequal variance model is used rather than letting the
# model be chosen via the BIC. This ensures that the calculation mirrors the
# maximum likelihood fitting and Bayesian inference done by baydem.
BchronDensityFast_modified <-
function(ages,ageSds,calCurves,pathToCalCurves=system.file('data',package='Bchron'),dfs=rep(100,length(ages)),samples=2000,G=30) {

if(length(ages)!=length(ageSds)) stop("ages and 1-sigma errors must be same length")
if(length(ages)!=length(calCurves)) stop("ages and Calibration curves must be same length")

# Calibrate ages
x = BchronCalibrate(ages=ages,ageSds=ageSds,calCurves=calCurves,pathToCalCurves=pathToCalCurves,dfs=rep(100,length(ages)))

# Get number of dates
n = length(x)

# Get a huge load of samples from the posteriors here
thetaBig = vector(length=n*samples)
for(i in 1:n) thetaBig[((i-1)*samples+1):(i*samples)] = sample(x[[i]]$ageGrid,size=samples,prob=x[[i]]$densities,replace=TRUE)

# Now run mclust
mclustOutput = mclust::densityMclust(data = thetaBig,
                                     G = G,
                                     modelNames = 'V')

output = list(out=mclustOutput,calAges=x)
class(output) = 'BchronDensityRunFast'
return(output)
}

# A wrapper function to do the mixture fit with BchronDensityFast_modified
do_BchronDensityFast_modified_fit <- function(trc_m,sig_trc_m,tau,K) {
  N <- length(trc_m)
  bchronFit <- BchronDensityFast_modified(ages=round(trc_m),
                                          ageSds=round(sig_trc_m),
                                          calCurves=rep('intcal20',N),G=K)
  f <- rep(0,length(tau))
  for(k in 1:K) {
    f <- f + bchronFit[[1]]$parameters$pro[k]*
      dnorm(1950-tau,
            bchronFit[[1]]$parameters$mean[k],
            sqrt(bchronFit[[1]]$parameters$variance$sigmasq[k]))
  }

  return(list(tau=tau,f=f,bchronFit=bchronFit))
}

# ------------------------------------------------------------------------------
# (1) Generate simulated data for N=10000 observations using the function
# baydem::simulate_rc_data, then subset the simulation to create simulations
# with N=100 and 1000 observations. If a save file with the simulation data
# has already been created, load the simulation data rather than creating
# the data again.
# ------------------------------------------------------------------------------

# Create the "target", simulation parameter vector, th_sim. The simulation
# distribution is a truncated Gaussian mixture with ordering pi1, pi2, mu1,
# mu2, sig1, sig2, and for which the truncation limits are tau_min = AD 600 and
# tau_max = AD 1300. The truncation has negligible effect since the probability
# density is exceedingly small at the truncation boundaries.
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

# Minimum calendar date (years AD)
tau_min <- 600
# Maximum calendar date (years AD)
tau_max <- 1300
# Spacing of calendar dates for the sampling grid
dtau <- 1

# Use the intcal20 calibration curve
calib_curve <- "intcal20"
calib_df <- load_calib_curve(calib_curve)

# Create the simulation specification. The error speficiation is that the error
# for the fraction modern value of each sample is uniformily drawn from the
# interval 0.0021 to 0.0028. Use a random number seed for reproducibility.
sim_spec <- list(model_spec=
                   list(density_type = "trunc_gauss_mix",
                        th=c(th_sim,tau_min,tau_max),
                        error_spec=list(type="unif_fm",min=.0021,max=.0028),
                        is_AD=T),
                 N=10000,
                 calib_curve=calib_curve,
                 seed=93004)

sim_file <- file.path("outputs","sim10000.rds")

if (!file.exists(sim_file)) {
  sim10000 <- simulate_rc_data(sim_spec)
  saveRDS(sim10000, sim_file)
} else {
  sim10000 <- readRDS(sim_file)
}

# Write out information about the range and mean of uncertainties in radiocarbon
# years for the section Model Validation: Simulations.
yaml::write_yaml(
  list(sig_trc_min=plyr::round_any(min(sim10000$data$rc_meas$sig_trc_m),.1),
       sig_trc_mean=plyr::round_any(mean(sim10000$data$rc_meas$sig_trc_m),.1),
       sig_trc_max=plyr::round_any(max(sim10000$data$rc_meas$sig_trc_m),.1)),
  "sig_trc_summary.yaml"
)

# (2) For each of N=100, N=1000, and N=10000, write the KDE Plot model to file
# (checking first whether a save file already exists for each case). These files
# can be run in OxCal's web interface to yield the curves that are plotted in
# create_simulation_plots.R. The OxCal web inteface cannot successfully run the
# N=10000 model, but it is written out anyway.
Nvect <- c(100,1000,10000)
# Ensure reproducibility using pre-calculated KDE file checksums
known_checksums <- c("98a9baf001c016ef29e6840ffbefa36d",
                     "0efc6390fd86a9cef41aeb2dff53b54a",
                     "325fcc0f96403ed92e02b12d054f06ed")
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- file.path("outputs",paste0("KDE_input_",N,".txt"))
  if (!file.exists(save_file)) {
    trc_m     <- sim10000$data$rc_meas$trc_m    [1:N]
    sig_trc_m <- sim10000$data$rc_meas$sig_trc_m[1:N]
    write_KDE_Plot_model(trc_m,sig_trc_m,save_file)
  }
  testthat::expect_equal(
     as.character(tools::md5sum(save_file)),
     known_checksums[m_N]
  )
}

# ------------------------------------------------------------------------------
# (2) For each of N=100, N=1000, and N=10000, do a maximum likelihood fit to
# data (checking first whether a save file already exists for each maximum
# likelihood fit). Use random number seeds for reproducibility. If necessary
# (because fewer cores are available or memory is limited), the number of cores
# used for parallel computation can be reduced from 10.
# ------------------------------------------------------------------------------
max_lik_fit_seeds <- c(739037,217298,971210)
num_cores <- 10

max_lik_fit_list <- list()
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- file.path("outputs",paste0("max_lik_fit",N,".rds"))
  seed <- max_lik_fit_seeds[m_N]
  if (!file.exists(save_file)) {
    phi_m <- sim10000$data$rc_meas$phi_m[1:N]
    sig_m <- sim10000$data$rc_meas$sig_m[1:N]
    max_lik_fit_list[[m_N]] <- fit_trunc_gauss_mix(2,
                                                   phi_m,
                                                   sig_m,
                                                   tau_min,
                                                   tau_max,
                                                   dtau,
                                                   calib_df,
                                                   num_restarts=100,
                                                   maxfeval=40000,
                                                   num_cores=num_cores,
                                                   input_seed=seed)
    saveRDS(max_lik_fit_list[[m_N]],save_file)
  } else {
    max_lik_fit_list[[m_N]] <- readRDS(save_file)
    expect_equal(
      max_lik_fit_list[[m_N]]$base_seed,
      seed
    )
  }
}

# ------------------------------------------------------------------------------
# (3) For each of N=100, N=1000, and N=10000, do a Bchron mixture fit (checking
# first whether a save file already exists for each Bchron mixture fit). Use
# random number seeds for reproducibility.
# ------------------------------------------------------------------------------
bchron_fit_seeds <- c(736877,860006,375035)
bchron_fit_list <- list()
tau_plot <- seq(tau_min,tau_max,by=dtau)
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- file.path("outputs",paste0("bchron_fit",N,".rds"))
  seed <- bchron_fit_seeds[m_N]
  if (!file.exists(save_file)) {
    set.seed(bchron_fit_seeds[m_N])
    trc_m     <- sim10000$data$rc_meas$trc_m    [1:N]
    sig_trc_m <- sim10000$data$rc_meas$sig_trc_m[1:N]
    bchron_fit_list[[m_N]] <- do_BchronDensityFast_modified_fit(trc_m,
                                                                sig_trc_m,
                                                                tau_plot,
                                                                2)
    saveRDS(bchron_fit_list[[m_N]],save_file)
  } else {
    bchron_fit_list[[m_N]] <- readRDS(save_file)
  }
}

# ------------------------------------------------------------------------------
# (4) Do Bayesian sampling for each N. As with previous steps, use random number
# seeds for reproducibility and load simulations from a save file if one exists.
# Also call summarize_bayesian_inference for each simulation.
# ------------------------------------------------------------------------------
# Set the hyperparameters
hp <- list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100
    alpha_r = (10 - 1) / 50,
    # The spacing for the Riemann sum (i.e., for the measurement matrix)
    dtau = 1
)

# Define the density model
density_model <- list(type="trunc_gauss_mix",
                      tau_min=tau_min,
                      tau_max=tau_max,
                      K=2)

# Save hp and density_model to file
saveRDS(hp, file.path("outputs","hp.rds"))
saveRDS(density_model, file.path("outputs","density_model.rds"))

# Do the inference
bayesian_soln_list <- list()
stan_seeds <- c(433582,774538,979639)
control <- list(samps_per_chain=4500,
                warmup = 2000)
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- file.path("outputs",paste0("bayesian_soln",N,".rds"))
  seed <- stan_seeds[m_N]
  if (!file.exists(save_file)) {
    phi_m <- sim10000$data$rc_meas$phi_m[1:N]
    sig_m <- sim10000$data$rc_meas$sig_m[1:N]
    rc_meas <- list(phi_m=phi_m,sig_m=sig_m)
    bayesian_soln <- sample_theta(rc_meas,
                         density_model,
                         hp,
                         calib_df,
                         th0=max_lik_fit_list[[m_N]]$th,
                         stan_seed=stan_seeds[m_N],
                         control=control)

    bayesian_soln_list[[m_N]] <- bayesian_soln
    saveRDS(bayesian_soln_list[[m_N]],save_file)
  } else {
    bayesian_soln_list[[m_N]] <- readRDS(save_file)
    expect_equal(
      bayesian_soln_list[[m_N]]$final_stan_seed,
      seed
    )
  }
}

# Calculate summary measures for the inference
bayesian_summ_list <- list()
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- file.path("outputs",paste0("bayesian_summ",N,".rds"))
  if (!file.exists(save_file)) {
    phi_m <- sim10000$data$rc_meas$phi_m[1:N]
    sig_m <- sim10000$data$rc_meas$sig_m[1:N]
    rc_meas <- list(phi_m=phi_m,sig_m=sig_m)
    bayesian_summ <- summarize_bayesian_inference(bayesian_soln_list[[m_N]],
                                                  rc_meas,
                                                  density_model,
                                                  calib_df,
                                                  dtau=1,
                                                  th_sim=th_sim)

    bayesian_summ_list[[m_N]] <- bayesian_summ
    saveRDS(bayesian_summ_list[[m_N]],save_file)
  } else {
    bayesian_summ_list[[m_N]] <- readRDS(save_file)
  }
}
