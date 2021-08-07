rm(list = ls())

library(baydem)
library(magrittr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mesorad_file <- file.path("outputs", "mesorad_final.csv")
if (!file.exists(mesorad_file)) {
  stop("outputs/mesorad_final.csv is missing. Run process_mesorad_dataset.R.")
}
mesorad <- read.csv(mesorad_file, stringsAsFactors=FALSE, check.names=FALSE)

# Check that the sums of the following two Mesorad columns have not changed
# (processing should be deterministic since a random number seed is used):
# "Conventional 14C age (BP)" and "Error".
testthat::expect_equal(
  sum(mesorad[,"Conventional 14C age (BP)"]),
  2126606
)

testthat::expect_equal(
  sum(mesorad[,"Error"]),
  58821
)

# Add columns to mesorad with the uncalibrated years BP and corresponding error
# using the baydem naming conventions.
mesorad$trc_m <- mesorad[,"Conventional 14C age (BP)"]
mesorad$sig_trc_m <- mesorad[,"Error"]

# Calculate the fraction modern and associated uncertainty for all data
mesorad$phi_m <- exp(-mesorad$trc_m/8033)
mesorad$sig_m <- mesorad$sig_trc_m * mesorad$phi_m / 8033

data_dir <- "outputs"

hp <-
  list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100
    alpha_r = (10 - 1) / 100,
    # Spacing for the measurement matrix (years)
    dtau = 1
  )
saveRDS(hp, file.path("outputs","maya_hp.rds"))

# Do Bayesian inference with 2 through 14 mixtures:
Kvect <- 2:14
num_models <- length(Kvect)
# Define random number seeds for (1) initializing the Bayesian inference and
# (2) the sampling in Stan. The following seeds were created via the following
# set of commands:

base_seed <- 516797 # from random.org, between 1 and 1,000,000
set.seed(base_seed)
seed_mat <- matrix(sample.int(1000000,2*num_models),ncol=2)
init_seed_vect <- seed_mat[,1]
stan_seed_vect <- seed_mat[,2]
# This yields the following values:
#> print(init_seed_vect)
# [1] 936488 407207 972265 494435 801442  13600 192234 501397 538238 168588
#[11]  73910 520818 699346
#> print(stan_seed_vect)
# [1] 452234 273402 243017 355225 732891 861065 658936 479337 262233 514932
#[11] 921114 493516 747145

# To reduce memory usage, use a for loop here to do Bayesian inference rather
# than using the standard pipeline functions. To allow interupted runs to be
# continued, only do the inference if a save file does not exist.

# Use the intcal20 calibration curve
calib_curve <- "intcal20"
calib_df <- load_calib_curve(calib_curve)

# Use only the Tikal observations
rc_meas <- mesorad[mesorad$Site == "Tikal",]

# Define the base density model (without K)
density_model0 <- list(type="trunc_gauss_mix",
                       tau_min=-1100,
                       tau_max=1900)
# Save a list containing the variables used for the inference in the following
# for loop
tikal_inference_inputs <- list(Kvect=Kvect,
                               density_model0=density_model0,
                               rc_meas=rc_meas,
                               hp=hp,
                               calib_df=calib_df,
                               init_seed_vect=init_seed_vect,
                               stan_seed_vect=stan_seed_vect)
saveRDS(tikal_inference_inputs, file.path("outputs","tikal_inference_inputs.rds"))
for (m_K in 1:length(Kvect)) {
  K <- Kvect[m_K]
  save_file <- file.path("outputs", paste0("tikal_K",K,".rds"))
  if(!file.exists(save_file)) {
    density_model <- density_model0
    density_model$K <- K
    t0 <- Sys.time() # start time
    bayesian_solution <- sample_theta(rc_meas,
                                      density_model,
                                      hp,
                                      calib_df,
                                      init_seed=init_seed_vect[m_K],
                                      stan_seed=stan_seed_vect[m_K])
    t1 <- Sys.time() # end time
    run_time_sec <- as.numeric(difftime(t1,t0,units="secs"))
    # Calculate PSIS-LOO CV
    log_lik_mat <- rstan::extract(bayesian_solution$fit,"logh")[[1]]
    loo_analysis <- loo::loo(log_lik_mat)
    loo_value <- loo_analysis["estimates"][[1]]["elpd_loo","Estimate"]
    # Save the results to file
    saveRDS(list(bayesian_solution=bayesian_solution,
                 run_time_sec=run_time_sec,
                 loo_value=loo_value),
            save_file)
  }
}