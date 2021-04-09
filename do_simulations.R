library(baydem)
library(magrittr)
library(Bchron)
library(testthat)

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

# Load the result of a KDE fit done externally using the KDE web inteface.
load_KDE_Plot_fit <- function(file_name) {
  data_frame <- read.table(file_name,sep="\t")
  tau <- as.vector(data_frame[,1])
  f   <- as.vector(data_frame[,2])
  dtau <- unique(diff(tau))
  if (length(dtau) != 1) {
    stop("Spacing of tau is not even")
  }
  f <- f/sum(f)/dtau
  return(list(tau=tau,f=f))
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
# with N=100 and 1000 observations.
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

sim10000 <- simulate_rc_data(sim_spec)

# ------------------------------------------------------------------------------
# (2) For each of N=100, N=1000, and N=10000, do a maximum likelihood fit to
# data (checking first whether a save file already exists for each maximum
# likelihood fit). Use random number seeds for reproducibility. If necessary
# (because fewer cores are available or memory is limited), the number of cores
# used for parallel computation can be reduced from 10.
# ------------------------------------------------------------------------------
max_lik_fit_seeds <- c(739037,217298,971210)
num_cores <- 10

Nvect <- c(100,1000,10000)
max_lik_fit_list <- list()
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- paste0("max_lik_fit",N,".rds")
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
# (3) Create plots for non-Bayesian "fits" for each N. Curves included are:
#
# (a) The target density
# (b) The maximum likelihood fit
# (c) The summed density
# (d) A Bchron mixture fit
# (e) An OxCal KDE (for N=100 and N=1000)
#
# This yields the following publication result(s):
# Fig1_non_bayesian_fits.pdf
# ------------------------------------------------------------------------------

# Write the KDE plots
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- paste0("KDE_input_",N,".txt")
  trc_m <- sim10000$data$rc_meas$trc_m[1:N]
  sig_trc_m <- sim10000$data$rc_meas$sig_trc_m[1:N]
  write_KDE_Plot_model(trc_m,sig_trc_m,save_file)
}

# Show non-Bayesian fits
num_plots <- 3
tau_plot <- seq(tau_min,tau_max,by=dtau)
f_sim <- calc_gauss_mix_pdf(th_sim,tau_plot,tau_min=tau_min,tau_max=tau_max)

M100 <- calc_meas_matrix(tau_plot,
                         sim10000$data$rc_meas$phi_m[1:100],
                         sim10000$data$rc_meas$sig_m[1:100],
                         calib_df
                         )
M100 <- M100 / replicate(length(tau_plot),rowSums(M100)*dtau)
f_spdf100 <- colMeans(M100)
f_ml100 <- calc_gauss_mix_pdf(max_lik_fit_list[[1]]$th,
                              tau_plot,
                              tau_min=tau_min,
                              tau_max=tau_max)
kde100 <- load_KDE_Plot_fit("posterior_from_oxcal_100.txt")
bchron100 <-
  do_BchronDensityFast_modified_fit(sim10000$data$rc_meas$trc_m[1:100],
                                    sim10000$data$rc_meas$sig_trc_m[1:100],
                                    tau_plot,
                                    2)
M1000 <- calc_meas_matrix(tau_plot,
                          sim10000$data$rc_meas$phi_m[1:1000],
                          sim10000$data$rc_meas$sig_m[1:1000],
                          calib_df
                          )
M1000 <- M1000 / replicate(length(tau_plot),rowSums(M1000)*dtau)
f_spdf1000 <- colMeans(M1000)
f_ml1000 <- calc_gauss_mix_pdf(max_lik_fit_list[[2]]$th,
                               tau_plot,
                               tau_min=tau_min,
                               tau_max=tau_max)
kde1000 <- load_KDE_Plot_fit("posterior_from_oxcal_1000.txt")
bchron1000 <-
  do_BchronDensityFast_modified_fit(sim10000$data$rc_meas$trc_m[1:1000],
                                    sim10000$data$rc_meas$sig_trc_m[1:1000],
                                    tau_plot,
                                    2)
M10000 <- calc_meas_matrix(tau_plot,
                           sim10000$data$rc_meas$phi_m,
                           sim10000$data$rc_meas$sig_m,
                           calib_df
                           )
M10000 <- M10000 / replicate(length(tau_plot),rowSums(M10000)*dtau)
f_spdf10000 <- colMeans(M10000)
f_ml10000 <- calc_gauss_mix_pdf(max_lik_fit_list[[3]]$th,
                                tau_plot,
                                tau_min=tau_min,
                                tau_max=tau_max)
bchron10000 <-
  do_BchronDensityFast_modified_fit(sim10000$data$rc_meas$trc_m,
                                    sim10000$data$rc_meas$sig_trc_m,
                                    tau_plot,
                                    2)

pdf("Fig2_non_bayesian_fits.pdf",width=5,height=2.5*num_plots)

  par(
    mfrow = c(num_plots, 1),
    xaxs = "i", # No padding for x-axis
    yaxs = "i", # No padding for y-axis
    # outer margins with ordering bottom, left, top, right:
    oma = c(4, 2, 2, 2),
    # plot margins with ordering bottom, left, top, right:
    mar = c(2, 4, 0, 0)
  )

  # N=100
  plot(tau_plot,
       f_sim,
       xlim=c(tau_min,tau_max),
       ylim=c(0,0.01),
       xlab = "",
       ylab = "Density",
       xaxt = "n",
       yaxt = "n",
       col="blue",
       type="l",
       lwd=2
      )
  lines(tau_plot,f_spdf100,col="black",lwd=2)
  lines(tau_plot,f_ml100,col="red",lwd=2)
  lines(kde100$tau,kde100$f,col="grey",lwd=2)
  lines(bchron100$tau,bchron100$f,col="grey",lwd=2,lty=3)
  # Add a label inicating the value of N
  text(
    labels = "N = 100",
    x = 600,
    y = 0.009,
    pos = 4,
    cex = 2
  )
  # Add a legend (only to the top plot)
  legend("topright",
       legend = c("Target","Max. Lik.","SPD","Bchron Mix.","OxCal KDE"),
       lty = c(1,1,1,3,1),
       col = c("Blue","Red","Black","Grey","Grey"),
       lwd = 2)

  # N=1000
  plot(tau_plot,
       f_sim,
       xlim=c(tau_min,tau_max),
       ylim=c(0,0.01),
       xlab = "",
       ylab = "Density",
       xaxt = "n",
       yaxt = "n",
       col="blue",
       type="l",
       lwd=2
      )
  lines(tau_plot,f_spdf1000,col="black",lwd=2)
  lines(tau_plot,f_ml1000,col="red",lwd=2)
  lines(kde1000$tau,kde1000$f,col="grey",lwd=2)
  lines(bchron1000$tau,bchron1000$f,col="grey",lwd=2,lty=3)
  # Add a label inicating the value of N
  text(
    labels = "N = 1000",
    x = 600,
    y = 0.009,
    pos = 4,
    cex = 2
  )

  # N=10000
  plot(tau_plot,
       f_sim,
       xlim=c(tau_min,tau_max),
       ylim=c(0,0.01),
       xlab = "",
       ylab = "Density",
       xaxt = "n",
       yaxt = "n",
       col="blue",
       type="l",
       lwd=2
      )
  lines(tau_plot,f_spdf10000,col="black",lwd=2)
  lines(tau_plot,f_ml10000,col="red",lwd=2)
  lines(bchron10000$tau,bchron10000$f,col="grey",lwd=2,lty=3)
  # Add a label inicating the value of N
  text(
    labels = "N = 10000",
    x = 600,
    y = 0.009,
    pos = 4,
    cex = 2
  )

  axis(side = 1)
  mtext("Calendar Date [AD]", side = 1, line = 2.5, cex = 0.75)

dev.off()

# ------------------------------------------------------------------------------
# (4) Do Bayesian sampling for each N
# ------------------------------------------------------------------------------
# Set the hyperparameters
hp <- list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 3,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100
    alpha_r = (3 - 1) / 100,
    # The spacing for the Riemann sum (i.e., for the measurement matrix)
    dtau = 1
)

# Define the density model
density_model <- list(type="trunc_gauss_mix",
                      tau_min=tau_min,
                      tau_max=tau_max,
                      K=2)

bayesian_soln_list <- list()
stan_seeds <- c(433582,774538,979639)
control <- list(samps_per_chain=4500,
                warmup = 2000)
for(m_N in 1:length(Nvect)) {
  N <- Nvect[m_N]
  save_file <- paste0("bayesian_soln",N,".rds")
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
