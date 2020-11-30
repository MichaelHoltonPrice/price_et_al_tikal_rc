library(baydem)
library(magrittr)
library(Bchron)

## Generate simulated data.

# The simulation distribution: a Gaussian mixture with ordering pi1, pi2, mu1,
# mu2, sig1, sig2
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

# Set the hyperparameters
hp <-
  list(
    # Class of fit (Gaussian mixture)
    fitType = "gaussmix",
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 3,
    # The gamma distribution rate parameter for sigma, yielding a mode of 50
    alpha_r = (3 - 1) / 50,
    # Minimum calendar date (years BC/AD)
    taumin = 600,
    # Maximum calendar date (years BC/AD)
    taumax = 1300,
    # Spacing for the measurement matrix (years)
    dtau = 1,
    # Number of mixtures
    K = 2
  )

# Either load or generate a sample of 10,000 simulated observations
sim_data_file <- 'sim_data.rds'
if(file.exists(sim_data_file)) {
  sim_data <- readRDS(sim_data_file)
} else {
  # For reproducibility, set the random number seed (from random.org between
  # 0 and 1,000,000)
  set.seed(999959)
  
  # The following piping first samples for the true calendar dates of the
  # samples using bd_sample_gauss_mix, then, for each calendar date, simulates
  # the radiocarbon measurement process by calling bd_draw_rc_meas_using_date.
  sim_data <-
  tibble::tibble(
    date_AD = baydem::bd_sample_gauss_mix(
      N = 10000,
      th = th_sim,
      taumin = hp$taumin,
      taumax = hp$taumax
    )
  ) %>%
  dplyr::bind_cols(
    .,
    baydem::bd_draw_rc_meas_using_date(
      t_e = .$date_AD,

      # Load the calibration data frame by calling bd_load_calib_curve
      calibDf = bd_load_calib_curve("intcal20"),

      # For simulating radiocarbon measurements, a draw is made for the standard
      # deviation of the fraction modern from a uniform density on the interval 0.0021
      # to 0.0028. This is specified via the list errorSpec
      errorSpec = list(
        min = .0021,
        max = .0028
      ),
      isAD = T
    ) %>%
      tibble::as_tibble()
  )

  saveRDS(sim_data,sim_data_file)
}

# If necessary, do the inference for each case of N=100, N=1000, and N=10000.
# This involves :
#
# (1) Generating the problem
# (2) Doing the Bayesian inference (sampling from the posterior distribution)
# (3) Running some standard analyses.
#
# Where possible (if the pertinent calculation has already been done), results
# are loaded from file, with checks made on the validity of saved results. 

# To ensure reproducibility, explicitly set the random number seeds used to
# initialize the sampling and by Stan (seeds from random.org between 1 and
# 1,000,000).
Nvect        <- c(   100,  1000, 10000) # Number of samples for each run
initSeedVect <- c(945814,369420,791667) # Initialization seed (see bd_do_inference)
stanSeedVect <- c(306928,973096, 25406) # Stan seed (see bd_do_inference)

# Iterate over "runs", storing the inference result in the list sim_soln and the
# analyses that use the solutions in the list sim_anal.
sim_soln <- list()
sim_anal <- list()

for(r in 1:length(Nvect)) {
  N <- Nvect[r] # Number of samples for this "run"

  # If necessary, build the problem and do the inference
  sim_soln_file <- paste0('sim_soln_',N,'.rds')
  if(file.exists(sim_soln_file)) {
    # Since the inference file exists, load it from file
    sim_soln[[length(sim_soln) + 1]] <- readRDS(sim_soln_file)

    # Check the dimensions and random number seeds for this simulation "run"
    if(length(sim_soln[[length(sim_soln)]]$prob$phi_m) != N) {
      stop('Wrong length for phi_m in solution read from file')
    }
    
    if(sim_soln[[length(sim_soln)]]$prob$control$initSeed != initSeedVect[r]) {
      stop('Wrong initSeed in solution read from file')
    }
    
    if(sim_soln[[length(sim_soln)]]$prob$control$stanSeed != stanSeedVect[r]) {
      stop('Wrong stanSeed in solution read from file')
    }

  } else {
    # Since the inference file does not exist, do the inference and save it to
    # file

    # Specify the problem as a list in the format expected by bd_do_inference
    prob <-
        list(
          phi_m = sim_data$phi_m[1:N],
          sig_m = sim_data$sig_m[1:N],
          hp = hp,
          calibDf = bd_load_calib_curve("intcal20"),
          # Define the control parameters for the call to Stan. Use 4500 total MCMC
          # samples, of which 2000 are warmup samples. Since four chains are used, this
          # yields 4*(4500-2000) = 10,000 total samples.
          control = list(
            sampsPerChain = 4500,
            warmup = 2000,
            initSeed = initSeedVect[r],
            stanSeed = stanSeedVect[r]
          )
        )
  
     sim_soln[[length(sim_soln) + 1]] <- baydem::bd_do_inference(prob)
     saveRDS(sim_soln[[length(sim_soln)]],sim_soln_file)
  }

  # If necessary, call bd_analyze_soln do some standard analyses on the
  # posterior samples.
  sim_anal_file <- paste0('sim_anal_',N,'.rds')
  if(file.exists(sim_anal_file)) {
    # Since the analysis file exists, load it from file
    sim_anal[[length(sim_anal) + 1]] <- readRDS(sim_anal_file)
  } else {
    # Since the analysis file does not exist, do the analysis and save it to
    # file
    sim_anal[[length(sim_anal) + 1]] <- baydem::bd_analyze_soln(sim_soln[[length(sim_soln)]],th_sim=th_sim)
    saveRDS(sim_anal[[length(sim_anal)]],sim_anal_file)
  }
}

#### Make Simulated Plots ####

set.seed(416755)

nplots <- length(sim_soln) + 1 # four plots total

# The following function is based on:
# https://github.com/andrewcparnell/Bchron/blob/master/R/BchronDensityFast.R
#
# The only difference is that in the call to mclust::densityMclust a one
# dimensional variable/unequal variance model is used rather than letting the
# model be chosen via the BIC. This ensures that this calculation most closely
# mirrors the inference done with baydem.
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

do_BchronDensityFast_modified_fit <- function(sim_data,N,tau,K) {
  bchronFit <- BchronDensityFast_modified(ages=round(sim_data$trc_m[1:N]),ageSds=round(sim_data$sig_trc_m[1:N]),calCurves=rep('intcal20',N),G=K)
  f_bc <- rep(0,length(tau))
  for(k in 1:K) {
    f_bc <- f_bc + bchronFit[[1]]$parameters$pro[k]*dnorm(1950-tau,bchronFit[[1]]$parameters$mean[k],sqrt(bchronFit[[1]]$parameters$variance$sigmasq[k]))
  }

  return(list(tau=tau,f_bc=f_bc,bchronFit=bchronFit))
}

# If necessary, do the BchronDensityFast_modified fits
tau <- sim_anal[[1]]$tau

sim_bc <- list()
for(r in 1:length(Nvect)) {
  N <- Nvect[r]
  sim_bc_file <- paste0('sim_bc_',N,'.rds')
  if(file.exists(sim_bc_file)) {
    # Since the Bchron file exists, load it from file
    sim_bc[[length(sim_bc) + 1]] <- readRDS(sim_bc_file)
  } else {
    # Since the Bchron file does not exist, do the fit and save it to file
    sim_bc[[length(sim_bc) + 1]] <- do_BchronDensityFast_modified_fit(sim_data,Nvect[r],tau,2)
    saveRDS(sim_bc[[length(sim_bc)]],sim_bc_file)
  }
}

# Generate a 4 x 1 graph figure summarizing the simulation results
pdf("Fig2_sim_inference.pdf", width = 5, height = 2.5 * nplots)

par(
  mfrow = c(nplots, 1),
  xaxs = "i", # No padding for x-axis
  yaxs = "i", # No padding for y-axis
  # outer margins with ordering bottom, left, top, right:
  oma = c(4, 2, 2, 2),
  # plot margins with ordering bottom, left, top, right:
  mar = c(2, 4, 0, 0)
  # Don't add data if it falls outside plot window
  # xpd = F
)

# (1) Add the calibration curve (first plot)
par(mar = c(0, 4, 0, 0))
bd_vis_calib_curve(min(tau),
  max(tau),
  sim_soln[[1]]$prob$calibDf,
  xlab = "",
  ylab = "Fraction Modern",
  xaxt = "n",
  invertCol = "gray80"
)
box()

# (2) Add the density plots (remaining plots)

par(mar = c(0, 4, 0, 0))
for(r in 1:length(Nvect)) {
  # Make a blank plot
  bd_make_blank_density_plot(sim_anal[[r]],
    ylim = c(0, 0.01),
    xlab = "",
    ylab = "Density",
    xaxt = "n",
    yaxt = "n"
  )

  # Add the shaded quantiles
  bd_add_shaded_quantiles(sim_anal[[r]],
    col = "gray80"
  )

  # Add the summed probability density
  bd_plot_summed_density(sim_anal[[r]],
    lwd = 2,
    add = T,
    col = "black"
  )

  # Add the Bchron fit
  lines(tau,sim_bc[[r]]$f_bc,
    lwd = 2,
    col = "black",
    lty=3
  )

  # Add solid 50% quantile
  bd_plot_50_percent_quantile(sim_anal[[r]],
    lwd = 2,
    add = T,
    col = "red"
  )

  # Plot the known, target distribution
  bd_plot_known_sim_density(sim_anal[[r]],
    lwd = 2,
    add = T,
    col = "blue"
  )

  # Add a label inicating the value of N for each subplot
  text(
    labels = paste0("N = ", Nvect[r]),
    x = 600,
    y = 0.009,
    pos = 4,
    cex = 2
  )
   
  # Add y-axis tickmarks at four locations
  axis(
    side = 2,
    at = c(0, 0.002, 0.004, 0.006, 0.008)
  )

}

# Add an x-axis label to the bottom plot
axis(side = 1)
mtext("Calendar Date [AD]", side = 1, line = 2.5, cex = 0.75)

dev.off()
