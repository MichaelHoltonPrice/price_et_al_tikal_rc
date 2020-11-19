library(baydem)
library(magrittr)
library(Bchron)

## Generate simulated data.
# Sample sizes of the simulations
sim_samp <- c(10, 100, 1000, 10000)

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
    # The gamma distribution rate parameter for sigma, yielding a mode of 100 or 300, respectively
    alpha_r = list(
      (3 - 1) / 100,
      (3 - 1) / 300
    ),
    # Minimum calendar date (years BC/AD)
    taumin = 600,
    # Maximum calendar date (years BC/AD)
    taumax = 1300,
    # Spacing for the measurement matrix (years)
    dtau = 1,
    # Number of mixtures
    K = 2
  ) %>%
  purrr::cross()

# Generate the largest sample;
# We'll subset this in the function below, simulating
# the process of retrieving more radiocarbon dates
# Set the random number seed (seed from random.org)
set.seed(806372)
sim_dates <-
  tibble::tibble(
    date_AD = baydem::bd_sample_gauss_mix(
      # A really large number of samples from which to draw the
      # test datasets, in case someone wants to run more than 10,000
      N = 100000,
      th = th_sim,
      taumin = hp[[1]]$taumin,
      taumax = hp[[1]]$taumax
    )
  ) %>%
  dplyr::bind_cols(
    .,
    baydem::bd_draw_rc_meas_using_date(
      t_e = .$date_AD,

      # Load the calibration data frame by calling bd_load_calib_curve
      calibDf = bd_load_calib_curve("intcal13"),

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

sim_inference <-
  sim_samp %>%
  magrittr::set_names(., .) %>%
  as.list() %>%
  purrr::cross2(hp) %>%
  purrr::map(function(x) {
    x %<>%
      magrittr::set_names(c("n", "hp"))

    x$sim_dates <- sim_dates[1:x$n, ]

    x
  })

# If simulated data exists for runs, load it. If are missing for any run, generate
# new data for all runs.
# Doing inference involves three steps:
#
# (1) Generate the problem
# (2) Do the Bayesian sampling
# (3) Run some standard analyses
out_file <- here::here("sim_inference.rds")

if (
  !file.exists(out_file) ||
    !identical(
      out_file %>%
        readr::read_rds() %>%
        purrr::map(magrittr::extract, c("n", "hp", "sim_dates")),
      sim_inference
    )
) {
  if (file.exists(out_file)) {
    saved_results <-
      out_file %>%
      readr::read_rds() %>%
      purrr::map(magrittr::extract, c("n", "hp", "sim_dates"))

    sim_inference %<>%
      setdiff(saved_results)
  }

  sim_inference %<>%
    purrr::map(
      function(x) {
        prob <-
          list(
            phi_m = x$sim_dates$phi_m,
            sig_m = x$sim_dates$sig_m,
            hp = x$hp,
            calibDf = bd_load_calib_curve("intcal13"),
            # Define the control parameters for the call to Stan. Use 4500 total MCMC
            # samples, of which 2000 are warmup samples. Since four chains are used, this
            # yields 4*(4500-2000) = 10,000 total samples.
            control = list(
              sampsPerChain = 4500,
              warmup = 2000
            )
          )

        soln <-
          baydem::bd_do_inference(prob)

        anal <-
          baydem::bd_analyze_soln(
            soln = soln,
            th_sim = th_sim
          )

        x$sim_output <-
          tibble::lst(
            prob,
            soln,
            anal
          )

        return(x)
      }
    )

  if (file.exists(out_file)) {
    saved_results <-
      out_file %>%
      readr::read_rds()

    # Discard whatever results were just calculated
    saved_results %<>%
      purrr::discard(saved_results %>%
        purrr::map(magrittr::extract, c("n", "hp", "sim_dates")) %>%
        magrittr::is_in(sim_inference %>%
          purrr::map(magrittr::extract, c("n", "hp", "sim_dates"))))

    sim_inference <- base::union(
      saved_results,
      sim_inference
    )
  }

  # Save the full result set
  sim_inference %>%
    readr::write_rds(out_file,
      compress = "gz"
    )
}

# Load the data.
sim_inference <-
  "sim_inference.rds" %>%
  here::here() %>%
  readr::read_rds()

#### Make Simulated Plots ####
sim_inference %<>%
  purrr::keep(function(x) x$hp$alpha_r == ((3 - 1) / 300)) %>%
  purrr::keep(function(x) x$n != 10)

set.seed(416755)

nplots <- length(sim_inference) + 1


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

do_BchronDensityFast_modified_fit <- function(sim,tau_bc,K) {
  N <- nrow(sim$sim_dates)
  bchronFit <- BchronDensityFast_modified(ages=round(sim$sim_dates$trc_m),ageSds=round(sim$sim_dates$sig_trc_m),calCurves=rep('intcal13',N),G=K)
  f_bc <- rep(0,length(tau_bc))
  for(k in 1:K) {
    f_bc <- f_bc + bchronFit[[1]]$parameters$pro[k]*dnorm(1950-tau_bc,bchronFit[[1]]$parameters$mean[k],sqrt(bchronFit[[1]]$parameters$variance$sigmasq[k]))
  }

  return(list(tau_bc=tau_bc,f_bc=f_bc,bchronFit=bchronFit))
}

simData <- readRDS('sim_inference.rds')

tau_bc <- sim_inference[[1]]$sim_output$anal$tau
bc100   <- do_BchronDensityFast_modified_fit(simData[[2]],tau_bc,2)
bc1000  <- do_BchronDensityFast_modified_fit(simData[[3]],tau_bc,2)
bc10000 <- do_BchronDensityFast_modified_fit(simData[[4]],tau_bc,2)

# Add the BchronDensityFast_modified fits to sim_inference
sim_inference[[1]]$tau_bc <- bc100$tau_bc
sim_inference[[1]]$f_bc   <- bc100$f_bc
sim_inference[[2]]$tau_bc <- bc1000$tau_bc
sim_inference[[2]]$f_bc   <- bc1000$f_bc
sim_inference[[3]]$tau_bc <- bc10000$tau_bc
sim_inference[[3]]$f_bc   <- bc10000$f_bc

# Generate a 4 x 1 graph figure summarizing the simulation results
pdf(here::here("Fig1_sim_inference.pdf"), width = 5, height = 2.5 * nplots)

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

# (1) Calibration curve
par(mar = c(0, 4, 0, 0))
bd_vis_calib_curve(min(sim_inference[[1]]$sim_output$anal$tau),
  max(sim_inference[[1]]$sim_output$anal$tau),
  sim_inference[[1]]$sim_output$prob$calibDf,
  xlab = "",
  ylab = "Fraction Modern",
  xaxt = "n",
  invertCol = "gray80"
)
box()

# (2-nplots) Density plots
sim_inference %>%
  purrr::walk(function(x) {
    par(mar = c(0, 4, 0, 0))

    bd_make_blank_density_plot(x$sim_output$anal,
      ylim = c(0, 0.01),
      xlab = "",
      ylab = "Density",
      xaxt = "n",
      yaxt = "n"
    )

    bd_add_shaded_quantiles(x$sim_output$anal,
      col = "gray80"
    )

    bd_plot_summed_density(x$sim_output$anal,
      lwd = 2,
      add = T,
      col = "black"
    )

    lines(x$tau_bc,x$f_bc,
      lwd = 2,
      col = "black",
      lty=3
    )

    bd_plot_50_percent_quantile(x$sim_output$anal,
      lwd = 2,
      add = T,
      col = "red"
    )

    bd_plot_known_sim_density(x$sim_output$anal,
      lwd = 2,
      add = T,
      col = "blue"
    )

    text(
      labels = paste0("n = ", x$n),
      x = 600,
      y = 0.009,
      pos = 4,
      cex = 2
    )

    axis(
      side = 2,
      at = c(0, 0.002, 0.004, 0.006, 0.008)
    )

    box()
  })

axis(side = 1)
mtext("Calendar Date [AD]", side = 1, line = 2.5, cex = 0.75)

dev.off()
