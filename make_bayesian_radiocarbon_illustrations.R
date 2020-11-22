rm(list = ls())

library(baydem)
library(Bchron)
library(rcarbon)

# There are many excellent tools to can be used to calibrate a calendar date.
# Before plotting the illustration plots, compare baydem's calibration against
# rcarbon's to ensure they give the same results.



# The date to be calibrated, about 833 AD
t_m   <- 1950 - 830 # The uncalibrated radiocarbon date in years before present
sig_t_m <- 20       # The uncertainty of the measurement (20 "years")

# Do the calibration in rcarbon
rcarbonCalib <- rcarbon::calibrate(x=t_m,errors=sig_t_m,calCurves='intcal13',F14C=T,eps=1e-6)

# Now, do the calibration in baydem using the same gridding for the calendar
# dates (tau)
tau <- 1950 - rcarbonCalib$grids[[1]]$calBP
dtau <- tau[2] - tau[1]
tau_min <- min(tau)
tau_mid <- mean(range(tau))
tau_max <- max(tau)

# The calibration dataframe, which has three columns:
# yearBP	uncalYearBP	uncalYearBPError
calibDf = baydem::bd_load_calib_curve("intcal13")

# The reference decay rate of carbon
kappa <- 8033

# The fraction modern and associated uncertainty
phi_m <- exp(-t_m/kappa)
sig_m <- sig_t_m * phi_m / kappa

# Calculate the likelihood using baydem
lik <- as.vector(bd_calc_meas_matrix(tau,phi_m,sig_m,calibDf,T,F))

# Calculate the posterior density assuming a uniform prior
f_prior_unif <- rep(1,length(tau)) / (tau_max-tau_min)

f_posterior_unif <- f_prior_unif * lik
f_posterior_unif <- f_posterior_unif / sum(f_posterior_unif*dtau)

# Do the actual check to ensure rcarbon and baydem give the same results
# (within an appropriate tolerance)
testthat::expect_equal(
  f_posterior_unif,
  rcarbonCalib$grids[[1]]$PrDens,
  tol=1e-4
)

# Calculate the posterior density assuming a triangular prior

f_prior_tria <- (tau_mid - tau_min)-abs(tau-tau_mid)
f_prior_tria <- f_prior_tria / sum(dtau*f_prior_tria)

f_posterior_tria <- f_prior_tria * lik
f_posterior_tria <- f_posterior_tria / sum(f_posterior_tria*dtau)



fileName <- here::here("Fig1_single_date_calibration.pdf")
### FIGURE 1: Calibration of a single date (using two priors)
pdf(fileName, width = 10, height = 6)
  plot(tau,f_prior_unif,col='black',lwd=3,type='l',xlab='Calendar Date [AD]',ylab='Density',lty=3,ylim=c(0,max(f_posterior_tria)))
  lines(tau,f_prior_tria,col='grey',lwd=3,lty=3)
  lines(tau,f_posterior_unif,col='black',lwd=3)
  lines(tau,f_posterior_tria,col='grey',lwd=3)
dev.off()
