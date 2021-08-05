rm(list = ls())
library(dplyr)
library(RColorBrewer)
library(baydem)
library(evd)
library(doParallel)
library(foreach)
library(MASS)
registerDoParallel(detectCores())

# Redirect print outputs to file
sink(file.path("outputs", "create_identif_results_gm_sink.txt"))

source("bayesian_radiocarbon_functions.R")

set.seed(75372) # from random.org between 1 and 1,000,000

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

hp <-
  list(
    # Class of fit (Gaussian mixture)
    fitType = "gaussmix",
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 50 
    alpha_r = (10 - 1) / 50,
    # Minimum calendar date (years BC/AD)
    tau_min = 600,
    # Maximum calendar date (years BC/AD)
    tau_max = 1300,
    # Spacing for the measurement matrix (years)
    dtau = 1,
    # Number of mixtures
    K = 2
  )

# Locations for calendar date grid (spacing of 1 year)
tau_vect <- seq(hp$tau_min, hp$tau_max, by = hp$dtau)
G <- length(tau_vect)

S <- 8
TH <- matrix(NA, S, 6)

for (s in 1:S) {
  TH[s, ] <- sample_gm(hp)
}
Fmat <- calc_gauss_mix_pdf_mat(TH, tau_vect)

tau_min <- hp$tau_min
tau_max <- hp$tau_max

calib_df <- load_calib_curve("intcal20")
tau_curve <- 1950 - calib_df$year_BP
phi_curve <- exp(-calib_df$uncal_year_BP / 8033)

# Calculate the calibration curve fraction modern at the locations of the calendar grid
phi_interp <- approx(tau_curve, phi_curve, tau_vect)
phi_interp <- phi_interp$y

# Use a measurement error for the fraction modern of 1e-3
meas_error <- 1e-3

# Determine the range of values for the fraction modern, then increase the
# range by four times the measurement error on both edges of the range
phi_min <- min(phi_interp)
phi_max <- max(phi_interp)
phi_min <- phi_min - meas_error * 4
phi_max <- phi_max + meas_error * 4
phi_vect <- seq(phi_min, phi_max, len = G * 4)

# Calculate the fraction modern without adding calibration uncertainty
M <- calc_meas_matrix(tau_vect,
                      phi_vect,
                      rep(meas_error,length(phi_vect)),
                      calib_df,
                      add_calib_unc=F)

# Plot some sample curves (including their fraction modern probablity density)
# along with a visualization of the calibration curve
equif_data <- assess_calib_curve_equif(calib_df)
can_invert <- equif_data$can_invert
inv_span_list <- equif_data$inv_span_list

phi_min_plot <- 0
phi_max_plot <- 0.08
col_vector <- mapply(brewer.pal, S, "Set1")
pdf(file.path("outputs","FigS2_gm_example.pdf"), width = 20, height = 12)
par(mfrow = c(3, 1))

vis_calib_curve(tau_min,
                tau_max,
                calib_df,
                xlab="Calendar Date [AD]",
                ylab="Fraction Modern",
                invert_col="gray80")


plot(NULL,
     type="n",
     xlim=c(tau_min,tau_max),
     ylim=c(0,max(Fmat)),
     xlab="Calendar Date [AD]",
     ylab="Probability Density")

for (ii in 1:length(inv_span_list)) {
  inv_reg <- inv_span_list[[ii]]
  if (between(inv_reg$tau_left, tau_min, tau_max) ||
    between(inv_reg$tau_right, tau_min, tau_max)) {
    rect(inv_reg$tau_left,
         0,
         inv_reg$tau_right,
         max(Fmat),
         border=NA,
         col="gray80")
  }
}

# Highlight regions 153 and 155 in all the graphs
rect(inv_span_list[[153]]$tau_right,
     max(Fmat) - .0004,
     inv_span_list[[155]]$tau_left,
     max(Fmat),
     border=NA,
     col="indianred1")

for (s in 1:S) {
  th_s <- TH[s, ]
  f_s <- calc_gauss_mix_pdf(th_s, tau_vect)
  lines(tau_vect, f_s, col = col_vector[s], lwd = 4)
  P <- calc_perturb_mat_gauss_mix(tau_vect, th_s, tau_min, tau_max)
  N <- MASS::Null(t(M %*% P))
  if (ncol(N) != 0) {
    stop(paste("Identifiability problem with sample", s))
  }
}

# Create matrix of fraction modern data for plotting
phi_pdf_mat <- matrix(NA, S, nrow(M))
for (s in 1:S) {
  phi_pdf_mat[s, ] <- M %*% as.matrix(Fmat[s, ])
}

pdf_min <- 0
pdf_max <- max(phi_pdf_mat)
plot(NULL,
     type="n",
     xlim=c(phi_min,phi_max),
     ylim=c(pdf_min, pdf_max),
     xlab="Fraction Modern",
     ylab="Probability Density")


for (ii in 1:length(inv_span_list)) {
  inv_reg <- inv_span_list[[ii]]
  if (between(inv_reg$phi_left, phi_min, phi_max) ||
    between(inv_reg$phi_right, phi_min, phi_max)) {
    rect(inv_reg$phi_left,
         pdf_min,
         inv_reg$phi_right,
         pdf_max,
         border=NA,
         col="gray80")
  }
}

rect(inv_span_list[[153]]$phi_right,
     max(phi_pdf_mat) - 5,
     inv_span_list[[155]]$phi_left,
     max(phi_pdf_mat),
     border=NA,
     col="indianred1")

for (s in 1:S) {
  lines(phi_vect, phi_pdf_mat[s, ], col = col_vector[s], lwd = 4)
}

dev.off()

# Check local identifiability of simulation parameter vector
if (!is_identified(th_sim, M, tau_vect, tau_min, tau_max)) {
  stop("Simulation parameter vector is not identified")
} else {
  print(paste0("Simulation parameter vector is identified, ",
               "with a measurement error of ", meas_error))
}

# Check local identifiability for a large number of random samples
S <- 100000

identified <- rep(F, S)
print(paste("Checking local identifiability for", S, "samples"))
TH_local <- matrix(NA, S, 6)
for (s in 1:S) {
  TH_local[s, ] <- sample_gm(hp)
}

identified <- foreach(s = 1:S,
                      .combine=cbind,
                      .packages=c('baydem','MASS')) %dopar% {
  output <- is_identified(TH_local[s, ], M, tau_vect, tau_min, tau_max)
}

num_bad <- sum(!identified)
print(paste0(num_bad,
             " (out of ", S, ") non-identifiable parameter vectors found"))

# If non-identifiable parameters are found, determine whether it is caused by P
# being non-identified. If not, throw an error.
if (sum(!identified) > 0) {
  ind_bad <- which(!identified)
  for (ii in 1:length(ind_bad)) {
    s <- ind_bad[ii]
    print(paste("Sample", s, "is not identified"))
    th_s <- TH_local[s, ]
    print(th_s)
    P <- calc_perturb_mat_gauss_mix(tau_vect, th_s, tau_min, tau_max)
    N_P <- MASS::Null(t(P))
    print(paste("The null size of P is", ncol(N_P)))
    if (ncol(N_P) == 0) {
      stop("Sample is not identified even though the null size of P is zero")
    }
  }
}

num_checks <- 100000
bad_loc_list <- list()
print(paste("Checking for non-identifiable pairs with 2 mixtures"))
start_time <- Sys.time()
print("----")
num_loc_phi <- 0
num_loc_f <- 0
num_loc_phi_and_f <- 0
num_loc_phi_not_f <- 0
num_pair <- 0
rel_tol <- 1e-6 # relative tolerance for checking fraction modern equality
for (cc in 1:num_checks) {
  th_a <- sample_gm(hp)
  th_b <- sample_gm(hp)
  f_a <- calc_gauss_mix_pdf(th_a, tau_vect, tau_min, tau_max)
  f_b <- calc_gauss_mix_pdf(th_b, tau_vect, tau_min, tau_max)

  phi_pdf_a <- M %*% f_a
  phi_pdf_b <- M %*% f_b
  equal_tol <- mean(c(phi_pdf_a, phi_pdf_b)) * rel_tol
  if (all(abs(phi_pdf_b - phi_pdf_a) <= equal_tol)) {
    num_pair <- num_pair + 1
    print("----")
    print(p)
    print(th_a)
    print(th_b)
  }

  # Check both parameterizations for local identifiability
  Pa <- calc_perturb_mat_gauss_mix(tau_vect, th_a, tau_min, tau_max)
  Pb <- calc_perturb_mat_gauss_mix(tau_vect, th_a, tau_min, tau_max)
  Na <- MASS::Null(t(M %*% Pa))
  Na_P <- MASS::Null(t(Pa))
  Nb <- MASS::Null(t(M %*% Pb))
  Nb_P <- MASS::Null(t(Pb))

  phiBad_a <- ncol(Na) != 0
  fBad_a <- ncol(Na_P) != 0
  if (phiBad_a) {
    num_loc_phi <- num_loc_phi + 1
  }

  if (fBad_a) {
    num_loc_f <- num_loc_f + 1
  }

  if (phiBad_a && fBad_a) {
    num_loc_phi_and_f <- num_loc_phi_and_f + 1
  }

  if (phiBad_a && (!fBad_a)) {
    num_loc_phi_not_f <- num_loc_phi_not_f + 1
    bad_loc_list[[length(badLocList) + 1]] <- th_a
  }

  phi_bad_b <- ncol(Nb) != 0
  f_bad_b <- ncol(Nb_P) != 0
  if (phi_bad_b) {
    num_loc_phi <- num_loc_phi + 1
  }

  if (f_bad_b) {
    num_loc_f <- num_loc_f + 1
  }

  if (phi_bad_b && f_bad_b) {
    num_loc_phi_and_f <- num_loc_phi_and_f + 1
  }

  if (phi_bad_b && (!f_bad_b)) {
    num_loc_phi_not_f <- num_loc_phi_not_f + 1
    bad_loc_list[[length(badLocList) + 1]] <- th_b
  }
  end_time <- Sys.time()
}
print(paste("Finished checking in", as.character(end_time - start_time)))
print(paste("Number of pairs checked is", as.character(num_checks)))
print(paste(as.character(num_pair),
            "observationally equivalent pairs found"))
print(paste(as.character(num_loc_phi_not_f),
            "parameterizations fail for phi but not f"))
print(paste(as.character(num_loc_f), "parameterizations fail for f"))
print(paste(as.character(num_loc_phi), "parameterizations fail for phi"))
print(paste(as.character(num_loc_phi_and_f),
            "parameterizations fail for phi and f"))

# End re-direct of print statements
sink()