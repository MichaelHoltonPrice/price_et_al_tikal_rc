# Visualize equifinality of exponentials in the vector rVect.
vis_exp_equif <- function(r_vect, tau_min, tau_max, calib_df, meas_error, G = 1000) {
  tau_curve <- 1950 - calib_df$year_BP
  phi_curve <- exp(-calib_df$uncal_year_BP / 8033)
  tau <- seq(tau_min, tau_max, len = G)
  phi_interp <- approx(tau_curve, phi_curve, tau)
  phi_interp <- phi_interp$y

  phi_min <- min(phi_interp)
  phi_max <- max(phi_interp)
  # Extend the range of possible phi measurements by 4 times the uncertainty
  phi_min <- phi_min - meas_error * 4
  phi_max <- phi_max + meas_error * 4
  phi_vect <- seq(phi_min, phi_max, len = G * 4)

  M <- calc_meas_matrix2(tau,
                         phi_vect,
                         rep(meas_error, length(phi_vect)),
                         calib_df)
  phi_pdf_mat <- matrix(NA, length(phi_vect), length(r_vect))
  for (j in 1:length(r_vect)) {
    f_j <- calc_exp_pdf(tau, r_vect[j], tau_min, tau_max)
    phi_pdf <- M %*% as.matrix(f_j)
    phi_pdf_mat[, j] <- phi_pdf
  }
  equi_list <- calc_calib_curve_equif_dates(calib_df)
  temp <- assess_calib_curve_equif(calib_df, equi_list)
  inv_span_list <- temp$inv_span_list

  pdf_min <- min(phi_pdf_mat)
  pdf_max <- max(phi_pdf_mat)
  # Add a blank plot
  plot(NULL,
       type="n",
       xlim=c(phi_min,phi_max),
       ylim=c(pdf_min,pdf_max),
       xlab="Fraction Modern",
       ylab="Probability Density")
  # Add grey bands to the phi plot where there is no equifinality in the
  # calibration curve
  for (ii in 1:length(inv_span_list)) {
    inv_reg <- inv_span_list[[ii]]
    if (dplyr::between(inv_reg$phi_left, phi_min, phi_max) ||
      dplyr::between(inv_reg$phi_right, phi_min, phi_max)) {
      rect(inv_reg$phi_left,
           pdf_min,
           inv_reg$phi_right,
           pdf_max,
           border=NA,
           col="gray80")
    }
  }
  # Plot density function for each r-value. Also check that there is no
  # identifiability problem. Although Null in package MASS is used, this
  # really amounts to showing that N, which is a column vector, is not the
  # zero vector.
  for (j in 1:length(r_vect)) {
    lines(phi_vect, phi_pdf_mat[, j], lwd = 2)
    P <- calc_perturb_mat_exp(tau, r_vect[j], tau_min, tau_max)
    N <- MASS::Null(t(M %*% P))
    if (ncol(N) != 0) {
      stop(paste("r = ",
                 as.character(r_vect[n]),
                 " is not locally identifiable", sep = ""))
    }
  }
}

# Calculate the probability density function for exponential growth / decay
# assuming a the distribution integrates to 1 on the interval taumin to
# taumax
calc_exp_pdf <- function(tau, r, tau_min, tau_max) {
  G <- length(tau)
  dTAU <- tau_max - tau_min # The inteval length
  if (r == 0) {
    f <- rep(1, length(tau)) / dTAU
  } else {
    f <- r * exp(r * (tau - tau_max)) / (1 - exp(-r * dTAU))
  }
  return(f)
}

# Calculate the perturbation matrix for the exponential growth / decay
# probablity density function
calc_perturb_mat_exp <- function(tau, r, tau_min, tau_max) {
  dTAU <- tau_max - tau_min # The inteval length
  if (r == 0) {
    S <- as.matrix(.5 + (tau - tau_max) / dTAU)
  } else {
    f <- r * exp(r * (tau - tau_max)) / (1 - exp(-r * dTAU))
    S <- f * (1 / r + tau - tau_max - dTAU * exp(-r * dTAU) / (1 - exp(-r * dTAU)))
    S <- as.matrix(S)
  }
  return(S)
}

calc_h <- function(th_reduced,M,tau,sig_min=0) {
  # If the parameter vector is invalid, return infinity
  if (!is_th_reduced_valid(th_reduced,sig_min)) {
    stop("th_reduced is invalid")
  }

  tau_min <- min(tau)
  tau_max <- max(tau)

  # Add an undersore to pi, the mixture proportions, since pi is 3.14... in R
  K <- (length(th_reduced) + 1)/3
  pi_ <- th_reduced[1:(K-1)]

  th <- c(1-sum(pi_),th_reduced)
  # Calculate v, the vector of densities
  v <- calc_gauss_mix_pdf(th,tau,tau_min,tau_max)
  h <- M %*% v
  return(h)
}

calc_gradient <- function(th_reduced,M,tau,sig_min=0) {
  tau_min <- min(tau)
  tau_max <- max(tau)
  K <- (length(th_reduced) + 1)/3
  th <- c(1-sum(th_reduced[1:(K-1)]),th_reduced)
  P <- calc_perturb_mat_gauss_mix(tau, th, tau_min, tau_max)

  h <- as.vector(calc_h(th_reduced,M,tau,sig_min))
  X1 <- M %*% P
  X2 <- replicate(ncol(P),1/h)
  X <- X1 * X2
  return(-as.vector(colSums(X)))
}

# Calculate the perturbation matrix for the (possibly) truncated Gaussian
# mixture distribution
calc_perturb_mat_gauss_mix <- function(tau, th, tau_min = NA, tau_max = NA) {
  K <- length(th) / 3 # Number of mixtures
  G <- length(tau) # Number of grid points
  P <- matrix(NA, G, 3 * K - 1) # perturbation matrix
  # eta is the normalization for (possible) truncation
  if (!is.na(tau_min)) {
    eta <- diff(calc_gauss_mix_pdf(th,
                                   c(tau_min,tau_max),
                                   type="cumulative"))
  } else {
    eta <- 1
  }

  # The first mixture proportion is fixed by the others:
  # pi_1 = pi_2 + pi_3 + ... + pi_K
  #
  # It is omitted from the perturbation matrix, and the preceding constraint
  # must be accounted for in calculating the perturbation matrix

  # Extract the parameters for mixture 1
  mu_1 <- th[K + 1]
  s_1 <- th[2 * K + 1]
  f_1 <- dnorm(tau, mu_1, s_1) # density for mixture 1

  # Iterate over mixtures to calculate the derivatives with respect to
  # pi_k, mu_k, and s_k.
  for (k in 1:K) {
    if (!is.na(tau_min)) {
      # The unnormalized density as a function of tau
      z <- calc_gauss_mix_pdf(th, tau)
    }
    # Extract the parameters for mixture k
    pi_k <- th[k]
    mu_k <- th[K + k]
    s_k <- th[2 * K + k]
    f_k <- dnorm(tau, mu_k, s_k) # density for mixture k
    if (k > 1) { # Since pi's sum to 1, omit pi_1 from P
      # pi_k term. f_1 enters the calculation since pi_1 is fixed by the other
      # mixture proportions (see above)
      P[, k - 1] <- (f_k - f_1) / eta

      # If a normalization is used, account for the partial derivative of eta
      if (!is.na(tau_min)) {
        P[, k - 1] <- P[, k - 1] - z / eta^2 *
          (pnorm((tau_max - mu_k) / s_k)
            - pnorm((tau_min - mu_k) / s_k)
            - pnorm((tau_max - mu_1) / s_1)
            + pnorm((tau_min - mu_1) / s_1))
      }
    }
    # mu_k and s_k terms
    P[, K + k - 1] <- f_k * pi_k * (tau - mu_k) / s_k^2 / eta # mu term
    P[, 2 * K + k - 1] <- f_k * pi_k *
      (-1 / s_k + (tau - mu_k)^2 / s_k^3) / eta # s term
    # If a normalization is used, account for the partial derivatives of eta
    if (!is.na(tau_min)) {
      P[, K + k - 1] <- P[, K + k - 1] + pi_k / s_k *
        (dnorm((tau_max - mu_k) / s_k) -
          dnorm((tau_min - mu_k) / s_k)) * z / eta^2
      P[, 2 * K + k - 1] <- P[, 2 * K + k - 1] + pi_k / s_k^2 *
        (dnorm((tau_max - mu_k) / s_k) * (tau_max - mu_k)
          - dnorm((tau_min - mu_k) / s_k) * (tau_min - mu_k)) * z / eta^2
    }
  }
  return(P)
}

# Determine if the sample is locally identified (only check for the density
# function of the fraction modern, which is equivalent to checking the
# size of the null space of M %*% P.
is_identified <- function(th, M, tau, taumin = NA, taumax = NA) {
  P <- calc_perturb_mat_gauss_mix(tau, th, taumin, taumax) # baydem package
  N <- Null(t(M %*% P)) # MASS package
  return(ncol(N) == 0)
}

# Sample for the parater vector of the Gaussian mixture using the
# hyperparameter list hp
sample_gm <- function(hp) {
  mu <- runif(hp$K, hp$tau_min, hp$tau_max)
  s <- rgamma(hp$K, shape = hp$alpha_s, rate = hp$alpha_r)
  pi_vect <- gtools::rdirichlet(1, rep(hp$alpha_d, hp$K))
  ths <- c(pi_vect, mu, s)
  return(ths)
}

# Sometimes ygrid is intended to be evenly spaced but, numerically, not quite
# so due to rounding issues. To align calculations with the assumption in the
# manuscript of evenly spaced tau-grids for the Riemann integration, use this
# measurement matrix function, which unlike calc_meas_matrix does not check
# whether the grid is regularly spaced.
calc_meas_matrix2 <- function(tau, phi_m, sig_m, calib_df, add_calib_unc = T) {
  # 	# tau is in AD
  #        if(!all(is.na(phiLim))) {
  #            phiMin <- phiLim[1]
  #            phiMax <- phiLim[2]
  # 	}



  # tau is in AD
  tau_BP <- 1950 - tau

  # extract the calibration curve variables and convert to fraction modern
  tau_curve <- rev(calib_df$year_BP)
  mu_c_curve <- exp(-rev(calib_df$uncal_year_BP) / 8033)
  sig_c_curve <- rev(calib_df$uncal_year_BP_error) * mu_c_curve / 8033

  # Interpolate curves at tau_BP to yield mu_c and sig_c
  mu_c <- stats::approx(tau_curve, mu_c_curve, tau_BP)
  mu_c <- mu_c$y
  sig_c <- stats::approx(tau_curve, sig_c_curve, tau_BP)
  sig_c <- sig_c$y

  PHI_m <- replicate(length(tau_BP), phi_m)
  SIG_m <- replicate(length(tau_BP), sig_m)

  MU_c <- t(replicate(length(phi_m), mu_c))
  if (add_calib_unc) {
    SIG_c <- t(replicate(length(sig_m), sig_c))
    SIG_sq <- SIG_m^2 + SIG_c^2
  } else {
    SIG_sq <- SIG_m^2
  }

  M <- exp(-(PHI_m - MU_c)^2 / (SIG_sq) / 2) / sqrt(SIG_sq) / sqrt(2 * pi)

  # Multiply by the the integration width
  dtau <- tau[2] - tau[1]
  M <- M * dtau
  return(M)
}
