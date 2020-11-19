# Visualize equifinality of exponentials in the vector rVect.
visExpEquif <- function(rVect, taumin, taumax, calibDf, measError, G = 1000) {
  tau_curve <- 1950 - calibDf$yearBP
  phi_curve <- exp(-calibDf$uncalYearBP / 8033)
  tau <- seq(taumin, taumax, len = G)
  phiInterp <- approx(tau_curve, phi_curve, tau)
  phiInterp <- phiInterp$y

  phiMin <- min(phiInterp)
  phiMax <- max(phiInterp)
  # Extend the range of possible phi measurements by 4 times the uncertainty
  phiMin <- phiMin - measError * 4
  phiMax <- phiMax + measError * 4
  phiVect <- seq(phiMin, phiMax, len = G * 4)
  dPhi <- phiVect[2] - phiVect[1]

  M <- calc_meas_matrix2(tau, phiVect, rep(measError, length(phiVect)), calibDf)
  phiPdfMat <- matrix(NA, length(phiVect), length(rVect))
  for (j in 1:length(rVect)) {
    f_j <- calcExpPdf(tau, rVect[j], taumin, taumax)
    phiPdf <- M %*% as.matrix(f_j)
    phiPdfMat[, j] <- phiPdf
  }
  equiList <- bd_calc_calib_curve_equif_dates(calibDf)
  temp <- bd_assess_calib_curve_equif(calibDf, equiList)
  canInvert <- temp$canInvert
  invSpanList <- temp$invSpanList

  pdfMin <- min(phiPdfMat)
  pdfMax <- max(phiPdfMat)
  # Add a blank plot
  plot(NULL, type = "n", xlim = c(phiMin, phiMax), ylim = c(pdfMin, pdfMax), xlab = "Fraction Modern", ylab = "Probability Density")
  # Add grey bands to the phi plot where there is no equifinality in the
  # calibration curve
  for (ii in 1:length(invSpanList)) {
    invReg <- invSpanList[[ii]]
    if (dplyr::between(invReg$phi_left, phiMin, phiMax) || dplyr::between(invReg$phi_right, phiMin, phiMax)) {
      rect(invReg$phi_left, pdfMin, invReg$phi_right, pdfMax, border = NA, col = "gray80")
    }
  }
  # Plot density function for each r-value. Also check that there is no
  # identifiability problem. Although Null in package MASS is used, this
  # really amounts to showing that N, which is a column vector, is not the
  # zero vector.
  for (j in 1:length(rVect)) {
    lines(phiVect, phiPdfMat[, j], lwd = 2)
    P <- calcPerturbMatExp(tau, phiVect, rVect[j], taumin, taumax)
    N <- MASS::Null(t(M %*% P))
    if (ncol(N) != 0) {
      stop(paste("r = ", as.character(rVect[n]), " is not locally identifiable", sep = ""))
    }
  }
}

# Calculate the probability density function for exponential growth / decay
# assuming a the distribution integrates to 1 on the interval taumin to
# taumax
calcExpPdf <- function(tau, r, taumin, taumax) {
  G <- length(tau)
  dTAU <- taumax - taumin # The inteval length
  if (r == 0) {
    f <- rep(1, length(tau)) / dTAU
  } else {
    f <- r * exp(r * (tau - taumax)) / (1 - exp(-r * dTAU))
  }
  return(f)
}

# Calculate the perturbation matrix for the exponential growth / decay
# probablity density function
calcPerturbMatExp <- function(tau, phiVect, r, taumin, taumax) {
  G <- length(tau)
  dTAU <- taumax - taumin # The inteval length
  if (r == 0) {
    S <- as.matrix(.5 + (tau - taumax) / dTAU)
  } else {
    f <- r * exp(r * (tau - taumax)) / (1 - exp(-r * dTAU))
    S <- f * (1 / r + tau - taumax - dTAU * exp(-r * dTAU) / (1 - exp(-r * dTAU)))
    S <- as.matrix(S)
  }
  return(S)
}

# Calculate the perturbation matrix for the (possibly) truncated Gaussian
# mixture distribution
calcPerturbMatGaussMix <- function(tau, th, taumin = NA, taumax = NA) {
  K <- length(th) / 3 # Number of mixtures
  G <- length(tau) # Number of grid points
  P <- matrix(NA, G, 3 * K - 1) # perturbation matrix
  # eta is the normalization for (possible) truncation
  if (!is.na(taumin)) {
    eta <- diff(bd_calc_gauss_mix_pdf(th, c(taumin, taumax), type = "cumulative"))
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
    if (!is.na(taumin)) {
      # The unnormalized density as a function of tau
      z <- bd_calc_gauss_mix_pdf(th, tau)
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
      if (!is.na(taumin)) {
        P[, k - 1] <- P[, k - 1] - z / eta^2 * (pnorm((taumax - mu_k) / s_k) - pnorm((taumin - mu_k) / s_k) - pnorm((taumax - mu_1) / s_1) + pnorm((taumin - mu_1) / s_1))
      }
    }
    # mu_k and s_k terms
    P[, K + k - 1] <- f_k * pi_k * (tau - mu_k) / s_k^2 / eta # mu term
    P[, 2 * K + k - 1] <- f_k * pi_k * (-1 / s_k + (tau - mu_k)^2 / s_k^3) / eta # s term
    # If a normalization is used, account for the partial derivatives of eta
    if (!is.na(taumin)) {
      P[, K + k - 1] <- P[, K + k - 1] + pi_k / s_k * (dnorm((taumax - mu_k) / s_k) - dnorm((taumin - mu_k) / s_k)) * z / eta^2
      P[, 2 * K + k - 1] <- P[, 2 * K + k - 1] + pi_k / s_k^2 * (dnorm((taumax - mu_k) / s_k) * (taumax - mu_k) - dnorm((taumin - mu_k) / s_k) * (taumin - mu_k)) * z / eta^2
    }
  }
  return(P)
}

# Determine if the sample is locally identified (only check for the density
# function of the fraction modern, which is equivalent to checking the
# size of the null space of M %*% P.
is_identified <- function(th, M, tau, taumin = NA, taumax = NA) {
  P <- calcPerturbMatGaussMix(tau, th, taumin, taumax)
  N <- MASS::Null(t(M %*% P))
  return(ncol(N) == 0)
}

# Sample for the parater vector of the Gaussian mixture using the
# hyperparameter list hp
sample_gm <- function(hp) {
  mu <- runif(hp$K, hp$taumin, hp$taumax)
  s <- rgamma(hp$K, shape = hp$alpha_s, rate = hp$alpha_r)
  piVect <- gtools::rdirichlet(1, rep(hp$alpha_d, hp$K))
  ths <- c(piVect, mu, s)
  return(ths)
}

# Sometimes ygrid is intended to be evenly spaced but, numerically, not quite
# so due to rounding issues. To align calculations with the assumption in the
# manuscript of evenly spaced tau-grids for the Riemann integration, use this
# measurement matrix function, which unlike bd_calc_meas_matrix does not check
# whether the grid is regularly spaced.
calc_meas_matrix2 <- function(tau, phi_m, sig_m, calibDf, addCalibUnc = T) {
  # 	# tau is in AD
  #        if(!all(is.na(phiLim))) {
  #            phiMin <- phiLim[1]
  #            phiMax <- phiLim[2]
  # 	}



  # tau is in AD
  tau_BP <- 1950 - tau

  # extract the calibration curve variables and convert to fraction modern
  tau_curve <- rev(calibDf$yearBP)
  mu_c_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_c_curve <- rev(calibDf$uncalYearBPError) * mu_c_curve / 8033

  # Interpolate curves at tau_BP to yield mu_c and sig_c
  mu_c <- stats::approx(tau_curve, mu_c_curve, tau_BP)
  mu_c <- mu_c$y
  sig_c <- stats::approx(tau_curve, sig_c_curve, tau_BP)
  sig_c <- sig_c$y

  PHI_m <- replicate(length(tau_BP), phi_m)
  SIG_m <- replicate(length(tau_BP), sig_m)

  MU_c <- t(replicate(length(phi_m), mu_c))
  if (addCalibUnc) {
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
