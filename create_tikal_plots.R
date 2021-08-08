library(baydem)
rm(list = ls())

# Read Tikal inference inputs
tikal_inputs_file <- file.path("outputs","tikal_inference_inputs.rds")
if (!file.exists(tikal_inputs_file)) {
  stop("outputs/tikal_inference_inputs.rds not found")
}

tikal_inputs <- readRDS(tikal_inputs_file)
# For code clarity, extract variables in tikal_inputs
Kvect          <- tikal_inputs$Kvect
density_model0 <- tikal_inputs$density_model0
rc_meas        <- tikal_inputs$rc_meas
hp             <- tikal_inputs$hp
calib_df       <- tikal_inputs$calib_df
init_seed_vect <- tikal_inputs$init_seed_vect
stan_seed_vect <- tikal_inputs$stan_seed_vect

# Create a vector of loo values
loo_vect <- rep(NA,length(Kvect))
for (m_K in 1:length(Kvect)) {
  K <- Kvect[m_K]
  save_file <- file.path("outputs", paste0("tikal_K",K,".rds"))
  if(!file.exists(save_file)) {
    stop(paste0("Missing save file for K=",K))
  }
  bayesian_soln <- readRDS(save_file)
  loo_vect[m_K] <- bayesian_soln$loo
}

### FIGURES S3: loo plot
# Plot the loo vs K
pdf(file.path("outputs","FigS3_tikal_loo.pdf"),width=8,height=6)
  plot(Kvect,loo_vect,xlab="K",ylab="loo")
dev.off()

# Identify the best model and write a summary of the loo
m_K_best <- which.max(loo_vect)
K_best <- Kvect[m_K_best]
loo_summary <- list(m_K_best=m_K_best,K_best=K_best,loo_vect=loo_vect)
yaml::write_yaml(loo_summary, file.path("outputs","loo_summary.yaml"))

# Do the Baysian summary for the best model
density_model <- density_model0
density_model$K <- K_best
save_file <- file.path("outputs", paste0("tikal_K",K_best,".rds"))
bayesian_soln <- readRDS(save_file)
bayesian_soln <- bayesian_soln$bayesian_solution
bayesian_summ <- summarize_bayesian_inference(bayesian_soln,rc_meas,density_model,calib_df,hp$dtau)

### FIGURE 2: Tikal densities
# Plot the Tikal reconstruction using the best Bayesian solution (per loo)
file_name <- file.path("outputs",paste0("Fig4_tikal_inference.pdf"))

pdf(file_name, width = 10, height = 10)

par(
  mfrow = c(2, 1),
  xaxs = "i", # No padding for x-axis
  yaxs = "i", # No padding for y-axis
  oma = c(4, 2, 2, 2),
  mar = c(0, 4, 0, 0)
)
# First: plot the calibration curve
vis_calib_curve(
  density_model$tau_min,
  density_model$tau_max,
  calib_df,
  xlab = "",
  ylab = "Fraction Modern",
  xaxt = "n",
  invert_col = "gray80"
)
box()

# Second: plot the Tikal reconstruction
make_blank_density_plot(bayesian_summ,
  xlab = "",
  xaxt = "n",
  yaxt = "n",
  xlim = c(-1100, 1900),
  ylim = c(0, 0.004)
)
plot_summed_density(bayesian_summ,
                    lwd = 3,
                    add = T,
                    col = "black")
plot_50_percent_quantile(bayesian_summ,
                         add = T,
                         lwd = 3,
                         col = "blue")
add_shaded_quantiles(bayesian_summ,
                     col = adjustcolor("blue", alpha.f = 0.25))
xat <- seq(-800, 1600, 200)
yat <- c(0, .001, .002, .003)
xlab <- xat
ylab <- yat
axis(side = 1, at = xat, labels = xlab)
axis(side = 2, at = yat, labels = ylab)
mtext("Year (AD)", side = 1, line = 2.5)
dev.off()

### FIGURE 3: Tikal Expert Comparison

calc_interval_densities <- function(dat) {
  # dat is matrix-like read in with read_excel
  # return value is a matrix with three columns:
  #
  # lower calendar date for the interval
  # upper calendar date for the interval
  # normalized density
  tau <- as.numeric(unlist(dat[, 1]))
  f <- as.numeric(unlist(dat[, 2]))
  numInt <- length(tau) / 2 # number of intervals
  flo <- f[2 * (1:numInt) - 1]
  fhi <- f[2 * (1:numInt)]
  if (!all(flo == fhi)) {
    stop("f values are not consistent")
  }
  if (!all(diff(tau)[seq(2, (numInt * 2 - 2), by = 2)] == 0)) {
    stop("tau values are not consistent")
  }
  # Interval durations
  taulo <- tau[seq(1, 2 * numInt, by = 2)]
  tauhi <- tau[seq(2, 2 * numInt, by = 2)]
  intDur <- tauhi - taulo
  return(cbind(taulo, tauhi, flo / sum(flo * intDur)))
}

add_interval_density_to_plot <- function(dat, ...) {
  if (ncol(dat) == 2) {
    lines(dat[, 1], dat[, 2], col = "black", lwd = 3)
  } else {
    for (i in 1:nrow(dat)) {
      if (i == 1) {
        lines(c(dat[i, 1], dat[i, 1]), c(0, dat[i, 3]), lwd = 3, ...)
      } else {
        lines(c(dat[i, 1], dat[i, 1]), c(dat[i - 1, 3], dat[i, 3]), lwd = 3, ...)
      }
      lines(c(dat[i, 1], dat[i, 2]), c(dat[i, 3], dat[i, 3]), lwd = 3, ...)
      if (i == nrow(dat)) {
        lines(c(dat[i, 2], dat[i, 2]), c(dat[i, 3], 0), lwd = 3, ...)
      }
    }
  }
}

calc_point_densities <- function(dat) {
  # dat is matrix-like read in with read_excel
  # return value is a matrix with the two columns:
  #
  # calendar date
  # normalized density
  tau <- as.numeric(unlist(dat[, 1]))
  f <- as.numeric(unlist(dat[, 2]))
  weightVect <- calc_trapez_weights(tau)
  f <- f / sum(f * weightVect)
  return(cbind(tau, f))
}

expert_recons <-
  c(
    "Haviland" = "Haviland2003",
    "Turner" = "Turner-broaderTikal",
    "Culbert" = "Culbert-central-adjusted",
    "Fry" = "Fry-centralTikal",
    "Santley" = "Santley-tikal"
  ) %>%
  purrr::map(~ readxl::read_excel("Tikal_Demography.xlsx",
    sheet = .x
  ))

# Extract the parameter matrices
TH <- extract_param(bayesian_soln$fit)

# Because the intervals used by experts differ, restrict the end-to-end
# Bayesian reconstruction to each expert interval seperately and normalize
# appropriately
expert_recons[c(
  "Haviland",
  "Turner",
  "Culbert"
)] %<>%
  purrr::map(function(x) {
    dens <- calc_interval_densities(x)
    tau <- seq(min(dens[, 1]), max(dens[, 2]), by = 1)
    fMat <- calc_gauss_mix_pdf_mat(TH,
      tau,
      tau_min = min(tau),
      tau_max = max(tau)
    )
    # Make the summary object explicitly to get the right integration limits
    summ <- list(tau = tau, probs = c(.025, .5, .975))
    summ$Qdens <- calc_quantiles(fMat, probs = summ$probs)

    tibble::lst(
      data = x,
      dens,
      tau,
      fMat,
      summ
    )
  })

expert_recons[c(
  "Fry",
  "Santley"
)] %<>%
  purrr::map(function(x) {
    dens <- calc_point_densities(x)
    tau <- seq(min(dens[, 1]), max(dens[, 1]), by = 1)
    fMat <- calc_gauss_mix_pdf_mat(TH,
      tau,
      tau_min = min(tau),
      tau_max = max(tau)
    )
    # Make the summary object explicitly to get the right integration limits
    summ <- list(tau = tau, probs = c(.025, .5, .975))
    summ$Qdens <- calc_quantiles(fMat, probs = summ$probs)

    tibble::lst(
      data = x,
      dens,
      tau,
      fMat,
      summ
    )
  })

# Now that preliminaries are out of the way, make the actual plot comparing our
# reconstruction to previous expert reconstructions
# TODO: consider extending density curve beyond current, expert ranges
file_name <- file.path("outputs","Fig5_tikal_prev_expert_comparison.pdf")
pdf(file_name, width = 6, height = 12)
par(
  mfrow = c(5, 1),
  xaxs = "i", # No padding for x-axis
  yaxs = "i", # No padding for y-axis
  oma = c(4, 2, 2, 2),
  mar = c(0, 4, 0, 0)
)

# Set ranges for plotting
fMax <- expert_recons %>%
  purrr::map_dbl(function(x) {
    max(x$summ$Qdens)
  }) %>%
  max()
tauRange <- expert_recons %>%
  purrr::map(function(x) {
    x$tau
  }) %>%
  unlist() %>%
  range()

# Set locations of tick marks
tauticks <- seq(-750, 1250, 250)
taulabs <- tauticks
taulabs[taulabs == 0] <- "-1/1"
fticks <- seq(0, 3.75, by = .5) / 1000

expert_recons %>%
  purrr::iwalk(function(x, i) {
    plot(NULL,
      xlim = tauRange,
      ylim = c(0, fMax),
      xaxt = "n",
      yaxt = "n",
      ylab = "Density x 1000"
    )
    text(
      x = -500,
      y = 0.0025,
      labels = i,
      cex = 2,
      adj = 0
    )
    axis(side = 1, at = tauticks, labels = taulabs)
    axis(side = 2, at = fticks, labels = fticks * 1000)
    add_interval_density_to_plot(x$dens, col = "black")
    plot_50_percent_quantile(x$summ, add = T, lwd = 3, col = "blue")
    add_shaded_quantiles(x$summ, col = adjustcolor("blue", alpha.f = 0.25))
  })

mtext("Year (AD)", side = 1, line = 2.5)
dev.off()

### FIGURE 4: Peak Population Histograms
# Generate a histogram plot of the peak population date for the Tikal
# reconstruction. Each value used to generate the histogram comes from a
# separate Bayesian sample.
# Extract the dates of the peak value
tpeak <- unlist(lapply(bayesian_summ$summ_list, function(s) {
  s$t_peak
}))

get_hist_breaks <- function(v,dv) {
  # For the input vector v, create histogram breaks with the spacing dv.
  v_min <- floor(min(v)/dv)*dv
  v_max <- ceiling(max(v)/dv)*dv
  v_breaks <- seq(v_min, v_max, by = dv)
  return(v_breaks)
}

# To improve interpretability of the histogram (by showing it on a smaller
# timespan) these samples below AD 630 are placed in a bin at AD 630 and the
# samples above 835 are placed in a bin at 835.
#
# TODO: update the following comment and test with the latest analysis data
# There are 76 such observations, most of which correspond to samples for which the
# Pre-Classic peak is very sharp.
#testthat::expect_equal(
#  sum(tpeak < 600),
#  76
#)

# Create a modified vector with the samples below 630 set to 627.5 and those
# above 835 set to 837.5
tpeak_cutoff_lo <- 630
tpeak_cutoff_hi <- 835
tpeak_modified <- tpeak
tpeak_modified[tpeak_modified < tpeak_cutoff_lo] <- tpeak_cutoff_lo - 2.5
tpeak_modified[tpeak_modified > tpeak_cutoff_hi] <- tpeak_cutoff_hi + 2.5

tpeak_breaks <- get_hist_breaks(tpeak_modified,5)
file_name <- file.path("outputs","Fig6_tikal_peak_population_histogram.pdf")

pdf(file_name, width = 10, height = 6)
  hist(tpeak_modified,
       breaks = tpeak_breaks,
       col = rgb(0, 0, 1, .5),
       freq = F)
  text(tpeak_cutoff_lo-2.5,0.0075,paste0("< AD ",tpeak_cutoff_lo),srt=90,cex=1)
  text(tpeak_cutoff_hi+2.5,0.0075,paste0("> AD ",tpeak_cutoff_hi),srt=90,cex=1)
dev.off()