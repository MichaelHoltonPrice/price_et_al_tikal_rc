library(baydem)
rm(list = ls())

tikal_file <- file.path("outputs","tikal.rds")
if (!file.exists(tikal_file)) {
  stop("Missing tikal.rds")
}
tikal <- readRDS(tikal_file)

all_file <- file.path("outputs","all.rds")
if (!file.exists(all_file)) {
  stop("Missing all.rds")
}
all <- readRDS(all_file)

### FIGURES S3 and S4: loo plots
# Plot the All and Tikal reconstructions together using the best result for each
# (per loo). These are Figures S3 and S4 in the supplement.
pdf(file.path("outputs","FigS3_tikal_loo.pdf"),width=8,height=6)
  plot(tikal$density_model$K,tikal$loo_vect,xlab="K",ylab="loo",main="Tikal")
dev.off()

pdf(file.path("outputs","FigS4_all_loo.pdf"),width=8,height=6)
  plot(all$density_model$K,all$loo_vect,xlab="K",ylab="loo",main="All")
dev.off()

### FIGURE 2: Maya/Tikal densities
# Plot the All and Tikal reconstructions together using the best result for each
# (per loo)
file_name <- file.path("outputs",paste0("Fig3_maya_inference.pdf"))

# TODO: consider fixing the y tick locations
pdf(file_name, width = 10, height = 10)
xat <- seq(-1000, 1800, 200)
xlab <- xat

par(
  mfrow = c(2, 1),
  xaxs = "i", # No padding for x-axis
  yaxs = "i", # No padding for y-axis
  oma = c(4, 2, 2, 2),
  mar = c(0, 4, 0, 0)
)
# First plot of Figure 2 [All]
make_blank_density_plot(all$bayesian_summary,
  xlab = "",
  xaxt = "n",
  xlim = c(-1100, 1900),
  ylim = c(0, 0.002)
)
plot_summed_density(all$bayesian_summary,
                    lwd = 3,
                    add = T,
                    col = "black")
plot_50_percent_quantile(all$bayesian_summary,
                         add = T,
                         lwd = 3,
                         col = "red")
add_shaded_quantiles(all$bayesian_summary,
                     col = adjustcolor("red", alpha.f = 0.25))
# Second plot of Figure 2 [Tikal]
make_blank_density_plot(tikal$bayesian_summary,
  xlab = "",
  xaxt = "n",
  xlim = c(-1100, 1900),
  ylim = c(0, 0.004)
)
plot_summed_density(tikal$bayesian_summary,
                    lwd = 3,
                    add = T,
                    col = "black")
plot_50_percent_quantile(tikal$bayesian_summary,
                         add = T,
                         lwd = 3,
                         col = "blue")
add_shaded_quantiles(tikal$bayesian_summary,
                     col = adjustcolor("blue", alpha.f = 0.25))
axis(side = 1, at = xat, labels = xlab)
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
TH_all <- extract_param(  all$bayesian_solutions[[  all$m_K_best]]$fit)
TH_tik <- extract_param(tikal$bayesian_solutions[[tikal$m_K_best]]$fit)

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
    fMat <- calc_gauss_mix_pdf_mat(TH_tik,
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
    fMat <- calc_gauss_mix_pdf_mat(TH_tik,
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
# TODO: fix cutoff of Bayesian densities
file_name <- file.path("outputs","Fig4_tikal_prev_expert_comparison.pdf")
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

# Generate a histogram plot of the peak population date for both All sites and
# Tikal. The underlying data for the histograms comes from each Bayesian sample.
# Extract the dates of the peak value
tpeak_all <- unlist(lapply(all$bayesian_summary$summ_list, function(s) {
  s$t_peak
}))
tpeak_tik <- unlist(lapply(tikal$bayesian_summary$summ_list, function(s) {
  s$t_peak
}))

get_hist_breaks <- function(v_all,v_tik,dv) {
  # For the input All and Tikal data, create histogram breaks with the spacing
  # dv.
  vmin <- floor(min(v_all,v_tik)/dv)*dv
  vmax <- ceiling(max(v_all,v_tik)/dv)*dv
  vBreaks <- seq(vmin, vmax, by = dv)
  return(vBreaks)
}

# To improve interpretability of the histogram (by showing it on a smaller
# timespan) these samples below AD 600 are placed in a bin at AD 600. There
# are 76 such observations, most of which correspond to samples for which the
# Pre-Classic peak is very sharp.
#testthat::expect_equal(
#  sum(tpeak_tik < 600),
#  76
#)

# Create a modified vector with the samples below 600 set to 597.5
tpeak_cutoff_lo <- 630
tpeak_cutoff_hi <- 835
tpeak_tik_modified <- tpeak_tik
tpeak_tik_modified[tpeak_tik_modified < tpeak_cutoff_lo] <- tpeak_cutoff_lo - 2.5
tpeak_tik_modified[tpeak_tik_modified > tpeak_cutoff_hi] <- tpeak_cutoff_hi + 2.5

tpeakBreaks <- get_hist_breaks(tpeak_all,tpeak_tik_modified,5)
file_name <- file.path("outputs","Fig5_maya_peak_population_histograms.pdf")

pdf(file_name, width = 10, height = 6)

hist(tpeak_all,
     breaks = tpeakBreaks,
     xlab = "Year (AD) of Peak Population",
     ylab = "Density",
     main = NULL,
     col = rgb(1, 0, 0, .5),
     freq = F)
hist(tpeak_tik_modified,
     breaks = tpeakBreaks,
     col = rgb(0, 0, 1, .5),
     add = T,
     freq = F)
text(tpeak_cutoff_lo-2.5,0.0075,paste0("< AD ",tpeak_cutoff_lo),srt=90,cex=1)
text(tpeak_cutoff_hi+2.5,0.0075,paste0("> AD ",tpeak_cutoff_hi),srt=90,cex=1)

dev.off()