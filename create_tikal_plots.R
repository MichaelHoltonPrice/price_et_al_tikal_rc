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

### FIGURE S3: loo plot
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

### FIGURE 4: Tikal densities
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

### FIGURE 5: Tikal Expert Comparison
calc_interval_densities <- function(recon, tau0, Qdens) {
  # recon (for reconstruction) is a matrix-like object read in with read_excel
  # This function normalizes the reconstruction to integrate to 1 on the
  # interval over which it is defined, and creates a second curve that is the
  # 50% quantile in the regions where the interval density is not defined (
  # normalized appropriately).
  #
  # Qdens is a matrix of quantiles at the calendar dates in tau0

  tau <- as.numeric(unlist(recon[, 1]))
  f   <- as.numeric(unlist(recon[, 2]))
  num_int <- length(tau) / 2 # number of intervals
  flo <- f[2 * (1:num_int) - 1]
  fhi <- f[2 * (1:num_int)]
  if (!all(flo == fhi)) {
    stop("f values are not consistent")
  }
  if (!all(diff(tau)[seq(2, (num_int * 2 - 2), by = 2)] == 0)) {
    stop("tau values are not consistent")
  }
  # Interval durations
  tau_lo <- tau[seq(1, 2 * num_int, by = 2)]
  tau_hi <- tau[seq(2, 2 * num_int, by = 2)]
  int_dur <- tau_hi - tau_lo # interval durations

  # Calculate the normalized density for the expert reconstruction
  f_expert <- flo / sum(flo * int_dur)

  # Create the normalized quantile curve on three intervals:
  # For the plotting, subset the quantile curve into three intervals:
  # (a) to the left     of the expert reconstruction
  # (b) On the interval of the expert reconstruction
  # (c) to the right of    the expert reconstruction
  ind_left   <- which(tau0 <= min(tau))
  ind_middle <- which( (tau0 >= min(tau)) & (max(tau) >= tau0) )
  ind_right  <- which(max(tau) <= tau0)

  # Normalize the quantiles such that the 50% quantile integrates to 1 on the
  # middle interval.
  q50  <- Qdens[2,]
  dtau0 <- tau0[2] - tau0[1] # assume tau0 is evenly spaced
  norm_factor <- sum(q50[ind_middle]) * dtau0
  Qdens <- Qdens / norm_factor

  f_max <- max(f_expert, Qdens)

  # Calculate the normalized density in the
  return(list(tau_lo=tau_lo,
              tau_hi=tau_hi,
              f_expert=f_expert,
              tau0=tau0,
              ind_left=ind_left,
              ind_middle=ind_middle,
              ind_right=ind_right,
              Qdens=Qdens,
              f_max=f_max))
}

add_interval_density_to_plot <- function(ido) {
  # ido stands for interval_density_object
  tau_lo <- ido$tau_lo
  tau_hi <- ido$tau_hi
  f      <- ido$f_expert
  for (i in 1:length(f)) {
    if (i == 1) {
      # For the first entry, plot an initial vertical line from 0
      lines(c(tau_lo[1], tau_lo[1]), c(0, f[1]), lwd=3, col="red")
    } else {
      # For other entries, plot a "regular" initial vertical line
      lines(c(tau_lo[i], tau_lo[i]), c(f[i - 1], f[i]), lwd=3, col="red")
    }
    # For all entries, plot a horizontal line
    lines(c(tau_lo[i], tau_hi[i]), c(f[i], f[i]), lwd=3, col="red")

    # For the final entry, plot a vertical line ending at 0
    if (i == length(f)) {
      lines(c(tau_hi[i], tau_hi[i]), c(f[i], 0), lwd=3, col="red")
    }
  }
}

add_quantiles_to_plot <- function(do) {
  # do stands for density_object (either for an interval or point density).
  #
  # Add the shaded quantiles for the middle region
  tau_middle <- do$tau0[do$ind_middle]
  qlo <- do$Qdens[1,do$ind_middle]
  qhi <- do$Qdens[3,do$ind_middle]
  graphics::polygon(c(tau_middle, rev(tau_middle)),
                    c(qlo, rev(qhi)),
                    border = NA,
                    xlab = NULL,
                    col = adjustcolor("blue", alpha.f = 0.25))

  # If necessary, add the left 50% quantile
  if (length(do$ind_left) > 0) {
    tau_left <- do$tau0[do$ind_left]
    q50_left <- do$Qdens[2,do$ind_left]
    lines(tau_left,q50_left,lwd=3,col="blue")
  }
  # Add the middle 50% quantile
  q50_middle <- do$Qdens[2,do$ind_middle]
  lines(tau_middle,q50_middle,lwd=3,col="blue")

  # If necessary, add the right 50% quantile
  if (length(do$ind_right) > 0) {
    tau_right <- do$tau0[do$ind_right]
    q50_right <- do$Qdens[2,do$ind_right]
    lines(tau_right,q50_right,lwd=3,col="blue")
  }
}

calc_point_densities <- function(recon, tau0, Qdens) {
  # recon (for reconstruction) is a matrix-like object read in with read_excel
  # This function normalizes the reconstruction to integrate to 1 on the
  # interval over which it is defined, and creates a second curve that is the
  # 50% quantile in the regions where the point density is not defined (
  # normalized appropriately).
  #
  # Qdens is a matrix of quantiles at the calendar dates in tau0


  tau <- as.numeric(unlist(recon[, 1]))
  f   <- as.numeric(unlist(recon[, 2]))

  weight_vect <- calc_trapez_weights(tau)
  f_expert <- f / sum(f * weight_vect)

  # Create the normalized quantile curve on three intervals:
  # For the plotting, subset the quantile curve into three intervals:
  # (a) to the left     of the expert reconstruction
  # (b) On the interval of the expert reconstruction
  # (c) to the right of    the expert reconstruction
  ind_left   <- which(tau0 <= min(tau))
  ind_middle <- which( (tau0 >= min(tau)) & (max(tau) >= tau0) )
  ind_right  <- which(max(tau) <= tau0)

  # Normalize the quantiles such that the 50% quantile integrates to 1 on the
  # middle interval.
  q50  <- Qdens[2,]
  dtau0 <- tau0[2] - tau0[1] # assume tau0 is evenly spaced
  norm_factor <- sum(q50[ind_middle]) * dtau0
  Qdens <- Qdens / norm_factor

  f_max <- max(f_expert, Qdens)

  # Calculate the normalized density in the
  return(list(tau=tau,
              f_expert=f_expert,
              tau0=tau0,
              ind_left=ind_left,
              ind_middle=ind_middle,
              ind_right=ind_right,
              Qdens=Qdens,
              f_max=f_max))
}

# Load the repert reconstructions. The result is a named list of tibbles for
# each expert reconstruction.
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

# Extract the baseline calendar dates and 50% quantile
tau0 <- bayesian_summ$tau
q50  <- bayesian_summ$Qdens[2,]

# The various expert reconstructions are defined on different intervals. In
# addition, some provide interval estimates (Haviland, Turner, and Culbert) while
# others provide point estimates (Fry and Stanley). Normalize the various
# curves appropriately, accounting for the difference between interval and
# point estimates. Also, calculate the normalized 50% quantile density on
# three intervals: to the left, middle, and right of the expert recontructions.
# Normalize these quantile curves to integrate to one on the middle interval.
interval_experts <- c("Haviland", "Turner", "Culbert")
interval_densities <- list()
for(n in 1:length(interval_experts)) {
  expert <- interval_experts[n]
  interval_densities[[n]] <- calc_interval_densities(expert_recons[[expert]],
                                                    tau0, bayesian_summ$Qdens)
}

point_experts <- c("Fry", "Santley")
point_densities <- list()
for(n in 1:length(point_experts)) {
  expert <- point_experts[n]
  point_densities[[n]] <- calc_point_densities(expert_recons[[expert]],
                                               tau0, bayesian_summ$Qdens)
}

# Now that preliminaries are out of the way, make the actual plot comparing our
# reconstruction to previous expert reconstructions
file_name <- file.path("outputs","Fig5_tikal_prev_expert_comparison.pdf")
# Set ranges for plotting
f_max <- max(
  unlist(lapply(interval_densities,function(ido){ido$f_max})),
  unlist(lapply(   point_densities,function(pdo){pdo$f_max}))
)

# Set locations of tick marks
tauticks <- seq(-750, 1250, 250)
taulabs <- tauticks
fticks <- seq(0, 3.75, by = .5) / 1000

pdf(file_name, width = 6, height = 12)
  par(
    mfrow = c(5, 1),
    xaxs = "i", # No padding for x-axis
    yaxs = "i", # No padding for y-axis
    oma = c(4, 2, 2, 2),
    mar = c(0, 4, 0, 0)
  )

  for (n in 1:length(interval_densities)) {
    plot(NULL,
      xlim = range(tau0),
      ylim = c(0, f_max),
      xaxt = "n",
      yaxt = "n",
      ylab = "Density"
    )
    add_quantiles_to_plot(interval_densities[[n]])
    add_interval_density_to_plot(interval_densities[[n]])
    text(
      x = -500,
      y = 0.0025,
      labels = interval_experts[n],
      cex = 2,
      adj = 0
    )
    axis(side = 1, at = tauticks, labels = taulabs)
    axis(side = 2, at = fticks, labels = fticks * 1000)
  }

  for (n in 1:length(point_densities)) {
    plot(NULL,
      xlim = range(tau0),
      ylim = c(0, f_max),
      xaxt = "n",
      yaxt = "n",
      ylab = "Density"
    )
    add_quantiles_to_plot(point_densities[[n]])
    lines(point_densities[[n]]$tau,point_densities[[n]]$f_expert, lwd=3,col="red")
    text(
      x = -500,
      y = 0.0025,
      labels = point_experts[n],
      cex = 2,
      adj = 0
    )
    axis(side = 1, at = tauticks, labels = taulabs)
    axis(side = 2, at = fticks, labels = fticks * 1000)
  }
  mtext("Year (AD)", side = 1, line = 2.5)
dev.off()

### FIGURE 6: Peak Population Histograms
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
# timespan) samples below AD 600 are placed in a bin at AD 600 and the samples
# above 810 are placed in a bin at 810. The samples below AD 600 range from
# -84 to 38 and those above 810 (there are two) are 1279 and 1318. These
# occur when the side peaks are especially narrow and the main peak is
# especially wide.

#Create a modified vector with the samples below 600 set to 597.5 and those
# above 810 set to 812.5. Check the counts for these bins.
tpeak_cutoff_lo <- 600
tpeak_cutoff_hi <- 810

ind_lo <- tpeak < tpeak_cutoff_lo
ind_hi <- tpeak > tpeak_cutoff_hi

num_lo <- sum(ind_lo)
num_hi <- sum(ind_hi)

range_lo <- range(tpeak[ind_lo])
range_hi <- range(tpeak[ind_hi])

testthat::expect_equal(
  num_lo,
  32
)
testthat::expect_equal(
  num_hi,
  2
)

tpeak_modified <- tpeak
tpeak_modified[tpeak_modified < tpeak_cutoff_lo] <- tpeak_cutoff_lo - 2.5
tpeak_modified[tpeak_modified > tpeak_cutoff_hi] <- tpeak_cutoff_hi + 2.5

# Write peak information to file
tpeak_summary <- list(tpeak_cutoff_lo=tpeak_cutoff_lo,
                      tpeak_cutoff_hi=tpeak_cutoff_hi,
                      tpeak=tpeak,
                      tpeak_modified=tpeak_modified,
                      num_lo=num_lo,
                      num_hi=num_hi,
                      range_lo=range_lo,
                      range_hi=range_hi)

yaml::write_yaml(tpeak_summary, file.path("outputs","tpeak_summary.yaml"))


tpeak_breaks <- get_hist_breaks(tpeak_modified,5)
file_name <- file.path("outputs","Fig6_tikal_peak_population_histogram.pdf")

pdf(file_name, width = 10, height = 6)
  hist(tpeak_modified,
       breaks = tpeak_breaks,
       col = rgb(0, 0, 1, .5),
       freq = F,
       main=NULL,
       xlab="Calendar Date [AD]")
  text(tpeak_cutoff_lo-2.5,0.0075,paste0("< AD ",tpeak_cutoff_lo),srt=90,cex=1)
  text(tpeak_cutoff_hi+2.5,0.0075,paste0("> AD ",tpeak_cutoff_hi),srt=90,cex=1)
dev.off()
