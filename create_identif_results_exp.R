library(baydem)

source(here::here("bayesian_radiocarbon_functions.R"))

set.seed(183450) # from random.org between 1 and 1,000,000

# Load the calibration curve and call the code to assess equifinality
calib_df <- load_calib_curve("intcal20")
tau_curve <- 1950 - calib_df$year_BP # AD
equi_list <- calc_calib_curve_equif_dates(calib_df)
phi_curve <- calc_calib_curve_frac_modern(calib_df)
equif_data <- assess_calib_curve_equif(calib_df)
canInvert <- equif_data$canInvert
invSpanList <- equif_data$invSpanList


tau_min2 <- 970
tau_max2 <- 1035
ind2 <- (tau_curve >= tau_min2) & (tau_curve <= tau_max2)
phi_min2 <- min(phi_curve[ind2])
phi_max2 <- max(phi_curve[ind2])
phi_vect2 <- unique(phi_curve[ind2])
tau_vect2 <- tau_curve[ind2]

for (kk in 1:length(equi_list)) {
  equi_entry <- equi_list[[kk]]
  if (equi_entry$ind_base %in% which(ind2)) {
    tau_vect2 <- unique(c(tau_vect2, equi_entry$tau_base))
    tau_vect2 <- unique(c(tau_vect2, equi_entry$tau_equi))
  }
}

tau_vect2 <- sort(tau_vect2)

pdf(here::here("FigS1_exp_example.pdf"), width = 20, height = 18)
# There are three rows of plots. The top row has a single, long plot. The next
# two rows rows have three plots each.
layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), 3, 3, byrow = TRUE))
tau_mina <- 700
tau_maxa <- 950
tau_minb <- 825
tau_maxb <- 875
tau_minc <- 900
tau_maxc <- 950
col0 <- "indianred1"
meas_error1 <- 1e-2
meas_error2 <- 1e-3

# Save named pairs of key information for reporting in the manuscript
outlog <- list()
outlog$`tau_mina` <- tau_mina
outlog$`tau_maxa` <- tau_maxa
outlog$`tau_minb` <- tau_minb
outlog$`tau_maxb` <- tau_maxb
outlog$`tau_minc` <- tau_minc
outlog$`tau_maxc` <- tau_maxc
outlog$`meas_error1` <- meas_error1
outlog$`meas_error2` <- meas_error2

kappa <- 8033 # Reference decay rate of carbon-14

# Calculate and save the uncertainty in uncalibrated calendar years for the two
# measurement error settings
meas_error_years1_700 <- kappa * meas_error1 / (exp(-(1950 - 700) / kappa))
meas_error_years2_700 <- kappa * meas_error2 / (exp(-(1950 - 700) / kappa))
meas_error_years1_950 <- kappa * meas_error1 / (exp(-(1950 - 950) / kappa))
meas_error_years2_950 <- kappa * meas_error2 / (exp(-(1950 - 950) / kappa))

outlog$meas_error_years1_700 <- meas_error_years1_700
outlog$meas_error_years2_700 <- meas_error_years2_700
outlog$meas_error_years1_950 <- meas_error_years1_950
outlog$meas_error_years2_950 <- meas_error_years2_950

# The vector of growth rates plot, -4% to 4% per annum
r_vect <- log(1 + seq(-.04, .04, by = .01))

# Show the calibration curve in the top row
vis_calib_curve(700,
                950,
                calib_df,
                xlab="Calendar Date [AD]",
                ylab="Fraction Modern",
                invert_col="gray80")
# Add red rectangles to the calibration curve plot to mark the time periods
# used to generate the fraction modern probability densities
rect(tau_mina, .865, tau_maxa, .866, border = NA, col = col0)
text(mean(c(tau_mina, tau_maxa)), .867, "Span 1", cex = 2)
rect(tau_minb, .855, tau_maxb, .856, border = NA, col = col0)
text(mean(c(tau_minb, tau_maxb)), .8535, "Span 2", cex = 2)
rect(tau_minc, .855, tau_maxc, .856, border = NA, col = col0)
text(mean(c(tau_minc, tau_maxc)), .8535, "Span 3", cex = 2)

# Add plots for the six cases (three timespans by two measurement error
# settings).
#
# If there is an identifiability problem for any cases, an error is thrown by
# visExpEquif.
vis_exp_equif(r_vect, tau_mina, tau_maxa, calib_df, meas_error1, 1000)
vis_exp_equif(r_vect, tau_minb, tau_maxb, calib_df, meas_error1, 1000)
vis_exp_equif(r_vect, tau_minc, tau_maxc, calib_df, meas_error1, 1000)
vis_exp_equif(r_vect, tau_mina, tau_maxa, calib_df, meas_error2, 1000)
vis_exp_equif(r_vect, tau_minb, tau_maxb, calib_df, meas_error2, 1000)
vis_exp_equif(r_vect, tau_minc, tau_maxc, calib_df, meas_error2, 1000)
dev.off()

tibble::tibble(
  Parameter = names(outlog),
  Value = unlist(outlog)
) %T>%
  readr::write_csv(here::here("SuppB_exp.csv"))
