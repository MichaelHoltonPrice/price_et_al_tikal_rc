rm(list = ls())

library(baydem)
library(magrittr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# We use locationless MesoRAD data, version 1.2
# Set a random number seed to make the Mesorad data generation reproducible
set.seed(276368)

# Begin generation of mesorad combined data (if necessary)
# Whitespace not updated
if (!file.exists(here::here("mesorad_combined.csv"))) {

  # Read the location-less data
  mesorad <-
    here::here("MesoRAD-v.1.2_no_locations.xlsx") %>%
    readxl::read_excel(sheet = "MesoRAD v 1.2 Dates")

  hygKeep <- c(
    NA,
    "Context and sample info not reported",
    "Duplicated lab number",
    "Lab number not reported",
    "Large measurement error",
    "Large measurement error, origin not reported",
    "Large measurement error, material pretreatment unclear",
    "Large measurement error, no contextual information reported",
    "Large measurement error; no contextual information reported",
    "No contextual information reported"
  )

  hygRej <- c(
    "Acid treatment only, anomalous date",
    "Anomalous date",
    "Conventional 14C yr not reported",
    "Date to early for context",
    "Date too early for context",
    "Date too early for context, reused beam",
    "Date too early for context, experimental dating technique",
    "Date too early for context; Experimental dating technique",
    "Date too late for context",
    "Date too late for context (modern)",
    "Date too late for context, pretreatment/purification unclear",
    "Date too late for context, Pretreatment/Purification unclear",
    "Date too late for context? (Ringle 2012)",
    "Excluded from sequence by investigator",
    "Large measurement error, 1975-76 date",
    "Large measurement error, date too early for context",
    "Large measurement errorr, date too early for context",
    "Large measurement error, date too late for context (modern)",
    "Large measurement error; date too late for context (modern)",
    "Old wood",
    "Outlier ignored by author",
    "outlier ignored by author",
    "Post-bomb (F14C > 1.0)",
    "Rejected by original researchers, 1975-76 date",
    "Reservoir affect unknown",
    "Reservoir correction unclear",
    "Since the Mayas didn't burn their dead, age of cremation should be early Olmeca or Tolteca",
    "Unconventional dating method"
  )

  # Test to make sure all hygiene categories are considered
  unrecognized <-
    mesorad %>%
    dplyr::filter(!(`Chronometric Hygiene/Issues with Dates` %in% c(hygKeep, hygRej))) %$%
    `Chronometric Hygiene/Issues with Dates` %>%
    unique() %>%
    sort()

  if (length(unrecognized) != 0) {
    stop(paste0("Unrecognized hygiene categories: ", unrecognized))
  }

  # Write table of hygiene data
  mesorad %>%
    dplyr::group_by(`Chronometric Hygiene/Issues with Dates`) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    readr::write_csv(here::here("log_mesorad_hygiene_counts.csv"))

  mesorad %<>%
    dplyr::filter(`Chronometric Hygiene/Issues with Dates` %in% hygKeep) %>%
    dplyr::mutate_at(
      .vars = dplyr::vars(
        `Conventional 14C age (BP)`,
        `Error`
      ),
      as.integer
    ) %>%
    dplyr::filter(
      !is.na(`Conventional 14C age (BP)`),
      !is.na(`Error`)
    ) %>%
    dplyr::arrange(Site, `Conventional 14C age (BP)`) %>%
    dplyr::mutate(
      `Duplicate/Replicate` =
        ifelse(!(duplicated(`Duplicate/Replicate`, incomparables = NA) |
          duplicated(`Duplicate/Replicate`, fromLast = TRUE, incomparables = NA)
        ),
        NA,
        `Duplicate/Replicate`
        )
    ) %T>%
    readr::write_csv(here::here("mesorad_filtered.csv")) %>%
    dplyr::select(
      Site,
      `Age BP` = `Conventional 14C age (BP)`,
      `Error` = `Error`,
      `Duplicate/Replicate`
    ) %>%
    dplyr::arrange(Site, `Age BP`)

  # Write a table of counts of dates by site, and the total number of dates
  mesorad %>%
    dplyr::group_by(Site, `Duplicate/Replicate`) %>%
    dplyr::count() %>%
    dplyr::mutate(no_replicates = ifelse(!is.na(`Duplicate/Replicate`), 1, n)) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      `N. dates` = sum(n),
      `N. dates with replicates combined` = sum(no_replicates)
    ) %>%
    dplyr::arrange(-`N. dates`) %T>%
    readr::write_csv(here::here("log_mesorad_filtered_site_counts.csv")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Site) %>%
    dplyr::summarise_all(sum)

  # Combine all duplicate dates
  combined_mesorad <-
    mesorad %>%
    dplyr::filter(!is.na(`Duplicate/Replicate`)) %>%
    dplyr::group_by(Site, `Duplicate/Replicate`) %>%
    dplyr::summarise(Combined = list(
      (function(age_bp, error) {
        samples <-
          ArchaeoChron::combination_Gauss(
            M = age_bp,
            s = error,
            studyPeriodMin = 0,
            studyPeriodMax = 10000,
            numberChains = 8,
            numberSample = 100000,
            thin = 10
          )
        return(
          tibble::tibble(
            `Age BP` = mean(unlist(samples)),
            `Error` = sd(unlist(samples))
          )
        )
      })(`Age BP`, `Error`)
    )) %>%
    dplyr::collect() %>%
    tidyr::unnest(Combined)

  # Now that the combined estimates are available, update the data
  mesorad %<>%
    dplyr::filter(is.na(`Duplicate/Replicate`)) %>%
    dplyr::bind_rows(combined_mesorad) %>%
    dplyr::arrange(Site, `Age BP`) %>%
    dplyr::select(-`Duplicate/Replicate`) %T>%
    readr::write_csv(here::here("mesorad_combined.csv"))
} # End generation of mesorad combined data (if necessary)


mesorad <-
  here::here("mesorad_combined.csv") %>%
  readr::read_csv()

runlog <- list()
runlog$`Total Sample Size` <- nrow(mesorad)

# Only include dates between 2850 and 200 BP (uncalibrated)
mesorad %<>%
  dplyr::filter(
    `Age BP` <= 2850,
    `Age BP` >= 200
  ) %T>%
  readr::write_csv(here::here("mesorad_final.csv"))

runlog$`Reduced Sample Size` <- nrow(mesorad)
runlog$`Tikal Sample Size` <- mesorad %>%
  dplyr::filter(Site == "Tikal") %>%
  nrow()

# Calculate the fraction modern and associated uncertainty for all data
mesorad %<>%
  dplyr::mutate(
    phi_m = exp(-`Age BP` / 8033),
    sig_m = `Error` * phi_m / 8033
  )

# Doing inference involves three steps:
#
# (1) Generate the problem
# (2) Do the Bayesian sampling
# (3) Run some standard analyses
#
# Since the inference may take quite some time, save the result of each of the
# ten runs to file and only do the inference if the file is missing. The
# ten runs are:
#	K	Site(s)
#	2	Tikal
#	4	Tikal
#	6	Tikal
#	8	Tikal
#	10	Tikal
#	2	All
#	4	All
#	6	All
#	8	All
#	10	All

# Specify the baseline hyperparameters. K is changed for different runs
hp0 <-
  list(
    # Class of fit (Gaussian mixture)
    fitType = "gaussmix",
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 500
    alpha_r = (10 - 1) / 500,
    # Minimum calendar date (years AD)
    taumin = -1100,
    # Maximum calendar date (years AD)
    taumax = 1900,
    # Spacing for the measurement matrix (years)
    dtau = 1
    # K is set to 2, 4, 6, 8, and 10 below
  )

# The calibration dataframe
calibDf = baydem::bd_load_calib_curve("intcal20")

# Define the control parameters for the call to Stan. Use 4500 total MCMC
# samples, of which 2000 are warmup samples. Since four chains are used, this
# yields 4*(4500-2000) = 10,000 total samples.
control0 = list(numChains = 4,sampsPerChain = 4500,warmup = 2000)


# To make runs reproducible, use the following random number seeds:
initSeed <- c(361591,403927,688927,  2917,204987,  86685,168649,132214,904995,328517)
stanSeed <- c(807472,264408,406443,778875, 57096, 257668,739622,726159,406443,605221)

# Iterate over number of mixtures
Kvect <- seq(2,10,by=2)
regions <- c('tik','all')

# Create a results list to store the problem, solution, and analysis
results <- list()
for(i in 1:length(Kvect)) {
  K <- Kvect[i]
  # Create the hyperparameter list, setting K
  hp <- hp0
  hp$K <- K
  # Iterate over regions
  for(j in 1:length(regions)) {
    runIndex <- i + (j-1)*length(Kvect)
    reg <- regions[j]
    fileName <- paste0('maya_inference_K',K,'_',reg,'.rds')
    if(!file.exists(fileName)) {
      # Do the inference
      # Create a sub-setting vector, ind (keep everything for all)
      if(reg == 'tik') {
        ind <- which(mesorad$Site == 'Tikal')
      } else {
        ind <- 1:nrow(mesorad)
      }
      # Set random number seeds in control
      control = control0
      control$initSeed <- initSeed[runIndex]
      control$stanSeed <- stanSeed[runIndex]
      # Create the list specifying the problem
      prob <- list(
                phi_m = mesorad$phi_m[ind],
                sig_m = mesorad$sig_m[ind],
                hp = hp,
                calibDf = calibDf,
                control = control
                )

      soln <- baydem::bd_do_inference(prob)
      anal <- baydem::bd_analyze_soln(soln)
      results[[runIndex]] <- list(prob=prob,soln=soln,anal=anal)
      readr::write_rds(results[[i + (j-1)*length(Kvect)]],fileName,compress='gz')
    } else {
      results[[i + (j-1)*length(Kvect)]] <- readr::read_rds(fileName)
    }
  }
}

# Double check the run consistency by checking the random number seeds
for(i in 1:length(Kvect)) {
  K <- Kvect[i]
  for(j in 1:length(regions)) {
    runIndex <- i + (j-1)*length(Kvect)
    reg <- regions[j]
    fileName <- paste0('maya_inference_K',K,'_',reg,'.rds')
    if(results[[runIndex]]$prob$control$initSeed != initSeed[runIndex]) {
      stop('Problem with run consistency for initSeed')
    }
    if(results[[runIndex]]$prob$control$stanSeed != stanSeed[runIndex]) {
      stop('Problem with run consistency for stanSeed')
    }
  }
}

library(ggplot2)

# Extract named "runs" from results list for code readability
out_tik_K2  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==2  && length(x$prob$phi_m)==sum(mesorad$Site=='Tikal')})))]]
out_tik_K4  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==4  && length(x$prob$phi_m)==sum(mesorad$Site=='Tikal')})))]]
out_tik_K6  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==6  && length(x$prob$phi_m)==sum(mesorad$Site=='Tikal')})))]]
out_tik_K8  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==8  && length(x$prob$phi_m)==sum(mesorad$Site=='Tikal')})))]]
out_tik_K10 <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==10 && length(x$prob$phi_m)==sum(mesorad$Site=='Tikal')})))]]
out_all_K2  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==2  && length(x$prob$phi_m)==nrow(mesorad)})))]]
out_all_K4  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==4  && length(x$prob$phi_m)==nrow(mesorad)})))]]
out_all_K6  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==6  && length(x$prob$phi_m)==nrow(mesorad)})))]]
out_all_K8  <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==8  && length(x$prob$phi_m)==nrow(mesorad)})))]]
out_all_K10 <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==10 && length(x$prob$phi_m)==nrow(mesorad)})))]]

### FIGURE 2: Maya/Tikal densities
# Plot the All and Tikal reconstructions together using the K=10 results
fileName <- here::here("Fig3_maya_inference_K10.pdf")

pdf(fileName, width = 10, height = 10)
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
bd_make_blank_density_plot(out_all_K10$anal,
  xlab = "",
  xaxt = "n",
  xlim = c(-1100, 1900),
  ylim = c(0, 0.002)
)
bd_plot_summed_density(out_all_K10$anal,lwd = 3,add = T,col = "black")
bd_plot_50_percent_quantile(out_all_K10$anal, add = T, lwd = 3, col = "red")
bd_add_shaded_quantiles(out_all_K10$anal, col = adjustcolor("red", alpha.f = 0.25))

# Second plot of Figure 2 [Tikal]
bd_make_blank_density_plot(out_tik_K10$anal,
  xlab = "",
  xaxt = "n",
  xlim = c(-1100, 1900),
  ylim = c(0, 0.0035)
)

bd_plot_summed_density(out_tik_K10$anal,lwd = 3,add = T,col = "black")
bd_plot_50_percent_quantile(out_tik_K10$anal, add = T, lwd = 3, col = "blue")
bd_add_shaded_quantiles(out_tik_K10$anal, col = adjustcolor("blue", alpha.f = 0.25))
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
  weightVect <- bd_calc_trapez_weights(tau)
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
  purrr::map(~ readxl::read_excel(here::here("Tikal_Demography.xlsx"),
    sheet = .x
  ))


# Extract the parameter matrices
TH_all <- bd_extract_param(out_all_K10$soln$fit)
TH_tik <- bd_extract_param(out_tik_K10$soln$fit)

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
    fMat <- bd_calc_gauss_mix_pdf_mat(TH_tik,
      tau,
      taumin = min(tau),
      taumax = max(tau)
    )
    # Make the analysis object explicitly to get the right integration limits
    anal <- list(tau = tau, probs = c(.025, .5, .975))
    anal$Qdens <- bd_calc_quantiles(fMat, probs = anal$probs)

    tibble::lst(
      data = x,
      dens,
      tau,
      fMat,
      anal
    )
  })

expert_recons[c(
  "Fry",
  "Santley"
)] %<>%
  purrr::map(function(x) {
    dens <- calc_point_densities(x)
    tau <- seq(min(dens[, 1]), max(dens[, 1]), by = 1)
    fMat <- bd_calc_gauss_mix_pdf_mat(TH_tik,
      tau,
      taumin = min(tau),
      taumax = max(tau)
    )
    # Make the analysis object explicitly to get the right integration limits
    anal <- list(tau = tau, probs = c(.025, .5, .975))
    anal$Qdens <- bd_calc_quantiles(fMat, probs = anal$probs)

    tibble::lst(
      data = x,
      dens,
      tau,
      fMat,
      anal
    )
  })


# Now that preliminaries are out of the way, make the actual plot comparing our
# reconstruction to previous expert reconstructions
fileName <- here::here("Fig4_tikal_prev_expert_comparison.pdf")
pdf(fileName, width = 6, height = 12)
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
    max(x$anal$Qdens)
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
    bd_plot_50_percent_quantile(x$anal, add = T, lwd = 3, col = "blue")
    bd_add_shaded_quantiles(x$anal, col = adjustcolor("blue", alpha.f = 0.25))
  })

mtext("Year (AD)", side = 1, line = 2.5)
dev.off()

### FIGURE 4: Peak Population Histograms

# Using the run with K=10 mixtures, generate a histogram plot of the peak
# population date for both All sites and Tikal. The underlying data for the
# histograms comes from each Bayesian sample.
#

fileName <- here::here("Fig5_maya_histograms.pdf")

pdf(fileName, width = 10, height = 6)

# Extract the dates of the peak value
tpeak_all <- unlist(lapply(out_all_K10$anal$summList, function(s) {
  s$tpeak
}))
tpeak_tik <- unlist(lapply(out_tik_K10$anal$summList, function(s) {
  s$tpeak
}))

get_hist_breaks <- function(v_all,v_tik,dv) {
  # For the input All and Tikal data, create histogram breaks with the spacing
  # dv.
  vmin <- floor(min(v_all,v_tik)/dv)*dv
  vmax <- ceiling(max(v_all,v_tik)/dv)*dv
  vBreaks <- seq(vmin, vmax, by = dv)
  return(vBreaks)
}

# For Tikal, there is exactly one sample for which the peak calendar date is
# ~AD 48 and for which there is a sharp, earlier peak and broader, later peak
# with a slightly higher value of the density for the sharp, earlier peak. In
# addition, there are 46 samples with peaks as low as AD 408. To improve
# interpretability of the histogram (by showing it on a smaller timespan) these
# 46 samples are placed in a single bin between AD 595 and 600, and the fact
# this is done is noted in the graph.

# Create a modified vector with the samples below 600 set to 597.5
tpeak_cutoff <- 600
tpeak_tik_modified <- tpeak_tik
tpeak_tik_modified[tpeak_tik_modified < tpeak_cutoff] <- tpeak_cutoff - 2.5

tpeakBreaks <- get_hist_breaks(tpeak_all,tpeak_tik_modified,5)
hist(tpeak_all, breaks = tpeakBreaks, xlab = "Year (AD) of Peak Population", ylab = "Density", main = NULL, col = rgb(1, 0, 0, .5), freq = F)
hist(tpeak_tik_modified, breaks = tpeakBreaks, col = rgb(0, 0, 1, .5), add = T, freq = F)
text(598.5,0.0075,paste0("< AD ",tpeak_cutoff),srt=90,cex=1)

dev.off()

### FIGURE S3: Summary plots of Maya data for all K

fileName <- here::here("FigS3_maya_inference_Kall.pdf")

pdf(fileName, width = 16, height = 5*length(Kvect))
xat <- seq(-1000, 1800, 200)
xlab <- xat

par(
  mfrow = c(length(Kvect), 2)
)

sitesToPlot <- c("All Sites", "Tikal")
counter <- 0
for (kk in 1:length(Kvect)) {
  for (ss in 1:length(sitesToPlot)) {
    if (sitesToPlot[ss] == "All Sites") {
      plotCol <- "red"
      # Subset of Maya data
      maya_sub <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==Kvect[kk]  && length(x$prob$phi_m)==nrow(mesorad)})))]]
    } else {
      plotCol <- "blue"
      # Subset of Maya data
      maya_sub <- results[[which(unlist(lapply(results,function(x){x$prob$hp$K==Kvect[kk]  && length(x$prob$phi_m)==sum(mesorad$Site=='Tikal')})))]]
    }
    counter <- counter + 1
    plotTitle <- paste0(sitesToPlot[ss], " [K = ", Kvect[kk], "]")

    bd_make_blank_density_plot(maya_sub$anal,
      xlab = "Calendar Date [AD]",
      ylab = "Density",
      xlim = c(-1100, 1900),
      ylim = c(0, 0.004),
      main = plotTitle
    )
    bd_plot_50_percent_quantile(maya_sub$anal, add = T, lwd = 3, col = plotCol)
    bd_add_shaded_quantiles(maya_sub$anal, col = adjustcolor(plotCol, alpha.f = 0.25))
  }
}
dev.off()
