rm(list = ls())

library(baydem)
library(magrittr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# If the "filtered" mesorad file does not exist, create it. Otherwise, load it.
mesorad_file <- file.path("outputs", "mesorad_final.csv")
if (!file.exists(mesorad_file)) {
  # Read the raw Mesorad v 1.2 data (location information has been removed)
  mesorad <- readxl::read_excel("MesoRAD-v.1.2_no_locations.xlsx",
                                sheet = "MesoRAD v 1.2 Dates")
  num_raw_observations <- nrow(mesorad)

  # Mesorad observatoins contain hygiene information. Certain hygiene
  # categories are acceptable for the Bayesian inference (hyg_keep), whereas
  # others should to be rejected for it (hyg_rej).
  hyg_keep <- c(
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

  hyg_rej <- c(
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
    dplyr::filter(!(`Chronometric Hygiene/Issues with Dates` %in% c(hyg_keep, hyg_rej))) %$%
    `Chronometric Hygiene/Issues with Dates` %>%
    unique() %>%
    sort()

  if (length(unrecognized) != 0) {
    stop(paste0("Unrecognized hygiene categories: ", unrecognized))
  }

  # (a) First filtration: observations with no 14C age or error
  keep <- rep(TRUE, nrow(mesorad))
  for (n in 1:nrow(mesorad)) {
    rc_year       <- mesorad[n,"Conventional 14C age (BP)"]
    rc_year_error <- mesorad[n,"Error"]
    keep[n] <- !is.na(rc_year) && !is.na(rc_year_error)
  }
  num_filtered_for_no_data <- sum(!keep)
  mesorad <- mesorad[keep,]

  # Write table of hygiene data, which are summarized in the publication
  mesorad %>%
    dplyr::group_by(`Chronometric Hygiene/Issues with Dates`) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    readr::write_csv(file.path("outputs","mesorad_hygiene_counts.csv"))

  # (b) Second filtration: observations with bad hygiene
  keep <- rep(TRUE, nrow(mesorad))
  for (n in 1:nrow(mesorad)) {
    hyg <- mesorad[n,"Chronometric Hygiene/Issues with Dates"]
    bad_hygiene <- FALSE
    if (!is.na(hyg)) {
      # Blank hygiene entries are NA
      if (hyg %in% hyg_rej) {
        # Do we reject for this hygiene category?
        bad_hygiene <- TRUE
      }
    }
    keep[n] <- !bad_hygiene
  }

  num_filtered_for_hygiene <- sum(!keep)
  mesorad <- mesorad[keep,]

  # Cast the tibble to a regular R data frame
  mesorad <- as.data.frame(mesorad)

  # (b) Duplicates and Replicates. For each Duplicate/Replicate group, choose
  # a random entry to represent the roup. Use a random number seed for
  # reproducibility.
  set.seed(276368)
  groups <- unique(na.omit(mesorad[,"Duplicate/Replicate"]))
  keep <- rep(TRUE,nrow(mesorad))

  for(group in groups) {
    group_ind <- which(group == mesorad[,"Duplicate/Replicate"])

    if(length(group_ind) == 0) {
      stop("Groups should have at least one entry")
    }
    if (length(group_ind) > 1) {
      # Randomly remove one entry
      ind_to_keep <- sample.int(length(group_ind), 1, replace=F)
      group_ind <- group_ind[-ind_to_keep]
      # Set the remaining observations to FALSE in keep
      keep[group_ind] <- FALSE
    }
  }
  num_filtered_for_dup_rep <- sum(!keep)
  mesorad <- mesorad[keep,]

  # Require observations to be between 2850 and 200 BP (uncalibrated)
  mesorad[,"Conventional 14C age (BP)"] <- as.numeric(mesorad[,"Conventional 14C age (BP)"])
  keep <- (mesorad[,"Conventional 14C age (BP)"] <= 2850) &
          (mesorad[,"Conventional 14C age (BP)"] >=  200)

  num_filtered_for_year <- sum(!keep)
  mesorad <- mesorad[keep,]

  filtration_log <- list(num_raw_observations=num_raw_observations,
                         num_filtered_for_no_data=num_filtered_for_no_data,
                         num_filtered_for_hygiene=num_filtered_for_hygiene,
                         num_filtered_for_dup_rep=num_filtered_for_dup_rep,
                         num_filtered_for_year=num_filtered_for_year)

  rlist::list.save(filtration_log,
                   file.path("outputs", "filtration_log.yaml"))

  write.csv(mesorad, mesorad_file, row.names=FALSE)
} else {
  # The "filtered" mesorad file already exists. Load it.
  mesorad <- read.csv(mesorad_file, stringsAsFactors=FALSE, check.names=FALSE)
}

# Before proceeding, check that the sums of the following two columns have not
# changes: "Conventional 14C age (BP)" and "Error". The checksum for the output
# .csv file is not necessarily consistent across setups.
testthat::expect_equal(
  sum(mesorad[,"Conventional 14C age (BP)"]),
  2126606
)

testthat::expect_equal(
  sum(mesorad[,"Error"]),
  58821
)
# Create a vector of site counts, sorted by size. Make it a list so it can be
# easily saved in a yaml file
site_counts_table <- rev(sort(table(mesorad$Site)))
site_counts <- list()
for (n in 1:length(site_counts_table)) {
  site_counts[[n]] <- as.numeric(site_counts_table[n])
}
names(site_counts) <- names(site_counts_table)
rlist::list.save(site_counts,
                 file.path("outputs", "site_counts.yaml"))


# Add columns to mesorad with the uncalibrated years BP and corresponding error
# using the baydem naming conventions.
mesorad$trc_m <- mesorad[,"Conventional 14C age (BP)"]
mesorad$sig_trc_m <- mesorad[,"Error"]

# Calculate the fraction modern and associated uncertainty for all data
mesorad$phi_m <- exp(-mesorad$trc_m/8033)
mesorad$sig_m <- mesorad$sig_trc_m * mesorad$phi_m / 8033

data_dir <- "outputs"

hp <-
  list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100
    alpha_r = (10 - 1) / 100,
    # Spacing for the measurement matrix (years)
    dtau = 1
  )

saveRDS(hp, file.path("outputs","maya_hp.rds"))

# If necessary, do the Bayesian inference for the Tikal observations
tikal_file <- file.path(data_dir, "tikal.rds")
if (!file.exists(tikal_file)) {
  set_rc_meas(data_dir,
              "tikal",
              mesorad[mesorad$Site == "Tikal",])
  set_density_model(data_dir, "tikal", list(type="trunc_gauss_mix",
                                            tau_min=-1100,
                                            tau_max=1900,
                                            K=2:6))
  set_calib_curve(data_dir, "tikal", "intcal20")
  do_bayesian_inference(data_dir,
                      "tikal",
                      hp,
                      input_seed=433653,
                      control=list(num_chains = 4,
                                   samps_per_chain = 4500,
                                   warmup = 2000)
                     )
}

# If necessary, do the Bayesian inference for All observations
all_file <- file.path(data_dir, "all.rds")
if (!file.exists(all_file)) {
  set_rc_meas(data_dir,
              "all",
              mesorad)

  set_density_model(data_dir, "all", list(type="trunc_gauss_mix",
                                          tau_min=-1100,
                                          tau_max=1900,
                                          K=2:6))

  set_calib_curve(data_dir, "all", "intcal20")

  do_bayesian_inference(data_dir,
                        "all",
                        hp,
                        input_seed=75539,
                        control=list(num_chains = 4,
                                     samps_per_chain = 4500,
                                     warmup = 2000)
                       )
}