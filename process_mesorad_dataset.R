rm(list = ls())
library(magrittr)

# Load the Mesorad dataset and apply a sequence of "filtrations". The
# filtrations are:
#
# (a) Remove observations with no 14C age or error
# (b) Remove observations with bad hygiene
# (c) For duplicates and replicates, randomly select a single observation to
#     represent the group.
#
# This yields a file "outputs/mesorad_final.csv" that contains the filtered
# dataset. In addition, the following three log files are created that
# summarize the filtration:
#
# outputs/filtration_log.yaml
# outputs/mesorad_hygiene_counts.csv
# outputs/site_counts.yaml
#
# While the entire Mesorad dataset is processed (i.e., all sites are included
# in mesorad_final), only the Tikal observations are used in our publication.

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

  # (c) Duplicates and Replicates. For each Duplicate/Replicate group, choose
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

# Check that the sums of the following two columns have not changes:
# "Conventional 14C age (BP)" and "Error" (processing should be deterministic
# since a random number seed is used).
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