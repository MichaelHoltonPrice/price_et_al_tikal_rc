# Create a vector of files to copy
file_names <- c()

# make_bayesian_radiocarbon_illustrations.R 
file_names <- c(file_names, "Fig1_single_date_calibration.pdf")

# create_identif_results_exp.R
file_names <- c(file_names, "FigS1_exp_example.pdf",
                            "SuppB_exp.csv")

# create_identif_results_gm.R
file_names <- c(file_names, "FigS1_exp_example.pdf",
                            "SuppB_exp.csv",
                            "sig_trc_summary.yaml",
                            "max_lik_fit100.rds",
                            "max_lik_fit1000.rds",
                            "max_lik_fit10000.rds",
                            "KDE_input_100.txt",
                            "KDE_input_1000.txt",
                            "KDE_input_10000.txt",
                            "Fig2_non_bayesian_fits.pdf",
                            "bayesian_soln100.rds",
                            "bayesian_soln1000.rds",
                            "bayesian_soln10000.rds",
                            "bayesian_summ100.rds",
                            "bayesian_summ1000.rds",
                            "bayesian_summ10000.rds",
                            "Fig3_bayesian_fits.pdf")

for (file_name in file_names) {
  if (!file.exists(file_name)) {
    stop(paste0(file_name, " does not exist"))
  }
}

file.copy(file_names, "/data")
