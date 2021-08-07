# Throw an error if /data/final_snapshot already exists.
snapshot_dir <- file.path("/data","final_snapshot")

if (dir.exists(snapshot_dir)) {
  stop("/data/final_snapshot already exists")
}

# Create the snapshot directory
dir.create(snapshot_dir)

# Copy the results folder
copyDirectory("outputs", snapshot_dir)

# Copy the results folder
copyDirectory("outputs", snapshot_dir)

# Copy the source code and input data files
files_to_copy <- c("MesoRAD-v.1.2_no_locations.xlsx",
                   "Tikal_Demography.xlsx",
                   "bayesian_radiocarbon_functions.R",
                   "create_Fig1.R",
                   "create_identif_results_exp.R",
                   "create_identif_results_gm.R",
                   "do_simulations.R",
                   "create_simulation_plots.R",
                   "preprocess_mesorad_dataset.R",
                   "do_tikal_inference.R",
                   "create_tikal_plots.R")
file.copy(file_to_copy,snapshot_dir)