This repository contains the data and code for our paper:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2021). *End-to-end Bayesian analysis for summarizing sets of radiocarbon dates*. In review.

# Installation
As with the R package baydem on which this analysis heavily relies, there are two options for installing necessary dependencies and running the analyses:

(1) Install on an existing computer with necessary dependencies
(2) Build a Docker image using the Dockerfile

While Option 1 is more conventional, Option 2 is more reliable since we have, intermittently, encountered difficulties using rstan, especially on Windows. Docker is a cross-platform solution that maximizes scientific reproducibility and, we believe, should be more widely used in scientific analyses.

## Option 1: Install on an existing computer with necessary dependencies
Follow the directions here to install baydem and its dependencies:

https://github.com/eehh-stanford/baydem

Install these additional packages in R:

```R
install.packages(c("rcarbon", "here", "readr", "gtools", 
                   "rlist", "readxl", "R.utils"))
```

## Option 2: Build a Docker image using the baydem Dockerfile
To force all docker material to be (re)downloaded prior to creating the Docker image -- a step you should be certain you want to take -- use: "docker system prune -a"

```console
git clone https://github.com/eehh-stanford/baydem
cd baydem
docker build -t michaelholtonprice/baydem .
```

Start a container. Use the -v tag to mirror a directory for passing files between the host machine and the Docker container. The directory to the left of the semicolon is for the host machine and the directory to the right of the semicolon is for the Docker container. The path for the host machine will need to be modified for your situation.

```console
docker run --name baydem -itv //c/mirrored_baydem_data:/data michaelholtonprice/baydem
```

Clone this github repository and install additional dependencies:

```console
git clone https://github.com/MichaelHoltonPrice/price_et_al_tikal_rc
cd price_et_al_tikal_rc
R
install.packages(c("rcarbon", "here", "readr", "gtools", 
                   "rlist", "readxl", "R.utils"))
```

# Run the analysis code

## Contents
This repository contains the following files (see below for why the first two input files are in the outputs folder):

-   Input data files:
    -   outputs/posterior_from_oxcal_100.txt
    -   outputs/posterior_from_oxcal_1000.txt
    -   MesoRAD-v.1.2_no_locations.xlsx
    -   Tikal_Demography.xlsx
-   R files:
    -   bayesian_radiocarbon_functions.R
    -   create_Fig1.R
    -   create_identif_results_exp.R
    -   create_identif_results_gm.R
    -   do_simulations.R
    -   create_simulation_plots.R
    -   preprocess_mesorad_dataset.R
    -   do_tikal_inference.R
    -   create_tikal_plots.R
-   License and README files:
    -   LICENSE.md
    -   README.md

## Run all analysis scripts
Rather than running each individual script as described in the remainder of this README, all the scripts can be run with the following command:
```R
source("run_all_analysis_scripts.R")
```

## Create the Bayesian inference publication plot
If necessary, set the R working directory to the directory with the files (e.g., using setwd). Then run the following script in R:

```R
source("create_Fig1.R")
```

This will generate the following file:
-   outputs/Fig1_single_date_calibration.pdf

## Create the exponential identifiability results
Run the following script in R:

```R
source("create_identif_results_exp.R")
```

This will generate the following files:
-   outputs/FigS1_exp_example.pdf
-   outputs/SuppB_exp.csv

## Create the Gaussian mixture (gm) identifiability results

Run the following script in R:

```R
source("create_identif_results_gm.R")
```

This will generate the following files:
-    outputs/FigS2_gm_example.pdf
-    outputs/create_identif_results_gm_sink.R

create_identif_results_gm_sink.R contains the print statements from running the script, which indicates whether any identifiability issues were identified (none were).

## Create the simulation results
The simulation results are created in three steps: (1) run "do_simulations.R", (2) use the Oxcal web interface, https://c14.arch.ox.ac.uk/oxcal/, to create the KDE results. (3) Run "create_simulation_plots.R".

More on the Oxcal web interface. To our knowledge, there is no good alternative to using the Oxcal web interface to do the KDE fits. To ensure full reproducibility of our analysis, which we have strived for, we therefore split the simulation script into two scripts. The first script, "do_simulations.R", creates files that are correctly specified Oxcal models (e.g., KDE_input_100.txt) that can be processed using the Oxcal web interface. Furthermore, we use a checksum for each input file to ensure that our pipeline for generating these input files always yields the same input files.

The only output files we track in this repository are the results from using the posterior Oxcal web interface,posterior_from_oxcal_100.txt and posterior_from_oxcal_1000.txt. However, we have archived our final publications results (that is, the full contents of the outputs directory) on the Digital Archaeological Record (tdar).

TODO: add our finals outputs to tdar.

(1) Run the following script in R:

```R
source("do_simulations.R")
```

This will generate the following files:
-    outputs/sim10000.rds
-    outputs/sig_trc_summary.yaml
-    outputs/KDE_input_100.txt
-    outputs/KDE_input_1000.txt
-    outputs/KDE_input_10000.txt
-    outputs/max_lik_fit100.rds
-    outputs/max_lik_fit1000.rds
-    outputs/max_lik_fit10000.rds
-    outputs/bchron_fit100.rds
-    outputs/bchron_fit1000.rds
-    outputs/bchron_fit10000.rds
-    outputs/hp.rds
-    outputs/density_model.rds
-    outputs/bayesian_soln100.rds
-    outputs/bayesian_soln1000.rds
-    outputs/bayesian_soln10000.rds
-    outputs/bayesian_summ100.rds
-    outputs/bayesian_summ1000.rds
-    outputs/bayesian_summ10000.rds

(2) Use the Oxcal web interface to create the following to files using the models in KDE_input_100.txt and KDE_input_1000.txt:

-    outputs/posterior_from_oxcal_100.txt
-    outputs/posterior_from_oxcal_1000.txt

We were unable to successfully run the N=10000 KDE model. Oxcal cannot handle that many observations. The files posterior_input_from_oxcal_100.txt and posterior_input_from_oxcal_1000.txt were obtained by fitting a KDE_Plot model, loading the model in the Oxcal web interface, choosing "Raw data" in the drop-down menu, and copying the Modelled Data: Posterior column.

(3) Run the following script in R:

```R
source("create_simulation_plots.R")
```

This will generate the following files:
-    outputs/Fig2_non_bayesian_fits.pdf
-    outputs/Fig3_bayesian_fits.pdf

## Create the Tikal results
The Tikal analysis requires running three scripts: preprocess_mesorad_dataset.R, do_tikal_inference.R, create_tikal_plots.R. The first script does necessary preprocessing (e.g., applying hygiene rules). The second script calls baydem::sample_theta to do Bayesian inference for truncated Gaussian mixture models with the number of mixtures ranging from K=2 through K=14. The third script creates the publication plots.

Run the following script in R to preprocess the Mesorad data used for the inference in the next step:
```R
source("preprocess_mesorad_dataset.R")
```
This will generate the following files:
-    outputs/filtration_log.yaml
-    outputs/mesorad_hygiene_counts.csv
-    outputs/mesorad_final.csv
-    outputs/site_counts.yaml
-    outputs/maya_hp.rds

Run the following script in R to do the Bayesian inference:

```R
source("do_tikal_inference.R")
```

This will create save files for each value of K (the number of mixture components) from 2 through 14, as well as an .rds file with inputs for the inference:
-    tikal_inference_inputs.rds
-    outputs/tikal_K2.rds
-    ...
-    outputs/tikal_K14.rds

Run the following script in R to create the publication plots:

```R
source("create_tikal_plots.R")
```

This will generate the following files:
-    outputs/loo_summary.yaml
-    outputs/FigS3_tikal_loo.pdf
-    outputs/Fig4_tikal_inference.pdf
-    outputs/Fig5_tikal_prev_expert_comparison.pdf
-    outputs/tpeak_summary.yaml
-    outputs/Fig6_tikal_peak_population_histogram.pdf

## If running in a Docker container, copy results to /data
Run the following script in R to copy the results, plus some additional files (the source code and inputs data files). This script is not included in run_all_analysis_scripts.R).

```R
source("create_final_snapshot.R")
``` 

This will populate the folder /data/final_files with the results, source code, and input data files. We have archived this folder in tDAR (the Digital Archaeological Record) with our final publication results.