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
                 "rlist", "readxl", "tools"))
```

## Option 2: Build a Docker image using the baydem Dockerfile

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
                 "rlist", "readxl", "tools"))
```

# Run the analysis code

## Contents
This repository contains the following files (see below for why the first two input files are in the outputs folder):

-   Input data files:
    -   outputs/posterior_from_oxcal_100.txt
    -   outputs/posterior_from_oxcal_1000.txt
    -   MesoRAD-v.1.2\_no_locations.xlsx
    -   Tikal_Demography.xlsx
-   R files:
    -   bayesian_radiocarbon_functions.R
    -   create_Fig1.R
    -   create_identif_results_exp.R
    -   create_identif_results_gm.R
    -   do_simulations.R
    -   create_simulation_plots.R
    -   do_maya_inference.R
-   License and README files:
    -   LICENSE.md
    -   README.md

## Create the Bayesian inference publication plot
If necessary, set the R working directory to the directory with the files (e.g., using setwd). Then run the following script in R:

```R
source("create_Fig1.R")
```

This will generate the following file:
-   outputs/Fig1\_single\_date\_calibration.pdf

## Create the exponential identifiability results
Run the following script in R:

```R
source("create_identif_results_exp.R")
```

This will generate the following files:
-   outputs/FigS1\_exp\_example.pdf
-   outputs/SuppB\_exp.csv

## Create the Gaussian mixture (gm) identifiability results

Run the following script in R:

```R
source("create_identif_results_gm.R")
```

This will generate the following files:
-    outputs/FigS2_gm_example.pdf
-    outputs/create_identif_results_gm_sink.R

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

We were unable to successfully run the N=10000 KDE model. Oxcal cannot handle that many observations.
TODO: explain precisely where the posterior data come from in the Oxcal results.

(3) Run the following script in R:

```R
source("create_simulation_plots.R")
```

This will generate the following files:
-    outputs/Fig2_non_bayesian_fits.pdf
-    outputs/Fig3_bayesian_fits.pdf

## Create the Maya results

Run the following script in R to do the Bayesian inference (the plots are created using the ensuing script):

```R
source("do_Maya_inference.R")
```

This will generate the following files:
-    outputs/filtration_log.yaml
-    outputs/mesorad_hygiene_counts.csv
-    outputs/mesorad_final.csv
-    outputs/site_counts.yaml
-    outputs/maya_hp.rds
-    outputs/tikal.rds
-    outputs/all.rds

TODO: add plotting script once it is finalized

## If running in a Docker container, copy results files
If necessary, exit R:
```R
q()
y
```

Copy the outputs folder to the mirrored /data directory

```console
cp outputs /data
```