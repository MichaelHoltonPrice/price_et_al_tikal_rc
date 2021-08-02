This repository contains the data and code for our paper:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2021). *End-to-end Bayesian analysis for summarizing sets of radiocarbon dates*. In review.

<!-- Our pre-print is online here: -->
<!-- > Authors, (YYYY). _End-to-end Bayesian analysis for summarizing sets of radiocarbon dates_. Name of journal/book, Accessed 01 Dec 2021. Online at <https://doi.org/xxx/xxx> -->
### How to cite

Please cite this compendium as:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2021). *Compendium of R code and data for End-to-end Bayesian analysis for summarizing sets of radiocarbon dates*. Accessed 01 Dec 2021.

# Installation
As with the R package baydem on which this analysis heavily relies, there are two options for installing necessary dependencies and running the analyses:

(1) Install on an existing computer with necessary dependencies
(2) Build a Docker image using the Dockerfile

## Option 1: Install on an existing computer with necessary dependencies
Follow the directions here to install baydem and its dependencies:

https://github.com/eehh-stanford/baydem

Install these additional packages in R:

```R
install.packages("rcarbon, here, readr, gtools")
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
install.packages("rcarbon, here, readr, gtools")
```

# Run the analysis code

## Contents
This repository contains the following files:

-   Input data files:
    -   posterior_from_oxcal_100.txt
    -   posterior_from_oxcal_1000.txt
    -   MesoRAD-v.1.2\_no\_locations.xlsx
    -   Tikal\_Demography.xlsx
-   R files:
    -   bayesian\_radiocarbon\_functions.R
    -   make\_bayesian\_radiocarbon\_illustrations.R
    -   create\_identif\_results\_exp.R
    -   create\_identif\_results\_gm.R
    -   do_simulations.R
    -   do\_maya\_inference.R
-   License and README files:
    -   LICENSE.md
    -   README.md

## Create the Bayesian inference publication plot
If necessary, set the R working directory to the directory with the files (e.g., using setwd). Then run the following script in R:

```R
source("make_bayesian_radiocarbon_illustrations.R")
```

This will generate the following file:
-   Fig1\_single\_date\_calibration.pdf

## Create the exponential identifiability results
Run the following script in R:

```R
source("create_identif_results_exp.R")
```

This will generate the following files:
-   FigS1\_exp\_example.pdf
-   SuppB\_exp.csv

## Create the Gaussian mixture (gm) identifiability results

Run the following script in R:

```R
source("create_identif_results_gm.R")
```

This will generate the following file:
-    FigS2\_gm\_example.pdf

## Create the simulation results

Run the following script in R:

```R
source("do_simulations.R")
```

This will generate the following files:
-    sig\_trc\_summary.yaml
-    max\_lik\_fit100.rds
     max\_lik\_fit1000.rds
-    max\_lik\_fit10000.rds
-    KDE\_input\_100.txt
-    KDE\_input\_1000.txt
-    KDE\_input\_10000.txt
-    Fig2\_non\_bayesian\_fits.pdf
-    bayesian\_soln100.rds
-    bayesian\_soln1000.rds
-    bayesian\_soln10000.rds
-    bayesian\_summ100.rds
-    bayesian\_summ1000.rds
-    bayesian\_summ10000.rds
-    Fig3\_bayesian\_fits.pdf

## Create the Maya results

Run the following script in R:

```R
source("do_Maya_inference.R")
```

TODO: add Maya files

TODO: add and describe script to copy all results files to /data if using a Docker container