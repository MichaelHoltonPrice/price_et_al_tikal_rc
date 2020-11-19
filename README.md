
<!-- README.md is generated from README.Rmd. Please edit that file -->
price\_et\_al\_tikal\_rc
========================

This repository contains the data and code for our paper:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2020). *End-to-end Bayesian analysis of 14C dates as an alternative to summed probability densities*. In review.

<!-- Our pre-print is online here: -->
<!-- > Authors, (YYYY). _End-to-end Bayesian analysis of 14C dates as an alternative to summed probability densities_. Name of journal/book, Accessed 18 Nov 2020. Online at <https://doi.org/xxx/xxx> -->
### How to cite

Please cite this compendium as:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2020). *Compendium of R code and data for End-to-end Bayesian analysis of 14C dates as an alternative to summed probability densities*. Accessed 18 Nov 2020.

Getting the code and input data
-------------------------------

The easiest way to get the code and input data is to clone this github repository. For example, enter the following at the command line to clone the repository, enter the newly-created directory, and list its contents:

``` console
git clone https://github.com/MichaelHoltonPrice/price_et_al_tikal_rc
cd price_et_al_tikal_rc
ls
```

Requirements
------------

This research compendium has been developed using the statistical programming language R. To work with the compendium, you will need the [R software](https://cloud.r-project.org/) itself and optionally [RStudio Desktop](https://rstudio.com/products/rstudio/download/).

To install the R packages that the code depends on enter the following in *R*:

``` r
# If necessary, install devtools
install.packages("devtools")

# Install baydem. On Linux, this may first require running the following in a terminal window:
# sudo apt-get install libv8-dev
devtools::install_github('eehh-stanford/baydem')

# If necessary, install additional package dependencies:
install.packages('Bchron')
install.packages('dplyr')
install.packages('here')
install.packages('magrittr')
install.packages('purrr')
install.packages('readr')
install.packages('readxl')
install.packages('tibble')
```

When baydem is installed, so is Rstan, which may have further depenencies. If so, see:

<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>

Contents
--------

-   Input data files:
    -   MesoRAD-v.1.1\_FINAL\_no\_locations.xlsx
    -   Tikal\_Demography.xlsx
-   R files:
    -   bayesian\_radiocarbon\_functions.R
    -   create\_identif\_results\_exp.R
    -   create\_identif\_results\_gm.R
    -   do\_maya\_inference.R
    -   do\_sim\_inference.R
-   License and README files:
    -   LICENSE.md
    -   README.md
    -   README.Rmd

Running the analyses
====================

Simulation
----------

If necessary, set the R working directory to the directory with the files (e.g., using setwd). Then type the following in R:

``` r
source('do_sim_inference.R')
```

This should take a few hours to finish. Once complete, the following new files are created:

-   sim\_inference.rds
-   Fig1\_sim\_inference.pdf

sim\_inference.rds stores the results of the Bayesian inference for N=10, 100, 1000,and 10000 simulated radiocarbon samples and two choices for the paramaterization of the prior (a total of eight cases). Fig1\_sim\_inference.pdf is Figure 1 in the article.

Maya results
------------

If necessary, set the R working directory to the directory with the files (e.g., using setwd). Then type the following in R:

``` r
source('do_maya_inference.R')
```

This should take about two days to finish. Once complete, the following new files are created:

-   Mesorad files
    -   log\_mesorad\_hygiene\_counts.csv
    -   mesorad\_filtered.csv
    -   log\_mesorad\_filtered\_site\_counts.csv
    -   mesorad\_combined.csv
    -   mesorad\_final.csv
-   Inference results
    -   maya\_inference\_K2\_tik.rds
    -   maya\_inference\_K4\_tik.rds
    -   maya\_inference\_K6\_tik.rds
    -   maya\_inference\_K8\_tik.rds
    -   maya\_inference\_K10\_tik.rds
    -   maya\_inference\_K2\_all.rds
    -   maya\_inference\_K4\_all.rds
    -   maya\_inference\_K6\_all.rds
    -   maya\_inference\_K8\_all.rds
    -   maya\_inference\_K10\_all.rds
-   Figures and count data
    -   Fig2\_maya\_inference\_K10.pdf
    -   Fig3\_tikal\_prev\_expert\_comparison.pdf
    -   Fig4\_maya\_histograms.pdf
    -   FigS3\_maya\_inference\_Kall.pdf
    -   FigS4\_maya\_inference\_K2\_and\_K10\_with\_rc\_curve.pdf
    -   supp\_count\_data.csv

The first set of files provides the Mesorad data at various stages of processing. The second set of files (maya\_inference\_K) provides the inference results for Tikal/All sites with K=2 to 10 numbers of mixtures components. The third set of files is the figures and count data for the main manuscript and supplement.

Identifiability results (supplement)
------------------------------------

If necessary, set the R working directory to the directory with the files (e.g., using setwd). There are two scripts to run, one for the exponential example and one for the Gaussian mixture example. For the exponential example type:

``` r
source('create_identif_results_exp.R')
```

This should take only a few minutes to run. Once complete, the following new files are created:

-   FigS1\_exp\_example.pdf
-   SuppB\_exp.csv

FigS1\_exp\_example.pdf is Figure S1 in the supplement. SuppB\_exp.csv provides summary information for the example, notably the mean measurement uncertainties reported in the supplement.

For the Gaussian mixture example type:

``` r
source('create_identif_results_gm.R')
```

This will take a few hours to run. Once complete, the following new file is created:

-   FigS2\_gm\_example.pdf

FigS2\_exp\_example.pdf is Figure S2 in the supplement. Information about the identifiability checks is printed out in R as the script runs.
