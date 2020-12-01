
<!-- README.md is generated from README.Rmd. Please edit that file -->
price\_et\_al\_tikal\_rc
========================

This repository contains the data and code for our paper:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2020). *End-to-end Bayesian analysis for summarizing sets of radiocarbon dates*. In review.

<!-- Our pre-print is online here: -->
<!-- > Authors, (YYYY). _End-to-end Bayesian analysis for summarizing sets of radiocarbon dates_. Name of journal/book, Accessed 01 Dec 2020. Online at <https://doi.org/xxx/xxx> -->
### How to cite

Please cite this compendium as:

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert, and J.H. Jones, (2020). *Compendium of R code and data for End-to-end Bayesian analysis for summarizing sets of radiocarbon dates*. Accessed 01 Dec 2020.

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
install.packages('rcarbon')
install.packages('tibble')
```

When baydem is installed, so is Rstan, which may have further depenencies. If so, see:

<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>

Contents
--------

-   Input data files:
    -   MesoRAD-v.1.2\_no\_locations.xlsx
    -   Tikal\_Demography.xlsx
-   R files:
    -   bayesian\_radiocarbon\_functions.R
    -   create\_identif\_results\_exp.R
    -   create\_identif\_results\_gm.R
    -   do\_maya\_inference.R
    -   do\_sim\_inference.R
    -   make\_bayesian\_radiocarbon\_illustrations.R
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

This should take a few hours to finish. Once complete, the following ten files are created:

-   Fig1\_sim\_inference.pdf
-   sim\_data.rds
-   sim\_soln\_100.rds
-   sim\_soln\_1000.rds
-   sim\_soln\_10000.rds
-   sim\_anal\_100.rds
-   sim\_anal\_1000.rds
-   sim\_anal\_10000.rds
-   sim\_bc\_100.rds
-   sim\_bc\_1000.rds
-   sim\_bc\_10000.rds

Fig1\_sim\_inference.pdf is for Figure 1 in the article. sim\_data.rds stores simulated samples. For N=100, N=1000, and N=10000, sim\_soln\*, sim\_anal\*, and sim\_bc\* store, respectively, the results (samples) of Bayesian inference, analyses based on the Bayesian inference, and results of the BchronDensityFast fits.

Maya results
------------

If necessary, set the R working directory to the directory with the files (e.g., using setwd). Then type the following in R:

``` r
source('do_maya_inference.R')
```

This may take a day or more to finish. Once complete, the following new files are created:

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
    -   Fig3\_maya\_inference\_K10.pdf
    -   Fig4\_tikal\_prev\_expert\_comparison.pdf
    -   Fig5\_maya\_histograms.pdf
    -   FigS3\_maya\_inference\_Kall.pdf

The first set of files provides the Mesorad data at various stages of processing. The second set of files (maya\_inference\_K) provides the inference results for Tikal/All sites with K=2 to 10 numbers of mixtures components. The third set of files is the figures the main manuscript and supplement.

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

Bayesian radiocarbon illustrations
----------------------------------

If necessary, set the R working directory to the directory with the files (e.g., using setwd). Type:

``` r
source('make_bayesian_radiocarbon_illustrations.R')
```

This should take only a few seconds to run. Once complete, the following new file is created:

-   FigS1\_exp\_example.pdf

Fig1\_single\_date\_calibration.pdf is Figure 1 in the main text.
