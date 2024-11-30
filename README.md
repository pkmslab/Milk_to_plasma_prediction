
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MP.prediction

<!-- badges: start -->
<!-- badges: end -->

The goal of MP.prediction is to predict the milk-to-plasma ratio (MP).
The function utilizes the name of a drug, searches the CHEMBL database
for parameters, and returns MPs for an array of different literature
methods (Phase distribution, Koshimichi et al., Logarithmically
transformed phase distribution, and Meskin and Lien) and also for an
extreme gradient boosting (XGBoost machine) learning model.
Alternatively, the user can input “custom”, create their own dataset,
and run the methodswith the desired parameters.

## Installation

You can install the development version of MP.prediction from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("akashp21/Milk_to_plasma_prediction")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MP.prediction)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
