
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MP.prediction

<!-- badges: start -->
<!-- badges: end -->

## Introduction

The goal of MP.prediction is to predict the milk-to-plasma ratio (MP)
using XGBoost. This function allows any user to use the developed
trained XGBoost model to make predictions. The package was designed to
take the reference name of any drug listed on the ChEMBL database as an
input and return the corresponding MPs. The function within the R
package has three features. The first feature is a scrape of the online
database ChEMBL, where the package extracts data from the application
programming interface (API) to obtain the physicochemical properties
and, if available, the fraction unbound in the plasma. The second
feature then presents the data to the user in the form of an editable
table prior to conducting the MP calculations, so that the user can
verify the values being input into all of the literature-based equations
and the XGBoost model. At this stage, the user can edit any of the
values belonging to each parameter. The final feature is making all the
MP predictions for the literature-based equations and the XGBoost model.

However, the user also has the option to bypass the data scrape and
input their own values for the parameters should the data not currently
be available on ChEMBL or discrepancies between values exist. By simply
inputting “custom” instead of a drug name into the function, the user
will be presented with a blank table to input parameters that they would
like to use. The function will then proceed to predict MPs with all the
literature-based equations and the XGBoost ML model according to their
own preferences.

## Installation

You can install the development version of MP.prediction from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pkmslab/Milk_to_plasma_prediction")
```

## Example

Below are basic examples which show you how to use the function and the
output you will obtain:

``` r
# Load the MP.prediction library
suppressPackageStartupMessages(library(MP.prediction))
```

Here is how you can call the function to determine the MPs for aspirin
and for a series of drugs. Each feature nested within the function will
produce a time stamp. This allows the user to see each feature run in
real-time.

``` r
## Example code
MP_prediction_function("aspirin")
MP_prediction_function(c("dextromethorphan", "acetaminophen", "ibuprofen"))
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
