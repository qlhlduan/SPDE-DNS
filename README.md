# SPDE-Term Structure

This repository contains the code for the paper "Forecasting the Term Structure of Interest Rates with SPDE approach" by Qihao Duan, Alexandre B. Simas, David Bolin, and Raphael Huser.

## Overview

The paper proposes an extension of the Dynamic Nelson-Siegel (DNS) model for forecasting the term structure of interest rates by incorporating a stochastic partial differential equation (SPDE) component to model the residual dependence across time and maturity. The model is implemented using the R-INLA and inlabru packages.

## Repository Structure

- `data_processing.R`: Downloads and processes the Fama-Bliss US Treasury zero-coupon bond yields.
- `models_DNS.R`: Defines the DNS model and the residual-adjusted models (stationary, non-stationary, anisotropic, spatio-temporal) and the joint estimation of lambda.
- `model_fitting.R`: Fits the models to the data.
- `forecasting.R`: Defines functions for out-of-sample forecasting and evaluation.
- `economic_value.R`: Implements the economic value analysis (portfolio optimization and performance fee).
- `requirements.R`: Lists the required R packages and installs them if missing.

## Data

The data is the Fama-Bliss US Treasury zero-coupon bond yields from January 1985 to December 2000. The data is downloaded from [https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt](https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt).

## Usage

1. Run `requirements.R` to install the required packages.
2. Run `data_processing.R` to download and process the data.
3. Run `model_fitting.R` to fit the models (this may take a long time).
4. Use the functions in `forecasting.R` and `economic_value.R` to reproduce the forecasting and economic value results.

## Note

This code is provided as is and may require adjustments to run on your system. The models are complex and computationally intensive.

## References

- Diebold, F. X., & Li, C. (2006). Forecasting the term structure of government bond yields. Journal of econometrics, 130(2), 337-364.
- Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(4), 423-498.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the royal statistical society: Series b (statistical methodology), 71(2), 319-392.
