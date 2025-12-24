# SPDE-Term Structure

This repository contains the code for the paper "Forecasting the Term Structure of Interest Rates with SPDE approach" by Qihao Duan, Alexandre B. Simas, David Bolin, and Raphaël Huser.

## Overview

The paper proposes an extension of the Dynamic Nelson-Siegel (DNS) model for forecasting the term structure of interest rates by incorporating a stochastic partial differential equation (SPDE) component to model the residual dependence across time and maturity. The model is implemented using the R-INLA and inlabru packages.

## Repository Structure

- `data_processing.R`: Downloads and processes Fama-Bliss US Treasury zero-coupon bond yields from CRSP.
- `nelson.R`: Implements the Dynamic Nelson-Siegel (DNS) model specification and state-space representation.
- `model_fitting.R`: Estimates DNS models with fixed lambda parameter using Kalman filtering.
- `model_fitting_joint.R`: Jointly estimates lambda parameter along with other model parameters using two alternative prior specifications.
- `prediction.R`: Generates out-of-sample yield curve forecasts and computes prediction errors.
- `crossval_metrics.R`: Calculates cross-validation metrics for model comparison.
- `economic_value.R`: Performs economic value analysis through portfolio optimization and calculates performance fees.
- `requirements.R`: Checks for and installs all required R packages and dependencies.
- `results/`: Contains saved model estimates, prediction outputs from fitted models.

## Data

The data is the Fama-Bliss US Treasury zero-coupon bond yields from January 1985 to December 2000. The data is downloaded from [https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt](https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt).

## Usage

1. Install dependencies by running requirements.R.
2. Execute data_processing.R to download and prepare the yield curve data.
3. Fit the models by running model_fitting.R and model_fitting_joint.R (note: estimation may take considerable time depending on your system).
4. Generate predictions using prediction.R and evaluate model performance with crossval_metrics.R.
5. Reproduce the economic value analysis by executing the functions in economic_value.R.

## Note

This code is provided as is and may require adjustments to run on your system. The models are complex and computationally intensive.

## References

- Diebold, F. X., & Li, C. (2006). Forecasting the term structure of government bond yields. Journal of econometrics, 130(2), 337-364.
- Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(4), 423-498.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the royal statistical society: Series b (statistical methodology), 71(2), 319-392.
