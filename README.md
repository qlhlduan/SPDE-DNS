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

- Arbia, G., & Di Marcantonio, M. (2015). Forecasting interest rates using geostatistical techniques. Econometrics, 3(4), 733-760.
- Bolin, D., & Kirchner, K. (2020). The rational SPDE approach for Gaussian random fields with general smoothness. Journal of Computational and Graphical Statistics, 29(2), 274-285.
- Bolin, D., Simas, A. B., & Xiong, Z. (2024). Covariance-based rational approximations of fractional SPDEs for computationally efficient Bayesian inference. Journal of Computational and Graphical Statistics, 33(1), 64-74.
- Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of proper scoring rules. Statistical Science, 38(1), 140-159.
- Campbell, J. Y., & Shiller, R. J. (1991). Yield spreads and interest rate movements: A bird's eye view. The Review of Economic Studies, 58(3), 495-514.
- Carriero, A., Kapetanios, G., & Marcellino, M. (2012). Forecasting government bond yields with large Bayesian vector autoregressions. Journal of Banking & Finance, 36(7), 2026-2047.
- Congedi, A., De Iaco, S., & Posa, D. (2025). A term structure geostatistical model with correlated residuals: A comparative analysis. Spatial Statistics, 67, 100886.
- Diebold, F. X., & Li, C. (2006). Forecasting the term structure of government bond yields. Journal of Econometrics, 130(2), 337-364.
- Diebold, F. X., Rudebusch, G. D., & Aruoba, S. B. (2006). The macroeconomy and the yield curve: A dynamic latent factor approach. Journal of Econometrics, 131(1-2), 309-338.
- Fama, E. F., & Bliss, R. R. (1987). The information in the term structure. The American Economic Review, 77(4), 680-692.
- Fuglstad, G.-A., Lindgren, F., Simpson, D., & Rue, H. (2015). Exploring a new class of non-stationary spatial Gaussian random fields with varying local anisotropy. Statistica Sinica, 115-133.
- Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules, prediction, and estimation. Journal of the American Statistical Association, 102(477), 359-378.
- Koopman, S. J., Mallee, M. I., & Van der Wel, M. (2010). Analyzing the term structure of interest rates using the dynamic Nelson-Siegel model with time-varying parameters. Journal of Business & Economic Statistics, 28(3), 329-343.
- Laurini, M. P., & Hotta, L. K. (2014). Forecasting the term structure of interest rates using integrated nested Laplace approximations. Journal of Forecasting, 33(3), 214-230.
- Lindgren, F., Bachl, F., Illian, J., Suen, M. H., Rue, H., & Seaton, A. E. (2024). inlabru: Software for fitting latent Gaussian models with non-linear predictors. arXiv preprint arXiv:2407.00791.
- Lindgren, F., Bakka, H., Bolin, D., Krainski, E., & Rue, H. (2024). A diffusion-based spatio-temporal extension of Gaussian Matérn fields. SORT-Statistics and Operations Research Transactions, 48(1), 3-66.
- Lindgren, F., & Rue, H. (2015). Bayesian spatial modelling with R-INLA. Journal of Statistical Software, 63(19).
- Markowitz, H. (1952). Portfolio selection. The Journal of Finance, 7(1), 77-91.
- Nelson, C. R., & Siegel, A. F. (1987). Parsimonious modeling of yield curves. The Journal of Business, 60(4), 473-489.
- Olafsdottir, H. K., Rootzén, H., & Bolin, D. (2024). Locally tail-scale invariant scoring rules for evaluation of extreme value forecasts. International Journal of Forecasting, 40(4), 1701-1720.
- Opschoor, D., & van der Wel, M. (2025). A smooth shadow-rate dynamic Nelson-Siegel model for yields at the zero lower bound. Journal of Business & Economic Statistics, 43(2), 298-311.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society Series B: Statistical Methodology, 71(2), 319-392.
- Seeger, M. (2004). Gaussian processes for machine learning. International Journal of Neural Systems, 14(02), 69-106.
- Stein, M. L. (1999). Interpolation of spatial data: Some theory for kriging. Springer, New York.
- Valente, F., & Laurini, M. (2024). Bayesian inference for long memory term structure models. Journal of Statistical Computation and Simulation, 94(8), 1735-1759.
- Van Niekerk, J., & Rue, H. (2024). Low-rank variational Bayes correction to the Laplace method. Journal of Machine Learning Research, 25(62), 1-25.
- Vasicek, O. (1977). An equilibrium characterization of the term structure. Journal of Financial Economics, 5(2), 177-188.
- West, K. D., Edison, H. J., & Cho, D. (1993). A utility-based comparison of some models of exchange rate volatility. Journal of International Economics, 35(1-2), 23-45.
- Xiang, J., & Zhu, X. (2013). A regime-switching Nelson-Siegel term structure model and interest rate forecasts. Journal of Financial Econometrics, 11(3), 522-555.
- Zhuang, K. (2021). Statistical properties of random field representation of the yield curve. PhD thesis, State University of New York at Stony Brook.
