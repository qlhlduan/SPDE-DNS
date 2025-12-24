library(INLA)
library(rSPDE)
library(fmesher)
library(inlabru)
library(Matrix)
library(sf)
library(fields)
source("nelson.R")

# Load and prepare data
yield <- read.csv("yield.csv")
data <- yield[yield$Date >= "1985-01-01" & yield$Date <= "2000-12-31",]
data <- data[,-c(1,2)]
rownames(data) <- 1:nrow(data)

# Set up parameters
maturity <- c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)
logmaturity <- log(maturity) * 24
dates <- 1:nrow(data)

# Prepare data for modeling
ix <- rep(dates, length(maturity))
iy <- rep(maturity, each = length(dates))
iylog <- rep(logmaturity, each = length(dates))

ym <- as.vector(t(t(as.matrix(data))))

data_joint <- data.frame(Y = ym, dates = ix, maturity = iylog, maturity_original = iy)

# Set up spatio-temporal field
s <- seq(from = 0, to = 120, length.out = 200)
t <- seq(from = 0, to = 192, length.out = 300)
model_st <- rspde.spacetime(space_loc = s, time_loc = t, alpha = 1, beta = 1, drift = FALSE)

# Define Nelson-Siegel factor loading functions with flexible lambda
f2_flex <- function(lambda, mat) {
  z <- lambda * mat
  z <- pmax(z, 1e-6)  # Avoid division by zero
  return((1 - exp(-z)) / z)
}

f3_flex <- function(lambda, mat) {
  f2_flex(lambda, mat) - exp(-lambda * mat)
}

# Define prior for lambda
lambda_gamma <- function(u) {
  shape <- 4       # shape parameter
  rate  <- shape / 0.068  # mean parameter
  qgamma(pnorm(u), shape = shape, rate = rate)
  }

lambda_lognormal <- function(u) {
  meanlog <- -2.707 # location parameter 
  sdlog <- 0.190 # log-scale spread 
  qlnorm(pnorm(u), meanlog = meanlog, sdlog = sdlog)
  }

# Define components with flexible lambda
components <- ~ lambda_internal(1, model = "linear", mean.linear = 0, prec.linear = 1) +
  beta1(dates, model = "ar1") +
  beta2(dates, model = "ar1") +
  beta3(dates, model = "ar1") +
  field(list(space = maturity, time = dates), model = model_st)

# Define formula with spatio-temporal field
formula_gamma <- Y ~ Intercept + beta1 + 
  f2_flex(lambda_gamma(lambda_internal), maturity_original) * beta2 + f2_flex(lambda_gamma(lambda_internal), maturity_original) + 
  f3_flex(lambda_gamma(lambda_internal), maturity_original) * beta3 + f3_flex(lambda_gamma(lambda_internal), maturity_original) + 
  field

formula_lognormal <- Y ~ Intercept + beta1 + 
  f2_flex(lambda_lognormal(lambda_internal), maturity_original) * beta2 + f2_flex(lambda_lognormal(lambda_internal), maturity_original) + 
  f3_flex(lambda_lognormal(lambda_internal), maturity_original) * beta3 + f3_flex(lambda_lognormal(lambda_internal), maturity_original) + 
  field

# Fit models
fit_gamma <- bru(components = components, formula = formula_gamma, family = "gaussian", data = data_joint)
fit_lognormal <- bru(components = components, formula = formula_lognormal, family = "gaussian", data = data_joint)

#Save models
save(fit_gamma,file="Fits/fit_gamma.RData")
save(fit_lognormal,file="Fits/fit_lognormal.RData")
