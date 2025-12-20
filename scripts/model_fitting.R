library(INLA) 
library(rSPDE)
library(fmesher)
library(inlabru)
library(Matrix)
library(sf)
library(fields)
source("nelson.R")

yield <- read.csv("yield.csv")
data <- yield[yield$Date >= "1985-01-01" & yield$Date <= "2000-12-31",]
data <- data[,-c(1,2)]
rownames(data) <- 1:nrow(data)

lambda <- 0.0609
maturity <- c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)
logmaturity <- log(maturity) * 24
dates <- 1:nrow(data)

f1 <- seq(1, 1, length = length(maturity))
f2 <- (1 - exp(-lambda * maturity)) / (lambda * maturity)
f3 <- ((1 - exp(-lambda * maturity)) / (lambda * maturity) - exp(-lambda * maturity))

F1 <- rep(f1, each = length(dates))
F2 <- rep(f2, each = length(dates))
F3 <- rep(f3, each = length(dates))

ix <- rep(dates, length(maturity))
iy <- rep(maturity, each = length(dates))
iylog <- rep(logmaturity, each = length(dates))

ym <- as.vector(t(t(as.matrix(data))))

data_ar <- data.frame(Y = ym, index_x = ix, index_y = iy, F1f = F1, F2f = F2, F3f = F3)
data_ngeo <- data.frame(Y = ym, index_x = ix, index_y = iy, dates = ix, maturity = iylog, F1f = F1, F2f = F2, F3f = F3)
data_geo <- st_as_sf(data_ngeo, coords = c("dates","maturity"))

ns <- build_nelsonsiegel_model(dates, maturity, lambda = 0.0609)
register_nelsonsiegel_inlabru()

#Bayesian DNS model
formula_ar1 <- Y ~ -1 + nelsonsiegel(cbind(index_x, index_y), model = ns) + F1f(F1f) + F2f(F2f) + F3f(F3f)

#Stationary SPDE model
data_locs = data.frame(dates = ix, maturity = iylog)
locs = st_as_sf(data_locs, coords = c("dates", "maturity"))
mesh <- fm_mesh_2d(loc = locs, max.edge = c(2, 4), offset = c(2, 4), cutoff = 1.5)
model_sf <- rspde.matern(mesh = mesh, parameterization = "spde")
formula_sf <- Y ~ -1 + Intercept(1) + nelsonsiegel(cbind(index_x, index_y), model = ns) + field(geometry, model = model_sf)

#Non-stationary SPDE model
sigma = cbind(0, 1, 0, (mesh$loc[,1]/192), (mesh$loc[,2]/120), 0, 0)
range = cbind(0, 0, 1, 0, 0, (mesh$loc[,1]/192), (mesh$loc[,2]/120))
model_nsf <- rspde.matern(mesh = mesh, B.sigma = sigma, B.range = range, parameterization = "matern")
formula_nsf <- Y ~ -1 + Intercept(1) + nelsonsiegel(cbind(index_x, index_y), model = ns) + field(geometry, model = model_nsf)

#Nonseparable Spatio-Temporal model
s <- seq(from = 0, to = 120, length.out = 200)
t <- seq(from = 0, to = 192, length.out = 300)
model_st <- rspde.spacetime(space_loc = s, time_loc = t, alpha = 1, beta = 1, drift = FALSE)
formula_st <- Y ~ -1 + Intercept(1) + nelsonsiegel(cbind(index_x, index_y), model = ns) + field(list(space = maturity, time = dates), model = model_st)

#Anisotropic model
model_aniso <- rspde.anistropic2d(mesh = mesh)
formula_aniso <- Y ~ -1 + Intercept(1) + nelsonsiegel(cbind(index_x, index_y), model = ns) + field(cbind(maturity, dates), model = model_aniso)

#Fit models
fit_ar1 <- bru(formula_ar1,  data = data_ar)
fit_sf <- bru(formula_sf,  data = data_geo)
fit_nsf <- bru(formula_nsf,  data = data_geo)
fit_st <- bru(formula_st, data = data_ngeo)
fit_aniso <- bru(formula_aniso, data = data_ngeo)


#Save models
saveRDS(list(fixed = fit_ar1$summary.fixed, hyper = fit_ar1$summary.hyperpar), 
        file = "Fits/fit_ar1_summary.rds")

saveRDS(list(fixed = fit_sf$summary.fixed, hyper = fit_sf$summary$hyperpar), 
        file = "Fits/fit_sf_summary.rds")

saveRDS(list(fixed = fit_nsf$summary.fixed, hyper = fit_nsf$summary.hyperpar), 
        file = "Fits/fit_nsf_summary.rds")

saveRDS(list(fixed = fit_st$summary.fixed, hyper = fit_st$summary$hyperpar), 
        file = "Fits/fit_st_summary.rds")

saveRDS(list(fixed = fit_aniso$summary.fixed, hyper = fit_aniso$summary.hyperpar), 
        file = "Fits/fit_aniso_summary.rds")

