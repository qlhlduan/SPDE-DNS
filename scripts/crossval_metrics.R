library(dplyr)
library(parallel)

# Read data
res <- read.csv("BDNS.csv")
yield <- read.csv("yield.csv")

# Filter and prepare yield data
data <- yield[yield$Date >= "1985-01-01" & yield$Date <= "2000-12-31",]
data <- data[,-c(1,2)]
rownames(data) <- 1:nrow(data)
target <- data[,c(1, 4, 10, 12, 17)]

# Create target_new (actual observations)
target_new <- do.call(rbind, lapply(120:180, function(i) {
  new_row <- c(as.numeric(target[i + 1, ]),
               as.numeric(target[i + 6, ]),
               as.numeric(target[i + 12, ]))
  return(new_row)
}))

target_new <- as.data.frame(target_new)
colnames(target_new) <- paste0("V", 1:15)

# Create threshold_new (1.05 times the previous observation)
threshold_new <- do.call(rbind, lapply(119:179, function(i) {
  new_row <- round(1.05 * c(as.numeric(target[i + 1, ]),
                            as.numeric(target[i + 6, ]),
                            as.numeric(target[i + 12, ])),3)
  return(new_row)
}))

threshold_new <- as.data.frame(threshold_new)
colnames(threshold_new) <- paste0("V", 1:15)

# Add group identifier to forecasts
res$group <- rep(1:61, each = 2000)

# Compute mean forecasts by group
res_grouped <- res %>%
  mutate(group = group) %>%  
  group_by(group) %>%
  summarise(across(everything(), mean))

res_result <- res_grouped %>% select(-group)

# RMSE function
rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}

# Calculate RMSE for each variable
rmse_values <- sapply(seq_along(target_new), function(i) {
  rmse(target_new[[i]], res_result[[i]])
})

rmse_output <- data.frame(Variable = names(target_new), RMSE = rmse_values)

# CRPS function
compute_crps <- function(forecast, y_obs) {
  N <- length(forecast)
  term1 <- mean(abs(forecast - y_obs))
  term2 <- mean(abs(outer(forecast, forecast, "-")))
  term1 - 0.5 * term2
}

# Calculate average CRPS
avg_crps <- sapply(names(target_new), function(col) {
  crps_by_group <- sapply(1:61, function(i) {
    y_obs <- target_new[i, col]                   
    forecast_samples <- res[res$group == i, col]  
    compute_crps(forecast_samples, y_obs)
  })
  mean(crps_by_group)
})

# SCRPS function
compute_scrps <- function(forecast, y_obs) {
  N <- length(forecast)
  term1 <- mean(abs(forecast - y_obs))
  term2 <- mean(abs(outer(forecast, forecast, "-")))
  term1 / term2 + 0.5 * log(term2)
}

# Calculate average SCRPS
avg_scrps <- sapply(names(target_new), function(col) {
  scrps_by_group <- sapply(1:61, function(i) {
    y_obs <- target_new[i, col]                   
    forecast_samples <- res[res$group == i, col]  
    compute_scrps(forecast_samples, y_obs)
  })
  mean(scrps_by_group)
})

crps_output <- data.frame(Variable = names(target_new), 
                          CRPS = avg_crps, 
                          SCRPS = avg_scrps)

# Helper function for weighted scores
compute_integral_indicator <- function(y1, y2, u) {
  lower <- min(y1, y2)
  upper <- max(y1, y2)
  return(max(0, upper - max(u, lower)))
}

# WCRPS function
compute_wcrps <- function(forecast, y_obs, threshold) {
  N <- length(forecast)
  term1 <- mean(sapply(forecast, compute_integral_indicator, y2 = y_obs, u = threshold))
  term2 <- mean(outer(forecast, forecast, Vectorize(function(y1, y2) compute_integral_indicator(y1, y2, threshold))))
  return(term1 - 0.5 * term2)
}

# WSCRPS function
compute_wscrps <- function(forecast, y_obs, threshold) {
  N <- length(forecast)
  term1 <- mean(abs(sapply(forecast, compute_integral_indicator, y2 = y_obs, u = threshold)))
  term2 <- mean(abs(outer(forecast, forecast, Vectorize(function(y1, y2) compute_integral_indicator(y1, y2, threshold)))))
  return(term1 / term2 + 0.5 * log(term2))
}

# Set number of cores for parallel processing
num_cores <- detectCores() - 1

# Calculate average WCRPS (parallel)
start_time <- proc.time()
avg_wcrps <- mclapply(names(target_new), function(col) {
  wcrps_by_group <- sapply(1:61, function(i) {
    y_obs <- target_new[i, col]                   
    forecast_samples <- res[res$group == i, col]
    threshold <- threshold_new[i, col]
    compute_wcrps(forecast_samples, y_obs, threshold)
  })
  mean(wcrps_by_group)
}, mc.cores = num_cores)
cat("WCRPS computation time:\n")
print(proc.time() - start_time)

# Calculate average WSCRPS (parallel)
start_time <- proc.time()
avg_wscrps <- mclapply(names(target_new), function(col) {
  wscrps_by_group <- sapply(1:61, function(i) {
    y_obs <- target_new[i, col]                   
    forecast_samples <- res[res$group == i, col]
    threshold <- threshold_new[i, col]
    compute_wscrps(forecast_samples, y_obs, threshold)
  })
  mean(wscrps_by_group)
}, mc.cores = num_cores)
cat("WSCRPS computation time:\n")
print(proc.time() - start_time)

wcrps_output <- data.frame(Variable = names(target_new), 
                           WCRPS = unlist(avg_wcrps), 
                           WSCRPS = unlist(avg_wscrps))
