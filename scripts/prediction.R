#load models
load("Fits/fit_bdns.RData")

#store the model for cross-validation
models <- list(model_name = fit_bdns)

# Cross-validation setup
index <- c(120:180)
mat <- c(0, 3, 9, 11, 16)
out <- c(1, 6, 12)

train_test_index <- list()
for (i in index) {
  temp <- list(
    train = list(rep(c(1:i), 17) + rep(192 * c(0:16), each = i)), 
    test = list(rep(i + 192 * mat, 3) + rep(out, each = 5))
  )
  train_test_index <- append(train_test_index, list(temp))
}

# Perform cross-validation
cv_result <- cross_validation(
  models_combined,
  n_samples = 1000,
  return_scores_folds = TRUE,
  include_best = TRUE,
  train_test_indexes = train_test_index,
  return_train_test = TRUE,
  return_post_samples = TRUE,
  return_true_test_values = TRUE,
  parallelize_RP = TRUE,
  n_cores_RP = parallel::detectCores() - 1,
  true_CV = TRUE,
  save_settings = TRUE,
  print = TRUE
)

# Save results
res <- matrix(unlist(cv_result$post_samples$model_name), 
                       ncol = 15, byrow = TRUE)
write.csv(res, "BDNS.csv", row.names = FALSE)
