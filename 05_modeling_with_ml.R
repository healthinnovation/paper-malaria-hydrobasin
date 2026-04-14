library(randomForest)
library(xgboost)
library(e1071)
library(caret)
library(kernlab)   # required for svmRadial in caret
library(tidyverse)

# 1. Reading PCA database ------------------------------------------------------
dist.pca <- read_csv("output/pca/db.dist.pca.csv")
hy6.pca  <- read_csv("output/pca/db.hy6.pca.csv")
hy7.pca  <- read_csv("output/pca/db.hy7.pca.csv")
hya.pca  <- read_csv("output/pca/db.hya.pca.csv")

# Metric functions -------------------------------------------------------------
rse_n <- function(actual, predicted, p = 1) {
  residuals <- actual - predicted
  n <- length(actual)
  output <- sqrt(sum(residuals^2) / (n - p))
  return(output)
}

rmse_n <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

# ==============================================================================
# HYPERPARAMETER GRIDS — defined once, shared across all spatial units
# ==============================================================================

# --- Random Forest ---
# mtry explored via tuneLength = 10 + search = "random" (caret)
# ntree fixed at 500: Breiman (2001); Hastie, Tibshirani & Friedman (2009)
rf.control <- trainControl(method = "cv", number = 5, search = "random",
                           verboseIter = TRUE)
rf.ntree <- 500

# --- XGBoost ---
# Manual 5-fold CV — required for xgboost >= 2.0 (caret's xgbTree incompatible)
# Ranges: Chen & Guestrin (2016); disease-modelling benchmarks
set.seed(2025)
xgb.grid <- expand.grid(
  nrounds   = c(50, 100, 200, 300),
  eta       = c(0.01, 0.05, 0.1, 0.2, 0.3),
  max_depth = c(3, 4, 6, 8)
) |> slice_sample(n = 20)

# XGBoost CV tuning — xgboost >= 2.0 API (x/y instead of data/params)
xgb_cv_tune <- function(x_train, y_train, grid, nfolds = 5, seed = 2025) {
  set.seed(seed)
  folds <- createFolds(y_train, k = nfolds, list = TRUE)
  
  results <- map_dfr(seq_len(nrow(grid)), function(i) {
    params <- grid[i, ]
    
    fold_rmse <- map_dbl(folds, function(idx) {
      x_tr <- as.matrix(x_train[-idx, ])
      x_va <- as.matrix(x_train[idx, ])
      y_tr <- y_train[-idx]
      y_va <- y_train[idx]
      
      model <- xgboost(
        x                = x_tr,
        y                = y_tr,
        nrounds          = params$nrounds,
        eta              = params$eta,
        max_depth        = params$max_depth,
        gamma            = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample        = 1,
        objective        = "reg:squarederror",
        verbose          = 0
      )
      preds <- predict(model, x_va)
      sqrt(mean((y_va - preds)^2))
    })
    
    tibble(
      nrounds   = params$nrounds,
      eta       = params$eta,
      max_depth = params$max_depth,
      cv_rmse   = mean(fold_rmse)
    )
  })
  
  results |> slice_min(cv_rmse, n = 1, with_ties = FALSE)
}

# --- SVM (RBF kernel) ---
# Log-scale ranges: Hsu, Chang & Lin (2003) A Practical Guide to SVM
set.seed(2025)
svm.grid <- expand.grid(
  C     = 2^seq(-3, 7, by = 1),
  sigma = 2^seq(-7, 1, by = 1)
) |> slice_sample(n = 20)

svm.control <- trainControl(method = "cv", number = 5, search = "grid",
                            verboseIter = TRUE)

# ==============================================================================

# 2. Split database ------------------------------------------------------------
set.seed(2025)

## 2.1 Districts
index.fal <- createDataPartition(dist.pca$fal, p = 0.8, list = FALSE)
index.viv <- createDataPartition(dist.pca$viv, p = 0.8, list = FALSE)
dist.train.data.fal <- dist.pca[index.fal, ]
dist.test.data.fal  <- dist.pca[-index.fal, ]
dist.train.data.viv <- dist.pca[index.viv, ]
dist.test.data.viv  <- dist.pca[-index.viv, ]

## 2.2 Hydrobasin level 6
index.fal <- createDataPartition(hy6.pca$fal, p = 0.8, list = FALSE)
index.viv <- createDataPartition(hy6.pca$viv, p = 0.8, list = FALSE)
hy6.train.data.fal <- hy6.pca[index.fal, ]
hy6.test.data.fal  <- hy6.pca[-index.fal, ]
hy6.train.data.viv <- hy6.pca[index.viv, ]
hy6.test.data.viv  <- hy6.pca[-index.viv, ]

## 2.3 Hydrobasin level 7
index.fal <- createDataPartition(hy7.pca$fal, p = 0.8, list = FALSE)
index.viv <- createDataPartition(hy7.pca$viv, p = 0.8, list = FALSE)
hy7.train.data.fal <- hy7.pca[index.fal, ]
hy7.test.data.fal  <- hy7.pca[-index.fal, ]
hy7.train.data.viv <- hy7.pca[index.viv, ]
hy7.test.data.viv  <- hy7.pca[-index.viv, ]

## 2.4 Hydrobasin ANA
index.fal <- createDataPartition(hya.pca$fal, p = 0.8, list = FALSE)
index.viv <- createDataPartition(hya.pca$viv, p = 0.8, list = FALSE)
hya.train.data.fal <- hya.pca[index.fal, ]
hya.test.data.fal  <- hya.pca[-index.fal, ]
hya.train.data.viv <- hya.pca[index.viv, ]
hya.test.data.viv  <- hya.pca[-index.viv, ]

# ==============================================================================
# 3. Districts
# ==============================================================================
variables <- c(
  "year", "month", "pob",
  "pc1_clim", "pc2_clim", "pc3_clim", "pc4_clim", "pc5_clim",
  "pc1_hidr", "pc2_hidr", "pc3_hidr", "pc4_hidr",
  "pc1_veg",  "pc2_veg",  "pc3_veg"
)

x.train.fal <- dist.train.data.fal[, variables]
x.test.fal  <- dist.test.data.fal[, variables]
y.train.fal <- dist.train.data.fal$fal
y.test.fal  <- dist.test.data.fal$fal

x.train.viv <- dist.train.data.viv[, variables]
x.test.viv  <- dist.test.data.viv[, variables]
y.train.viv <- dist.train.data.viv$viv
y.test.viv  <- dist.test.data.viv$viv

## 3.1 Random Forest -----------------------------------------------------------
set.seed(2025)
model.rf.fal.cv <- train(
  x = x.train.fal, y = y.train.fal,
  method     = "rf",
  trControl  = rf.control,
  tuneLength = 10,
  ntree      = rf.ntree)

model.rf.viv.cv <- train(
  x = x.train.viv, y = y.train.viv,
  method     = "rf",
  trControl  = rf.control,
  tuneLength = 10,
  ntree      = rf.ntree)

predictions.rf.fal <- predict(model.rf.fal.cv, x.test.fal) %>% round()
predictions.rf.viv <- predict(model.rf.viv.cv, x.test.viv) %>% round()

rse.rf.fal  <- rse_n(actual = y.test.fal, predicted = predictions.rf.fal) %>% round(2)
rmse.rf.fal <- rmse_n(actual = y.test.fal, predicted = predictions.rf.fal) %>% round(2)
rse.rf.viv  <- rse_n(actual = y.test.viv, predicted = predictions.rf.viv) %>% round(2)
rmse.rf.viv <- rmse_n(actual = y.test.viv, predicted = predictions.rf.viv) %>% round(2)

dist.pca.rf <- dist.pca |>
  mutate(
    predicted_fal_rf = as.vector(predict(model.rf.fal.cv, dist.pca) %>% round()),
    predicted_viv_rf = as.vector(predict(model.rf.viv.cv, dist.pca) %>% round()),
    residual_fal_rf  = fal - predicted_fal_rf,
    residual_viv_rf  = viv - predicted_viv_rf
  )

## 3.2 XGBoost -----------------------------------------------------------------
# Step 1: find best hyperparameters via manual 5-fold CV
best.xgb.fal <- xgb_cv_tune(x.train.fal, y.train.fal, xgb.grid)
best.xgb.viv <- xgb_cv_tune(x.train.viv, y.train.viv, xgb.grid)

# Step 2: refit on full training set with best parameters
model.xgb.fal <- xgboost(
  x                = as.matrix(x.train.fal),
  y                = y.train.fal,
  nrounds          = best.xgb.fal$nrounds,
  eta              = best.xgb.fal$eta,
  max_depth        = best.xgb.fal$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

model.xgb.viv <- xgboost(
  x                = as.matrix(x.train.viv),
  y                = y.train.viv,
  nrounds          = best.xgb.viv$nrounds,
  eta              = best.xgb.viv$eta,
  max_depth        = best.xgb.viv$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

predictions.xgb.fal <- predict(model.xgb.fal, as.matrix(x.test.fal)) %>% round()
predictions.xgb.viv <- predict(model.xgb.viv, as.matrix(x.test.viv)) %>% round()

rmse.xgb.fal <- rmse_n(actual = y.test.fal, predicted = predictions.xgb.fal) %>% round(2)
rse.xgb.fal  <- rse_n(actual = y.test.fal, predicted = predictions.xgb.fal) %>% round(2)
rmse.xgb.viv <- rmse_n(actual = y.test.viv, predicted = predictions.xgb.viv) %>% round(2)
rse.xgb.viv  <- rse_n(actual = y.test.viv, predicted = predictions.xgb.viv) %>% round(2)

dist.pca.xgb <- dist.pca.rf |>
  mutate(
    predicted_fal_xgb = as.vector(predict(model.xgb.fal,
                                          as.matrix(dist.pca[, variables])) %>% round()),
    predicted_viv_xgb = as.vector(predict(model.xgb.viv,
                                          as.matrix(dist.pca[, variables])) %>% round()),
    residual_fal_xgb  = fal - predicted_fal_xgb,
    residual_viv_xgb  = viv - predicted_viv_xgb)

## 3.3 SVM ---------------------------------------------------------------------
set.seed(2025)
model.svm.fal <- train(
  x = x.train.fal, y = y.train.fal,
  method    = "svmRadial",
  trControl = svm.control,
  tuneGrid  = svm.grid)

model.svm.viv <- train(
  x = x.train.viv, y = y.train.viv,
  method    = "svmRadial",
  trControl = svm.control,
  tuneGrid  = svm.grid)

predictions.svm.fal <- predict(model.svm.fal, x.test.fal) %>% round()
predictions.svm.viv <- predict(model.svm.viv, x.test.viv) %>% round()

rmse.svm.fal <- rmse_n(actual = y.test.fal, predicted = predictions.svm.fal) %>% round(2)
rse.svm.fal  <- rse_n(actual = y.test.fal, predicted = predictions.svm.fal) %>% round(2)
rmse.svm.viv <- rmse_n(actual = y.test.viv, predicted = predictions.svm.viv) %>% round(2)
rse.svm.viv  <- rse_n(actual = y.test.viv, predicted = predictions.svm.viv) %>% round(2)

dist.pca.svm <- dist.pca.xgb |>
  mutate(
    predicted_fal_svm = as.vector(predict(model.svm.fal,
                                          dist.pca[, variables]) %>% round()),
    predicted_viv_svm = as.vector(predict(model.svm.viv,
                                          dist.pca[, variables]) %>% round()),
    residual_fal_svm  = fal - predicted_fal_svm,
    residual_viv_svm  = viv - predicted_viv_svm)

## 3.4 Best hyperparameters ----------------------------------------------------
best.params.dist <- bind_rows(
  tibble(unit = "districts", model = "RF", disease = "P.Falciparum",
         param = c("mtry", "ntree"),
         value = c(model.rf.fal.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "districts", model = "RF", disease = "P.Vivax",
         param = c("mtry", "ntree"),
         value = c(model.rf.viv.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "districts", model = "XGBoost", disease = "P.Falciparum",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.fal$nrounds, best.xgb.fal$eta,
                   best.xgb.fal$max_depth)),
  tibble(unit = "districts", model = "XGBoost", disease = "P.Vivax",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.viv$nrounds, best.xgb.viv$eta,
                   best.xgb.viv$max_depth)),
  tibble(unit = "districts", model = "SVM", disease = "P.Falciparum",
         param = c("C", "sigma"),
         value = c(model.svm.fal$bestTune$C, model.svm.fal$bestTune$sigma)),
  tibble(unit = "districts", model = "SVM", disease = "P.Vivax",
         param = c("C", "sigma"),
         value = c(model.svm.viv$bestTune$C, model.svm.viv$bestTune$sigma))
)

## 3.5 Error metrics -----------------------------------------------------------
met_rf  <- data.frame(modelo = "rf",  disease = c("P.Falciparum", "P.Vivax"),
                      rmse = c(rmse.rf.fal,  rmse.rf.viv),
                      rse  = c(rse.rf.fal,   rse.rf.viv))
met_xgb <- data.frame(modelo = "xgb", disease = c("P.Falciparum", "P.Vivax"),
                      rmse = c(rmse.xgb.fal, rmse.xgb.viv),
                      rse  = c(rse.xgb.fal,  rse.xgb.viv))
met_svm <- data.frame(modelo = "svm", disease = c("P.Falciparum", "P.Vivax"),
                      rmse = c(rmse.svm.fal, rmse.svm.viv),
                      rse  = c(rse.svm.fal,  rse.svm.viv))

met_total <- bind_rows(met_rf, met_xgb, met_svm) |>
  mutate(region = "districts")

if (!dir.exists("output/metrics")) dir.create("output/metrics")
write_csv(met_total,        "output/metrics/error_metrics_districts.csv")
write_csv(dist.pca.svm,     "output/metrics/error_metrics_total_districts.csv")
write_csv(best.params.dist, "output/metrics/best_hyperparameters_districts.csv")


# ==============================================================================
# 4. Hydrobasin level 6
# ==============================================================================
variables <- c(
  "year", "month", "pob",
  "pc1_clim", "pc2_clim", "pc3_clim", "pc4_clim",
  "pc1_hidr", "pc2_hidr", "pc3_hidr", "pc4_hidr",
  "pc1_veg",  "pc2_veg",  "pc3_veg"
)

x.train.fal <- hy6.train.data.fal[, variables]
x.test.fal  <- hy6.test.data.fal[, variables]
y.train.fal <- hy6.train.data.fal$fal
y.test.fal  <- hy6.test.data.fal$fal

x.train.viv <- hy6.train.data.viv[, variables]
x.test.viv  <- hy6.test.data.viv[, variables]
y.train.viv <- hy6.train.data.viv$viv
y.test.viv  <- hy6.test.data.viv$viv

## 4.1 Random Forest -----------------------------------------------------------
set.seed(2025)
model.rf.fal.cv <- train(
  x = x.train.fal, y = y.train.fal,
  method = "rf", trControl = rf.control,
  tuneLength = 10, ntree = rf.ntree)

model.rf.viv.cv <- train(
  x = x.train.viv, y = y.train.viv,
  method = "rf", trControl = rf.control,
  tuneLength = 10, ntree = rf.ntree)

predictions.rf.fal <- predict(model.rf.fal.cv, x.test.fal) %>% round()
predictions.rf.viv <- predict(model.rf.viv.cv, x.test.viv) %>% round()

rse.rf.fal  <- rse_n(y.test.fal, predictions.rf.fal) %>% round(2)
rmse.rf.fal <- rmse_n(y.test.fal, predictions.rf.fal) %>% round(2)
rse.rf.viv  <- rse_n(y.test.viv, predictions.rf.viv) %>% round(2)
rmse.rf.viv <- rmse_n(y.test.viv, predictions.rf.viv) %>% round(2)

hy6.pca.rf <- hy6.pca |>
  mutate(
    predicted_fal_rf = as.vector(predict(model.rf.fal.cv, hy6.pca) %>% round()),
    predicted_viv_rf = as.vector(predict(model.rf.viv.cv, hy6.pca) %>% round()),
    residual_fal_rf  = fal - predicted_fal_rf,
    residual_viv_rf  = viv - predicted_viv_rf)

## 4.2 XGBoost -----------------------------------------------------------------
best.xgb.fal <- xgb_cv_tune(x.train.fal, y.train.fal, xgb.grid)
best.xgb.viv <- xgb_cv_tune(x.train.viv, y.train.viv, xgb.grid)

model.xgb.fal <- xgboost(
  x                = as.matrix(x.train.fal),
  y                = y.train.fal,
  nrounds          = best.xgb.fal$nrounds,
  eta              = best.xgb.fal$eta,
  max_depth        = best.xgb.fal$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

model.xgb.viv <- xgboost(
  x                = as.matrix(x.train.viv),
  y                = y.train.viv,
  nrounds          = best.xgb.viv$nrounds,
  eta              = best.xgb.viv$eta,
  max_depth        = best.xgb.viv$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

predictions.xgb.fal <- predict(model.xgb.fal, as.matrix(x.test.fal)) %>% round()
predictions.xgb.viv <- predict(model.xgb.viv, as.matrix(x.test.viv)) %>% round()

rmse.xgb.fal <- rmse_n(y.test.fal, predictions.xgb.fal) %>% round(2)
rse.xgb.fal  <- rse_n(y.test.fal, predictions.xgb.fal) %>% round(2)
rmse.xgb.viv <- rmse_n(y.test.viv, predictions.xgb.viv) %>% round(2)
rse.xgb.viv  <- rse_n(y.test.viv, predictions.xgb.viv) %>% round(2)

hy6.pca.xgb <- hy6.pca.rf |>
  mutate(
    predicted_fal_xgb = as.vector(predict(model.xgb.fal,
                                          as.matrix(hy6.pca[, variables])) %>% round()),
    predicted_viv_xgb = as.vector(predict(model.xgb.viv,
                                          as.matrix(hy6.pca[, variables])) %>% round()),
    residual_fal_xgb  = fal - predicted_fal_xgb,
    residual_viv_xgb  = viv - predicted_viv_xgb)

## 4.3 SVM ---------------------------------------------------------------------
set.seed(2025)
model.svm.fal <- train(
  x = x.train.fal, y = y.train.fal,
  method = "svmRadial", trControl = svm.control, tuneGrid = svm.grid)

model.svm.viv <- train(
  x = x.train.viv, y = y.train.viv,
  method = "svmRadial", trControl = svm.control, tuneGrid = svm.grid)

predictions.svm.fal <- predict(model.svm.fal, x.test.fal) %>% round()
predictions.svm.viv <- predict(model.svm.viv, x.test.viv) %>% round()

rmse.svm.fal <- rmse_n(y.test.fal, predictions.svm.fal) %>% round(2)
rse.svm.fal  <- rse_n(y.test.fal, predictions.svm.fal) %>% round(2)
rmse.svm.viv <- rmse_n(y.test.viv, predictions.svm.viv) %>% round(2)
rse.svm.viv  <- rse_n(y.test.viv, predictions.svm.viv) %>% round(2)

hy6.pca.svm <- hy6.pca.xgb |>
  mutate(
    predicted_fal_svm = as.vector(predict(model.svm.fal,
                                          hy6.pca[, variables]) %>% round()),
    predicted_viv_svm = as.vector(predict(model.svm.viv,
                                          hy6.pca[, variables]) %>% round()),
    residual_fal_svm  = fal - predicted_fal_svm,
    residual_viv_svm  = viv - predicted_viv_svm)

## 4.4 Best hyperparameters ----------------------------------------------------
best.params.hy6 <- bind_rows(
  tibble(unit = "hydrobasin_l6", model = "RF", disease = "P.Falciparum",
         param = c("mtry", "ntree"),
         value = c(model.rf.fal.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "hydrobasin_l6", model = "RF", disease = "P.Vivax",
         param = c("mtry", "ntree"),
         value = c(model.rf.viv.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "hydrobasin_l6", model = "XGBoost", disease = "P.Falciparum",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.fal$nrounds, best.xgb.fal$eta,
                   best.xgb.fal$max_depth)),
  tibble(unit = "hydrobasin_l6", model = "XGBoost", disease = "P.Vivax",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.viv$nrounds, best.xgb.viv$eta,
                   best.xgb.viv$max_depth)),
  tibble(unit = "hydrobasin_l6", model = "SVM", disease = "P.Falciparum",
         param = c("C", "sigma"),
         value = c(model.svm.fal$bestTune$C, model.svm.fal$bestTune$sigma)),
  tibble(unit = "hydrobasin_l6", model = "SVM", disease = "P.Vivax",
         param = c("C", "sigma"),
         value = c(model.svm.viv$bestTune$C, model.svm.viv$bestTune$sigma))
)

## 4.5 Error metrics -----------------------------------------------------------
met_total <- bind_rows(
  data.frame(modelo = "rf",  disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.rf.fal,  rmse.rf.viv),
             rse  = c(rse.rf.fal,   rse.rf.viv)),
  data.frame(modelo = "xgb", disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.xgb.fal, rmse.xgb.viv),
             rse  = c(rse.xgb.fal,  rse.xgb.viv)),
  data.frame(modelo = "svm", disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.svm.fal, rmse.svm.viv),
             rse  = c(rse.svm.fal,  rse.svm.viv))
) |> mutate(region = "hydrobasin level 6")

write_csv(met_total,       "output/metrics/error_metrics_hy6.csv")
write_csv(hy6.pca.svm,     "output/metrics/error_metrics_total_hydrobasinlevel06.csv")
write_csv(best.params.hy6, "output/metrics/best_hyperparameters_hy6.csv")


# ==============================================================================
# 5. Hydrobasin level 7
# ==============================================================================
variables <- c(
  "year", "month", "pob",
  "pc1_clim", "pc2_clim", "pc3_clim", "pc4_clim", "pc5_clim",
  "pc1_hidr", "pc2_hidr", "pc3_hidr", "pc4_hidr",
  "pc1_veg",  "pc2_veg"
)

x.train.fal <- hy7.train.data.fal[, variables]
x.test.fal  <- hy7.test.data.fal[, variables]
y.train.fal <- hy7.train.data.fal$fal
y.test.fal  <- hy7.test.data.fal$fal

x.train.viv <- hy7.train.data.viv[, variables]
x.test.viv  <- hy7.test.data.viv[, variables]
y.train.viv <- hy7.train.data.viv$viv
y.test.viv  <- hy7.test.data.viv$viv

## 5.1 Random Forest -----------------------------------------------------------
set.seed(2025)
model.rf.fal.cv <- train(
  x = x.train.fal, y = y.train.fal,
  method = "rf", trControl = rf.control,
  tuneLength = 10, ntree = rf.ntree)

model.rf.viv.cv <- train(
  x = x.train.viv, y = y.train.viv,
  method = "rf", trControl = rf.control,
  tuneLength = 10, ntree = rf.ntree)

predictions.rf.fal <- predict(model.rf.fal.cv, x.test.fal) %>% round()
predictions.rf.viv <- predict(model.rf.viv.cv, x.test.viv) %>% round()

rse.rf.fal  <- rse_n(y.test.fal, predictions.rf.fal) %>% round(2)
rmse.rf.fal <- rmse_n(y.test.fal, predictions.rf.fal) %>% round(2)
rse.rf.viv  <- rse_n(y.test.viv, predictions.rf.viv) %>% round(2)
rmse.rf.viv <- rmse_n(y.test.viv, predictions.rf.viv) %>% round(2)

hy7.pca.rf <- hy7.pca |>
  mutate(
    predicted_fal_rf = as.vector(predict(model.rf.fal.cv, hy7.pca) %>% round()),
    predicted_viv_rf = as.vector(predict(model.rf.viv.cv, hy7.pca) %>% round()),
    residual_fal_rf  = fal - predicted_fal_rf,
    residual_viv_rf  = viv - predicted_viv_rf)

## 5.2 XGBoost -----------------------------------------------------------------
best.xgb.fal <- xgb_cv_tune(x.train.fal, y.train.fal, xgb.grid)
best.xgb.viv <- xgb_cv_tune(x.train.viv, y.train.viv, xgb.grid)

model.xgb.fal <- xgboost(
  x                = as.matrix(x.train.fal),
  y                = y.train.fal,
  nrounds          = best.xgb.fal$nrounds,
  eta              = best.xgb.fal$eta,
  max_depth        = best.xgb.fal$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

model.xgb.viv <- xgboost(
  x                = as.matrix(x.train.viv),
  y                = y.train.viv,
  nrounds          = best.xgb.viv$nrounds,
  eta              = best.xgb.viv$eta,
  max_depth        = best.xgb.viv$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

predictions.xgb.fal <- predict(model.xgb.fal, as.matrix(x.test.fal)) %>% round()
predictions.xgb.viv <- predict(model.xgb.viv, as.matrix(x.test.viv)) %>% round()

rmse.xgb.fal <- rmse_n(y.test.fal, predictions.xgb.fal) %>% round(2)
rse.xgb.fal  <- rse_n(y.test.fal, predictions.xgb.fal) %>% round(2)
rmse.xgb.viv <- rmse_n(y.test.viv, predictions.xgb.viv) %>% round(2)
rse.xgb.viv  <- rse_n(y.test.viv, predictions.xgb.viv) %>% round(2)

hy7.pca.xgb <- hy7.pca.rf |>
  mutate(
    predicted_fal_xgb = as.vector(predict(model.xgb.fal,
                                          as.matrix(hy7.pca[, variables])) %>% round()),
    predicted_viv_xgb = as.vector(predict(model.xgb.viv,
                                          as.matrix(hy7.pca[, variables])) %>% round()),
    residual_fal_xgb  = fal - predicted_fal_xgb,
    residual_viv_xgb  = viv - predicted_viv_xgb)

## 5.3 SVM ---------------------------------------------------------------------
set.seed(2025)
model.svm.fal <- train(
  x = x.train.fal, y = y.train.fal,
  method = "svmRadial", trControl = svm.control, tuneGrid = svm.grid)

model.svm.viv <- train(
  x = x.train.viv, y = y.train.viv,
  method = "svmRadial", trControl = svm.control, tuneGrid = svm.grid)

predictions.svm.fal <- predict(model.svm.fal, x.test.fal) %>% round()
predictions.svm.viv <- predict(model.svm.viv, x.test.viv) %>% round()

rmse.svm.fal <- rmse_n(y.test.fal, predictions.svm.fal) %>% round(2)
rse.svm.fal  <- rse_n(y.test.fal, predictions.svm.fal) %>% round(2)
rmse.svm.viv <- rmse_n(y.test.viv, predictions.svm.viv) %>% round(2)
rse.svm.viv  <- rse_n(y.test.viv, predictions.svm.viv) %>% round(2)

hy7.pca.svm <- hy7.pca.xgb |>
  mutate(
    predicted_fal_svm = as.vector(predict(model.svm.fal,
                                          hy7.pca[, variables]) %>% round()),
    predicted_viv_svm = as.vector(predict(model.svm.viv,
                                          hy7.pca[, variables]) %>% round()),
    residual_fal_svm  = fal - predicted_fal_svm,
    residual_viv_svm  = viv - predicted_viv_svm)

## 5.4 Best hyperparameters ----------------------------------------------------
best.params.hy7 <- bind_rows(
  tibble(unit = "hydrobasin_l7", model = "RF", disease = "P.Falciparum",
         param = c("mtry", "ntree"),
         value = c(model.rf.fal.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "hydrobasin_l7", model = "RF", disease = "P.Vivax",
         param = c("mtry", "ntree"),
         value = c(model.rf.viv.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "hydrobasin_l7", model = "XGBoost", disease = "P.Falciparum",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.fal$nrounds, best.xgb.fal$eta,
                   best.xgb.fal$max_depth)),
  tibble(unit = "hydrobasin_l7", model = "XGBoost", disease = "P.Vivax",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.viv$nrounds, best.xgb.viv$eta,
                   best.xgb.viv$max_depth)),
  tibble(unit = "hydrobasin_l7", model = "SVM", disease = "P.Falciparum",
         param = c("C", "sigma"),
         value = c(model.svm.fal$bestTune$C, model.svm.fal$bestTune$sigma)),
  tibble(unit = "hydrobasin_l7", model = "SVM", disease = "P.Vivax",
         param = c("C", "sigma"),
         value = c(model.svm.viv$bestTune$C, model.svm.viv$bestTune$sigma))
)

## 5.5 Error metrics -----------------------------------------------------------
met_total <- bind_rows(
  data.frame(modelo = "rf",  disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.rf.fal,  rmse.rf.viv),
             rse  = c(rse.rf.fal,   rse.rf.viv)),
  data.frame(modelo = "xgb", disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.xgb.fal, rmse.xgb.viv),
             rse  = c(rse.xgb.fal,  rse.xgb.viv)),
  data.frame(modelo = "svm", disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.svm.fal, rmse.svm.viv),
             rse  = c(rse.svm.fal,  rse.svm.viv))
) |> mutate(region = "hydrobasin level 7")

write_csv(met_total,       "output/metrics/error_metrics_hy7.csv")
write_csv(hy7.pca.svm,     "output/metrics/error_metrics_total_hydrobasinlevel7.csv")
write_csv(best.params.hy7, "output/metrics/best_hyperparameters_hy7.csv")


# ==============================================================================
# 6. Hydrobasin ANA
# ==============================================================================
variables <- c(
  "year", "month", "pob",
  "pc1_clim", "pc2_clim", "pc3_clim", "pc4_clim", "pc5_clim",
  "pc1_hidr", "pc2_hidr", "pc3_hidr", "pc4_hidr",
  "pc1_veg",  "pc2_veg"
)

x.train.fal <- hya.train.data.fal[, variables]
x.test.fal  <- hya.test.data.fal[, variables]
y.train.fal <- hya.train.data.fal$fal
y.test.fal  <- hya.test.data.fal$fal

x.train.viv <- hya.train.data.viv[, variables]
x.test.viv  <- hya.test.data.viv[, variables]
y.train.viv <- hya.train.data.viv$viv
y.test.viv  <- hya.test.data.viv$viv

## 6.1 Random Forest -----------------------------------------------------------
set.seed(2025)
model.rf.fal.cv <- train(
  x = x.train.fal, y = y.train.fal,
  method = "rf", trControl = rf.control,
  tuneLength = 10, ntree = rf.ntree)

model.rf.viv.cv <- train(
  x = x.train.viv, y = y.train.viv,
  method = "rf", trControl = rf.control,
  tuneLength = 10, ntree = rf.ntree)

predictions.rf.fal <- predict(model.rf.fal.cv, x.test.fal) %>% round()
predictions.rf.viv <- predict(model.rf.viv.cv, x.test.viv) %>% round()

rse.rf.fal  <- rse_n(y.test.fal, predictions.rf.fal) %>% round(2)
rmse.rf.fal <- rmse_n(y.test.fal, predictions.rf.fal) %>% round(2)
rse.rf.viv  <- rse_n(y.test.viv, predictions.rf.viv) %>% round(2)
rmse.rf.viv <- rmse_n(y.test.viv, predictions.rf.viv) %>% round(2)

hya.pca.rf <- hya.pca |>
  mutate(
    predicted_fal_rf = as.vector(predict(model.rf.fal.cv, hya.pca) %>% round()),
    predicted_viv_rf = as.vector(predict(model.rf.viv.cv, hya.pca) %>% round()),
    residual_fal_rf  = fal - predicted_fal_rf,
    residual_viv_rf  = viv - predicted_viv_rf)

## 6.2 XGBoost -----------------------------------------------------------------
best.xgb.fal <- xgb_cv_tune(x.train.fal, y.train.fal, xgb.grid)
best.xgb.viv <- xgb_cv_tune(x.train.viv, y.train.viv, xgb.grid)

model.xgb.fal <- xgboost(
  x                = as.matrix(x.train.fal),
  y                = y.train.fal,
  nrounds          = best.xgb.fal$nrounds,
  eta              = best.xgb.fal$eta,
  max_depth        = best.xgb.fal$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

model.xgb.viv <- xgboost(
  x                = as.matrix(x.train.viv),
  y                = y.train.viv,
  nrounds          = best.xgb.viv$nrounds,
  eta              = best.xgb.viv$eta,
  max_depth        = best.xgb.viv$max_depth,
  gamma            = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample        = 1,
  objective        = "reg:squarederror",
  verbose          = 0)

predictions.xgb.fal <- predict(model.xgb.fal, as.matrix(x.test.fal)) %>% round()
predictions.xgb.viv <- predict(model.xgb.viv, as.matrix(x.test.viv)) %>% round()

rmse.xgb.fal <- rmse_n(y.test.fal, predictions.xgb.fal) %>% round(2)
rse.xgb.fal  <- rse_n(y.test.fal, predictions.xgb.fal) %>% round(2)
rmse.xgb.viv <- rmse_n(y.test.viv, predictions.xgb.viv) %>% round(2)
rse.xgb.viv  <- rse_n(y.test.viv, predictions.xgb.viv) %>% round(2)

hya.pca.xgb <- hya.pca.rf |>
  mutate(
    predicted_fal_xgb = as.vector(predict(model.xgb.fal,
                                          as.matrix(hya.pca[, variables])) %>% round()),
    predicted_viv_xgb = as.vector(predict(model.xgb.viv,
                                          as.matrix(hya.pca[, variables])) %>% round()),
    residual_fal_xgb  = fal - predicted_fal_xgb,
    residual_viv_xgb  = viv - predicted_viv_xgb)

## 6.3 SVM ---------------------------------------------------------------------
set.seed(2025)
model.svm.fal <- train(
  x = x.train.fal, y = y.train.fal,
  method = "svmRadial", trControl = svm.control, tuneGrid = svm.grid)

model.svm.viv <- train(
  x = x.train.viv, y = y.train.viv,
  method = "svmRadial", trControl = svm.control, tuneGrid = svm.grid)

predictions.svm.fal <- predict(model.svm.fal, x.test.fal) %>% round()
predictions.svm.viv <- predict(model.svm.viv, x.test.viv) %>% round()

rmse.svm.fal <- rmse_n(y.test.fal, predictions.svm.fal) %>% round(2)
rse.svm.fal  <- rse_n(y.test.fal, predictions.svm.fal) %>% round(2)
rmse.svm.viv <- rmse_n(y.test.viv, predictions.svm.viv) %>% round(2)
rse.svm.viv  <- rse_n(y.test.viv, predictions.svm.viv) %>% round(2)

hya.pca.svm <- hya.pca.xgb |>
  mutate(
    predicted_fal_svm = as.vector(predict(model.svm.fal,
                                          hya.pca[, variables]) %>% round()),
    predicted_viv_svm = as.vector(predict(model.svm.viv,
                                          hya.pca[, variables]) %>% round()),
    residual_fal_svm  = fal - predicted_fal_svm,
    residual_viv_svm  = viv - predicted_viv_svm)

## 6.4 Best hyperparameters ----------------------------------------------------
best.params.hya <- bind_rows(
  tibble(unit = "hydrobasin_ana", model = "RF", disease = "P.Falciparum",
         param = c("mtry", "ntree"),
         value = c(model.rf.fal.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "hydrobasin_ana", model = "RF", disease = "P.Vivax",
         param = c("mtry", "ntree"),
         value = c(model.rf.viv.cv$bestTune$mtry, rf.ntree)),
  tibble(unit = "hydrobasin_ana", model = "XGBoost", disease = "P.Falciparum",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.fal$nrounds, best.xgb.fal$eta,
                   best.xgb.fal$max_depth)),
  tibble(unit = "hydrobasin_ana", model = "XGBoost", disease = "P.Vivax",
         param = c("nrounds", "eta", "max_depth"),
         value = c(best.xgb.viv$nrounds, best.xgb.viv$eta,
                   best.xgb.viv$max_depth)),
  tibble(unit = "hydrobasin_ana", model = "SVM", disease = "P.Falciparum",
         param = c("C", "sigma"),
         value = c(model.svm.fal$bestTune$C, model.svm.fal$bestTune$sigma)),
  tibble(unit = "hydrobasin_ana", model = "SVM", disease = "P.Vivax",
         param = c("C", "sigma"),
         value = c(model.svm.viv$bestTune$C, model.svm.viv$bestTune$sigma))
)

## 6.5 Error metrics -----------------------------------------------------------
met_total <- bind_rows(
  data.frame(modelo = "rf",  disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.rf.fal,  rmse.rf.viv),
             rse  = c(rse.rf.fal,   rse.rf.viv)),
  data.frame(modelo = "xgb", disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.xgb.fal, rmse.xgb.viv),
             rse  = c(rse.xgb.fal,  rse.xgb.viv)),
  data.frame(modelo = "svm", disease = c("P.Falciparum", "P.Vivax"),
             rmse = c(rmse.svm.fal, rmse.svm.viv),
             rse  = c(rse.svm.fal,  rse.svm.viv))
) |> mutate(region = "hydrobasin ANA")

write_csv(met_total,       "output/metrics/error_metrics_hya.csv")
write_csv(hya.pca.svm,     "output/metrics/error_metrics_total_hydrobasinana.csv")
write_csv(best.params.hya, "output/metrics/best_hyperparameters_hya.csv")


# ==============================================================================
# 7. Consolidate all hyperparameter tables
# ==============================================================================
best.params.all <- bind_rows(
  best.params.dist,
  best.params.hy6,
  best.params.hy7,
  best.params.hya
)
write_csv(best.params.all, "output/metrics/best_hyperparameters_all.csv")
message(">>> Done. Key output: output/metrics/best_hyperparameters_all.csv")