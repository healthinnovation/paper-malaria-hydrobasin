library(tidyverse)
library(sf)
library(leaps)
library(MASS)
library(dplyr)

sf_use_s2(use_s2 = F)

# 1. Hydrobasin level 6 ---------------------------------------------------
sp.hy06 <- st_read(
  'rawdata/geometry.gpkg',
  layer = 'hydrobasin_06')|>
  mutate(new_id = 1:13) |> 
  mutate(hydroname = paste0(hydroname,'.06')) 

orden <- sp.hy06$hydroname

ddbb_hy06 <- st_read("data_for_model/hy6.gpkg") |> 
  mutate(hydroname = paste0(hydroname,'.06')) |> 
  mutate(new_id = as.integer(factor(hydroname,levels = orden,labels = 1:13))) %>%
  st_drop_geometry()


# Step forward Poisson ----------------------------------------------------

# Evaluando equidispersion
summary(ddbb_hy06$fal) # Mean 42.3
var(ddbb_hy06$fal) # Var 6908.8

summary(ddbb_hy06$viv) # Mean 160.6
var(ddbb_hy06$viv) # Var 59290.1

# No se cumple el supuesto de equidispersión, se usará Binomial Negativo

fit_fal_nb_null <- glm.nb(fal ~ offset(log(pob)), data = ddbb_hy06)
fit_viv_nb_null <- glm.nb(viv ~ offset(log(pob)), data = ddbb_hy06)

fit_step_fal <- stepAIC(
  object = fit_fal_nb_null, 
  scope = list(lower = fal ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_fal)

fit_step_viv <- stepAIC(
  object = fit_viv_nb_null, 
  scope = list(lower = viv ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_viv)


# Modelos finales cuenca Hy6: SFW Hy06------------------------------------------

final_fal_nb_hy06 <- glm.nb(fal ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                         t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                       data = ddbb_hy06)

summary(final_fal_nb_hy06)

pred_fal_hy06 <- predict (
  final_fal_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_fal_hy06 <- ddbb_hy06$fal
rmse_fal_hy06 <- sqrt(mean((obs_fal_hy06 - pred_fal_hy06)^2))
n     <- length(obs_fal_hy06)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_fal_hy06 - pred_fal_hy06

rse_fal_hy06 <- sqrt( sum(resid^2) / (n - p) )



final_viv_nb_hy06 <- glm.nb(viv ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                         t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                       data = ddbb_hy06)

summary(final_viv_nb_hy06)

pred_viv_hy06 <- predict (
  final_viv_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_viv_hy06 <- ddbb_hy06$viv
rmse_viv_hy06 <- sqrt(mean((obs_viv_hy06 - pred_viv_hy06)^2))
n     <- length(obs_viv_hy06)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_viv_hy06 - pred_viv_hy06

rse_viv_hy06 <- sqrt( sum(resid^2) / (n - p) )

#Modelos finales cuenca Hy6: SFW Hy07-------------------------------------------

final_fal_nb_hy06 <- glm.nb(fal ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_hy06)

summary(final_fal_nb_hy06)

pred_fal_hy06 <- predict (
  final_fal_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_fal_hy06 <- ddbb_hy06$fal
rmse_fal_hy06 <- sqrt(mean((obs_fal_hy06 - pred_fal_hy06)^2))
n     <- length(obs_fal_hy06)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_fal_hy06 - pred_fal_hy06

rse_fal_hy06 <- sqrt( sum(resid^2) / (n - p) )



final_viv_nb_hy06 <- glm.nb(viv ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_hy06)

summary(final_viv_nb_hy06)

pred_viv_hy06 <- predict (
  final_viv_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_viv_hy06 <- ddbb_hy06$viv
rmse_viv_hy06 <- sqrt(mean((obs_viv_hy06 - pred_viv_hy06)^2))
n     <- length(obs_viv_hy06)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_viv_hy06 - pred_viv_hy06

rse_viv_hy06 <- sqrt( sum(resid^2) / (n - p) )


#Modelos finales cuenca Hy6: SFW district---------------------------------------

final_fal_nb_hy06 <- glm.nb(fal ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy06)

summary(final_fal_nb_hy06)

pred_fal_hy06 <- predict (
  final_fal_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_fal_hy06 <- ddbb_hy06$fal
rmse_fal_hy06 <- sqrt(mean((obs_fal_hy06 - pred_fal_hy06)^2))
n     <- length(obs_fal_hy06)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_fal_hy06 - pred_fal_hy06

rse_fal_hy06 <- sqrt( sum(resid^2) / (n - p) )



final_viv_nb_hy06 <- glm.nb(viv ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy06)

summary(final_viv_nb_hy06)

pred_viv_hy06 <- predict (
  final_viv_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_viv_hy06 <- ddbb_hy06$viv
rmse_viv_hy06 <- sqrt(mean((obs_viv_hy06 - pred_viv_hy06)^2))
n     <- length(obs_viv_hy06)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_viv_hy06 - pred_viv_hy06

rse_viv_hy06 <- sqrt( sum(resid^2) / (n - p) )

#Modelos finales cuenca Hy6: SFW ana -------------------------------------------

final_fal_nb_hy06 <- glm.nb(fal ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy06)

summary(final_fal_nb_hy06)

pred_fal_hy06 <- predict (
  final_fal_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_fal_hy06 <- ddbb_hy06$fal
rmse_fal_hy06 <- sqrt(mean((obs_fal_hy06 - pred_fal_hy06)^2))
n     <- length(obs_fal_hy06)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_fal_hy06 - pred_fal_hy06

rse_fal_hy06 <- sqrt( sum(resid^2) / (n - p) )



final_viv_nb_hy06 <- glm.nb(viv ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy06)

summary(final_viv_nb_hy06)

pred_viv_hy06 <- predict (
  final_viv_nb_hy06,
  newdata = ddbb_hy06,
  type = "response"
)
obs_viv_hy06 <- ddbb_hy06$viv
rmse_viv_hy06 <- sqrt(mean((obs_viv_hy06 - pred_viv_hy06)^2))
n     <- length(obs_viv_hy06)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy06))      # cantidad de parámetros estimados
resid <- obs_viv_hy06 - pred_viv_hy06

rse_viv_hy06 <- sqrt( sum(resid^2) / (n - p) )


# 2. Hydrobasin level 7 ---------------------------------------------------

sp.hy07 <- st_read(
  'rawdata/geometry.gpkg',
  layer = 'hydrobasin_07')|>
  mutate(new_id = 1:76) |> 
  mutate(hydroname = paste0(hydroname,'.07')) 

orden <- sp.hy07$hydroname

ddbb_hy07 <- st_read("data_for_model/hy7.gpkg") |> 
  mutate(hydroname = paste0(hydroname,'.07')) |> 
  mutate(new_id = as.integer(factor(hydroname,levels = orden,labels = 1:76))) %>%
  st_drop_geometry() %>%
  drop_na()

# Step forward Poisson ----------------------------------------------------

# Evaluando equidispersion
summary(ddbb_hy07$fal) # Mean 7.2
var(ddbb_hy07$fal) # Var 607.3

summary(ddbb_hy07$viv) # Mean 27.5
var(ddbb_hy07$viv) # Var 6349.6

# No se cumple el supuesto de equidispersión, se usará Binomial Negativo

fit_fal_nb_null <- glm.nb(fal ~ offset(log(pob)), data = ddbb_hy07)
fit_viv_nb_null <- glm.nb(viv ~ offset(log(pob)), data = ddbb_hy07)

fit_step_fal <- stepAIC(
  object = fit_fal_nb_null, 
  scope = list(lower = fal ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_fal)

fit_step_viv <- stepAIC(
  object = fit_viv_nb_null, 
  scope = list(lower = viv ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_viv)

# Modelo final cuenca Hy7: SFW Hy06 --------------------------------------------

final_fal_nb_hy07 <- glm.nb(fal ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                              t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                       data = ddbb_hy07)

summary(final_fal_nb_hy07)

pred_fal_hy07 <- predict (
  final_fal_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_fal_hy07 <- ddbb_hy07$fal
rmse_fal_hy07 <- sqrt(mean((obs_fal_hy07 - pred_fal_hy07)^2))
n     <- length(obs_fal_hy07)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_fal_hy07 - pred_fal_hy07

rse_fal_hy07 <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_hy07 <- glm.nb(viv ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                              t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                       data = ddbb_hy07)

summary(final_viv_nb_hy07)

pred_viv_hy07 <- predict (
  final_viv_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_viv_hy07 <- ddbb_hy07$viv
rmse_viv_hy07 <- sqrt(mean((obs_viv_hy07 - pred_viv_hy07)^2))
n     <- length(obs_viv_hy07)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_viv_hy07 - pred_viv_hy07

rse_viv_hy07<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca Hy7: SFW Hy07 --------------------------------------------

final_fal_nb_hy07 <- glm.nb(fal ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_hy07)

summary(final_fal_nb_hy07)

pred_fal_hy07 <- predict (
  final_fal_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_fal_hy07 <- ddbb_hy07$fal
rmse_fal_hy07 <- sqrt(mean((obs_fal_hy07 - pred_fal_hy07)^2))
n     <- length(obs_fal_hy07)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_fal_hy07 - pred_fal_hy07

rse_fal_hy07 <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_hy07 <- glm.nb(viv ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_hy07)

summary(final_viv_nb_hy07)

pred_viv_hy07 <- predict (
  final_viv_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_viv_hy07 <- ddbb_hy07$viv
rmse_viv_hy07 <- sqrt(mean((obs_viv_hy07 - pred_viv_hy07)^2))
n     <- length(obs_viv_hy07)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_viv_hy07 - pred_viv_hy07

rse_viv_hy07<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca Hy7: SFW dist --------------------------------------------

final_fal_nb_hy07 <- glm.nb(fal ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy07)

summary(final_fal_nb_hy07)

pred_fal_hy07 <- predict (
  final_fal_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_fal_hy07 <- ddbb_hy07$fal
rmse_fal_hy07 <- sqrt(mean((obs_fal_hy07 - pred_fal_hy07)^2))
n     <- length(obs_fal_hy07)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_fal_hy07 - pred_fal_hy07

rse_fal_hy07 <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_hy07 <- glm.nb(viv ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy07)

summary(final_viv_nb_hy07)

pred_viv_hy07 <- predict (
  final_viv_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_viv_hy07 <- ddbb_hy07$viv
rmse_viv_hy07 <- sqrt(mean((obs_viv_hy07 - pred_viv_hy07)^2))
n     <- length(obs_viv_hy07)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_viv_hy07 - pred_viv_hy07

rse_viv_hy07<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca Hy7: SFW ana  --------------------------------------------

final_fal_nb_hy07 <- glm.nb(fal ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy07)

summary(final_fal_nb_hy07)

pred_fal_hy07 <- predict (
  final_fal_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_fal_hy07 <- ddbb_hy07$fal
rmse_fal_hy07 <- sqrt(mean((obs_fal_hy07 - pred_fal_hy07)^2))
n     <- length(obs_fal_hy07)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_fal_hy07 - pred_fal_hy07

rse_fal_hy07 <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_hy07 <- glm.nb(viv ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_hy07)

summary(final_viv_nb_hy07)

pred_viv_hy07 <- predict (
  final_viv_nb_hy07,
  newdata = ddbb_hy07,
  type = "response"
)
obs_viv_hy07 <- ddbb_hy07$viv
rmse_viv_hy07 <- sqrt(mean((obs_viv_hy07 - pred_viv_hy07)^2))
n     <- length(obs_viv_hy07)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_hy07))      # cantidad de parámetros estimados
resid <- obs_viv_hy07 - pred_viv_hy07

rse_viv_hy07<- sqrt( sum(resid^2) / (n - p) )

# 3. District --------------------------------------------------------------

sp.dis <- st_read(
  'rawdata/geometry.gpkg',
  layer = 'districts')|>
  mutate(new_id = 1:53)

orden <- sp.dis$hydroname

ddbb_dist <- st_read("data_for_model/dist.gpkg") |> 
  st_drop_geometry() %>%
  drop_na()

# Step forward Poisson ----------------------------------------------------

# Evaluando equidispersion
summary(ddbb_dist$fal) # Mean 10.4
var(ddbb_dist$fal) # Var 1210.8

summary(ddbb_dist$viv) # Mean 39.4
var(ddbb_dist$viv) # Var 7811.3

# No se cumple el supuesto de equidispersión, se usará Binomial Negativo

fit_fal_nb_null <- glm.nb(fal ~ offset(log(pob)), data = ddbb_dist)
fit_viv_nb_null <- glm.nb(viv ~ offset(log(pob)), data = ddbb_dist)

fit_step_fal <- stepAIC(
  object = fit_fal_nb_null, 
  scope = list(lower = fal ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_fal)

fit_step_viv <- stepAIC(
  object = fit_viv_nb_null, 
  scope = list(lower = viv ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_viv)

# Modelo final cuenca Distritos ------------------------------------------------

final_fal_nb <- glm.nb(fal ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                         evi + ro + lhveg + offset(log(pob)), 
                       data = ddbb_dist)

summary(final_fal_nb)

pred_fal <- predict (
  final_fal_nb,
  newdata = ddbb_hy06,
  type = "response"
)
obs_fal <- ddbb_hy06$fal
rmse_fal <- sqrt(mean((obs_fal - pred_fal)^2))
n     <- length(obs_fal)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb))      # cantidad de parámetros estimados
resid <- obs_fal - pred_fal

rse_fal <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb <- glm.nb(viv ~ lwveg + so + tmax + shum + t2m + lstd + ndvi +
                         evi + ro + lhveg + offset(log(pob)), 
                       data = ddbb_dist)

summary(final_viv_nb)

pred_viv <- predict (
  final_viv_nb,
  newdata = ddbb_hy06,
  type = "response"
)
obs_viv <- ddbb_hy06$viv
rmse_viv <- sqrt(mean((obs_viv - pred_viv)^2))
n     <- length(obs_viv)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb))      # cantidad de parámetros estimados
resid <- obs_viv - pred_viv

rse_viv <- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca dist: SFW Hy06 --------------------------------------------

final_fal_nb_dist <- glm.nb(fal ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                              t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_fal_nb_dist)

pred_fal_dist <- predict (
  final_fal_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_fal_dist <- ddbb_dist$fal
rmse_fal_dist <- sqrt(mean((obs_fal_dist - pred_fal_dist)^2))
n     <- length(obs_fal_dist)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_dist))      # cantidad de parámetros estimados
resid <- obs_fal_dist - pred_fal_dist

rse_fal_dist <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_dist <- glm.nb(viv ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                              t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_viv_nb_dist)

pred_viv_dist <- predict (
  final_viv_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_viv_dist <- ddbb_dist$viv
rmse_viv_dist <- sqrt(mean((obs_viv_dist - pred_viv_dist)^2))
n     <- length(obs_viv_dist)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_dist))      # cantidad de parámetros estimados
resid <- obs_viv_dist - pred_viv_dist

rse_viv_dist<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca dist: SFW hy07 --------------------------------------------

final_fal_nb_dist <- glm.nb(fal ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_fal_nb_dist)

pred_fal_dist <- predict (
  final_fal_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_fal_dist <- ddbb_dist$fal
rmse_fal_dist <- sqrt(mean((obs_fal_dist - pred_fal_dist)^2))
n     <- length(obs_fal_dist)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_dist))      # cantidad de parámetros estimados
resid <- obs_fal_dist - pred_fal_dist

rse_fal_dist <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_dist <- glm.nb(viv ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_viv_nb_dist)

pred_viv_dist <- predict (
  final_viv_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_viv_dist <- ddbb_dist$viv
rmse_viv_dist <- sqrt(mean((obs_viv_dist - pred_viv_dist)^2))
n     <- length(obs_viv_dist)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_dist))      # cantidad de parámetros estimados
resid <- obs_viv_dist - pred_viv_dist

rse_viv_dist<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca dist: SFW dist --------------------------------------------

final_fal_nb_dist <- glm.nb(fal ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_fal_nb_dist)

pred_fal_dist <- predict (
  final_fal_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_fal_dist <- ddbb_dist$fal
rmse_fal_dist <- sqrt(mean((obs_fal_dist - pred_fal_dist)^2))
n     <- length(obs_fal_dist)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_dist))      # cantidad de parámetros estimados
resid <- obs_fal_dist - pred_fal_dist

rse_fal_dist <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_dist <- glm.nb(viv ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_viv_nb_dist)

pred_viv_dist <- predict (
  final_viv_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_viv_dist <- ddbb_dist$viv
rmse_viv_dist <- sqrt(mean((obs_viv_dist - pred_viv_dist)^2))
n     <- length(obs_viv_dist)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_dist))      # cantidad de parámetros estimados
resid <- obs_viv_dist - pred_viv_dist

rse_viv_dist<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca dist: SFW ana  --------------------------------------------

final_fal_nb_dist <- glm.nb(fal ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_fal_nb_dist)

pred_fal_dist <- predict (
  final_fal_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_fal_dist <- ddbb_dist$fal
rmse_fal_dist <- sqrt(mean((obs_fal_dist - pred_fal_dist)^2))
n     <- length(obs_fal_dist)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_dist))      # cantidad de parámetros estimados
resid <- obs_fal_dist - pred_fal_dist

rse_fal_dist <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_dist <- glm.nb(viv ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_dist)

summary(final_viv_nb_dist)

pred_viv_dist <- predict (
  final_viv_nb_dist,
  newdata = ddbb_dist,
  type = "response"
)
obs_viv_dist <- ddbb_dist$viv
rmse_viv_dist <- sqrt(mean((obs_viv_dist - pred_viv_dist)^2))
n     <- length(obs_viv_dist)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_dist))      # cantidad de parámetros estimados
resid <- obs_viv_dist - pred_viv_dist

rse_viv_dist<- sqrt( sum(resid^2) / (n - p) )


# 4. Ana -----------------------------------------------------------------

sp.ana <- st_read(
  'rawdata/geometry.gpkg',
  layer = 'hydroana')|>
  mutate(new_id = 1:33)

orden <- sp.ana$hydroname

ddbb_ana <- st_read("data_for_model/dist.gpkg") |> 
  st_drop_geometry() %>%
  drop_na()

# Step forward Poisson ----------------------------------------------------

# Evaluando equidispersion
summary(ddbb_ana$fal) # Mean 10.4
var(ddbb_ana$fal) # Var 1210.8

summary(ddbb_ana$viv) # Mean 39.4
var(ddbb_ana$viv) # Var 7811.3

# No se cumple el supuesto de equidispersión, se usará Binomial Negativo

fit_fal_nb_null <- glm.nb(fal ~ offset(log(pob)), data = ddbb_ana)
fit_viv_nb_null <- glm.nb(viv ~ offset(log(pob)), data = ddbb_ana)

fit_step_fal <- stepAIC(
  object = fit_fal_nb_null, 
  scope = list(lower = fal ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_fal)

fit_step_viv <- stepAIC(
  object = fit_viv_nb_null, 
  scope = list(lower = viv ~ offset(log(pob)), upper = fal ~ ndvi + evi + fpar + lwveg + lhveg + pr + ro + etp + shum + so + 
                 tmax + tmin + lstd + lstn + t2m + offset(log(pob))),
  direction = "forward"
)

summary(fit_step_viv)

# Modelo final cuenca Ana  -----------------------------------------------------

final_fal_nb <- glm.nb(fal ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                         evi + ro + lhveg + offset(log(pob)), 
                       data = ddbb_ana)

summary(final_fal_nb)

pred_fal <- predict (
  final_fal_nb,
  newdata = ddbb_hy06,
  type = "response"
)
obs_fal <- ddbb_hy06$fal
rmse_fal <- sqrt(mean((obs_fal - pred_fal)^2))
n     <- length(obs_fal)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb))      # cantidad de parámetros estimados
resid <- obs_fal - pred_fal

rse_fal <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb <- glm.nb(viv ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi +
                         evi + ro + lhveg + offset(log(pob)), 
                       data = ddbb_ana)

summary(final_viv_nb)

pred_viv <- predict (
  final_viv_nb,
  newdata = ddbb_hy06,
  type = "response"
)
obs_viv <- ddbb_hy06$viv
rmse_viv <- sqrt(mean((obs_viv - pred_viv)^2))
n     <- length(obs_viv)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb))      # cantidad de parámetros estimados
resid <- obs_viv - pred_viv

rse_viv <- sqrt( sum(resid^2) / (n - p) )


# Modelo final cuenca ana: SFW Hy06 --------------------------------------------

final_fal_nb_ana <- glm.nb(fal ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                              t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_fal_nb_ana)

pred_fal_ana <- predict (
  final_fal_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_fal_ana <- ddbb_ana$fal
rmse_fal_ana <- sqrt(mean((obs_fal_ana - pred_fal_ana)^2))
n     <- length(obs_fal_ana)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_ana))      # cantidad de parámetros estimados
resid <- obs_fal_ana - pred_fal_ana

rse_fal_ana <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_ana <- glm.nb(viv ~ tmax + lwveg + so + fpar + shum + tmin + pr + 
                              t2m + lhveg + etp + evi + lstn + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_viv_nb_ana)

pred_viv_ana <- predict (
  final_viv_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_viv_ana <- ddbb_ana$viv
rmse_viv_ana <- sqrt(mean((obs_viv_ana - pred_viv_ana)^2))
n     <- length(obs_viv_ana)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_ana))      # cantidad de parámetros estimados
resid <- obs_viv_ana - pred_viv_ana

rse_viv_ana<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca ana: SFW hy07 --------------------------------------------

final_fal_nb_ana <- glm.nb(fal ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_fal_nb_ana)

pred_fal_ana <- predict (
  final_fal_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_fal_ana <- ddbb_ana$fal
rmse_fal_ana <- sqrt(mean((obs_fal_ana - pred_fal_ana)^2))
n     <- length(obs_fal_ana)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_ana))      # cantidad de parámetros estimados
resid <- obs_fal_ana - pred_fal_ana

rse_fal_ana <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_ana <- glm.nb(viv ~ lwveg + so + shum + t2m + ndvi + etp + lstd + 
                              fpar + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_viv_nb_ana)

pred_viv_ana <- predict (
  final_viv_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_viv_ana <- ddbb_ana$viv
rmse_viv_ana <- sqrt(mean((obs_viv_ana - pred_viv_ana)^2))
n     <- length(obs_viv_ana)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_ana))      # cantidad de parámetros estimados
resid <- obs_viv_ana - pred_viv_ana

rse_viv_ana<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca ana: SFW dist --------------------------------------------

final_fal_nb_ana <- glm.nb(fal ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_fal_nb_ana)

pred_fal_ana <- predict (
  final_fal_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_fal_ana <- ddbb_ana$fal
rmse_fal_ana <- sqrt(mean((obs_fal_ana - pred_fal_ana)^2))
n     <- length(obs_fal_ana)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_ana))      # cantidad de parámetros estimados
resid <- obs_fal_ana - pred_fal_ana

rse_fal_ana <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_ana <- glm.nb(viv ~ lwveg + so + tmax + shum + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_viv_nb_ana)

pred_viv_ana <- predict (
  final_viv_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_viv_ana <- ddbb_ana$viv
rmse_viv_ana <- sqrt(mean((obs_viv_ana - pred_viv_ana)^2))
n     <- length(obs_viv_ana)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_ana))      # cantidad de parámetros estimados
resid <- obs_viv_ana - pred_viv_ana

rse_viv_ana<- sqrt( sum(resid^2) / (n - p) )

# Modelo final cuenca ana: SFW ana  --------------------------------------------

final_fal_nb_ana <- glm.nb(fal ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_fal_nb_ana)

pred_fal_ana <- predict (
  final_fal_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_fal_ana <- ddbb_ana$fal
rmse_fal_ana <- sqrt(mean((obs_fal_ana - pred_fal_ana)^2))
n     <- length(obs_fal_ana)                 # cantidad de observaciones
p     <- length(coef(final_fal_nb_ana))      # cantidad de parámetros estimados
resid <- obs_fal_ana - pred_fal_ana

rse_fal_ana <- sqrt( sum(resid^2) / (n - p) )

final_viv_nb_ana <- glm.nb(viv ~ lwveg + so + tmax + shum + etp + t2m + lstd + ndvi + 
                              evi + ro + lhveg + offset(log(pob)), 
                            data = ddbb_ana)

summary(final_viv_nb_ana)

pred_viv_ana <- predict (
  final_viv_nb_ana,
  newdata = ddbb_ana,
  type = "response"
)
obs_viv_ana <- ddbb_ana$viv
rmse_viv_ana <- sqrt(mean((obs_viv_ana - pred_viv_ana)^2))
n     <- length(obs_viv_ana)                 # cantidad de observaciones
p     <- length(coef(final_viv_nb_ana))      # cantidad de parámetros estimados
resid <- obs_viv_ana - pred_viv_ana

rse_viv_ana<- sqrt( sum(resid^2) / (n - p) )
