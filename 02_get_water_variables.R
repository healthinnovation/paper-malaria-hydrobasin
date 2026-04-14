library(rgee)
library(sf)
library(tidyverse)
library(mapview)
library(formattable)
ee_Initialize(quiet = TRUE)

# 1. Reading spatial data -------------------------------------------------
districts <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'districts') |> 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

hydro_ana <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'hydroana') |> 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

hydro_06 <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'hydrobasin_06') |> 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

hydro_07 <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'hydrobasin_07') |> 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

# 2. Spatial data to earth engine objects ---------------------------------
ee_districts <- districts |> 
  select(codigo) |> 
  sf_as_ee(quiet = TRUE) 

ee_hydro_ana <- hydro_ana |> 
  select(hydroname) |> 
  sf_as_ee(quiet = TRUE)

ee_hydro_06  <- hydro_06  |> 
  select(hydroname) |> 
  sf_as_ee(quiet = TRUE)

ee_hydro_07  <- hydro_07  |> 
  select(hydroname) |> 
  sf_as_ee(quiet = TRUE)

## 2.1 Parameters
start_date <- 2009
end_date <- 2022

## 2.2 Water variables
ee_aet  <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('aet')

ee_def <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('def')

ee_pdsi <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('pdsi')

ee_pet <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('pet')

ee_pr <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('pr')

ee_ro <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('ro')

ee_soil <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('soil')

## 2.3 Water functions 
ee_r_aet <- function(){
  ee_reducer_ic <- ee_aet$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_def <- function(){
  ee_reducer_ic <- ee_def$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_pdsi <- function(){
  ee_reducer_ic <- ee_pdsi$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.01)
  return(ee_reducer_ic)
}

ee_r_pet <- function(){
  ee_reducer_ic <- ee_pet$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_pr <- function(){
  ee_reducer_ic <- ee_pr$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_ro <- function(){
  ee_reducer_ic <- ee_ro$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_soil <- function(){
  ee_reducer_ic <- ee_soil$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

# 3. Extracting water variables in the districts --------------------------

## 3.1 Actual evapotranspiration
ee_r_aet() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_aet

## 3.2 Climate water deficit
ee_r_def() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_def

## 3.3 Palmer Drought Severity Index
ee_r_pdsi() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_pdsi

## 3.4 Reference evapotranspiration
ee_r_pet() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_pet

## 3.5 Precipitation accumulation
ee_r_pr() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_pr

## 3.6 Runoff
ee_r_ro() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_ro

## 3.7 Soil Mousture
ee_r_soil() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_soil


## 3.8 Dataset final
districts |> 
  left_join(
    y = districts_aet,
    by = 'codigo') |> 
  left_join(
    y = districts_def,
    by = 'codigo') |>
  left_join(
    y = districts_pdsi,
    by = 'codigo') |> 
  left_join(
    y = districts_pet,
    by = 'codigo') |>
  left_join(
    y = districts_pr,
    by = 'codigo') |>
  left_join(
    y = districts_ro,
    by = 'codigo') |>
  left_join(
    y = districts_soil,
    by = 'codigo') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_aet:X202212_soil,
    names_to = 'variables',
    values_to = 'valor') |>
  separate(
    col = variables,
    into = c('fecha','variable'),
    sep = '_') |> 
  mutate(
    fecha = str_extract(pattern = '\\d+',fecha),
    fecha = case_when(nchar(fecha) == 6 ~ paste0(fecha,'01'),TRUE ~ fecha),
    fecha = ymd(fecha),
    year = year(fecha), 
    month = month(fecha),
    variable = str_to_lower(variable)) |> 
  relocate(c('variable','valor'),.after = month) -> districts_water

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/water')){dir.create('output/water')}
write_csv(districts_water,'output/water/districts_water.csv')

# 4. Extracting water variables in Hydroana -------------------------------

## 4.1 Actual evapotranspiration
ee_r_aet() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_aet

## 4.2 Climate water deficit
ee_r_def() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_def

## 4.3 Palmer Drought Severity Index
ee_r_pdsi() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_pdsi

## 4.4 Reference evapotranspiration
ee_r_pet() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_pet

## 4.5 Precipitation accumulation
ee_r_pr() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_pr

## 4.6 Runoff
ee_r_ro() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_ro
 
## 4.7 Soil Moiusture
ee_r_soil() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_soil

## 4.8 Dataset final
hydro_ana |> 
  left_join(
    y = hydroana_aet,
    by = 'hydroname') |> 
  left_join(
    y = hydroana_def,
    by = 'hydroname') |>
  left_join(
    y = hydroana_pdsi, 
    by = 'hydroname') |> 
  left_join(
    y = hydroana_pet,
    by = 'hydroname') |>
  left_join(
    y = hydroana_pr,
    by = 'hydroname') |>
  left_join(
    y = hydroana_ro,
    by = 'hydroname') |>
  left_join(
    y = hydroana_soil,
    by = 'hydroname') |>
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_aet:X202212_soil,
    names_to = 'variables',
    values_to = 'valor') |>
  separate(
    col = variables,
    into = c('fecha','variable'),
    sep = '_') |> 
  mutate(
    fecha = str_extract(pattern = '\\d+',fecha),
    fecha = case_when(nchar(fecha) == 6 ~ paste0(fecha,'01'),TRUE ~ fecha),
    fecha = ymd(fecha),
    year = year(fecha), 
    month = month(fecha),
    variable = str_to_lower(variable)) |> 
  relocate(c('variable','valor'),.after = month) -> hydroana_water

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/water')){dir.create('output/water')}
write_csv(hydroana_water,'output/water/hydroana_water.csv')

# 5. Extracting water variables in Hydrobasin level 6 ---------------------
## 5.1 Actual evapotranspiration
ee_r_aet() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_aet

## 5.2 Climate water deficit
ee_r_def() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_def

## 5.3 Palmer Drought Severity Index
ee_r_pdsi() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_pdsi

## 5.4 Reference evapotranspiration
ee_r_pet() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_pet

## 5.5 Precipitation accumulation
ee_r_pr() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_pr

## 5.6 Runoff
ee_r_ro() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_ro

## 5.7 Soil Moisture
ee_r_soil() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_soil

## 5.8 Dataset final
hydro_06 |> 
  left_join(
    y = hydro06_aet,
    by = 'hydroname') |> 
  left_join(
    y = hydro06_def,
    by = 'hydroname') |>
  left_join(
    y = hydro06_pdsi,
    by = 'hydroname') |> 
  left_join(
    y = hydro06_pet,
    by = 'hydroname') |>
  left_join(
    y = hydro06_pr,
    by = 'hydroname') |>
  left_join(
    y = hydro06_ro,
    by = 'hydroname') |>
  left_join(
    y = hydro06_soil,
    by = 'hydroname') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_aet:X202212_soil,
    names_to = 'variables',
    values_to = 'valor') |>
  separate(
    col = variables,
    into = c('fecha','variable'),
    sep = '_') |> 
  mutate(
    fecha = str_extract(pattern = '\\d+',fecha),
    fecha = case_when(nchar(fecha) == 6 ~ paste0(fecha,'01'),TRUE ~ fecha),
    fecha = ymd(fecha),
    year = year(fecha), 
    month = month(fecha),
    variable = str_to_lower(variable)) |> 
  relocate(c('variable','valor'),.after = month) -> hydro06_water

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/water')){dir.create('output/water')}
write_csv(hydro06_water,'output/water/hydro06_water.csv')

# 6. Extracting water variables in hydrobasin level 7 ---------------------
## 6.1 Actual evapotranspiration
ee_r_aet() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_aet

## 6.2 Climate water deficit
ee_r_def() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_def

## 6.3 Palmer Drought Severity Index
ee_r_pdsi() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_pdsi

## 6.4 Reference evapotranspiration
ee_r_pet() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_pet

## 6.5 Precipitation accumulation
ee_r_pr() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_pr

## 6.6 Runoff
ee_r_ro() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_ro

## 6.6 Soil Moisture
ee_r_soil() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_soil

## 6. Dataset final
hydro_07 |> 
  left_join(
    y = hydro07_aet,
    by = 'hydroname') |> 
  left_join(
    y = hydro07_def,
    by = 'hydroname') |>
  left_join(
    y = hydro07_pdsi,
    by = 'hydroname') |> 
  left_join(
    y = hydro07_pet,
    by = 'hydroname') |>
  left_join(
    y = hydro07_pr,
    by = 'hydroname') |>
  left_join(
    y = hydro07_ro,
    by = 'hydroname') |>
  left_join(
    y = hydro07_soil,
    by = 'hydroname') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_aet:X202212_soil,
    names_to = 'variables',
    values_to = 'valor') |>
  separate(
    col = variables,
    into = c('fecha','variable'),
    sep = '_') |> 
  mutate(
    fecha = str_extract(pattern = '\\d+',fecha),
    fecha = case_when(nchar(fecha) == 6 ~ paste0(fecha,'01'),TRUE ~ fecha),
    fecha = ymd(fecha),
    year = year(fecha), 
    month = month(fecha),
    variable = str_to_lower(variable)) |> 
  relocate(c('variable','valor'),.after = month) -> hydro07_water

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/water')){dir.create('output/water')}
write_csv(hydro07_water,'output/water/hydro07_water.csv')