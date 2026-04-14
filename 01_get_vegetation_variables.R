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

# 2. spatial data to earth engine object ----------------------------------
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

## 2.2 Vegetation variables
ee_leaf_area_index_high  <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_high_vegetation')

ee_leaf_area_index_high_min <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_high_vegetation_min')

ee_leaf_area_index_high_max <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_high_vegetation_max')

ee_leaf_area_index_low  <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_low_vegetation')

ee_leaf_area_index_low_min <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_low_vegetation_min')

ee_leaf_area_index_low_max <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_low_vegetation_max')

ee_evi <- ee$ImageCollection('MODIS/061/MOD13A3') |> 
  ee$ImageCollection$select('EVI')

ee_ndvi <- ee$ImageCollection('MODIS/061/MOD13A3') |> 
  ee$ImageCollection$select('NDVI')

## 2.3 Vegetation variables in R functions
ee_r_leaf_area_index_high <- function(){
  ee_reducer_ic <- ee_leaf_area_index_high$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_high_max <- function(){
  ee_reducer_ic <- ee_leaf_area_index_high_max$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_high_min <- function(){
  ee_reducer_ic <- ee_leaf_area_index_high_min$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_low <- function(){
  ee_reducer_ic <- ee_leaf_area_index_low$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_low_max <- function(){
  ee_reducer_ic <- ee_leaf_area_index_low_max$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_low_min <- function(){
  ee_reducer_ic <- ee_leaf_area_index_low_min$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_evi <- function(){
  ee_reducer_ic <- ee_evi$filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.0001)
  return(ee_reducer_ic)
}

ee_r_ndvi <- function(){
  ee_reducer_ic <- ee_ndvi$filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.0001)
  return(ee_reducer_ic)
}

# 3. Extracting vegetation variables in the Districts ---------------------

## 3.1 Leaf area index high vegetation
ee_r_leaf_area_index_high() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_high

names(districts_leaf_area_index_high) <- sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_high)))

## 3.2 Leaf area index high vegetation min
ee_r_leaf_area_index_high_min() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_high_min

names(districts_leaf_area_index_high_min) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_high_min)))

## 3.3 Leaf area index high vegetation max
ee_r_leaf_area_index_high_max() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_high_max

names(districts_leaf_area_index_high_max) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_high_max)))

## 3.4 Leaf area index low vegetation
ee_r_leaf_area_index_low() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_low

names(districts_leaf_area_index_low) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_low)))

## 3.5 Leaf area index low vegetation min
ee_r_leaf_area_index_low_min() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_low_min

names(districts_leaf_area_index_low_min) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_low_min)))

## 3.6 Leaf area index low vegetation max
ee_r_leaf_area_index_low_max() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_low_max

names(districts_leaf_area_index_low_max) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_low_max)))

## 3.7 EVI
ee_r_evi() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> districts_evi

names(districts_evi) <-  sub(
  sprintf("(.{%d})", 9),
  "\\1_",
  gsub('_','',names(districts_evi)))

## 3.8 NDVI
ee_r_ndvi() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> districts_ndvi

names(districts_ndvi) <-  sub(
  sprintf("(.{%d})", 9),
  "\\1_",
  gsub('_','',names(districts_ndvi)))

## 3.9 Dataset final
districts |> 
  left_join(
    y = districts_leaf_area_index_high,
    by = 'codigo') |> 
  left_join(
    y = districts_leaf_area_index_high_min,
    by = 'codigo') |>
  left_join(
    y = districts_leaf_area_index_high_max,
    by = 'codigo') |> 
  left_join(
    y = districts_leaf_area_index_low,
    by = 'codigo') |>
  left_join(
    y = districts_leaf_area_index_low_min,
    by = 'codigo') |>
  left_join(
    y = districts_leaf_area_index_low_max,
    by = 'codigo') |>
  left_join(
    y = districts_evi,
    by = 'codigo') |>
  left_join(
    y = districts_ndvi,
    by = 'codigo') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_leafareaindexhighvegetation:X20221201_NDVI,
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
  relocate(c('variable','valor'),.after = month) |> 
  mutate(variable = str_to_lower(variable)) -> districts_vegetation

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/vegetation')){dir.create('output/vegetation')}
write_csv(districts_vegetation,'output/vegetation/districts_vegetation.csv')

# 4. Extracting vegetation variables in Hydroana --------------------------
## 4.1 Leaf area index high vegetation
ee_r_leaf_area_index_high() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_leaf_area_index_high

names(hydroana_leaf_area_index_high) <- c(
  names(hydroana_leaf_area_index_high)[1],
  sub(
    sprintf("(.{%d})", 7),"\\1_",
    gsub('_','',names(hydroana_leaf_area_index_high)[-1])
    )
  )

## 4.2 Leaf area index high vegetation min
ee_r_leaf_area_index_high_min() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_leaf_area_index_high_min

names(hydroana_leaf_area_index_high_min) <-  c(
  names(hydroana_leaf_area_index_high_min)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydroana_leaf_area_index_high_min)[-1])
    )
  )

## 4.3 Leaf area index high vegetation max
ee_r_leaf_area_index_high_max() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_leaf_area_index_high_max

names(hydroana_leaf_area_index_high_max) <-  c(
  names(hydroana_leaf_area_index_high_max)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydroana_leaf_area_index_high_max)[-1])
    )
  )

## 4.4 Leaf area index low vegetation
ee_r_leaf_area_index_low() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_leaf_area_index_low

names(hydroana_leaf_area_index_low) <-  c(
  names(hydroana_leaf_area_index_low)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydroana_leaf_area_index_low)[-1])
    )
  )

## 4.5 Leaf area index low vegetation min
ee_r_leaf_area_index_low_min() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_leaf_area_index_low_min

names(hydroana_leaf_area_index_low_min) <-  c(
  names(hydroana_leaf_area_index_low_min)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydroana_leaf_area_index_low_min)[-1])
    )
  )

## 4.6 Leaf area index low vegetation max
ee_r_leaf_area_index_low_max() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydroana_leaf_area_index_low_max

names(hydroana_leaf_area_index_low_max) <-  c(
  names(hydroana_leaf_area_index_low_max)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydroana_leaf_area_index_low_max)[-1])
    )
  )

## 4.7 EVI
ee_r_evi() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> hydroana_evi

names(hydroana_evi) <-  c(
  names(hydroana_evi)[1],
  sub(
    sprintf("(.{%d})", 9),
    "\\1_",
    gsub('_','',names(hydroana_evi)[-1])
    )
  )

## 4.8 NDVI
ee_r_ndvi() |> 
  ee_extract(
    y = ee_hydro_ana,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> hydroana_ndvi

names(hydroana_ndvi) <- c(
  names(hydroana_ndvi)[1],
  sub(
    sprintf("(.{%d})", 9),
    "\\1_",
    gsub('_','',names(hydroana_ndvi)[-1])
    )
  ) 

## 4.9 Dataset final
hydro_ana |> 
  left_join(
    y = hydroana_leaf_area_index_high, 
    by = 'hydroname') |> 
  left_join(
    y = hydroana_leaf_area_index_high_min, 
    by = 'hydroname') |>
  left_join(
    y = hydroana_leaf_area_index_high_max, 
    by = 'hydroname') |> 
  left_join(
    y = hydroana_leaf_area_index_low,
    by = 'hydroname') |>
  left_join(
    y = hydroana_leaf_area_index_low_min,
    by = 'hydroname') |>
  left_join(
    y = hydroana_leaf_area_index_low_max,
    by = 'hydroname') |>
  left_join(
    y = hydroana_evi,
    by = 'hydroname') |>
  left_join(
    y = hydroana_ndvi,
    by = 'hydroname') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols =  X200901_leafareaindexhighvegetation:X20221201_NDVI,
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
  relocate(c('variable','valor'),.after = month) |> 
  mutate(variable = str_to_lower(variable)) -> hydroana_vegetation

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/vegetation')){dir.create('output/vegetation')}
write_csv(hydroana_vegetation,'output/vegetation/hydroana_vegetation.csv')

# 5. Extracting vegetation variables in Hydrobasin level 6 ----------------
## 5.1 Leaf area index high vegetation
ee_r_leaf_area_index_high() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_leaf_area_index_high

names(hydro06_leaf_area_index_high) <- c(
  names(hydro06_leaf_area_index_high)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydro06_leaf_area_index_high)[-1])
    )
  )

## 5.2 Leaf area index high vegetation min
ee_r_leaf_area_index_high_min() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_leaf_area_index_high_min

names(hydro06_leaf_area_index_high_min) <-  c(
  names(hydro06_leaf_area_index_high_min)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    sub('_','',names(hydro06_leaf_area_index_high_min)[-1])
    )
  )

## 5.3 Leaf area index high vegetation max
ee_r_leaf_area_index_high_max() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_leaf_area_index_high_max
names(hydro06_leaf_area_index_high_max) <-  c(names(hydro06_leaf_area_index_high_max)[1],sub(sprintf("(.{%d})", 7), "\\1_", gsub('_','',names(hydro06_leaf_area_index_high_max)[-1])))

## 5.4 Leaf area index low vegetation
ee_r_leaf_area_index_low() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_leaf_area_index_low
names(hydro06_leaf_area_index_low) <-  c(names(hydro06_leaf_area_index_low)[1],sub(sprintf("(.{%d})", 7), "\\1_", gsub('_','',names(hydro06_leaf_area_index_low)[-1])))

## 5.5 Leaf area index low vegetation min
ee_r_leaf_area_index_low_min() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_leaf_area_index_low_min
names(hydro06_leaf_area_index_low_min) <-  c(names(hydro06_leaf_area_index_low_min)[1],sub(sprintf("(.{%d})", 7), "\\1_", gsub('_','',names(hydro06_leaf_area_index_low_min)[-1])))
## 5.6 Leaf area index low vegetation max
ee_r_leaf_area_index_low_max() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro06_leaf_area_index_low_max
names(hydro06_leaf_area_index_low_max) <-  c(names(hydro06_leaf_area_index_low_max)[1],sub(sprintf("(.{%d})", 7), "\\1_", gsub('_','',names(hydro06_leaf_area_index_low_max)[-1])))

## 5.7 EVI
ee_r_evi() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> hydro06_evi
names(hydro06_evi) <-  c(names(hydro06_evi)[1],sub(sprintf("(.{%d})", 9), "\\1_", gsub('_','',names(hydro06_evi)[-1])))

## 5.8 NDVI
ee_r_ndvi() |> 
  ee_extract(
    y = ee_hydro_06,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> hydro06_ndvi
names(hydro06_ndvi) <- c(names(hydro06_ndvi)[1],sub(sprintf("(.{%d})", 9), "\\1_", gsub('_','',names(hydro06_ndvi)[-1]))) 

## 5.9 Dataset final
hydro_06 |> 
  left_join(y = hydro06_leaf_area_index_high, by = 'hydroname') |> 
  left_join(y = hydro06_leaf_area_index_high_min, by = 'hydroname') |>
  left_join(y = hydro06_leaf_area_index_high_max, by = 'hydroname') |> 
  left_join(y = hydro06_leaf_area_index_low,by = 'hydroname') |>
  left_join(y = hydro06_leaf_area_index_low_min,by = 'hydroname') |>
  left_join(y = hydro06_leaf_area_index_low_max,by = 'hydroname') |>
  left_join(y = hydro06_evi,by = 'hydroname') |>
  left_join(y = hydro06_ndvi,by = 'hydroname') |> 
  st_drop_geometry() |> 
  pivot_longer(cols =  X200901_leafareaindexhighvegetation:X20221201_NDVI,names_to = 'variables',values_to = 'valor') |>
  separate(col = variables,into = c('fecha','variable'),sep = '_') |> 
  mutate(
    fecha = str_extract(pattern = '\\d+',fecha),
    fecha = case_when(nchar(fecha) == 6 ~ paste0(fecha,'01'),TRUE ~ fecha),
    fecha = ymd(fecha),
    year = year(fecha), 
    month = month(fecha),
    variable = str_to_lower(variable)) |> 
  relocate(c('variable','valor'),.after = month) |> 
  mutate(variable = str_to_lower(variable))-> hydro06_vegetation

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/vegetation')){dir.create('output/vegetation')}
write_csv(hydro06_vegetation,'output/vegetation/hydro06_vegetation.csv')

# 6. Extracting vegetation variables in Hydrobasin level 7 ----------------

## 6.1 Leaf area index high vegetation
ee_r_leaf_area_index_high() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_leaf_area_index_high

names(hydro07_leaf_area_index_high) <- c(
  names(hydro07_leaf_area_index_high)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydro07_leaf_area_index_high)[-1])
    )
  )

## 6.2 Leaf area index high vegetation min
ee_r_leaf_area_index_high_min() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_leaf_area_index_high_min

names(hydro07_leaf_area_index_high_min) <-  c(
  names(hydro07_leaf_area_index_high_min)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydro07_leaf_area_index_high_min)[-1])
    )
  )

## 6.3 Leaf area index high vegetation max
ee_r_leaf_area_index_high_max() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_leaf_area_index_high_max

names(hydro07_leaf_area_index_high_max) <-  c(
  names(hydro07_leaf_area_index_high_max)[1],
  sub(sprintf("(.{%d})", 7),
      "\\1_",
      gsub('_','',names(hydro07_leaf_area_index_high_max)[-1])
      )
  )

## 6.4 Leaf area index low vegetation
ee_r_leaf_area_index_low() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_leaf_area_index_low
names(hydro07_leaf_area_index_low) <-  c(
  names(hydro07_leaf_area_index_low)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydro07_leaf_area_index_low)[-1])
    )
  )

## 6.5 Leaf area index low vegetation min
ee_r_leaf_area_index_low_min() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_leaf_area_index_low_min

names(hydro07_leaf_area_index_low_min) <-  c(
  names(hydro07_leaf_area_index_low_min)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydro07_leaf_area_index_low_min)[-1])
    )
  )

## 6.6 Leaf area index low vegetation max
ee_r_leaf_area_index_low_max() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> hydro07_leaf_area_index_low_max

names(hydro07_leaf_area_index_low_max) <-  c(
  names(hydro07_leaf_area_index_low_max)[1],
  sub(
    sprintf("(.{%d})", 7),
    "\\1_",
    gsub('_','',names(hydro07_leaf_area_index_low_max)[-1])
    )
  )

## 6.7 EVI
ee_r_evi() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> hydro07_evi

names(hydro07_evi) <-  c(
  names(hydro07_evi)[1],
  sub(
    sprintf("(.{%d})", 9),
    "\\1_",
    gsub('_','',names(hydro07_evi)[-1])
    )
  )

## 6.8 NDVI
ee_r_ndvi() |> 
  ee_extract(
    y = ee_hydro_07,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> hydro07_ndvi

names(hydro07_ndvi) <- c(
  names(hydro07_ndvi)[1],
  sub(
    sprintf("(.{%d})", 9),
    "\\1_",
    gsub('_','',names(hydro07_ndvi)[-1])
    )
  ) 

## 6.9 Dataset final
hydro_07 |> 
  left_join(
    y = hydro07_leaf_area_index_high, 
    by = 'hydroname') |> 
  left_join(
    y = hydro07_leaf_area_index_high_min,
    by = 'hydroname') |>
  left_join(
    y = hydro07_leaf_area_index_high_max, 
    by = 'hydroname') |> 
  left_join(
    y = hydro07_leaf_area_index_low,
    by = 'hydroname') |>
  left_join(
    y = hydro07_leaf_area_index_low_min,
    by = 'hydroname') |>
  left_join(
    y = hydro07_leaf_area_index_low_max,
    by = 'hydroname') |>
  left_join(
    y = hydro07_evi,
    by = 'hydroname') |>
  left_join(
    y = hydro07_ndvi,
    by = 'hydroname') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols =X200901_leafareaindexhighvegetation:X20221201_NDVI,
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
  relocate(c('variable','valor'),.after = month) -> hydro07_vegetation

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/vegetation')){dir.create('output/vegetation')}
write_csv(hydro07_vegetation,'output/vegetation/hydro07_vegetation.csv')