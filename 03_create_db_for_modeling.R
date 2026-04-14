library(tidyverse)
library(sf)
library(qgisprocess)
library(lubridate)
sf_use_s2(use_s2 = F)

# 1. Reading spatial data -------------------------------------------------
dist  <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'districts')

hy06  <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'hydrobasin_06')

hy07  <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'hydrobasin_07')

hyana <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'hydroana') |> 
  mutate(hydroname = str_remove(hydroname,pattern = "name"))

malaria <- read_rds(file = 'sources/rawdata/malaria_basins_edited_v2.rds') 

# 2. Covariables ----------------------------------------------------------
cov.dist.clim <- read_csv("output/climate/districts_climate.csv")
cov.hy6.clim <- read_csv("output/climate/hydro06_climate.csv")
cov.hy7.clim <- read_csv("output/climate/hydro07_climate.csv")
cov.hya.clim <- read_csv("output/climate/hydroana_climate.csv")

cov.dist.wat <- read_csv("output/water/districts_water.csv")
cov.hy6.wat <- read_csv("output/water/hydro06_water.csv")
cov.hy7.wat <- read_csv("output/water/hydro07_water.csv")
cov.hya.wat <- read_csv("output/water/hydroana_water.csv")

cov.dist.veg <- read_csv("output/vegetation/districts_vegetation.csv")
cov.hy6.veg <- read_csv("output/vegetation/hydro06_vegetation.csv")
cov.hy7.veg <- read_csv("output/vegetation/hydro07_vegetation.csv")
cov.hya.veg <- read_csv("output/vegetation/hydroana_vegetation.csv")

cov.dist <- bind_rows(cov.dist.clim, cov.dist.wat, cov.dist.veg)
cov.hy6 <- bind_rows(cov.hy6.clim, cov.hy6.wat, cov.hy6.veg)
cov.hy7 <- bind_rows(cov.hy7.clim, cov.hy7.wat, cov.hy7.veg)
cov.hya <- bind_rows(cov.hya.clim, cov.hya.wat, cov.hya.veg)

cov.dist.var <- cov.dist |> 
  select(ubigeo,year, month,variable,valor) |> 
  filter(year %in% 2009:2018)

cov.hy6.var  <- cov.hy6 |> 
  select(hydroname,year, month,variable,valor)|> 
  filter(year %in% 2009:2018)

cov.hy7.var  <- cov.hy7 |> 
  select(hydroname,year, month,variable,valor)|> 
  filter(year %in% 2009:2018)

cov.hya.var  <- cov.hya |> 
  select(hydroname,year, month,variable,valor) |> 
  mutate(hydroname = str_remove(hydroname,pattern = "name")) |> 
  filter(year %in% 2009:2018)

cov.dist.var <- cov.dist.var|> 
  pivot_wider(names_from = variable, values_from = valor)

cov.hy6.var <- cov.hy6.var  |> 
  pivot_wider(names_from = variable, values_from = valor)

cov.hy7.var <- cov.hy7.var  |>
  pivot_wider(names_from = variable, values_from = valor)

cov.hya.var <- cov.hya.var  |>
  pivot_wider(names_from = variable, values_from = valor)
  
# 3. Geoprocessing --------------------------------------------------------
point_to_polygon <- qgis_function(algorithm = "native:joinattributesbylocation")

## Districts
m.dist <- point_to_polygon(
  INPUT = dist,
  PREDICATE = 0,
  JOIN = malaria,
  OUTPUT = qgis_tmp_vector(ext = ".csv")) |> 
  qgis_extract_output(name = "OUTPUT") |> 
  read_csv()  |> 
  select(ubigeo,year,month,viv,fal,nrohab) |> 
  group_by(ubigeo,year,month) |> 
  summarise(
    viv  = sum(viv),
    fal  = sum(fal),
    pob  = round(mean(nrohab))) 

ddbb.dist <- m.dist |> 
  inner_join(cov.dist.var,by = c("ubigeo","year","month"))

## Hydrobasin level 06
m.hy6 <- point_to_polygon(
  INPUT = hy06,
  PREDICATE = 0,
  JOIN = malaria,
  OUTPUT = qgis_tmp_vector(ext = ".csv")) |> 
  qgis_extract_output(name = "OUTPUT") |> 
  read_csv()  |> 
  select(hydroname,year,month,viv,fal,nrohab) |> 
  group_by(hydroname,year,month) |> 
  summarise(
    viv  = sum(viv),
    fal  = sum(fal),
    pob  = round(mean(nrohab))) 

ddbb.hy6 <- m.hy6 |> 
  inner_join(cov.hy6.var,by = c("hydroname","year","month"))


## Hydrobasin level 07
m.hy7 <- point_to_polygon(
  INPUT = hy07,
  PREDICATE = 0,
  JOIN = malaria,
  OUTPUT = qgis_tmp_vector(ext = ".csv")) |> 
  qgis_extract_output(name = "OUTPUT") |> 
  read_csv()  |> 
  select(hydroname,year,month,viv,fal,nrohab) |> 
  group_by(hydroname,year,month) |> 
  summarise(
    viv  = sum(viv),
    fal  = sum(fal),
    pob  = round(mean(nrohab))) 

ddbb.hy7 <- m.hy7 |> 
  inner_join(cov.hy7.var,by = c("hydroname","year","month"))

## Hydrobasin ANA
m.hya <- point_to_polygon(
  INPUT = hyana,
  PREDICATE = 0,
  JOIN = malaria,
  OUTPUT = qgis_tmp_vector(ext = ".csv")) |> 
  qgis_extract_output(name = "OUTPUT") |> 
  read_csv()  |> 
  select(hydroname,year,month,viv,fal,nrohab) |> 
  group_by(hydroname,year,month) |> 
  summarise(
    viv  = sum(viv),
    fal  = sum(fal),
    pob  = round(mean(nrohab)))

ddbb.hya <- m.hya |> 
  inner_join(cov.hya.var,by = c("hydroname","year","month"))

# Exporting target data ---------------------------------------------------
if(!dir.exists('output/ddbb')){dir.create('output/ddbb')}
write_csv(ddbb.dist,'output/ddbb/ddbb_dist.csv')
write_csv(ddbb.hy6, 'output/ddbb/ddbb_hy6.csv')
write_csv(ddbb.hy7, 'output/ddbb/ddbb_hy7.csv')
write_csv(ddbb.hya, 'output/ddbb/ddbb_hya.csv')