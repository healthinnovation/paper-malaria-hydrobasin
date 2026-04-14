library(tidyverse)
library(factoextra)
library(viridis)
library(patchwork)

# 0. Reading data ---------------------------------------------------------
dist <- read_csv("output/ddbb/ddbb_dist.csv") %>% mutate(pdsi = pdsi*-1, vs = vs*-1, def = def*-1)
hy6  <- read_csv("output/ddbb/ddbb_hy6.csv") %>% mutate(pdsi  = pdsi*-1, vs = vs*-1, def = def*-1)
hy7  <- read_csv("output/ddbb/ddbb_hy7.csv") %>% mutate(pdsi  = pdsi*-1, vs = vs*-1, def = def*-1)
hya  <- read_csv("output/ddbb/ddbb_hya.csv") %>% mutate(pdsi  = pdsi*-1, vs = vs*-1, def = def*-1)

# 1. Variables Groups -----------------------------------------------------
climaticos <- c(
  "tmmn",
  "tmmx",
  "vs",
  "temperature2m",
  "soiltemperature",
  "potentialevaporation.max",
  "potentialevaporation.min",
  "runoff",
  "totalprecipitation",
  "lstday"
  )

vegetacion <- c(
  "leafareaindexhighvegetation",
  "leafareaindexhighvegetationmin", 
  "leafareaindexhighvegetationmax", 
  "leafareaindexlowvegetation",
  "leafareaindexlowvegetationmin", 
  "leafareaindexlowvegetationmax",
  "evi",
  "ndvi"
)

hidrologicos <- c(
  "aet",
  "def",
  "pdsi",
  "pet",
  "pr",
  "ro",
  "soil"
  )

# 2. PCA analysis ---------------------------------------------------------
## 2.1 Districts
### Climate
dist.pca.clima <- dist |> 
  select(all_of(climaticos)) |> 
  scale() |> 
  prcomp(scale. = F)

dist.pca.clima |> summary()

### Hydrology
dist.pca.hidro <- dist |> 
  select(all_of(hidrologicos)) |> 
  scale() |> 
  prcomp(scale. = F)

dist.pca.hidro |> summary()

### Vegetation
dist.pca.veg <- dist |> select(all_of(vegetacion)) |> 
  scale() |> 
  prcomp(scale. = F)

dist.pca.veg |> summary()

### New database
db.dist.pca <-  dist |> 
  mutate(
    pc1_clim = as.vector(dist.pca.clima$x[,1]),
    pc2_clim = as.vector(dist.pca.clima$x[,2]),
    pc3_clim = as.vector(dist.pca.clima$x[,3]),
    pc4_clim = as.vector(dist.pca.clima$x[,4]),
    pc5_clim = as.vector(dist.pca.clima$x[,5]),
    pc1_hidr = as.vector(dist.pca.hidro$x[,1]),
    pc2_hidr = as.vector(dist.pca.hidro$x[,2]),
    pc3_hidr = as.vector(dist.pca.hidro$x[,3]),
    pc4_hidr = as.vector(dist.pca.hidro$x[,4]),
    pc1_veg = as.vector(dist.pca.veg$x[,1]),
    pc2_veg = as.vector(dist.pca.veg$x[,2]),
    pc3_veg = as.vector(dist.pca.veg$x[,3])) |> 
  select(1:6,pc1_clim:pc3_veg)

### Plot PCA
nombres_climate <- c(
  "Min Temperature",
  "Max Temperature",
  "Wind Speed (10m)",
  "Air Temperature (2m)",
  "Soil Temperature (0-7cm)",
  "Max Potential Evaporation",
  "Min Potential Evaporation",
  "Surface and Subsurface Runoff",
  "Precipitation (Rain & Snow)",
  "Land Surface Temperature")

rownames(dist.pca.clima$rotation) <- nombres_climate
p1 <- dist.pca.clima |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Climate components")

nombres_hidro <- c(
  "Actual Evapotranspiration",
  "Climate Water Deficit",
  "Palmer Drought Severity Index",
  "Reference Evapotranspiration",
  "Precipitation Accumulation",
  "Runoff",
  "Soil Moisture")

rownames(dist.pca.hidro$rotation) <- nombres_hidro
p2 <- dist.pca.hidro |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Hidrology components")

nombres_veg <- c(
  "Leaf Area Index (High Vegetation)",
  "Min Leaf Area Index (High Vegetation)",
  "Max Leaf Area Index (High Vegetation)",
  "Leaf Area Index (Low Vegetation)",
  "Min Leaf Area Index (Low Vegetation)",
  "Max Leaf Area Index (Low Vegetation)",
  "Enhanced Vegetation Index (EVI)",
  "Normalized Difference Vegetation Index (NDVI)")

rownames(dist.pca.veg$rotation) <- nombres_veg
p3 <- dist.pca.veg |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Vegetatiom components")

dis_pca <- p1|p2|p3

### Export PCA database
if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/pca')){dir.create('output/pca')}
write_csv(db.dist.pca,'output/pca/db.dist.pca.csv')

## 2.2 Hydrobasin level 6
### Climate
hy6.pca.clima <- hy6 |> 
  select(all_of(climaticos)) |> 
  scale() |> 
  prcomp(scale. = F)

hy6.pca.clima |> summary()

### Hydrology
hy6.pca.hidro <- hy6 |> 
  select(all_of(hidrologicos)) |> 
  scale() |> 
  prcomp(scale. = F)

hy6.pca.hidro |> summary()

### Vegetation
hy6.pca.veg <- hy6 |> select(all_of(vegetacion)) |> 
  scale() |> 
  prcomp(scale. = F)

hy6.pca.veg |> summary()

### New database
db.hy6.pca <-  hy6 |> 
  mutate(
    pc1_clim = as.vector(hy6.pca.clima$x[,1]),
    pc2_clim = as.vector(hy6.pca.clima$x[,2]),
    pc3_clim = as.vector(hy6.pca.clima$x[,3]),
    pc4_clim = as.vector(hy6.pca.clima$x[,4]),
    pc1_hidr = as.vector(hy6.pca.hidro$x[,1]),
    pc2_hidr = as.vector(hy6.pca.hidro$x[,2]),
    pc3_hidr = as.vector(hy6.pca.hidro$x[,3]),
    pc4_hidr = as.vector(hy6.pca.hidro$x[,4]),
    pc1_veg = as.vector(hy6.pca.veg$x[,1]),
    pc2_veg = as.vector(hy6.pca.veg$x[,2]),
    pc3_veg = as.vector(hy6.pca.veg$x[,3])) |> 
  select(1:6,pc1_clim:pc3_veg)

### PCA plots
rownames(hy6.pca.clima$rotation) <-  nombres_climate
m1 <- hy6.pca.clima |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    lenged.title = "Contribution") + 
  labs(title = "Climate components")

rownames(hy6.pca.hidro$rotation) <-  nombres_hidro
m2 <- hy6.pca.hidro |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Hidrology components")

rownames(hy6.pca.veg$rotation) <-  nombres_veg
m3 <- hy6.pca.veg |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Vegetatiom components")

hy6_pca <- m1|m2|m3

### Export PCA database
if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/pca')){dir.create('output/pca')}
write_csv(db.hy6.pca,'output/pca/db.hy6.pca.csv')

## 2.3 Hydrobasin level 7
### Climate
hy7.pca.clima <- hy7 |> 
  select(all_of(climaticos)) |> 
  scale() |> 
  prcomp(scale. = F)

hy7.pca.clima |> summary()

### Hydrology
hy7.pca.hidro <- hy7 |> 
  select(all_of(hidrologicos)) |> 
  scale() |> 
  prcomp(scale. = F)

hy7.pca.hidro |> summary()

### Vegetation
hy7.pca.veg <- hy7 |> select(all_of(vegetacion)) |> 
  scale() |> 
  prcomp(scale. = F)

hy7.pca.veg |> summary()

### New database
db.hy7.pca <-  hy7 |> 
  mutate(
    pc1_clim = as.vector(hy7.pca.clima$x[,1]),
    pc2_clim = as.vector(hy7.pca.clima$x[,2]),
    pc3_clim = as.vector(hy7.pca.clima$x[,3]),
    pc4_clim = as.vector(hy7.pca.clima$x[,4]),
    pc5_clim = as.vector(hy7.pca.clima$x[,5]),
    pc1_hidr = as.vector(hy7.pca.hidro$x[,1]),
    pc2_hidr = as.vector(hy7.pca.hidro$x[,2]),
    pc3_hidr = as.vector(hy7.pca.hidro$x[,3]),
    pc4_hidr = as.vector(hy7.pca.hidro$x[,4]),
    pc1_veg = as.vector(hy7.pca.veg$x[,1]),
    pc2_veg = as.vector(hy7.pca.veg$x[,2])) |> 
  select(1:6,pc1_clim:pc2_veg)

### PCA plots
rownames(hy7.pca.clima$rotation) <-  nombres_climate
n1 <- hy7.pca.clima |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Climate components")

rownames(hy7.pca.hidro$rotation) <-  nombres_hidro
n2 <- hy7.pca.hidro |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Hidrology components")

rownames(hy7.pca.veg$rotation) <-  nombres_veg
n3 <- hy7.pca.veg |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Vegetatiom components")

hy7_pca <- n1|n2|n3

### Export PCA database
if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/pca')){dir.create('output/pca')}
write_csv(db.hy7.pca,'output/pca/db.hy7.pca.csv')

## 2.4 Hydrobasin ANA
### Climate
hya.pca.clima <- hya |> 
  select(all_of(climaticos)) |> 
  scale() |> 
  prcomp(scale. = F)

hya.pca.clima |> summary()

### Hydrology
hya.pca.hidro <- hya |> 
  select(all_of(hidrologicos)) |> 
  scale() |> 
  prcomp(scale. = F)

hya.pca.hidro |> summary()

### Vegetation
hya.pca.veg <- hya |> select(all_of(vegetacion)) |> 
  scale() |> 
  prcomp(scale. = F)

hya.pca.veg |> summary()

### New database
db.hya.pca <-  hya |> 
  mutate(
    pc1_clim = as.vector(hya.pca.clima$x[,1]),
    pc2_clim = as.vector(hya.pca.clima$x[,2]),
    pc3_clim = as.vector(hya.pca.clima$x[,3]),
    pc4_clim = as.vector(hya.pca.clima$x[,4]),
    pc5_clim = as.vector(hya.pca.clima$x[,5]),
    pc1_hidr = as.vector(hya.pca.hidro$x[,1]),
    pc2_hidr = as.vector(hya.pca.hidro$x[,2]),
    pc3_hidr = as.vector(hya.pca.hidro$x[,3]),
    pc4_hidr = as.vector(hya.pca.hidro$x[,4]),
    pc1_veg = as.vector(hya.pca.veg$x[,1]),
    pc2_veg = as.vector(hya.pca.veg$x[,2])) |> 
  select(1:6,pc1_clim:pc2_veg)

### PCA plots
rownames(hya.pca.clima$rotation) <-  nombres_climate
q1 <- hya.pca.clima |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Climate components")

rownames(hya.pca.hidro$rotation) <-  nombres_hidro
q2 <- hya.pca.hidro |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Hidrology components")

rownames(hya.pca.veg$rotation) <-  nombres_veg
q3 <- hya.pca.veg |> 
  fviz_pca_var(
    col.var = "contrib",
    gradient.cols = rocket(n = 8,direction = -1),
    repel = TRUE,
    legend.title = "Contribution") + 
  labs(title = "Vegetatiom components")

ana_pca <- q1|q2|q3

### Export PCA database
if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/pca')){dir.create('output/pca')}
write_csv(db.hya.pca,'output/pca/db.hya.pca.csv')

# 3. Final Plot -----------------------------------------------------------
dis_pca/hy6_pca/hy7_pca/ana_pca

ggsave(filename = "output/graphics/pca.svg",width = 30,height = 25)
