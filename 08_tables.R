library(tidyverse)
library(gt)
library(ggcorrplot)
library(patchwork)

# 1. Output of error metrics table  ---------------------------------------
if(!dir.exists('output/tables')){dir.create('output/tables')}
fun_round <- function(x){return(round(x,2))}
met_viv <-  list.files("output/metrics/",pattern = "*.csv",full.names = T) |> 
  lapply(read_csv) |> 
  map_df(.f = as.data.frame) |>
  arrange(-desc(rmse)) |> 
  filter(disease %in% "P.Vivax") |> 
  relocate(region,.before = modelo) |> 
  mutate_if(is.numeric,fun_round)

met_viv |> 
  gt() |> 
  gtsave("output/tables/ml_error_metrics_vivax.docx")
# 2. Covariables Statistic ------------------------------------------------
climatic <- c(
  "tmmn",
  "tmmx",
  "vs",
  "temperature2m",
  "soiltemperature",
  "lstday",
  "totalprecipitation",
  "pr",
  "pet")

hydrologic <- c(
  "runoff",
  "ro",
  "soil",
  "potentialevaporation.max",
  "potentialevaporation.min")

vegetation <- c(
  "evi",
  "ndvi",
  "leafareaindexhighvegetation",
  "leafareaindexhighvegetationmin", 
  "leafareaindexhighvegetationmax",
  "leafareaindexlowvegetation",
  "leafareaindexlowvegetationmin",
  "leafareaindexlowvegetationmax")

# Variables names 
variable_mapping <- c(
  "tmmn" = "Minimum temperature",
  "tmmx" = "Maximum temperature",
  "vs" = "Wind speed at 10 m",
  "temperature2m" = "Temperature of air at 2m",
  "soiltemperature" = "Temperature of the soil in layer (0-7cm)",
  "potentialevaporation.max" = "Potential evaporation max",
  "potentialevaporation.min" = "Potential evaporation min",
  "totalprecipitation" = "Some water from rainfall",
  "pr" = "Precipitation accumulation",
  "lstday" = "Land Surface Temperature",
  "leafareaindexhighvegetation" = "Leaf area index high vegetation",
  "leafareaindexhighvegetationmin" = "Leaf area index high vegetation min",
  "leafareaindexhighvegetationmax" = "Leaf area index high vegetation max",
  "leafareaindexlowvegetation" = "Leaf area index low vegetation",
  "leafareaindexlowvegetationmin" = "Leaf area index low vegetation min",
  "leafareaindexlowvegetationmax" = "Leaf area index low vegetation max",
  "evi" = "Enhanced Vegetation Index (EVI)",
  "ndvi" = "Normalized Difference Vegetation Index (NDVI)",
  "aet" = "Actual evapotranspiration",
  "def" = "Climate water deficit",
  "pdsi" = "Palmer Drought Severity Index",
  "pet" = "Reference evapotranspiration",
  "ro" = "Runoff",
  "runoff" = "Runoff accumulation",
  "soil" = "Soil Moisture"
)

dis <- read_csv("output/ddbb/ddbb_dist.csv")|> 
  select(tmmn:ndvi) |> 
  summarise(
    across(
      tmmn:ndvi,
      list(media = mean, sd = sd),
      na.rm = TRUE)
    )

# Convert wide format to long format
dis_long <- dis |> 
  pivot_longer(cols = everything(), 
               names_to = c("Variable", ".value"), 
               names_sep = "_")

# Assign group labels
dis_long <- dis_long |> 
  mutate(Group = case_when(
    Variable %in% climatic ~ "Climatic",
    Variable %in% hydrologic ~ "Hydrologic",
    Variable %in% vegetation ~ "Vegetation",
    TRUE ~ "Unknown"
  ))

# Arrange for better visualization
dis_long <- dis_long |> arrange(Group, Variable) |> 
  mutate_if(is.numeric,fun_round) |> 
  relocate(Group,.before = Variable) |> 
  mutate(region = "Districts")

hy6 <- read_csv("output/ddbb/ddbb_hy6.csv") |>
  select(tmmn:ndvi) |> 
  summarise(
    across(
      tmmn:ndvi,
      list(media = mean, sd = sd),
      na.rm = TRUE)
  )

# Convert wide format to long format
hy6_long <- hy6 |> 
  pivot_longer(cols = everything(), 
               names_to = c("Variable", ".value"), 
               names_sep = "_")

# Assign group labels
hy6_long <- hy6_long |> 
  mutate(Group = case_when(
    Variable %in% climatic ~ "Climatic",
    Variable %in% hydrologic ~ "Hydrologic",
    Variable %in% vegetation ~ "Vegetation",
    TRUE ~ "Unknown"
  ))

# Arrange for better visualization
hy6_long <- hy6_long |> arrange(Group, Variable) |> 
  mutate_if(is.numeric,fun_round) |> 
  relocate(Group,.before = Variable) |> 
  mutate(region = "Hydrobasin level 6")


hy7 <- read_csv("output/ddbb/ddbb_hy7.csv") |> 
  select(tmmn:ndvi) |> 
  summarise(
    across(
      tmmn:ndvi,
      list(media = mean, sd = sd),
      na.rm = TRUE)
  )

# Convert wide format to long format
hy7_long <- hy7 |> 
  pivot_longer(cols = everything(), 
               names_to = c("Variable", ".value"), 
               names_sep = "_")

# Assign group labels
hy7_long <- hy7_long |> 
  mutate(Group = case_when(
    Variable %in% climatic ~ "Climatic",
    Variable %in% hydrologic ~ "Hydrologic",
    Variable %in% vegetation ~ "Vegetation",
    TRUE ~ "Unknown"
  ))

# Arrange for better visualization
hy7_long <- hy7_long |> arrange(Group, Variable) |> 
  mutate_if(is.numeric,fun_round) |> 
  relocate(Group,.before = Variable) |> 
  mutate(region = "Hydrobasin level 7")

hya <- read_csv("output/ddbb/ddbb_hya.csv") |> 
  select(tmmn:ndvi) |> 
  summarise(
    across(
      tmmn:ndvi,
      list(media = mean, sd = sd),
      na.rm = TRUE)
  )

# Convert wide format to long format
hya_long <- hya |> 
  pivot_longer(cols = everything(), 
               names_to = c("Variable", ".value"), 
               names_sep = "_")

# Assign group labels
hya_long <- hya_long |> 
  mutate(Group = case_when(
    Variable %in% climatic ~ "Climatic",
    Variable %in% hydrologic ~ "Hydrologic",
    Variable %in% vegetation ~ "Vegetation",
    TRUE ~ "Unknown"
  ))

# Arrange for better visualization
hya_long <- hya_long |> arrange(Group, Variable) |> 
  mutate_if(is.numeric,fun_round) |> 
  relocate(Group,.before = Variable) |> 
  mutate(region = "Hydrobasin ANA")

## Final table 
table_stats <-  bind_rows(dis_long,hy6_long,hy7_long,hya_long) |> 
  relocate(region,.after = Group) |> 
  mutate(Variable = recode(Variable, !!!variable_mapping))

table_stats |> 
  gt() |> 
  gtsave("output/tables/stats_variables.docx")


# dis_top <- table_stats |> 
#   group_by(Group,region) |> 
#   slice_max(order_by = media, n = 3) |> 
#   arrange(Group, desc(media))


# Correlation matrix ------------------------------------------------------
abbreviations <- c(
  "tmmn" = "tmin",
  "tmmx" = "tmax",
  "vs" = "wind_10m",
  "temperature2m" = "t2m",
  "soiltemperature" = "tsoil",
  "potentialevaporation.max" = "pet_max",
  "potentialevaporation.min" = "pet_min",
  "totalprecipitation" = "rainfall",
  "pr" = "prec",
  "lstday" = "lst",
  "leafareaindexhighvegetation" = "lai_high",
  "leafareaindexhighvegetationmin" = "lai_high_min",
  "leafareaindexhighvegetationmax" = "lai_high_max",
  "leafareaindexlowvegetation" = "lai_low",
  "leafareaindexlowvegetationmin" = "lai_low_min",
  "leafareaindexlowvegetationmax" = "lai_low_max",
  "evi" = "evi",
  "ndvi" = "ndvi",
  "aet" = "aet",
  "def" = "cwd",
  "pdsi" = "pdsi",
  "pet" = "ref_et",
  "ro" = "runoff",
  "runoff" = "runoff_acc",
  "soil" = "soil_m"
)

dis <- read_csv("output/ddbb/ddbb_dist.csv") |> 
  pivot_longer(cols = tmmn:ndvi,names_to = "variable",values_to = "valor") |> 
  mutate(variable = recode(variable, !!!abbreviations)) |> 
  pivot_wider(names_from = variable,values_from = valor) |> 
  select(tmin:ndvi) |> 
  scale() |> 
  cor() |> 
  ggcorrplot(outline.color = "white") +
  scale_fill_viridis_c(option = "rocket",direction = -1) + 
  labs(title = "Correlation Matrix (Districts)") + 
  theme(legend.position = "bottom")

hy6 <- read_csv("output/ddbb/ddbb_hy6.csv") |> 
  pivot_longer(cols = tmmn:ndvi,names_to = "variable",values_to = "valor") |> 
  mutate(variable = recode(variable, !!!abbreviations)) |> 
  pivot_wider(names_from = variable,values_from = valor) |> 
  select(tmin:ndvi) |> 
  scale() |> 
  cor() |> 
  ggcorrplot(outline.color = "white") +
  scale_fill_viridis_c(option = "rocket",direction = -1) + 
  labs(title = "Correlation Matrix (Hydrobasin level 6)") + 
  theme(legend.position = "bottom")

hy7 <- read_csv("output/ddbb/ddbb_hy7.csv") |> 
  pivot_longer(cols = tmmn:ndvi,names_to = "variable",values_to = "valor") |> 
  mutate(variable = recode(variable, !!!abbreviations)) |> 
  pivot_wider(names_from = variable,values_from = valor) |> 
  select(tmin:ndvi) |> 
  scale() |> 
  cor() |> 
  ggcorrplot(outline.color = "white") +
  scale_fill_viridis_c(option = "rocket",direction = -1) +
  labs(title = "Correlation Matrix (Hydrobasin level 7)") + 
  theme(legend.position = "bottom")

hya <- read_csv("output/ddbb/ddbb_hya.csv") |> 
  pivot_longer(cols = tmmn:ndvi,names_to = "variable",values_to = "valor") |> 
  mutate(variable = recode(variable, !!!abbreviations)) |> 
  pivot_wider(names_from = variable,values_from = valor) |> 
  select(tmin:ndvi) |> 
  scale() |> 
  cor() |> 
  ggcorrplot(outline.color = "white") +
  scale_fill_viridis_c(option = "rocket",direction = -1) + 
  labs(title = "Correlation Matrix (Hydrobasin ANA)") + 
  theme(legend.position = "bottom")


(dis|hy6)/(hy7|hya) 
if(!dir.exists('output/graphics')){dir.create('output/graphics')}
ggsave("output/graphics/correlation_matrix.png",plot = last_plot(),width = 15,height = 15,dpi = 300,bg = "white")
