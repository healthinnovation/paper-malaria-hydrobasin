library(sf)
library(tidyverse)
library(viridis)
library(extrafont)
library(cptcity)
library(ggtext)

minmax <- function(x){
  lyr <- (x - (min(x)))/(max(x) - (min(x)))
  return(lyr)
}

# 1. Districts ------------------------------------------------------------
variables <- c(
  "ubigeo",
  "year",
  "month",
  "residual_fal_rf",
  "residual_viv_rf",
  "residual_fal_xgb",
  "residual_viv_xgb",
  "residual_fal_svm",
  "residual_viv_svm"
)


mis_etiquetas <- c(
  "residual_fal_rf"  = "*P.Falciparum* RF",
  "residual_fal_svm" = "*P.Falciparum* SVM",
  "residual_fal_xgb" = "*P.Falciparum* XGB",
  "residual_viv_rf"  = "*P.Vivax* RF",
  "residual_viv_svm" = "*P.Vivax* SVM",
  "residual_viv_xgb" = "*P.Vivax* XGB"
)


dist <- st_read(
  dsn = "sources/rawdata/geometry.gpkg",
  layer = "districts")

error_dist <- read_csv("output/metrics/error_metrics_total_districts.csv") |> 
  select(all_of(variables)) |> 
  mutate(ubigeo = as.character(ifelse(test = nchar(as.character(ubigeo))< 5, paste0('0',ubigeo),no = ubigeo))) |> 
  group_by(ubigeo) |> 
  summarise(across(
    .cols = starts_with("resi"),
    .fns = mean,   
    .names = "{.col}"
  ))

error_dist <- error_dist |> 
  mutate(across(
    .cols = starts_with("resi"),
    .fns = minmax,   
    .names = "{.col}"
  ))

dist_sf <- dist |> 
  inner_join(error_dist,by = "ubigeo") |> 
  pivot_longer(cols = residual_fal_rf:residual_viv_svm)

m1 <-  dist_sf |> 
  ggplot() + 
  geom_sf(aes(fill = value), color = "black", linewidth =0.1) + 
  facet_wrap(name~.,labeller = as_labeller(mis_etiquetas)) + 
  scale_fill_gradientn(name = "Residual values",colours = cpt(pal = "jjg_polarity_RdYlGn",rev = 1)) + 

  theme(
    legend.position = "none",
    legend.title = element_markdown(face = "italic", size = 12,hjust = 0.5),
    strip.text = element_markdown(face = "italic",family = "Arvo"),
    panel.border = element_rect(colour = "gray",fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(colour = "gray",linetype = "dashed"),
    axis.text = element_text(family = "Arvo",size = 6)
  )  

# Hydrobasin level 06 -----------------------------------------------------
hy06 <- st_read("sources/rawdata/geometry.gpkg",layer = "hydrobasin_06")

variables <- c(
  "hydroname",
  "year",
  "month",
  "residual_fal_rf",
  "residual_viv_rf",
  "residual_fal_xgb",
  "residual_viv_xgb",
  "residual_fal_svm",
  "residual_viv_svm"
)

hy6 <- st_read(
  dsn = "sources/rawdata/geometry.gpkg",
  layer = "hydrobasin_06")

error_hy6 <- read_csv("output/metrics/error_metrics_total_hydrobasinlevel06.csv") |> 
  select(all_of(variables)) |> 
  group_by(hydroname) |> 
  summarise(across(
    .cols = starts_with("resi"),
    .fns = mean,   
    .names = "{.col}"
  ))

error_hy6 <- error_hy6 |> 
  mutate(across(
    .cols = starts_with("resi"),
    .fns = minmax,   
    .names = "{.col}"
  ))

hy6_sf <- hy6 |> 
  inner_join(error_hy6,by = "hydroname") |> 
  pivot_longer(cols = residual_fal_rf:residual_viv_svm)


m2 <-  hy6_sf |> 
  ggplot() + 
  geom_sf(aes(fill = value), color = "black", linewidth =0.1) + 
  facet_wrap(name~.,labeller = as_labeller(mis_etiquetas)) + 
  scale_fill_gradientn(name = "Residual values",colours = cpt(pal = "jjg_polarity_RdYlGn",rev = 1)) + 

  theme(
    legend.position = "none",
    legend.title = element_markdown(face = "italic", size = 12,hjust = 0.5),
    strip.text = element_markdown(face = "italic",family = "Arvo"),
    panel.border = element_rect(colour = "gray",fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(colour = "gray",linetype = "dashed"),
    axis.text = element_text(family = "Arvo",size = 6)
  ) 

# Hydrobasin level 7 ------------------------------------------------------
bbox_ref <- st_bbox(hy06)  
hy7 <- st_read(
  dsn = "sources/rawdata/geometry.gpkg",
  layer = "hydrobasin_07")

error_hy7 <- read_csv("output/metrics/error_metrics_total_hydrobasinlevel7.csv") |> 
  select(all_of(variables)) |> 
  group_by(hydroname) |> 
  summarise(across(
    .cols = starts_with("resi"),
    .fns = mean,   
    .names = "{.col}"
  ))

error_hy7 <- error_hy7 |> 
  mutate(across(
    .cols = starts_with("resi"),
    .fns = minmax,   
    .names = "{.col}"
  ))

hy7_sf <- hy7 |> 
  inner_join(error_hy7,by = "hydroname") |> 
  pivot_longer(cols = residual_fal_rf:residual_viv_svm)

m3 <-  hy7_sf |> 
  ggplot() + 
  geom_sf(aes(fill = value), color = "black", linewidth =0.1) + 
  facet_wrap(name~.,labeller = as_labeller(mis_etiquetas)) + 
  scale_fill_gradientn(name = "Residual values",colours = cpt(pal = "jjg_polarity_RdYlGn",rev = 1)) + 

  theme(
    legend.position = "none",
    legend.title = element_markdown(face = "italic", size = 12,hjust = 0.5),
    strip.text = element_markdown(face = "italic",family = "Arvo"),
    panel.border = element_rect(colour = "gray",fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(colour = "gray",linetype = "dashed"),
    axis.text = element_text(family = "Arvo",size = 6),
    plot.margin = margin(0, 0, 0, 0),
  ) +
  coord_sf(
    xlim = c(bbox_ref["xmin"], bbox_ref["xmax"]),
    ylim = c(bbox_ref["ymin"], c(bbox_ref["ymax"])),
    expand = FALSE
  )


# Hydrobasin ANA ----------------------------------------------------------
hya <- st_read(
  dsn = "sources/rawdata/geometry.gpkg",
  layer = "hydroana") |> 
  mutate(hydroname = str_remove(string = hydroname,pattern = "name"))

error_hya <- read_csv("output/metrics/error_metrics_total_hydrobasinana.csv") |> 
  select(all_of(variables)) |> 
  group_by(hydroname) |> 
  summarise(across(
    .cols = starts_with("resi"),
    .fns = mean,   
    .names = "{.col}"
  ))

error_hya <- error_hya |> 
  mutate(across(
    .cols = starts_with("resi"),
    .fns = minmax,   
    .names = "{.col}"
  ))

hya_sf <- hya |> 
  inner_join(error_hya, by = "hydroname") |> 
  pivot_longer(cols = residual_fal_rf:residual_viv_svm)

m4 <-  hya_sf |> 
  ggplot() + 
  geom_sf(aes(fill = value), color = "black", linewidth =0.1) + 
  facet_wrap(name~.,labeller = as_labeller(mis_etiquetas)) + 
  scale_fill_gradientn(name = "Residual values",colours = cpt(pal = "jjg_polarity_RdYlGn",rev = 1)) + 
  theme(
    legend.position = "none",
    legend.title = element_markdown(face = "italic", size = 12,hjust = 0.5),
    strip.text = element_markdown(face = "italic",family = "Arvo"),
    panel.border = element_rect(colour = "gray",fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(colour = "gray",linetype = "dashed"),
    axis.text = element_text(family = "Arvo",size = 6)
  ) 

# Panel  ------------------------------------------------------------------
(m4|m1) / (m2| m3)
ggsave("output/graphics/error_maps_v3.svg", last_plot(),width = 12,height = 11)