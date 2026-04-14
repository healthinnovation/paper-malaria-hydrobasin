library(tidyverse)
library(sf)
library(patchwork)
library(viridis)
library(ggtext)
library(ggforce)
library(patchwork)
# HeatMap -----------------------------------------------------------------
data <- st_read("sources/rawdata/cdc_data_malaria_con_pob.gpkg") |> 
  st_drop_geometry() |> 
  group_by(anio,mes) |> 
  summarise(
    p.viv = sum(viv,na.rm = T),
    p.fal = sum(fal,na.rm = T)) |> 
  mutate(
    anio = factor(anio),
    mes = factor(mes,labels = month.abb)
  )

p1 <- ggplot()+ 
  geom_raster(data = data, aes(x = mes,y = anio, fill = p.fal)) +
  scale_fill_gradientn(
    name = "<i><b>P.Falciparum</b></i>",
    colours = rocket(n = 5,direction = -1)) + 
  labs( y = "", x = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  guides(fill = guide_colourbar(title.position = "top",barwidth = 15,barheight = .5)) + 
  theme(
    legend.title = element_markdown(hjust = 0.5,family = "Bitter",colour = "black"),
    legend.position = "top",
    axis.text = element_text(family = "Bitter",colour = "black"),
    legend.text = element_text(family = "Bitter",colour = "black")
  )

p2 <- ggplot()+ 
  geom_raster(data = data, aes(x = mes,y = anio, fill = p.viv))+
  scale_fill_gradientn(
    name = "<i><b>P.Vivax</b></i>",
    colours = mako(n = 5,direction = -1)) + 
  labs( y = "", x = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  guides(fill = guide_colourbar(title.position = "top",barwidth = 15,barheight = .5)) + 
  theme(
    legend.position = "top",
    axis.text = element_text(family = "Bitter",colour = "black"),
    legend.text = element_text(family = "Bitter",colour = "black"),
    legend.title = element_markdown(hjust = 0.5,family = "Bitter",colour = "black")
  )

(p <- p1 | p2)
# ggsave(
#   plot = last_plot(),
#   filename = 'output/graphics/heatmap.png',
#   dpi=900,
#   width = 25,
#   height = 15,
#   bg = 'white',
#   units='cm')

# Error circles -----------------------------------------------------------
## RMSE
data <- read_csv("output/metrics/ml_metrics_hydrobasin.csv") %>% 
  janitor::clean_names()

rmse_limits <- data %>% 
  filter(specie == "P. Falciparum") %>%
  summarise(
    rmse_min = min (rmse),
    rmse_max = max(rmse)
    )

rmse_breaks <-  seq(from = rmse_limits$rmse_min, to = rmse_limits$rmse_max, length.out = 5)

p1 <- data %>% 
  filter(specie == "P. Falciparum") %>% 
  group_by(geographic_boundaries) %>% 
  mutate(mean_rmse = mean(rmse)) %>%  
  ungroup() %>% 
  mutate(geographic_boundaries = factor(geographic_boundaries, 
                                        levels = unique(geographic_boundaries[order(mean_rmse)]))) %>%  
  ggplot(aes(x = machine_learning_model, y = geographic_boundaries, fill = rmse)) + 
  geom_tile(size = 2) +  
  coord_radial(inner.radius = 0.3, expand = FALSE) +  
  scale_fill_viridis("<i><b>P.Falciparum</b></i><br><br>RMSE(Error)",,option = "rocket", direction = -1, 
                     limits = c(rmse_limits$rmse_min,rmse_limits$rmse_max), breaks = rmse_breaks) +  # Fijar límites y cortes
  labs(y = "", x = "", fill = "RMSE (Error)") +  
  theme_minimal(base_family = "Bitter") +  
  theme(
    legend.position = "top",
    legend.title = element_markdown(hjust = 0.5,family = "Bitter",colour = "black"),
    axis.text.y = element_blank()
    ) + 
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = 15,
      barheight = .5
    )
  )


rmse_limits <- data %>% 
  filter(specie == "P. Vivax") %>%
  summarise(
    rmse_min = min (rmse),
    rmse_max = max(rmse)
  )

rmse_breaks <-  seq(from = rmse_limits$rmse_min, to = rmse_limits$rmse_max, length.out = 5)

p2 <- data %>% 
  filter(specie == "P. Vivax") %>% 
  group_by(geographic_boundaries) %>% 
  mutate(mean_rmse = mean(rmse)) %>%  
  ungroup() %>% 
  mutate(geographic_boundaries = factor(geographic_boundaries, 
                                        levels = unique(geographic_boundaries[order(mean_rmse)]))) %>%  
  ggplot(aes(x = machine_learning_model, y = geographic_boundaries, fill = rmse)) + 
  geom_tile(size = 2) +  
  coord_radial(inner.radius = 0.3, expand = FALSE) +  
  scale_fill_viridis("<i><b>P.Vivax</b></i><br><br>RMSE(Error)",option = "mako", direction = -1, 
                     limits = c(rmse_limits$rmse_min,rmse_limits$rmse_max), breaks = rmse_breaks) +  # Fijar límites y cortes
  labs(y = "", x = "", fill = "RMSE (Error)") +  
  theme_minimal() +  
  theme(
    legend.position = "top",
    legend.title = element_markdown(hjust = 0.5,family = "Bitter",colour = "black"),
    axis.text.y = element_blank()
    ) + 
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = 15,
      barheight = .5
    )
  )

# Panel Plot
p1 + p2
# ggsave(
#   filename = "output/graphics/rmse.svg",
#   plot = last_plot(),
#   width = 15,
#   height = 10
#   )

## RSE
data <- read_csv("output/metrics/ml_metrics_hydrobasin.csv") %>% 
  janitor::clean_names()

rse_limits <- data %>% 
  filter(specie == "P. Falciparum") %>%
  summarise(
    rse_min = min (rse),
    rse_max = max(rse)
  )

rse_breaks <-  seq(from = rse_limits$rse_min, to = rse_limits$rse_max, length.out = 5)

p1 <- data %>% 
  filter(specie == "P. Falciparum") %>% 
  group_by(geographic_boundaries) %>% 
  mutate(mean_rse = mean(rse)) %>%  
  ungroup() %>% 
  mutate(geographic_boundaries = factor(geographic_boundaries, 
                                        levels = unique(geographic_boundaries[order(mean_rse)]))) %>%  
  ggplot(aes(x = machine_learning_model, y = geographic_boundaries, fill = rse)) + 
  geom_tile(size = 2) +  
  coord_radial(inner.radius = 0.3, expand = FALSE) +  
  scale_fill_viridis("<i><b>P.Falciparum</b></i><br><br>rse(Error)",,option = "rocket", direction = -1, 
                     limits = c(rse_limits$rse_min,rse_limits$rse_max), breaks = rse_breaks) +  # Fijar límites y cortes
  labs(y = "", x = "", fill = "rse (Error)") +  
  theme_minimal(base_family = "Bitter") +  
  theme(
    legend.position = "top",
    legend.title = element_markdown(hjust = 0.5,family = "Bitter",colour = "black"),
    axis.text.y = element_blank()
    ) + 
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = 15,
      barheight = .5
    )
  )


rse_limits <- data %>% 
  filter(specie == "P. Vivax") %>%
  summarise(
    rse_min = min (rse),
    rse_max = max(rse)
  )

rse_breaks <-  seq(from = rse_limits$rse_min, to = rse_limits$rse_max, length.out = 5)

p2 <- data %>% 
  filter(specie == "P. Vivax") %>% 
  group_by(geographic_boundaries) %>% 
  mutate(mean_rse = mean(rse)) %>%  
  ungroup() %>% 
  mutate(geographic_boundaries = factor(geographic_boundaries, 
                                        levels = unique(geographic_boundaries[order(mean_rse)]))) %>%  
  ggplot(aes(x = machine_learning_model, y = geographic_boundaries, fill = rse)) + 
  geom_tile(size = 2) +  
  coord_radial(inner.radius = 0.3, expand = FALSE) +  
  scale_fill_viridis("<i><b>P.Vivax</b></i><br><br>rse(Error)",option = "mako", direction = -1, 
                     limits = c(rse_limits$rse_min,rse_limits$rse_max), breaks = rse_breaks) +  # Fijar límites y cortes
  labs(y = "", x = "", fill = "rse (Error)") +  
  theme_minimal() +  
  theme(
    legend.position = "top",
    legend.title = element_markdown(hjust = 0.5,family = "Bitter",colour = "black"),
    axis.text.y = element_blank()
  ) + 
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = 15,
      barheight = .5
    )
  )

# Panel Plot
p1 + p2
# ggsave(
#   filename = "output/graphics/rse.svg",
#   plot = last_plot(),
#   width = 15,
#   height = 10
#   )