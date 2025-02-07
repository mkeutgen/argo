# Load necessary libraries
library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
library(segmented)
library(dplyr)
library(conflicted)
library(pracma)
library(fs)

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(viridis)
library(sp)
library(spdep)
library(mgcv)
library(patchwork)

library(fitdistrplus)
library(poweRlaw)
library(ggplot2)
library(tidync)

# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load Qiu's data
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
setwd("/data/GLOBARGO/src/")
df_argo_clean <- read_csv("data/df_argo_loc.csv")
df_argo_clean$TIME %>% min(na.rm = T) # "2010-05-30"
df_argo_clean$TIME %>% max(na.rm = T) # "2024-05-02"

df_complete_clean <- read_csv("data/df_eddy_subduction_anom.csv")
df_carbon_clean <- read_csv("data/df_carbon_subduction_anom.csv")
df_eke <- read_table("data/Qiu2018_data_fig5.asc")

df_eke$LON_num %>% unique() %>% sort() # 3° res



# Produce maps of EKE : 
df_eke$LON <- df_eke$`lon(E)` 
df_eke$LAT <- df_eke$`lat(N)`
df_eke$EKE_BL <- df_eke$`eke(bl)(m/s)^2` %>%  as.numeric()
df_eke$EKE_UNBL <- df_eke$`eke(unbl)(m/s)^2` %>%  as.numeric()

df_eke <- df_eke %>% select(LON,LAT,EKE_BL,EKE_UNBL)



df_eke <- df_eke %>%
  mutate(
    LON = ifelse(as.numeric(LON) > 180, as.numeric(LON) - 360, as.numeric(LON))
  )





breaks_eke_bl <- c(1e-4,1e-3 ,2.5e-3,5e-3 ,7.5e-3,
            1e-2,2.5e-2,5e-2,7.5e-2,1e-1,1)

breaks_eke_unbl <- c(1e-4,1e-3 ,2e-3,3e-3 ,4e-3,5e-3,
            6e-3,7e-3,8e-3,9e-3,1e-2,2.5e-2,1e-1)



df_eke <- df_eke %>%
  mutate(
    EKE_BL_bin = cut(EKE_BL,
                     breaks = breaks_eke_bl,
                     include.lowest = TRUE),
    EKE_UNBL_bin = cut(EKE_UNBL,
                      breaks = breaks_eke_unbl,
                      include.lowest = TRUE),
    EKE_UNBL_bin_bl = cut(EKE_UNBL,
                          breaks = breaks_eke_unbl,
                          include.lowest = TRUE)
  )





# 1) Balanced EKE map

eke_balanced_map <- ggplot() +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = 53, ymax = 90),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -90, ymax = -67),
    fill = "lightblue",
    inherit.aes = FALSE
  )+
  geom_tile(
    data = df_eke,
    aes(x = as.numeric(LON), 
        y = as.numeric(LAT), 
        fill = EKE_BL_bin)
  ) +
  scale_fill_viridis_d() + 
  geom_sf(data = world, fill = "grey", color = "grey", inherit.aes = FALSE) +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  labs(
    title = "Balanced Eddy Kinetic Energy",
    x = "Longitude",
    y = "Latitude",
    fill = "EKE (BL)"
  ) +
  theme_minimal() 



# 2) Unbalanced EKE map
eke_unbalanced_map <- ggplot() +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = 53, ymax = 90),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -90, ymax = -67),
    fill = "lightblue",
    inherit.aes = FALSE
  )+
  geom_tile(
    data = df_eke,
    aes(x = as.numeric(LON), 
        y = as.numeric(LAT), 
        fill = EKE_UNBL_bin)
  ) +
  geom_contour(
    data = df_eke,
    aes(x = as.numeric(LON), 
        y = as.numeric(LAT), 
        z = EKE_UNBL),
    color = "white", alpha = 0.3
  ) +
  scale_fill_viridis_d(
    option = "viridis") +
  geom_sf(data = world, fill = "grey", color = "grey", inherit.aes = FALSE) +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  labs(
    title = "Unbalanced Eddy Kinetic Energy",
    x = "Longitude",
    y = "Latitude",
    fill = "EKE (BL)"
  ) +
  theme_minimal() 


# Finally display
print(eke_balanced_map)
print(eke_unbalanced_map)

# Smoothing with GAMs
# Convert LON, LAT to numeric and ensure positivity of EKE
df_eke <- df_eke %>%
  mutate(
    LON_num = as.numeric(LON),
    LAT_num = as.numeric(LAT),
    EKE_BL_pos = ifelse(EKE_BL <= 0, NA_real_, EKE_BL),
    EKE_UNBL_pos = ifelse(EKE_UNBL <= 0, NA_real_, EKE_UNBL),
    logEKE_BL   = log(EKE_BL_pos),
    logEKE_UNBL = log(EKE_UNBL_pos)
  ) %>%
  # filter out any rows with NA (which includes any EKE <= 0)
  filter(!is.na(LON_num), !is.na(LAT_num),
         !is.na(logEKE_BL), !is.na(logEKE_UNBL))



# Gaussian on log(EKE_BL)
gam_bl_gaussian_log <- gam(
  logEKE_BL ~ s(LAT_num,LON_num, bs = "sos", k = 600),
  data   = df_eke,
  family = gaussian(),
  method = "REML"
)


gam_bl_gaussian_log %>% summary() #R-sq.(adj) =  0.913   Deviance explained = 92.3%

gam_unbl_gaussian_log <- gam(
  logEKE_UNBL ~ s(LAT_num,LON_num, bs = "sos", k = 600),
  data   = df_eke,
  family = gaussian(),
  method = "REML"
)


# Suppose we define a function to make a simple lat-lon grid:
make_prediction_grid <- function(df, lon_step = 0.25, lat_step = 0.25) {
  lon_seq <- seq(min(df$LON_num, na.rm = TRUE),
                 max(df$LON_num, na.rm = TRUE),
                 by = lon_step)
  lat_seq <- seq(min(df$LAT_num, na.rm = TRUE),
                 max(df$LAT_num, na.rm = TRUE),
                 by = lat_step)
  expand.grid(LON_num = lon_seq, LAT_num = lat_seq)
}

pred_grid <- make_prediction_grid(df_eke)

pred_grid$EKE_bl_hat_gaussian_log <- exp(
  predict(gam_bl_gaussian_log, newdata = pred_grid, type = "link")
)

pred_grid$EKE_unbl_hat_gaussian_log <- exp(
  predict(gam_unbl_gaussian_log, newdata = pred_grid, type = "link")
)



pred_grid$EKE_bl_hat_binned <- cut(pred_grid$EKE_bl_hat_gaussian_log,
                                   breaks = breaks_eke_bl,
                                   include.lowest = TRUE)

pred_grid$EKE_unbl_hat_binned <- cut(pred_grid$EKE_unbl_hat_gaussian_log,
                                     breaks = breaks_eke_unbl,
                                     include.lowest = TRUE)


pred_grid$EKE_bl_hat_binned %>% levels()


pred_grid$EKE_unbl_hat_binned %>% levels()



# Custom function to format legend labels
format_eke_labels <- function(bins) {
  labels <- gsub("[\\(\\)\\[\\]]", "", bins)  # Remove all types of brackets
  labels <- gsub(",", " - ", labels)          # Replace commas with " - "
  return(labels)
}

pred_grid <- pred_grid %>% na.omit() 

# Plot with improved aesthetics
map_eke_bl <- ggplot(pred_grid, aes(x = LON_num, y = LAT_num, fill = EKE_bl_hat_binned)) +
  
  # Tile plot with discrete color bins
  geom_tile(alpha = 0.95) +  # Slight transparency for better visibility
  # Custom discrete color scale with legend improvements
  scale_fill_viridis_d(
    option = "viridis",  
    name = "EKE (m²/s²)", # Updated legend title
    guide = guide_legend(reverse = TRUE,na.translate = FALSE)  # Reverse legend order
  ) + 
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = 53, ymax = 90),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -90, ymax = -67),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  # Add the world map
  geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
  # Coordinate system (keep only coord_sf)
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  
  # Labels and theme improvements
  labs(
    title = "Balanced Eddy Kinetic Energy",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +  # Increased text size for readability
  theme(
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 12, face = "bold"),  # Improve legend title
    legend.text = element_text(size = 10),  # Improve legend labels
    panel.grid = element_blank()  # Remove background grid for cleaner visuals
  )

# Display the plot
print(map_eke_bl)

# Unbalanced EKE Map
map_eke_unbl <- ggplot(pred_grid, aes(x = LON_num, y = LAT_num, fill = EKE_unbl_hat_binned)) +
  
  # Tile plot with discrete color bins
  geom_tile(alpha = 0.95) +  # Slight transparency for better visibility
  # Custom discrete color scale with legend improvements
  scale_fill_viridis_d(
    option = "viridis",  
    name = "EKE (m²/s²)", # Updated legend title
    guide = guide_legend(reverse = TRUE,na.translate = FALSE)  # Reverse legend order
  ) + 
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = 53, ymax = 90),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -90, ymax = -67),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  # Contour lines for additional information
  #geom_contour(
  #  data = pred_grid,
  #  aes(x = LON_num, y = LAT_num, z = EKE_unbl_hat_gaussian_log),
  #  color = "white", 
  #  alpha = 0.5,       # More subtle contrast
  #  bins = 10,
  #  linewidth = 0.5    # Thin contour lines for better readability
  #) +
  # Add the world map
  geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
  # Coordinate system (keep only coord_sf)
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  
  # Labels and theme improvements
  labs(
    title = "Unbalanced Eddy Kinetic Energy",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +  # Increased text size for readability
  theme(
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 12, face = "bold"),  # Improve legend title
    legend.text = element_text(size = 10),  # Improve legend labels
    panel.grid = element_blank()  # Remove background grid for cleaner visuals
  )

ggsave(plot=map_eke_bl,filename = "figures/smoothed_map_eke_bl.png",width = 10,height = 10)
ggsave(plot=map_eke_unbl,filename = "figures/smoothed_map_eke_unbl.png",width = 10,height = 10)

combined_discrete_eke_map <- (eke_balanced_map + map_subd_full) / (eke_unbalanced_map + map_subd_full)


combined_smoothed_eke_map <- (map_eke_bl + map_subd_full) / (map_eke_unbl + map_subd_full)

ggsave("figures/combined_discrete_eke_map.png", combined_discrete_eke_map, width = 14, height = 10, dpi = 300)

ggsave("figures/combined_smoothed_eke_map.png", combined_smoothed_eke_map, width = 14, height = 10, dpi = 300)

###############
## CHLOROPHYLL DATA from https://www.oceancolour.org/browser/
# Or from https://data.jrc.ec.europa.eu/dataset/d6f9abd9-777c-4a0c-a5f7-669612f83307#dataaccess
################
library(ncdf4)

chloro_data <- tidync("/data/GLOBARGO/data/chlor_a_mean_climatology_2010_2024_9km.nc")
chloro_data

df_chloro <- chloro_data %>% 
  hyper_tibble() %>%
  rename(LAT = lat, LON = lon, MONTH = month, CHLA = chlor_a_mean) %>% 
  mutate(LAT = as.numeric(LAT),
         LON = as.numeric(LON))

df_chloro_jan <- df_chloro %>% filter(MONTH == 1)

df_chloro_jan %>% str()

ggplot(df_chloro_jan, aes(x = LON, y = LAT, fill = log(CHLA))) +
  
  # Tile plot with discrete color bins
  geom_tile(alpha = 0.95) + scale_fill_viridis()+ # Slight transparency for better visibility
  # Add the world map
  geom_sf(data = world, fill = "darkgrey", color = "white", inherit.aes = FALSE) +
  # Coordinate system (keep only coord_sf)
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  
  # Labels and theme improvements
  labs(
    title = "Chlorophyll Alpha",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +  # Increased text size for readability
  theme(
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 12, face = "bold"),  # Improve legend title
    legend.text = element_text(size = 10),  # Improve legend labels
    panel.grid = element_blank()  # Remove background grid for cleaner visuals
  )

# Alternative EU dataset
chloro_data_jan <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_01.nc")
chloro_data_jan

df_chloro_jan <- chloro_data_jan %>% 
  hyper_tibble() %>%
  rename(LAT = lat, LON = lon, CHLA = Chl_a) %>% 
  mutate(LAT = as.numeric(LAT),
         LON = as.numeric(LON))

df_chloro_jan$chla_exp <- 10^(df_chloro_jan$CHLA)

df_chloro_jan$chla_exp %>% summary()
ggplot(df_chloro_jan, aes(x = LON, y = LAT, fill = chla_exp)) +
  
  # Tile plot with discrete color bins
  geom_tile(alpha = 0.95) + scale_fill_viridis(limits = c(0, 1))+ # Slight transparency for better visibility
  # Add the world map
  geom_sf(data = world, fill = "darkgrey", color = "white", inherit.aes = FALSE) +
  # Coordinate system (keep only coord_sf)
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  
  # Labels and theme improvements
  labs(
    title = "Chlorophyll Alpha",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +  # Increased text size for readability
  theme(
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 12, face = "bold"),  # Improve legend title
    legend.text = element_text(size = 10),  # Improve legend labels
    panel.grid = element_blank()  # Remove background grid for cleaner visuals
  )
View(df_chlor)

chloro_data_jul <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_07.nc")
chloro_data_jul

df_chloro_jul <- chloro_data_jul %>% 
  hyper_tibble() %>%
  rename(LAT = lat, LON = lon, CHLA = Chl_a) %>% 
  mutate(LAT = as.numeric(LAT),
         LON = as.numeric(LON))

df_chloro_jul$chla_exp <- 10^(df_chloro_jul$CHLA)

df_chloro_jul$chla_exp %>% summary()
p <- ggplot(df_chloro_jul, aes(x = LON, y = LAT, fill = chla_exp)) +
  
  # Tile plot with discrete color bins
  geom_tile(alpha = 0.95) + scale_fill_viridis(limits = c(0, 1))+ # Slight transparency for better visibility
  # Add the world map
  geom_sf(data = world, fill = "darkgrey", color = "white", inherit.aes = FALSE) +
  # Coordinate system (keep only coord_sf)
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  
  # Labels and theme improvements
  labs(
    title = "Climatology of Chlorophyll Alpha (July)",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal(base_size = 14) +  # Increased text size for readability
  theme(
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 12, face = "bold"),  # Improve legend title
    legend.text = element_text(size = 10),  # Improve legend labels
    panel.grid = element_blank()  # Remove background grid for cleaner visuals
  )

ggsave(plot = p,"figures/chlorophyll_seawifs_chla_july.png")
