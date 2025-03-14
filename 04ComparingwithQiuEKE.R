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
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(viridis)
library(sp)
library(spdep)
library(mgcv)
library(patchwork)
library(scales)
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
  lon_seq <- seq(-180,180,
                 by = lon_step)
  lat_seq <- seq(-90,
                 90,
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

# Define global scales for probability contours (same as before)
global_contour_breaks <- c(0.05, 0.10, 0.15, 0.20, 0.25)
global_contour_labels <- as.character(round(global_contour_breaks * 100, 0))


# 1. Aggregate Argo float locations to 5° bins for the season
argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 180, by = 5),
  lat_bin_stipple = seq(-90, 90, by = 5)
) %>%
  as_tibble()

# Left join the existing argo_bins to the full grid and replace NAs with 0
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(count = ifelse(is.na(count), 0, count))


# 2. Identify undersampled areas (e.g., fewer than 5 profiles)
undersampled <- argo_bins_full %>% filter(count < 1)




undersampled_corners <- undersampled %>%
  rowwise() %>%
  mutate(corners = list(
    data.frame(
      LON = c(lon_bin_stipple,
              lon_bin_stipple + stipple_resolution,
              lon_bin_stipple,
              lon_bin_stipple + stipple_resolution),
      LAT = c(lat_bin_stipple,
              lat_bin_stipple,
              lat_bin_stipple + stipple_resolution,
              lat_bin_stipple + stipple_resolution)
    )
  )) %>%
  ungroup() %>%
  unnest(corners)

# # ----- Prepare the prediction grid for probability of subduction (yearly) -----
# # (Assumes pred_carb_year exists with columns: lon_bin, lat_bin, proportion)
# prediction_grid_year <- pred_carb_year %>%
#   mutate(
#     lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
#     lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
#   ) %>%
#   left_join(argo_bins_year, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
#   mutate(undersampled = ifelse(is.na(count), FALSE, TRUE))
# 
# # ----- Clean and prepare the EKE prediction grid (yearly) -----
# # (Assumes pred_grid contains the balanced EKE data with columns LON_num, LAT_num, and EKE_bl_hat_binned)
# pred_grid <- pred_grid %>% na.omit()
# 
# # Build undersampled corners
# argo_bins <- df_argo_clean %>%
#   mutate(
#     lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
#     lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
#   ) %>%
#   group_by(lon_bin_stipple, lat_bin_stipple) %>%
#   summarize(count = n(), .groups = "drop")
# 
# # 2. Identify undersampled areas (e.g., fewer than 5 profiles)
# undersampled <- argo_bins %>% filter(count < 5)
# 
# # 3. Generate corner points for each undersampled cell (4 corners per cell)
# undersampled_corners <- undersampled %>%
#   rowwise() %>%
#   mutate(corners = list(
#     data.frame(
#       LON = c(lon_bin_stipple,
#               lon_bin_stipple + stipple_resolution,
#               lon_bin_stipple,
#               lon_bin_stipple + stipple_resolution),
#       LAT = c(lat_bin_stipple,
#               lat_bin_stipple,
#               lat_bin_stipple + stipple_resolution,
#               lat_bin_stipple + stipple_resolution)
#     )
#   )) %>%
#   ungroup() %>%
#   unnest(corners)
# 
# ggplot()+geom_contour(data = pred_full_subd,
#              aes(x = lon_bin, y = lat_bin, z = proportion,
#                  color = factor(round(after_stat(level) * 100, 0))),
#              breaks = global_contour_breaks,
#              alpha = 1,
#              inherit.aes = FALSE,
#              show.legend = TRUE)
# 
# # ----- Build the Plot -----
# Manually defined labels for the balanced EKE legend.
eke_labels_bl <- c(
  "< 1e-4",
  "(1e-4, 1e-3)",
  "(1e-3, 2.5e-3)",
  "(2.5e-3, 5e-3)",
  "(5e-3, 7.5e-3)",
  "(7.5e-3, 1e-2)",
  "(1e-2, 2.5e-2)",
  "(2.5e-2, 5e-2)",
  "(5e-2, 7.5e-2)",
  "(7.5e-2, 1e-1)",
  "(1e-1, 1)"
)

# Manually defined labels for the unbalanced EKE legend.
eke_labels_unbl <- c(
  "< 1e-4",
  "(1e-4, 1e-3)",
  "(1e-3, 2e-3)",
  "(2e-3, 3e-3)",
  "(3e-3, 4e-3)",
  "(4e-3, 5e-3)",
  "(5e-3, 6e-3)",
  "(6e-3, 7e-3)",
  "(7e-3, 8e-3)",
  "(8e-3, 9e-3)",
  "(9e-3, 1e-2)",
  "(1e-2, 2.5e-2)",
  "(2.5e-2, 1e-1)"
)

# Define a common theme without panel borders
common_theme <- theme_bw(base_size = 25) +
  theme(
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

map_eke_bl <- ggplot() +
  # 2. Rectangles for polar regions (optional)
  geom_rect(aes(xmin = -180, xmax = 180, ymin = 50, ymax = 90),
            fill = "lightblue", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = -60),
            fill = "lightblue", inherit.aes = FALSE) +
  
  # 1. Tile plot for yearly EKE background
  geom_tile(data = pred_grid, 
            aes(x = LON_num, y = LAT_num, fill = EKE_bl_hat_binned),
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Balanced EKE\n (m²/s²)",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE),
    labels=eke_labels_bl,
  ) +
  # 3. Overlay probability of subduction contours in red
  new_scale_color() +
  geom_contour(data = pred_full_subd,
               aes(x = lon_bin, y = lat_bin, z = proportion,
                   color = factor(round(after_stat(level) * 100, 0))),
               breaks = global_contour_breaks,
               alpha = 0.75,
               linewidth = 1.5,
               inherit.aes = FALSE,
               show.legend = TRUE) +
  scale_color_brewer(
    palette = "Reds",
    name = "Probability of \nsubduction (%)",
    limits = global_contour_labels,
    breaks = global_contour_labels,
    labels = function(x) paste0(x, "%")
  ) +
  
  # 4. Add stippling for undersampled cells
  geom_point(data = undersampled_corners,
             aes(x = LON, y = LAT),
             color = "white", alpha = 0.6, size = 1.5, shape = 20) +
  
  # 5. Overlay the world map
  geom_sf(data = world,
          fill = "white", color = "white",
          inherit.aes = FALSE) +
  
  # 6. Coordinate system and labels
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = paste("Balanced Eddy Kinetic Energy & Subduction Probability"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank()
  )

# Display the plot
print(map_eke_bl)
# Unbalanced EKE Map

map_eke_unbl <- ggplot() +
  # 2. Rectangles for polar regions (optional)
  geom_rect(aes(xmin = -180, xmax = 180, ymin = 50, ymax = 90),
            fill = "lightblue", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = -60),
            fill = "lightblue", inherit.aes = FALSE) +
  
  # 1. Tile plot for yearly EKE background
  geom_tile(data = pred_grid, 
            aes(x = LON_num, y = LAT_num, fill = EKE_unbl_hat_binned),
            alpha = 0.75) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Unbalanced EKE\n (m²/s²)",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE),
    labels=eke_labels_unbl,
  ) +
  # 3. Overlay probability of subduction contours in red
  new_scale_color() +
  geom_contour(data = pred_full_subd,
               aes(x = lon_bin, y = lat_bin, z = proportion,
                   color = factor(round(after_stat(level) * 100, 0))),
               breaks = global_contour_breaks,
               alpha = 1,
               linewidth = 1.5,
               inherit.aes = FALSE,
               show.legend = TRUE) +
  scale_color_brewer(
    palette = "Reds",
    name = "Probability of \nsubduction (%)",
    limits = global_contour_labels,
    breaks = global_contour_labels,
    labels = function(x) paste0(x, "%")
  ) +
  
  # 4. Add stippling for undersampled cells
  geom_point(data = undersampled_corners,
             aes(x = LON, y = LAT),
             color = "white", alpha = 0.6, size = 1.5, shape = 20) +
  
  # 5. Overlay the world map
  geom_sf(data = world,
          fill = "white", color = "white",
          inherit.aes = FALSE) +
  
  # 6. Coordinate system and labels
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = paste("Unbalanced Eddy Kinetic Energy & Subduction Probability"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank()
  )



ggsave(plot=map_eke_bl,filename = "figures/smoothed_map_eke_bl.png",width = 18,height = 10,dpi = 300)
ggsave(plot=map_eke_unbl,filename = "figures/smoothed_map_eke_unbl.png",width = 18,height = 10,dpi = 300)

ombined_discrete_eke_map <- (eke_balanced_map + map_subd_full) / (eke_unbalanced_map + map_subd_full)


combined_smoothed_eke_map <- (map_eke_bl + map_subd_full) / (map_eke_unbl + map_subd_full)

ggsave("figures/combined_discrete_eke_map.png", combined_discrete_eke_map, width = 14, height = 10, dpi = 300)

ggsave("figures/combined_smoothed_eke_map.png", combined_smoothed_eke_map, width = 14, height = 10, dpi = 300)

###############
## CHLOROPHYLL DATA from https://www.oceancolour.org/browser/
# Or from https://data.jrc.ec.europa.eu/dataset/d6f9abd9-777c-4a0c-a5f7-669612f83307#dataaccess
################
library(ncdf4)

chloro_data <- tidync("/data/GLOBARGO/data/chlor_a_mean_climatology_2010_2024_9km(1).nc")
chloro_data
df_chloro <- chloro_data %>% 
  hyper_tibble() %>%
  rename(LATITUDE = lat, LONGITUDE = lon, MONTH = month, CHLA = chlor_a_mean) %>% 
  mutate(LAT = as.numeric(LATITUDE),
         LON = as.numeric(LONGITUDE))

df_chloro$MONTH <- df_chloro$MONTH %>% as.numeric()


add_month_region <- function(df) {
  df %>%
    mutate(
      month = month(MONTH),  # Extract month as a factor
      region = case_when(
        LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
        LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
        LATITUDE >= 0 & LATITUDE <= 30 ~ "Northern Tropics",                     # 0° to 30°N
        LATITUDE < 0 & LATITUDE >= -30 ~ "Southern Tropics",                     # 0° to 30°S
        LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
        TRUE ~ NA_character_  # Exclude undefined regions
      )
    )
}

df_chloro_with_regions <- add_month_region(df_chloro)

df_chloro_summary <- df_chloro_with_regions %>%
  group_by(region,month) %>%
  summarize(CHLA_med = median(CHLA,na.rm = T))

write_csv(df_chloro_summary,
          file = "/data/GLOBARGO/data/df_chloro_monthly_climatology_aqua_modis_2010_2024.csv")

# ----- Prepare Seasonal Chlorophyll Data -----
df_chloro_djf <- df_chloro %>% 
  filter(MONTH %in% c(12, 1, 2)) %>% 
  group_by(LAT, LON) %>% 
  summarise(CHLA = mean(CHLA,na.rm = T), .groups = "drop")

df_chloro_mam <- df_chloro %>% 
  filter(MONTH %in% c(3, 4, 5)) %>% 
  group_by(LAT, LON) %>% 
  summarise(CHLA = mean(CHLA,na.rm = T), .groups = "drop")

df_chloro_jja <- df_chloro %>% 
  filter(MONTH %in% c(6, 7, 8)) %>% 
  group_by(LAT, LON) %>% 
  summarise(CHLA = mean(CHLA,na.rm = T), .groups = "drop")

df_chloro_son <- df_chloro %>% 
  filter(MONTH %in% c(9, 10, 11)) %>% 
  group_by(LAT, LON) %>% 
  summarise(CHLA = mean(CHLA,na.rm = T), .groups = "drop")

write_csv(df_chloro_djf,"/data/GLOBARGO/data/df_chloro_djf_climatology_aquamodis.csv")
write_csv(df_chloro_mam,"/data/GLOBARGO/data/df_chloro_mam_climatology_aquamodis.csv")
write_csv(df_chloro_jja,"/data/GLOBARGO/data/df_chloro_jja_climatology_aquamodis.csv")
write_csv(df_chloro_son,"/data/GLOBARGO/data/df_chloro_son_climatology_aquamodis.csv")


# ----- Set Parameters for Stippling -----
stipple_resolution <- 5  # degrees
# Define global scales (compute these once, based on your full dataset)
global_chla_limits <- c(0.01, 85)  # Adjust as needed
global_contour_breaks <- c(0.01,0.05,0.10,0.15,0.20,0.25)  # Underlying proportion values
global_contour_labels <- as.character(round(global_contour_breaks * 100, 0))

plot_season <- function(df_chloro_season, pred_carb, argo_months, season_label) {
  
  # 1. Aggregate Argo float locations to 5° bins for the season
  argo_bins <- df_argo_clean %>%
    filter(month(TIME) %in% argo_months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  # 2. Identify undersampled areas (e.g., fewer than 5 profiles)
  undersampled <- argo_bins %>% filter(count < 5)
  
  # 3. Generate corner points for each undersampled cell (4 corners per cell)
  undersampled_corners <- undersampled %>%
    rowwise() %>%
    mutate(corners = list(
      data.frame(
        LON = c(lon_bin_stipple,
                lon_bin_stipple + stipple_resolution,
                lon_bin_stipple,
                lon_bin_stipple + stipple_resolution),
        LAT = c(lat_bin_stipple,
                lat_bin_stipple,
                lat_bin_stipple + stipple_resolution,
                lat_bin_stipple + stipple_resolution)
      )
    )) %>%
    ungroup() %>%
    unnest(corners)
  
  # 4. Prepare the prediction grid for the season
  prediction_grid_season <- pred_carb %>%
    mutate(
      lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
    ) %>%
    left_join(undersampled, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(undersampled = ifelse(is.na(count), FALSE, TRUE))
  
  # 5. Build the ggplot for this season
  p <- ggplot() +
    # Main tile layer for chlorophyll data
    geom_tile(data = df_chloro_season,
              aes(x = LON, y = LAT, fill = CHLA),
              alpha = 0.75) +
    scale_fill_viridis(
      trans = "log",
      limits = global_chla_limits,
      name = expression("Chlorophyll " * alpha ~ " (g/m"^3*")"),
      breaks = c(0.01, 0.1, 1, 10),
      labels = comma_format()
    ) +
    
    # New color scale for the contour layer
    new_scale_color() +
    geom_contour(data = prediction_grid_season,
                 aes(x = lon_bin, y = lat_bin, z = proportion,
                     color = factor(round(after_stat(level) * 100, 0))),
                 breaks = global_contour_breaks,
                 alpha = 1,
                 inherit.aes = FALSE,
                 show.legend = TRUE,
                 straight = TRUE) +
    scale_color_brewer(
      palette = "Reds",  # You can choose palettes like "Set1", "Dark2", "Paired", etc.
      name = "Probability of carbon subduction (%)",
      limits = global_contour_labels,
      breaks = global_contour_labels,
      labels = function(x) paste0(x, "%")
    )+
    
    # Add stippling: plot four corner points per undersampled grid cell
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT),
               color = "black", alpha = 0.6, size = 1.5, shape = 20) +
    
    # Add world map layer
    geom_sf(data = world,
            fill = "lightgrey", color = "lightgrey",
            inherit.aes = FALSE) +
    
    # Set coordinate system
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    
    # Labels and theme
    labs(
      title = paste("(", season_label, ")", sep = ""),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, "cm"),
      panel.grid = element_blank()
    )
  
  return(p)
}

plot_djf <- plot_season(df_chloro_djf, pred_djf_carb, c(12, 1, 2), "DJF")
plot_mam <- plot_season(df_chloro_mam, pred_mam_carb, c(3, 4, 5), "MAM")
plot_jja <- plot_season(df_chloro_jja, pred_jja_carb, c(6, 7, 8), "JJA")
plot_son <- plot_season(df_chloro_son, pred_son_carb, c(9, 10, 11), "SON")

# ----- Save Individual Seasonal Plots -----
ggsave("AQUA_Chloro_DJF.png", plot = plot_djf, width = 12, height = 10, dpi = 300)
ggsave("AQUA_Chloro_MAM.png", plot = plot_mam, width = 12, height = 10, dpi = 300)
ggsave("AQUA_Chloro_JJA.png", plot = plot_jja, width = 12, height = 10, dpi = 300)
ggsave("AQUA_Chloro_SON.png", plot = plot_son, width = 12, height = 10, dpi = 300)



# Alternative EU dataset
# Define the months you want to import.
# (Change to 1:12 if you have files for all 12 months.)
months <- 1:12

# Function to import and process a single month
read_chloro_month <- function(month) {
  # Construct the file name assuming two-digit month numbers
  file_path <- sprintf("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_%02d.nc", month)
  message("Reading file: ", file_path)
  
  df <- tidync(file_path) %>% 
    hyper_tibble() %>%
    rename(LAT = lat, LON = lon, CHLA = Chl_a) %>% 
    mutate(
      LAT = as.numeric(LAT),
      LON = as.numeric(LON),
      MONTH = month
    )
  
  # Convert from log10-scale to actual chlorophyll values
  df <- df %>% mutate(CHLA = 10^(CHLA))
  
  return(df)
}

# Import and combine the data for all months
df_chloro_all <- map_dfr(months, read_chloro_month)

# (Optional) Check the first few rows of the combined data
head(df_chloro_all)

df_chloro_all$LATITUDE <- df_chloro_all$LAT
df_chloro_all$LONGITUDE <- df_chloro_all$LON

# --- Assign Regions ---
# For this sanity check we use a simple classification:
#   - Southern Ocean: LAT < -40
#   - North Atlantic: points in a specific bounding box (e.g., LON between -100 and 20 and LAT between 0 and 60)
#   - Pacific Ocean: all remaining points
df_chloro_all <- df_chloro_all %>%
  mutate(
    region = case_when(
      LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
      LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
      LATITUDE >= 0 & LATITUDE <= 30 ~ "Northern Tropics",                     # 0° to 30°N
      LATITUDE < 0 & LATITUDE >= -30 ~ "Southern Tropics",                     # 0° to 30°S
      LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
      TRUE ~ NA_character_  # Exclude undefined regions
    )
  )
# --- Compute Summary Statistics ---
# Here we compute the median and mean chlorophyll per month and region.
df_chloro_summary <- df_chloro_all %>%
  group_by(region, MONTH) %>%
  summarise(
    CHLA_med = median(CHLA, na.rm = TRUE),
    CHLA_mean = mean(CHLA, na.rm = TRUE),
    .groups = "drop"
  )

# Check the summary dataset
print(df_chloro_summary)

write_csv(df_chloro_summary,"/data/GLOBARGO/data/df_chloro_monthly_climatology_seawifs.csv")

# --- Plot as a Sanity Check ---
# This plot shows the median chlorophyll per month for each region.
ggplot(df_chloro_summary, aes(x = MONTH, y = CHLA_mean)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  facet_wrap(~ region) +
  labs(
    title = "Mean Chlorophyll (Seawifs) per Month by Region",
    x = "Month",
    y = "Mean Chlorophyll (mg m^-3)"
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold", size = 12))



chloro_data_jan <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_01.nc")
chloro_data_feb <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_02.nc")
chloro_data_march <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_03.nc")
chloro_data_april <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_04.nc")
chloro_data_may <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_05.nc")
chloro_data_jun <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_06.nc")
chloro_data_jul <- tidync("/data/GLOBARGO/data/chlorophyll_monthly_seawifs/GMIS_V_CHLA_06.nc")


df_chloro_jan <- chloro_data_jan %>% 
  hyper_tibble() %>%
  rename(LAT = lat, LON = lon, CHLA = Chl_a) %>% 
  mutate(LAT = as.numeric(LAT),
         LON = as.numeric(LON),
         MONTH = 1)

df_chloro_jan$CHLA <- 10^(df_chloro_jan$CHLA)

df_chloro_jan %>% head()

df_chloro_feb <- chloro_data_feb %>% 
  hyper_tibble() %>%
  rename(LAT = lat, LON = lon, CHLA = Chl_a) %>% 
  mutate(LAT = as.numeric(LAT),
         LON = as.numeric(LON),
         MONTH = 2)


df_chloro_feb$chla <- 10^(df_chloro_feb$CHLA)

df_chloro_dec <- chloro_data_dec %>% 
  hyper_tibble() %>%
  rename(LAT = lat, LON = lon, CHLA = Chl_a) %>% 
  mutate(LAT = as.numeric(LAT),
         LON = as.numeric(LON),
         MONTH = 12)


df_chloro_dec$chla <- 10^(df_chloro_dec$CHLA)

df_chloro_djf2 <- dplyr::bind_rows(df_chloro_feb,df_chloro_jan,df_chloro_dec)
df_chloro_djf2 <- df_chloro_djf2 %>% group_by(LAT,LON) %>% summarise(chla = mean(chla,na.rm=T))
df_chloro_djf2$CHLA <- df_chloro_djf2$chla

plot_djf_v2 <- plot_season(df_chloro_djf2, pred_djf_carb, c(12, 1, 2), "DJF")
ggsave("AQUA_Chloro_DJF_seawifs.png", plot = plot_djf, width = 12, height = 10, dpi = 300)


# ----- Create and Save Combined Plot with Collected Legends -----
g <- ggarrange(plot_djf_v2,plot_mam,plot_jja,plot_son,
               nrow = 2,ncol= 2, common.legend = T, legend = "bottom")

ggsave("AQUA_Chloro_Combined_Seasons.png", plot = g, width = 20, height = 10, dpi = 300)
