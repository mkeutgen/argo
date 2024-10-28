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

conflicts_prefer(dplyr::filter)

# Step 1: Investigate the structure of the data
df_complete_clean <- read_csv(file = "data/df_eddy_subduction_anom.csv")
df_complete_clean %>% head()
df_argo_clean <- read_csv(file = "data/df_argo_loc.csv")
df_argo_clean %>% head()

# Step 2: Define months for the four distinct seasons: DJF, MAM, JJA, SON
djf_months <- c(12, 1, 2)   # December, January, February
mam_months <- c(3, 4, 5)    # March, April, May
jja_months <- c(6, 7, 8)    # June, July, August
son_months <- c(9, 10, 11)  # September, October, November

# Step 3: Filter data for each season based on the month of the year
df_argo_djf <- df_argo_clean %>% filter(month(TIME) %in% djf_months)
df_argo_mam <- df_argo_clean %>% filter(month(TIME) %in% mam_months)
df_argo_jja <- df_argo_clean %>% filter(month(TIME) %in% jja_months)
df_argo_son <- df_argo_clean %>% filter(month(TIME) %in% son_months)

df_complete_djf <- df_complete_clean %>% filter(month(TIME) %in% djf_months)
df_complete_mam <- df_complete_clean %>% filter(month(TIME) %in% mam_months)
df_complete_jja <- df_complete_clean %>% filter(month(TIME) %in% jja_months)
df_complete_son <- df_complete_clean %>% filter(month(TIME) %in% son_months)

# Step 4: Bin the data into longitude and latitude grids

# Define bin size for longitude and latitude
bin_size <- 5
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE)), ceiling(max(df_argo_clean$LONGITUDE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE)), ceiling(max(df_argo_clean$LATITUDE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Function to assign bins to Argo profiles and anomalies for each season
assign_bins <- function(df) {
  df %>%
    mutate(
      lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
      lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
    ) %>%
    mutate(
      lon_bin = as.numeric(as.character(lon_bin)),
      lat_bin = as.numeric(as.character(lat_bin))
    )
}

# Step 5: Assign bins to Argo profiles and anomalies
df_argo_djf <- assign_bins(df_argo_djf)
df_argo_mam <- assign_bins(df_argo_mam)
df_argo_jja <- assign_bins(df_argo_jja)
df_argo_son <- assign_bins(df_argo_son)

df_complete_djf <- assign_bins(df_complete_djf)
df_complete_mam <- assign_bins(df_complete_mam)
df_complete_jja <- assign_bins(df_complete_jja)
df_complete_son <- assign_bins(df_complete_son)

# Step 6: Compute counts of total profiles and subduction events for each season

compute_counts <- function(df_argo, df_complete, years = 14) {
  total_counts <- df_argo %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_total = n(), .groups = 'drop')
  
  anomaly_counts <- df_complete %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_anomaly = n(), .groups = 'drop')
  
  # Merge counts and compute frequency per month by dividing by number of years (14)
  merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin")) %>%
    mutate(
      count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
      frequency_per_month = count_anomaly / (years * 3)  # Divide by the total number of months in the season (3 months)
    ) %>%
    filter(count_total > 0) %>%
    filter(!is.na(lat_bin))
  
  return(merged_counts)
}

# Step 7: Compute frequency of subduction events per month for each season
frequency_djf <- compute_counts(df_argo_djf, df_complete_djf)
frequency_mam <- compute_counts(df_argo_mam, df_complete_mam)
frequency_jja <- compute_counts(df_argo_jja, df_complete_jja)
frequency_son <- compute_counts(df_argo_son, df_complete_son)

# Step 8: Fit Poisson GAM model for each season

fit_seasonal_gam <- function(data, season) {
  gam_model <- gam(
    count_anomaly ~ s(lat_bin, lon_bin, bs = "sos", k = 200),  # Spatial smoothing
    family = poisson(link = "log"),                           # Poisson for counts
    offset = log(count_total),                                # Offset for number of Argo profiles
    data = data,                                              # Data for the season
    method = "REML"
  )
  cat("\nModel for", season, "\n")
  summary(gam_model)
  return(gam_model)
}

# Fit models for each season
gam_djf <- fit_seasonal_gam(frequency_djf, "DJF")
gam_mam <- fit_seasonal_gam(frequency_mam, "MAM")
gam_jja <- fit_seasonal_gam(frequency_jja, "JJA")
gam_son <- fit_seasonal_gam(frequency_son, "SON")

# Step 9: Create prediction grids for each season

create_prediction_grid <- function(data) {
  lon_seq <- seq(min(data$lon_bin), max(data$lon_bin), by = 1)
  lat_seq <- seq(min(data$lat_bin), max(data$lat_bin), by = 1)
  expand.grid(lon_bin = lon_seq, lat_bin = lat_seq)
}

# Create prediction grids for each season
prediction_grid_djf <- create_prediction_grid(frequency_djf)
prediction_grid_mam <- create_prediction_grid(frequency_mam)
prediction_grid_jja <- create_prediction_grid(frequency_jja)
prediction_grid_son <- create_prediction_grid(frequency_son)

# Step 10: Predict the frequency of subduction events for each grid cell

prediction_grid_djf$frequency_per_month <- predict(gam_djf, newdata = prediction_grid_djf, type = "response")
prediction_grid_mam$frequency_per_month <- predict(gam_mam, newdata = prediction_grid_mam, type = "response")
prediction_grid_jja$frequency_per_month <- predict(gam_jja, newdata = prediction_grid_jja, type = "response")
prediction_grid_son$frequency_per_month <- predict(gam_son, newdata = prediction_grid_son, type = "response")

# Function to adjust frequency from 5x5 degree grid to 1x1 degree resolution
adjust_frequency_for_resolution <- function(prediction_grid) {
  prediction_grid <- prediction_grid %>%
    mutate(frequency_per_month = frequency_per_month / 25)  # Adjust by factor of 25
  return(prediction_grid)
}

# Adjust the predicted frequencies for each season
prediction_grid_djf <- adjust_frequency_for_resolution(prediction_grid_djf)
prediction_grid_mam <- adjust_frequency_for_resolution(prediction_grid_mam)
prediction_grid_jja <- adjust_frequency_for_resolution(prediction_grid_jja)
prediction_grid_son <- adjust_frequency_for_resolution(prediction_grid_son)

# Updated binned color scale for frequency maps (can keep the same scale)
binned_color_scale_gam <- scale_fill_viridis_b(
  name = "Events/Month",
  oob = scales::squish
)

# Function to plot the adjusted frequency map with updated scaling
plot_gam_frequency_map <- function(prediction_grid, title) {
  ggplot() +
    geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = frequency_per_month)) +
    geom_contour(data = prediction_grid, aes(x = lon_bin, y = lat_bin, z = frequency_per_month), color = "white", alpha = 0.3) +
    geom_sf(data = world, fill = "gray80", color = "gray80") +
    coord_sf(xlim = range(prediction_grid$lon_bin), ylim = range(prediction_grid$lat_bin), expand = FALSE) +
    binned_color_scale_gam +
    labs(title = title, x = "Longitude", y = "Latitude", fill = "Events/Month") +
    theme_minimal()
}

# Recreate the plots for each season with the adjusted frequency
gam_map_djf <- plot_gam_frequency_map(prediction_grid_djf, "Predicted Frequency of Subduction Events (DJF, Adjusted for 1x1)")
gam_map_mam <- plot_gam_frequency_map(prediction_grid_mam, "Predicted Frequency of Subduction Events (MAM, Adjusted for 1x1)")
gam_map_jja <- plot_gam_frequency_map(prediction_grid_jja, "Predicted Frequency of Subduction Events (JJA, Adjusted for 1x1)")
gam_map_son <- plot_gam_frequency_map(prediction_grid_son, "Predicted Frequency of Subduction Events (SON, Adjusted for 1x1)")

# Combine all maps into one layout with the adjusted frequency
combined_gam_maps <- gam_map_djf + gam_map_mam + gam_map_jja + gam_map_son + 
  plot_layout(ncol = 2, nrow = 2, guides = "collect") + 
  plot_annotation(title = "Predicted Frequency of Subduction Events (Adjusted for 1x1 Grid, Events/Month)")

# Save the updated combined plot
ggsave(filename = "figures/TimeSpaceVar/4SEASONS/combined_gam_frequency_maps_adjusted.png", plot = combined_gam_maps, width = 15, height = 12)
