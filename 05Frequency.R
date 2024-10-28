library(tidyverse)
library(sf)
library(ggplot2)
library(mgcv)

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

# first investigate structure of data
df_complete_clean <- read_csv(file = "data/df_eddy_subduction_anom.csv")
df_complete_clean %>% head()
df_argo_clean <- read_csv(file = "data/df_argo_loc.csv")
df_argo_clean %>% head()

# Define bin size for longitude and latitude
bin_size <- 5
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE)), ceiling(max(df_argo_clean$LONGITUDE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE)), ceiling(max(df_argo_clean$LATITUDE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Function to assign bins to longitude and latitude
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

# Apply binning to both df_argo_clean and df_complete_clean
df_argo_clean <- assign_bins(df_argo_clean)
df_complete_clean <- assign_bins(df_complete_clean)


# Define months for the four distinct periods: DJF, MAM, JJA, SON
djf_months <- c(12, 1, 2)   # December, January, February
mam_months <- c(3, 4, 5)    # March, April, May
jja_months <- c(6, 7, 8)    # June, July, August
son_months <- c(9, 10, 11)  # September, October, November

# Calculate total months for each season (using full years)
seasonal_months <- list(
  DJF = length(djf_months),
  MAM = length(mam_months),
  JJA = length(jja_months),
  SON = length(son_months)
)

# Function to compute event frequency for a specific season
compute_seasonal_frequency <- function(df_argo, df_complete, months_in_season) {
  df_complete_season <- df_complete %>%
    filter(month(TIME) %in% months_in_season)
  
  # Compute frequency of subduction events per grid cell per month for the season
  subduction_frequency <- df_complete_season %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(
      total_events = n(),
      frequency_per_month = total_events / (length(unique(year(TIME))) * length(months_in_season)),  # Events per month
      .groups = 'drop'
    )
  
  return(subduction_frequency)
}

# Compute seasonal frequencies for DJF, MAM, JJA, SON
frequency_djf <- compute_seasonal_frequency(df_argo_clean, df_complete_clean, djf_months)
frequency_mam <- compute_seasonal_frequency(df_argo_clean, df_complete_clean, mam_months)
frequency_jja <- compute_seasonal_frequency(df_argo_clean, df_complete_clean, jja_months)
frequency_son <- compute_seasonal_frequency(df_argo_clean, df_complete_clean, son_months)

# Convert to spatial objects for each season (transform to Robinson projection)
argo_with_frequency_djf_sf <- df_argo_clean %>%
  left_join(frequency_djf, by = c("lon_bin", "lat_bin")) %>%
  filter(!is.na(frequency_per_month)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

argo_with_frequency_mam_sf <- df_argo_clean %>%
  left_join(frequency_mam, by = c("lon_bin", "lat_bin")) %>%
  filter(!is.na(frequency_per_month)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

argo_with_frequency_jja_sf <- df_argo_clean %>%
  left_join(frequency_jja, by = c("lon_bin", "lat_bin")) %>%
  filter(!is.na(frequency_per_month)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

argo_with_frequency_son_sf <- df_argo_clean %>%
  left_join(frequency_son, by = c("lon_bin", "lat_bin")) %>%
  filter(!is.na(frequency_per_month)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

# Check and filter out non-finite values before plotting
frequency_djf <- frequency_djf %>% filter(is.finite(frequency_per_month))
frequency_mam <- frequency_mam %>% filter(is.finite(frequency_per_month))
frequency_jja <- frequency_jja %>% filter(is.finite(frequency_per_month))
frequency_son <- frequency_son %>% filter(is.finite(frequency_per_month))


# Define discrete bins without using log transformation
binned_color_scale <- scale_fill_viridis_b(
  name = "Events/Month",
  breaks = c(0.1, 0.5, 1, 2, 4, 8),  # Define breaks from 0.1 to 8
  limits = c(0.1, 8),  # Set limits to match the range
  labels = c("0.1", "0.5", "1", "2", "4", "8"),  # Keep labels in original values
  oob = scales::squish  # Squish values outside the limit into closest bin
)

# Function to plot frequency map with discrete bins
plot_frequency_tile_map_discrete <- function(frequency_data, title, output_file) {
  ggplot(frequency_data, aes(x = lon_bin, y = lat_bin)) +
    geom_tile(aes(fill = frequency_per_month)) +  # Directly use the frequency_per_month
    geom_contour(aes(z = frequency_per_month), color = "white", alpha = 0.1) +
    geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
    coord_sf(xlim = range(frequency_data$lon_bin), ylim = range(frequency_data$lat_bin)) +
    binned_color_scale +  # Apply the discrete color scale
    labs(title = title, x = "Longitude", y = "Latitude", fill = "Events/Month") +
    theme_minimal()
  
  # Save the plot
  ggsave(output_file, plot = last_plot(), width = 10, height = 6)
}

# Plot the frequency maps for DJF, MAM, JJA, SON using the discrete scale from 0.1 to 8
plot_frequency_tile_map_discrete(frequency_djf, "Subduction Event Frequency (DJF)", "figures/TimeSpaceVar/subduction_event_frequency_djf_discrete.png")
plot_frequency_tile_map_discrete(frequency_mam, "Subduction Event Frequency (MAM)", "figures/TimeSpaceVar/subduction_event_frequency_mam_discrete.png")
plot_frequency_tile_map_discrete(frequency_jja, "Subduction Event Frequency (JJA)", "figures/TimeSpaceVar/subduction_event_frequency_jja_discrete.png")
plot_frequency_tile_map_discrete(frequency_son, "Subduction Event Frequency (SON)", "figures/TimeSpaceVar/subduction_event_frequency_son_discrete.png")

