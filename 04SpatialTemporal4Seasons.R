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

# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

df_argo_clean <- read_csv("data/df_argo_loc.csv")
df_complete_clean <- read_csv("data/df_eddy_subduction_anom.csv")
df_carbon_clean <- read_csv("data/df_carbon_subduction_anom.csv")


# Define months for the four distinct periods: DJF, MAM, JJA, SON
djf_months <- c(12, 1, 2)   # December, January, February
mam_months <- c(3, 4, 5)    # March, April, May
jja_months <- c(6, 7, 8)    # June, July, August
son_months <- c(9, 10, 11)  # September, October, November

# Filter data for DJF, MAM, JJA, and SON
df_argo_djf <- df_argo_clean %>% filter(month(TIME) %in% djf_months)
df_argo_mam <- df_argo_clean %>% filter(month(TIME) %in% mam_months)
df_argo_jja <- df_argo_clean %>% filter(month(TIME) %in% jja_months)
df_argo_son <- df_argo_clean %>% filter(month(TIME) %in% son_months)

df_complete_djf <- df_complete_clean %>% filter(month(TIME) %in% djf_months)
df_complete_mam <- df_complete_clean %>% filter(month(TIME) %in% mam_months)
df_complete_jja <- df_complete_clean %>% filter(month(TIME) %in% jja_months)
df_complete_son <- df_complete_clean %>% filter(month(TIME) %in% son_months)

df_carbon_djf <- df_carbon_clean %>% filter(month(TIME) %in% djf_months)
df_carbon_mam <- df_carbon_clean %>% filter(month(TIME) %in% mam_months)
df_carbon_jja <- df_carbon_clean %>% filter(month(TIME) %in% jja_months)
df_carbon_son <- df_carbon_clean %>% filter(month(TIME) %in% son_months)




# Convert data frames to spatial objects and transform to Robinson projection
df_argo_djf_sf <- df_argo_djf %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_mam_sf <- df_argo_mam %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_jja_sf <- df_argo_jja %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_son_sf <- df_argo_son %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

# Load world map data and transform to Robinson projection
world <- ne_countries(scale = "medium", returnclass = "sf")
meridian <- -180
wld.new <- st_break_antimeridian(world, lon_0 = meridian)

wld.rob.sf <- st_transform(wld.new, paste("+proj=robin +lon_0=", meridian, "+k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# Function to plot Argo profiles for each season
plot_argo_profiles <- function(df_sf, title, output_file) {
  ggplot() +
    geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
    geom_sf(data = df_sf, aes(geometry = geometry), color = "blue", alpha = 0.1, size = 0.1) +
    coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "gray90")) -> p
  
  # Save the plot
  ggsave(filename = output_file, plot = p, width = 10, height = 6)
}

# Generate maps for each season
plot_argo_profiles(df_argo_djf_sf, "DJF Argo Profiles", "figures/TimeSpaceVar/map_djf_argo_prof.jpg")
plot_argo_profiles(df_argo_mam_sf, "MAM Argo Profiles", "figures/TimeSpaceVar/map_mam_argo_prof.jpg")
plot_argo_profiles(df_argo_jja_sf, "JJA Argo Profiles", "figures/TimeSpaceVar/map_jja_argo_prof.jpg")
plot_argo_profiles(df_argo_son_sf, "SON Argo Profiles", "figures/TimeSpaceVar/map_son_argo_prof.jpg")
# Function for KDE analysis
kde_analysis <- function(df, n_grid) {
  MASS::kde2d(df$LONGITUDE, df$LATITUDE, n = n_grid)
}

# Run KDE for each season
kde_djf <- kde_analysis(df_argo_djf, 180)
kde_mam <- kde_analysis(df_argo_mam, 180)
kde_jja <- kde_analysis(df_argo_jja, 180)
kde_son <- kde_analysis(df_argo_son, 180)

# Convert KDE to data frame for plotting
kde_to_dataframe <- function(kde_result) {
  expand.grid(x = kde_result$x, y = kde_result$y) %>%
    mutate(z = as.vector(kde_result$z))
}

# Convert KDE results to data frames
kde_djf_df <- kde_to_dataframe(kde_djf)
kde_mam_df <- kde_to_dataframe(kde_mam)
kde_jja_df <- kde_to_dataframe(kde_jja)
kde_son_df <- kde_to_dataframe(kde_son)


# Function to plot KDE results
plot_kde <- function(kde_df, title, output_file) {
  ggplot() +
    geom_tile(data = kde_df, aes(x = x, y = y, fill = z)) +
    geom_contour(data = kde_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
    coord_fixed() +
    geom_sf(data = world, fill = "gray80", color = "black") +
    scale_fill_viridis_c() +
    labs(title = title, x = "Longitude", y = "Latitude", fill = "Density") +
    theme_minimal() -> p
  
  # Save the plot
  ggsave(output_file, plot = p, width = 10, height = 6)
}

# Plot KDE maps for each season
plot_kde(kde_djf_df, "Density of DJF Argo Profiles", "figures/TimeSpaceVar/kde_djf_argo_profiles.png")
plot_kde(kde_mam_df, "Density of MAM Argo Profiles", "figures/TimeSpaceVar/kde_mam_argo_profiles.png")
plot_kde(kde_jja_df, "Density of JJA Argo Profiles", "figures/TimeSpaceVar/kde_jja_argo_profiles.png")
plot_kde(kde_son_df, "Density of SON Argo Profiles", "figures/TimeSpaceVar/kde_son_argo_profiles.png")

# Bin data
# Define bin size for longitude and latitude
bin_size <- 5
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE)), ceiling(max(df_argo_clean$LONGITUDE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE)), ceiling(max(df_argo_clean$LATITUDE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins to Argo profiles and anomalies for each season
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

# Assign bins to Argo profiles and anomalies
df_argo_djf <- assign_bins(df_argo_djf)
df_argo_mam <- assign_bins(df_argo_mam)
df_argo_jja <- assign_bins(df_argo_jja)
df_argo_son <- assign_bins(df_argo_son)

df_complete_djf <- assign_bins(df_complete_djf)
df_complete_mam <- assign_bins(df_complete_mam)
df_complete_jja <- assign_bins(df_complete_jja)
df_complete_son <- assign_bins(df_complete_son)


df_carbon_djf <- assign_bins(df_carbon_djf)
df_carbon_mam <- assign_bins(df_carbon_mam)
df_carbon_jja <- assign_bins(df_carbon_jja)
df_carbon_son <- assign_bins(df_carbon_son)


# Compute counts of total profiles and anomalies for each season
compute_counts <- function(df_argo, df_complete) {
  total_counts <- df_argo %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_total = n(), .groups = 'drop')
  
  anomaly_counts <- df_complete %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_anomaly = n(), .groups = 'drop')
  
  # Merge counts and compute proportions
  merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin")) %>%
    mutate(
      count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
      proportion = count_anomaly / count_total
    ) %>%
    filter(count_total > 0) %>%
    filter(proportion < 0.8) %>%  # Filter out high proportions (anomalies > 80%)
    filter(!is.na(lat_bin))
  
  return(merged_counts)
}

# Compute merged_counts of SUBDUCTION for each season
merged_counts_djf <- compute_counts(df_argo_djf, df_complete_djf)
merged_counts_mam <- compute_counts(df_argo_mam, df_complete_mam)
merged_counts_jja <- compute_counts(df_argo_jja, df_complete_jja)
merged_counts_son <- compute_counts(df_argo_son, df_complete_son)

# Compute merged_counts of CARBON SUBDUCTION for each season
merged_carbon_counts_djf <- compute_counts(df_argo_djf, df_complete_djf)
merged_carbon_counts_mam <- compute_counts(df_argo_mam, df_complete_mam)
merged_carbon_counts_jja <- compute_counts(df_argo_jja, df_complete_jja)
merged_carbon_counts_son <- compute_counts(df_argo_son, df_complete_son)


# Create proportion maps of SUBDUCTION (w/o carbon) without saving them to files directly
prop_map_djf <- ggplot(merged_counts_djf, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (DJF)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_mam <- ggplot(merged_counts_mam, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_mam$lon_bin), ylim = range(merged_counts_mam$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (MAM)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_jja <- ggplot(merged_counts_jja, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_jja$lon_bin), ylim = range(merged_counts_jja$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (JJA)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_son <- ggplot(merged_counts_son, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_son$lon_bin), ylim = range(merged_counts_son$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (SON)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

# Combine proportion maps into a single figure using patchwork
combined_proportion_maps <- prop_map_djf + prop_map_mam + prop_map_jja + prop_map_son + 
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(title = "Proportion of Subduction Events Across Seasons")

# Save combined figure
ggsave(filename = "figures/TimeSpaceVar/4SEASONS/combined_proportion_maps.png", plot = combined_proportion_maps, width = 15, height = 12)

# Create proportion maps of CARBON SUBDUCTION without saving them to files directly
prop_map_djf_carbon <- ggplot(merged_carbon_counts_djf, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (DJF)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_mam_carbon <- ggplot(merged_carbon_counts_mam, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_mam$lon_bin), ylim = range(merged_counts_mam$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (MAM)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_jja_carbon <- ggplot(merged_carbon_counts_jja, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_jja$lon_bin), ylim = range(merged_counts_jja$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (JJA)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_son_carbon <- ggplot(merged_carbon_counts_son, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_son$lon_bin), ylim = range(merged_counts_son$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (SON)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

# Combine proportion maps into a single figure using patchwork
combined_proportion_maps_carbon <- prop_map_djf_carbon + prop_map_mam_carbon +
  prop_map_jja_carbon + prop_map_son_carbon + 
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(title = "Proportion of Carbon Subduction Events Across Seasons")

# Save combined figure
ggsave(filename = "figures/TimeSpaceVar/4SEASONS/combined_proportion_maps_carbon.png", plot = combined_proportion_maps_carbon, width = 15, height = 12)

###############################################################################
# 1) GAM Fitting Function
###############################################################################
# This function expects a data frame (merged_counts) with:
#  - lon_bin, lat_bin (numeric)
#  - count_total
#  - count_anomaly
# and fits a GAM:   cbind(anomaly, total - anomaly) ~ s(lat, lon).
# 
# We specify k=300 or some other suitable value. Feel free to adjust.

fit_gam_season <- function(merged_counts, k_value = 300) {
  # Filter out any invalid rows if necessary
  df <- merged_counts %>%
    filter(!is.na(lon_bin), !is.na(lat_bin), count_total > 0)
  
  gam_model <- gam(
    cbind(count_anomaly, count_total - count_anomaly) ~
      s(lat_bin, lon_bin, bs = "sos", k = k_value),
    family = binomial(link = "logit"),
    data = df,
    method = "REML"
  )
  return(gam_model)
}

###############################################################################
# 2) Create Prediction Grid Function
###############################################################################
# We'll sample lat_bin, lon_bin in 1° increments (adjust as desired).
# Make sure the range covers the data domain for that season.

create_prediction_grid <- function(merged_counts, step = 1) {
  lon_seq <- seq(min(merged_counts$lon_bin, na.rm = TRUE),
                 max(merged_counts$lon_bin, na.rm = TRUE),
                 by = step)
  lat_seq <- seq(min(merged_counts$lat_bin, na.rm = TRUE),
                 max(merged_counts$lat_bin, na.rm = TRUE),
                 by = step)
  expand.grid(lon_bin = lon_seq, lat_bin = lat_seq)
}

###############################################################################
# 3) Predict Smoothed Proportions
###############################################################################

predict_gam <- function(gam_model, merged_counts, step = 1) {
  pred_grid <- create_prediction_grid(merged_counts, step)
  pred_grid$proportion <- predict(gam_model, newdata = pred_grid, type = "response")
  pred_grid <- pred_grid %>% filter(!is.na(proportion))
  return(pred_grid)
}

###############################################################################
# 4) Plot a GAM Result (Lat-Lon) With Discrete Color Bins
###############################################################################
# We'll define the scale in a separate function so we can re-use it across plots.

make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.1) {
  # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
  breaks_vec <- seq(prob_min, prob_max, by = binwidth)
  
  scale_fill_viridis_b(
    name = "Probability",
    breaks = breaks_vec,
    limits = c(prob_min, prob_max),
    oob = scales::squish
  )
}

plot_gam_map <- function(pred_grid, world_data, season_label, event_label, common_scale) {
  x_limits <- range(pred_grid$lon_bin, na.rm = TRUE)
  y_limits <- range(pred_grid$lat_bin, na.rm = TRUE)
  
  ggplot() +
    geom_tile(data = pred_grid,
              aes(x = lon_bin, y = lat_bin, fill = proportion)) +
    geom_contour(data = pred_grid,
                 aes(x = lon_bin, y = lat_bin, z = proportion),
                 color = "white", alpha = 0.3) +
    geom_sf(data = world_data, fill = "gray80", color = "gray40", inherit.aes = FALSE) +
    coord_sf(xlim = x_limits, ylim = y_limits, expand = FALSE, crs = st_crs(4326)) +
    labs(
      title = paste0("Estimated ", event_label, " Probability (", season_label, ")"),
      x = "Longitude", y = "Latitude"
    ) +
    common_scale +    # apply the discrete color scale
    theme_minimal()
}

###############################################################################
# 5) Build and Plot Subduction GAM for Each Season
###############################################################################

# Fit GAMs
gam_djf_subd <- fit_gam_season(merged_counts_djf,  k_value = 300)
gam_mam_subd <- fit_gam_season(merged_counts_mam,  k_value = 300)
gam_jja_subd <- fit_gam_season(merged_counts_jja,  k_value = 300)
gam_son_subd <- fit_gam_season(merged_counts_son,  k_value = 300)

# Predict
pred_djf_subd <- predict_gam(gam_djf_subd, merged_counts_djf, step = 1)
pred_mam_subd <- predict_gam(gam_mam_subd, merged_counts_mam, step = 1)
pred_jja_subd <- predict_gam(gam_jja_subd, merged_counts_jja, step = 1)
pred_son_subd <- predict_gam(gam_son_subd, merged_counts_son, step = 1)

# Determine global min/max for subduction proportions
subd_values <- c(pred_djf_subd$proportion, pred_mam_subd$proportion,
                 pred_jja_subd$proportion, pred_son_subd$proportion)
subd_min <- 0  # typically 0
subd_max <- max(subd_values, na.rm = TRUE)

# Build a discrete color scale common to all subduction maps
# e.g. 0.1 step. Adjust if desired.
subd_scale <- make_discrete_scale(subd_min, subd_max, binwidth = 0.05)

# Plot each season with the common discrete scale
map_subd_djf <- plot_gam_map(pred_djf_subd, world, "DJF", "Subduction", subd_scale)
map_subd_mam <- plot_gam_map(pred_mam_subd, world, "MAM", "Subduction", subd_scale)
map_subd_jja <- plot_gam_map(pred_jja_subd, world, "JJA", "Subduction", subd_scale)
map_subd_son <- plot_gam_map(pred_son_subd, world, "SON", "Subduction", subd_scale)

combined_subd <- (map_subd_djf + map_subd_mam) / (map_subd_jja + map_subd_son) +
  plot_annotation(title = "GAM-Estimated Probability of Subduction Events (Discrete Scale)")
combined_subd

ggsave("figures/TimeSpaceVar/4SEASONS/gam_subduction_discrete_4seasons.png",
        combined_subd, width = 14, height = 10)

###############################################################################
# 6) Build and Plot Carbon Subduction GAM for Each Season
###############################################################################

gam_djf_carb <- fit_gam_season(merged_carbon_counts_djf,  k_value = 300)
gam_mam_carb <- fit_gam_season(merged_carbon_counts_mam,  k_value = 300)
gam_jja_carb <- fit_gam_season(merged_carbon_counts_jja,  k_value = 300)
gam_son_carb <- fit_gam_season(merged_carbon_counts_son,  k_value = 300)

# Predict
pred_djf_carb <- predict_gam(gam_djf_carb, merged_carbon_counts_djf, step = 1)
pred_mam_carb <- predict_gam(gam_mam_carb, merged_carbon_counts_mam, step = 1)
pred_jja_carb <- predict_gam(gam_jja_carb, merged_carbon_counts_jja, step = 1)
pred_son_carb <- predict_gam(gam_son_carb, merged_carbon_counts_son, step = 1)

# Global min/max for carbon subduction
carb_values <- c(pred_djf_carb$proportion, pred_mam_carb$proportion,
                 pred_jja_carb$proportion, pred_son_carb$proportion)
carb_min <- 0
carb_max <- max(carb_values, na.rm = TRUE)

# Discrete color scale for carbon subduction
carb_scale <- make_discrete_scale(carb_min, carb_max, binwidth = 0.05)

map_carb_djf <- plot_gam_map(pred_djf_carb, world, "DJF", "Carbon Subduction", carb_scale)
map_carb_mam <- plot_gam_map(pred_mam_carb, world, "MAM", "Carbon Subduction", carb_scale)
map_carb_jja <- plot_gam_map(pred_jja_carb, world, "JJA", "Carbon Subduction", carb_scale)
map_carb_son <- plot_gam_map(pred_son_carb, world, "SON", "Carbon Subduction", carb_scale)

combined_carb <- (map_carb_djf + map_carb_mam) / (map_carb_jja + map_carb_son) +
  plot_annotation(title = "GAM-Estimated Probability of Carbon Subduction Events (Discrete Scale)")

 ggsave("figures/TimeSpaceVar/4SEASONS/gam_carbon_subduction_discrete_4seasons.png",
        combined_carb, width = 14, height = 10)

###############################################################################
# Done!
###############################################################################

# Summary:
# - We create 8 total plots: 4 for subduction, 4 for carbon subduction, 
#   each with a discrete color scale and a uniform range per dataset.
# - The scale is set by the global min (0) and max proportion among 
#   the four seasons, so they're all directly comparable. 
# - Adjust 'binwidth' in 'make_discrete_scale()' if you want finer or coarser bins.

###############################################################################
# Done!
###############################################################################

# - You now have 2 sets of maps (subduction vs. carbon subduction), each with 
#   4 seasons. 
# - Adjust 'k_value', 'step' in 'create_prediction_grid()', and color scales 
#   as needed for best results.
# - If you want a consistent color scale across all seasons, you can find the 
#   max proportion across all prediction grids and use e.g. 
#   scale_fill_viridis_c(limits = c(0, global_max)).
 
 
# Now let's do histograms faceted by region showing the probability of (carbon) subduction
# by months :
 
 ###############################################################################
 # Seasonal Cycle of Subduction (Carbon and Non-Carbon) by Region and Month
 ###############################################################################
 
 # Add Month and Region Information to Data
 add_month_region <- function(df) {
   df %>%
     mutate(
       month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
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
 
 df_argo_clean <- add_month_region(df_argo_clean)
 df_complete_clean <- add_month_region(df_complete_clean)
 df_carbon_clean <- add_month_region(df_carbon_clean)
 
 # Step 2: Compute Probabilities by Region and Month
 compute_monthly_probabilities <- function(df_argo, df_events) {
   # Total Argo profiles by region and month
   total_counts <- df_argo %>%
     group_by(region, month) %>%
     summarize(count_total = n(), .groups = "drop")
   
   # Subduction events by region and month
   event_counts <- df_events %>%
     group_by(region, month) %>%
     summarize(count_event = n(), .groups = "drop")
   
   # Merge and compute proportions
   merged <- full_join(total_counts, event_counts, by = c("region", "month")) %>%
     mutate(
       count_event = replace_na(count_event, 0),
       proportion = ifelse(count_total > 0, count_event / count_total, NA)
     )
   return(merged)
 }
 
 # Compute probabilities for subduction and carbon subduction
 monthly_probs_subduction <- compute_monthly_probabilities(df_argo_clean, df_complete_clean)
 monthly_probs_carbon <- compute_monthly_probabilities(df_argo_clean, df_carbon_clean)
 
 # Step 3: Plot Monthly Histograms Faceted by Region
 plot_monthly_histograms <- function(monthly_probs, title) {
   monthly_probs <- monthly_probs %>%
     filter(!is.na(region), !is.na(proportion), !is.na(month))  # Remove NA regions, proportions, and months
   
   ggplot(monthly_probs, aes(x = month, y = proportion, fill = region)) +
     geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
     facet_wrap(~region, ncol = 2,axes="all") +
     scale_fill_viridis_d() +
     labs(
       title = title,
       x = "Month",
       y = "Probability of Subduction",
       fill = "Region"
     ) +
     theme_minimal() +
     theme(
       strip.text = element_text(size = 12, face = "bold"),
       axis.text.x = element_text(angle = 45, hjust = 1)
     )
 }
 
 
# Plot subduction probabilities
 plot_subduction <- plot_monthly_histograms(monthly_probs_subduction, "Probability that an Argo Float captures a subduction event")
 
 # Plot carbon subduction probabilities
 plot_carbon_subduction <- plot_monthly_histograms(monthly_probs_carbon, "Probability that an Argo Float captures a subduction event with carbon")
 
 # Display the plots
 plot_subduction
 plot_carbon_subduction
 
 # Optional: Save the plots
 ggsave("figures/seasonal_subduction_by_region.png", plot = plot_subduction, width = 10, height = 8)
 ggsave("figures/seasonal_carbon_subduction_by_region.png", plot = plot_carbon_subduction, width = 10, height = 8)
 
 
