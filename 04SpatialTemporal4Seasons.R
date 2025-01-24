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
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

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
merged_carbon_counts_djf <- compute_counts(df_argo_djf, df_carbon_djf)
merged_carbon_counts_mam <- compute_counts(df_argo_mam, df_carbon_mam)
merged_carbon_counts_jja <- compute_counts(df_argo_jja, df_carbon_jja)
merged_carbon_counts_son <- compute_counts(df_argo_son, df_carbon_son)


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

###########################################################
### 7. Add stippling DJF with the SAME discrete subd_scale
############################################################

# We assume subd_scale is already defined, e.g.:
# subd_scale <- make_discrete_scale(subd_min, subd_max, binwidth = 0.05)

gam_resolution <- 1       # Resolution for GAM predictions
stipple_resolution <- 5   # Coarser resolution for stippling

# 1) Aggregate Argo float locations to 5° resolution for stippling
argo_bins <- df_argo_clean %>%
  filter(month(TIME) %in% c(12, 1, 2)) %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")  # Count profiles per 5° bin

# 2) Identify undersampled areas
undersampled <- argo_bins %>%
  filter(count < 5)  # or any other threshold you prefer

# 3) Grab DJF prediction grid from the existing subduction GAM results
#    pred_djf_subd is the data frame with columns: lon_bin, lat_bin, proportion
prediction_grid <- pred_djf_subd  

# Keep track of the coarser bins for stippling
prediction_grid <- prediction_grid %>%
  mutate(
    lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
  )

# 4) Merge GAM predictions with undersampled areas
prediction_grid <- prediction_grid %>%
  left_join(undersampled,
            by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(
    undersampled = ifelse(is.na(count), FALSE, TRUE)
  )

# 5) Plot with the shared 'subd_scale' instead of scale_fill_viridis_c()
gam_map_with_stippling_djf <- ggplot() +
  # Tile layer of subduction probability at 1° resolution
  geom_tile(data = prediction_grid,
            aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid,
               aes(x = lon_bin, y = lat_bin, z = proportion),
               color = "white", alpha = 0.3) +
  # Add stippling where undersampled == TRUE
  geom_point(
    data = filter(prediction_grid, undersampled == TRUE),
    aes(x = lon_bin_stipple, y = lat_bin_stipple),
    color = "white", alpha = 0.6, size = 0.25, shape = 20
  ) +
  # World map for context
  geom_sf(data = world, fill = "black", color = "black", inherit.aes = FALSE) +
  coord_sf(expand = FALSE) +
  labs(
    title = "GAM-Estimated Probability of Subduction Events (DJF + Stippling)",
    x = "Longitude",
    y = "Latitude",
    fill = "Probability"
  ) +
  # Here is where we reuse the COMMON discrete scale
  subd_scale +
  theme_minimal()

# Display or save the final figure
gam_map_with_stippling_djf
# ggsave("figures/TimeSpaceVar/4SEASONS/gam_subduction_djf_stipple_sameScale.png",
#        gam_map_with_stippling_djf, width = 10, height = 7)



###########################################################
### 7. Add stippling DJF ########################################
#############################################################
# Define resolutions
gam_resolution <- 1   # Resolution for GAM predictions
stipple_resolution <- 5  # Coarser resolution for undersampling shading

# Step 1: Aggregate Argo float locations to 10° resolution
argo_bins <- df_argo_clean %>%
  filter(month(TIME) %in% c(12, 1, 2))  %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")  # Count profiles per 10° bin

# Step 2: Identify undersampled areas
undersampled <- argo_bins %>%
  filter(count < 5)  # Define threshold for undersampling

# Step 3: Prepare GAM prediction data at 1° resolution
# Assuming `prediction_grid` is your GAM prediction output (from predict_gam)
prediction_grid <- pred_djf_subd  # Example for DJF GAM predictions

# Keep GAM grid at 1° resolution
prediction_grid <- prediction_grid %>%
  mutate(
    lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
  )

# Step 4: Merge GAM predictions with undersampled areas
# Add a TRUE/FALSE column indicating undersampled cells at 10° resolution
prediction_grid <- prediction_grid %>%
  left_join(undersampled, by = c("lon_bin_stipple" = "lon_bin_stipple",
                                 "lat_bin_stipple" = "lat_bin_stipple")) %>%
  mutate(
    undersampled = ifelse(is.na(count), FALSE, TRUE)  # Mark as TRUE if undersampled
  )

# Step 5: Plot GAM map with stippling for undersampled areas
gam_map_with_stippling_djf <- ggplot() +
  # GAM predictions as a tile layer at 1° resolution
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid, aes(x = lon_bin, y = lat_bin, z = proportion),
               color = "white", alpha = 0.3) +
  # Add stippling for unwhitedersampled areas at 10° resolution
  geom_point(data = filter(prediction_grid, undersampled == TRUE),
             aes(x = lon_bin_stipple, y = lat_bin_stipple),
             color = "white", alpha = 0.6, size = 0.25, shape = 20) +
  # World map for context
  geom_sf(data = world, fill = "black", color = "black", inherit.aes = FALSE) +
  coord_sf(expand = FALSE) +
  labs(
    title = "GAM-Estimated Probability of Subduction Events (DJF)",
    x = "Longitude",
    y = "Latitude",
    fill = "Probability"
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

# Display the plot
print(gam_map_with_stippling_djf)

#####
##  stippling MAM
#####
# Define resolutions
gam_resolution <- 1   # Resolution for GAM predictions
stipple_resolution <- 5  # Coarser resolution for undersampling shading

# Step 1: Aggregate Argo float locations to 10° resolution
argo_bins <- df_argo_clean %>%
  filter(month(TIME) %in% c(3, 4, 5))  %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")  # Count profiles per 10° bin

# Step 2: Identify undersampled areas
undersampled <- argo_bins %>%
  filter(count < 5)  # Define threshold for undersampling

# Step 3: Prepare GAM prediction data at 1° resolution
# Assuming `prediction_grid` is your GAM prediction output (from predict_gam)
prediction_grid <- pred_mam_subd  # Example for DJF GAM predictions

# Keep GAM grid at 1° resolution
prediction_grid <- prediction_grid %>%
  mutate(
    lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
  )

# Step 4: Merge GAM predictions with undersampled areas
# Add a TRUE/FALSE column indicating undersampled cells at 10° resolution
prediction_grid <- prediction_grid %>%
  left_join(undersampled, by = c("lon_bin_stipple" = "lon_bin_stipple",
                                 "lat_bin_stipple" = "lat_bin_stipple")) %>%
  mutate(
    undersampled = ifelse(is.na(count), FALSE, TRUE)  # Mark as TRUE if undersampled
  )

# Step 5: Plot GAM map with stippling for undersampled areas
gam_map_with_stippling_mam <- ggplot() +
  # GAM predictions as a tile layer at 1° resolution
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid, aes(x = lon_bin, y = lat_bin, z = proportion),
               color = "white", alpha = 0.3) +
  # Add stippling for unwhitedersampled areas at 10° resolution
  geom_point(data = filter(prediction_grid, undersampled == TRUE),
             aes(x = lon_bin_stipple, y = lat_bin_stipple),
             color = "white", alpha = 0.6, size = 0.25, shape = 20) +
  # World map for context
  geom_sf(data = world, fill = "black", color = "black", inherit.aes = FALSE) +
  coord_sf(expand = FALSE) +
  labs(
    title = "GAM-Estimated Probability of Subduction Events (MAM)",
    x = "Longitude",
    y = "Latitude",
    fill = "Probability"
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

# Display the plot
print(gam_map_with_stippling_mam)


# Define resolutions
gam_resolution <- 1   # Resolution for GAM predictions
stipple_resolution <- 5  # Coarser resolution for undersampling shading

# Step 1: Aggregate Argo float locations to 10° resolution
argo_bins <- df_argo_clean %>%
  filter(month(TIME) %in% c(6, 7, 8))  %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")  # Count profiles per 10° bin

# Step 2: Identify undersampled areas
undersampled <- argo_bins %>%
  filter(count < 5)  # Define threshold for undersampling

# Step 3: Prepare GAM prediction data at 1° resolution
# Assuming `prediction_grid` is your GAM prediction output (from predict_gam)
prediction_grid <- pred_jja_subd  # Example for DJF GAM predictions

# Keep GAM grid at 1° resolution
prediction_grid <- prediction_grid %>%
  mutate(
    lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
  )

# Step 4: Merge GAM predictions with undersampled areas
# Add a TRUE/FALSE column indicating undersampled cells at 10° resolution
prediction_grid <- prediction_grid %>%
  left_join(undersampled, by = c("lon_bin_stipple" = "lon_bin_stipple",
                                 "lat_bin_stipple" = "lat_bin_stipple")) %>%
  mutate(
    undersampled = ifelse(is.na(count), FALSE, TRUE)  # Mark as TRUE if undersampled
  )

# Step 5: Plot GAM map with stippling for undersampled areas
gam_map_with_stippling_jja <- ggplot() +
  # GAM predictions as a tile layer at 1° resolution
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid, aes(x = lon_bin, y = lat_bin, z = proportion),
               color = "white", alpha = 0.3) +
  # Add stippling for unwhitedersampled areas at 10° resolution
  geom_point(data = filter(prediction_grid, undersampled == TRUE),
             aes(x = lon_bin_stipple, y = lat_bin_stipple),
             color = "white", alpha = 0.6, size = 0.25, shape = 20) +
  # World map for context
  geom_sf(data = world, fill = "black", color = "black", inherit.aes = FALSE) +
  coord_sf(expand = FALSE) +
  labs(
    title = "GAM-Estimated Probability of Subduction Events (JJA)",
    x = "Longitude",
    y = "Latitude",
    fill = "Probability"
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

# Display the plot
print(gam_map_with_stippling_jja)


# SON #

# Define resolutions
gam_resolution <- 1   # Resolution for GAM predictions
stipple_resolution <- 5  # Coarser resolution for undersampling shading

# Step 1: Aggregate Argo float locations to 10° resolution
argo_bins <- df_argo_clean %>%
  filter(month(TIME) %in% c(9, 10, 11))  %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")  # Count profiles per 10° bin

# Step 2: Identify undersampled areas
undersampled <- argo_bins %>%
  filter(count < 5)  # Define threshold for undersampling

# Step 3: Prepare GAM prediction data at 1° resolution
# Assuming `prediction_grid` is your GAM prediction output (from predict_gam)
prediction_grid <- pred_son_subd  

# Keep GAM grid at 1° resolution
prediction_grid <- prediction_grid %>%
  mutate(
    lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
  )

# Step 4: Merge GAM predictions with undersampled areas
# Add a TRUE/FALSE column indicating undersampled cells at 10° resolution
prediction_grid <- prediction_grid %>%
  left_join(undersampled, by = c("lon_bin_stipple" = "lon_bin_stipple",
                                 "lat_bin_stipple" = "lat_bin_stipple")) %>%
  mutate(
    undersampled = ifelse(is.na(count), FALSE, TRUE)  # Mark as TRUE if undersampled
  )

# Step 5: Plot GAM map with stippling for undersampled areas
gam_map_with_stippling_son <- ggplot() +
  # GAM predictions as a tile layer at 1° resolution
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid, aes(x = lon_bin, y = lat_bin, z = proportion),
               color = "white", alpha = 0.3) +
  # Add stippling for unwhitedersampled areas at 10° resolution
  geom_point(data = filter(prediction_grid, undersampled == TRUE),
             aes(x = lon_bin_stipple, y = lat_bin_stipple),
             color = "white", alpha = 0.6, size = 0.25, shape = 20) +
  # World map for context
  geom_sf(data = world, fill = "black", color = "black", inherit.aes = FALSE) +
  coord_sf(expand = FALSE) +
  labs(
    title = "GAM-Estimated Probability of Subduction Events (JJA)",
    x = "Longitude",
    y = "Latitude",
    fill = "Probability"
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

# Display the plot
print(gam_map_with_stippling_son)



###############################################################################
# Done!
###############################################################################

# Summary:
# - We create 8 total plots: 4 for subduction, 4 for carbon subduction,
#   each with a discrete color scale and a uniform range per dataset.
# - The scale is set by the global min (0) and max proportion among
#   the four seasons, so they're all directly comparable.
# - Adjust 'binwidth' in 'make_discrete_scale()' if I want finer or coarser bins.

###############################################################################
# Done!
###############################################################################

# Now let's do histograms faceted by region showing the probability of (carbon) subduction
# by months :

###############################################################################
# Seasonal Cycle of Subduction (Carbon and Non-Carbon) by Region and Month
###############################################################################

# Add Month and Region Information to Data
# Splitting Northern Tropics and Southern Tropics as a sanity check (should not be seasonal there)
# add_month_region <- function(df) {
#   df %>%
#     mutate(
#       month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
#       region = case_when(
#         LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
#         LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
#         LATITUDE >= 0 & LATITUDE <= 30 ~ "Northern Tropics",                     # 0° to 30°N
#         LATITUDE < 0 & LATITUDE >= -30 ~ "Southern Tropics",                     # 0° to 30°S
#         LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
#         TRUE ~ NA_character_  # Exclude undefined regions
#       )
#     )
# }


add_month_region <- function(df) {
  df %>%
    mutate(
      month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
      region = case_when(
        LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
        LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
        LATITUDE <= 30 & LATITUDE >= -30 ~ "Tropics",                     # 0° to 30°S
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
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 16),
      title = element_text(size=15)
    )+  
    guides(fill = "none")
}


# Plot subduction probabilities
plot_subduction <- plot_monthly_histograms(monthly_probs_subduction, "Probability that an Argo Float captures a subduction event")

# Plot carbon subduction probabilities
plot_carbon_subduction <- plot_monthly_histograms(monthly_probs_carbon, "Probability that an Argo Float captures a subduction event with carbon")


# Save the plots
ggsave("figures/seasonal_subduction_by_region.png", plot = plot_subduction, width = 10, height = 8)
ggsave("figures/seasonal_carbon_subduction_by_region.png", plot = plot_carbon_subduction, width = 10, height = 8)

### Add a statistical test to bring evidence against H0 that the monthly frequency of subduction and of carbon subduction
# sampled from a uniform distribution to test whether subduction is or not seasonal
# Remove NA months from both datasets
monthly_probs_subduction <- monthly_probs_subduction %>% filter(!is.na(month))
monthly_probs_carbon <- monthly_probs_carbon %>% filter(!is.na(month))

# Function to test uniformity for each region
test_seasonality_proportions <- function(monthly_probs) {
  monthly_probs %>%
    group_by(region) %>%
    summarize(
      chisq_stat = chisq.test(
        x = proportion * count_total,  # Observed weighted counts
        p = rep(1 / 12, 12),           # Expected uniform proportions
        rescale.p = TRUE               # Rescale expected to match total
      )$statistic,
      p_value = chisq.test(
        x = proportion * count_total,
        p = rep(1 / 12, 12),
        rescale.p = TRUE
      )$p.value,
      .groups = "drop"
    )
}


t <- monthly_probs_carbon %>% filter(region == "North Atlantic")
chisq.test(t$count_event)

# Test seasonality for subduction
seasonality_test_subduction <- test_seasonality_proportions(monthly_probs_subduction)

# Test for carbon subduction
seasonality_test_carbon <- test_seasonality_proportions(monthly_probs_carbon)

# Display results
print("Subduction Seasonality Test:")
print(seasonality_test_subduction)

print("Carbon Subduction Seasonality Test:")
print(seasonality_test_carbon)



################################################""
#################################################"

# Gather all subduction probability values across the four seasons
common_subd_values <- c(
  pred_djf_subd$proportion,
  pred_mam_subd$proportion,
  pred_jja_subd$proportion,
  pred_son_subd$proportion
)

common_subd_min <- 0   # Typically 0
common_subd_max <- max(common_subd_values, na.rm = TRUE)

# Option 1: Using a discrete binned scale
#           (If you already have a function `make_discrete_scale`, use it here.)
make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.05) {
  breaks_vec <- seq(prob_min, prob_max, by = binwidth)
  scale_fill_viridis_b(
    name = "Probability",
    breaks = breaks_vec,
    limits = c(prob_min, prob_max),
    oob = scales::squish
  )
}
subd_scale <- make_discrete_scale(common_subd_min, common_subd_max, binwidth = 0.05)

# OR Option 2: A continuous scale — e.g.:
# subd_scale <- scale_fill_viridis_c(limits = c(common_subd_min, common_subd_max),
#                                    oob = scales::squish)
#
# Either way, you'll have a single color scale object, e.g. "subd_scale".


library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

create_stippled_map <- function(df_argo_clean,
                                pred_grid,          # e.g. pred_djf_subd
                                months,             # c(12,1,2) for DJF
                                resolution = 5,     # coarser resolution for stippling
                                threshold = 5,      # <5 profiles = undersampled
                                season_label = "DJF",
                                color_scale = subd_scale,
                                suppress = "none") { # should the legend be suppressed, "none" if you want to combine maps
  
  # 1) Coarse bin Argo data for "undersampled" identification
  argo_bins <- df_argo_clean %>%
    filter(lubridate::month(TIME) %in% months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / resolution) * resolution,
      lat_bin_stipple = floor(LATITUDE / resolution) * resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  # 2) Identify undersampled
  undersampled <- argo_bins %>%
    filter(count < threshold)
  
  # 3) Prepare the prediction data for mapping
  prediction_grid <- pred_grid %>%
    mutate(
      lon_bin_stipple = floor(lon_bin / resolution) * resolution,
      lat_bin_stipple = floor(lat_bin / resolution) * resolution
    )
  
  # 4) Merge with undersampled info
  prediction_grid <- prediction_grid %>%
    left_join(undersampled,
              by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(undersampled = ifelse(is.na(count), FALSE, TRUE))
  
  # 5) Plot
  p <- ggplot() +
    geom_tile(data = prediction_grid,
              aes(x = lon_bin, y = lat_bin, fill = proportion)) +
    geom_contour(data = prediction_grid,
                 aes(x = lon_bin, y = lat_bin, z = proportion),
                 color = "white", alpha = 0.3) +
    geom_point(
      data = filter(prediction_grid, undersampled == TRUE),
      aes(x = lon_bin_stipple, y = lat_bin_stipple),
      color = "white", alpha = 0.6, size = 0.25, shape = 20
    ) +
    geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    labs(
      title = paste("(", season_label, ")", sep=""),
      x = "Longitude",
      y = "Latitude",
      fill = "Probability"
    ) +
    # Use the *common* scale
    color_scale +  
    guides(fill = suppress)+
    theme_minimal()
  
  return(p)
}


# Make sure you have your 'world' sf object
# e.g. world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")


color_legend <- ggplot() +
  geom_tile(data = data.frame(x = c(1), y = c(1), z = c(1)),
            aes(x = x, y = y, fill = z)) +
  subd_scale +
  labs(fill = "Probability") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.width = unit(4, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


# Generate stippled maps for all seasons
gam_map_with_stippling_djf <- create_stippled_map(
  df_argo_clean, pred_djf_subd, months = c(12, 1, 2),
  resolution = 5, threshold = 5, season_label = "DJF",
  color_scale = subd_scale
)

gam_map_with_stippling_djf_with_legend <- create_stippled_map(
  df_argo_clean, pred_djf_subd, months = c(12, 1, 2),
  resolution = 5, threshold = 5, season_label = "DJF",
  color_scale = subd_scale, suppress = "legend"
)+labs(title = "Probability that a profile contains a subduction event (DJF)")

ggsave(plot = gam_map_with_stippling_djf_with_legend,
       filename = "figures/gam_map_with_stippling_djf_with_legend.png")

gam_map_with_stippling_mam <- create_stippled_map(
  df_argo_clean, pred_mam_subd, months = c(3, 4, 5),
  resolution = 5, threshold = 5, season_label = "MAM",
  color_scale = subd_scale
)

gam_map_with_stippling_jja <- create_stippled_map(
  df_argo_clean, pred_jja_subd, months = c(6, 7, 8),
  resolution = 5, threshold = 5, season_label = "JJA",
  color_scale = subd_scale
)

gam_map_with_stippling_son <- create_stippled_map(
  df_argo_clean, pred_son_subd, months = c(9, 10, 11),
  resolution = 5, threshold = 5, season_label = "SON",
  color_scale = subd_scale
)

# Combine maps into a 2x2 grid
combined_stippling_maps <- (gam_map_with_stippling_djf + gam_map_with_stippling_mam) /
  (gam_map_with_stippling_jja + gam_map_with_stippling_son)

# Add the unified legend at the bottom
final_plot <- combined_stippling_maps +
  plot_layout(guides = "collect") &
  plot_annotation(title = "Subduction Probability Across Seasons") &
  theme(
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),  # Make the legend wider
    legend.key.height = unit(0.5, "cm"),  # Adjust height for clarity
    legend.text = element_text(size = 8),
    legend.text.position = "bottom", # Adjust font size
    legend.title = element_text(size = 12)  # Adjust title size
  )

# save
ggsave("figures/TimeSpaceVar/4SEASONS/gam_subduction_combined_stippled.png",
       final_plot, width = 7, height = 5)



create_stippled_carbon_map <- function(df_argo_clean,
                                       pred_grid,
                                       months,
                                       resolution = 5,
                                       threshold = 5,
                                       season_label = "DJF",
                                       color_scale = carb_scale) {
  # 1) Coarse bin Argo data (this is for identifying undersampling)
  argo_bins <- df_argo_clean %>%
    filter(lubridate::month(TIME) %in% months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / resolution) * resolution,
      lat_bin_stipple = floor(LATITUDE / resolution) * resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  # 2) Identify undersampled bins
  undersampled <- argo_bins %>%
    filter(count < threshold)
  
  # 3) Prepare the prediction grid for mapping
  #    by matching it to the same stipple resolution
  prediction_grid <- pred_grid %>%
    mutate(
      lon_bin_stipple = floor(lon_bin / resolution) * resolution,
      lat_bin_stipple = floor(lat_bin / resolution) * resolution
    ) %>%
    left_join(undersampled, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(undersampled = ifelse(is.na(count), FALSE, TRUE))
  
  # 4) Plot
  p <- ggplot() +
    # Probability as tiles
    geom_tile(data = prediction_grid,
              aes(x = lon_bin, y = lat_bin, fill = proportion)) +
    # Probability contours
    geom_contour(data = prediction_grid,
                 aes(x = lon_bin, y = lat_bin, z = proportion),
                 color = "white", alpha = 0.3) +
    # Stippling for undersampled bins
    geom_point(
      data = filter(prediction_grid, undersampled == TRUE),
      aes(x = lon_bin_stipple, y = lat_bin_stipple),
      color = "white", alpha = 0.6, size = 0.25, shape = 20
    ) +
    # World map for context
    geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    labs(
      title = paste("(", season_label, ")", sep = ""),
      x = "Longitude",
      y = "Latitude"
    ) +
    color_scale +
    theme_minimal()
  
  return(p)
}

gam_map_carb_djf_stippled <- create_stippled_carbon_map(
  df_argo_clean,
  pred_djf_carb,
  months = c(12,1,2),
  resolution = 5,
  threshold = 5,
  season_label = "DJF",
  color_scale = carb_scale
)

gam_map_carb_mam_stippled <- create_stippled_carbon_map(
  df_argo_clean,
  pred_mam_carb,
  months = c(3,4,5),
  resolution = 5,
  threshold = 5,
  season_label = "MAM",
  color_scale = carb_scale
)

gam_map_carb_jja_stippled <- create_stippled_carbon_map(
  df_argo_clean,
  pred_jja_carb,
  months = c(6,7,8),
  resolution = 5,
  threshold = 5,
  season_label = "JJA",
  color_scale = carb_scale
)

gam_map_carb_son_stippled <- create_stippled_carbon_map(
  df_argo_clean,
  pred_son_carb,
  months = c(9,10,11),
  resolution = 5,
  threshold = 5,
  season_label = "SON",
  color_scale = carb_scale
)

color_legend_carb <- ggplot(data = data.frame(x = c(1), y = c(1), z = c(0.5))) +
  geom_tile(aes(x = x, y = y, fill = z)) +
  carb_scale +
  labs(fill = "Probability") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width  = unit(2, "cm"),
    legend.key.height = unit(0.4, "cm"),
    axis.title        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    panel.grid        = element_blank()
  )

# 2×2 layout of the four seasonal maps
combined_carb_maps <- (gam_map_carb_djf_stippled + gam_map_carb_mam_stippled) /
  (gam_map_carb_jja_stippled + gam_map_carb_son_stippled)

# Add a main title (optional)
final_carb_plot <- combined_carb_maps +
  plot_layout(guides = "collect") &
  plot_annotation(title = "Carbon Subduction Probability Across Seasons") &
  theme(
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),  # Make the legend wider
    legend.key.height = unit(0.5, "cm"),  # Adjust height for clarity
    legend.text = element_text(size = 8),
    legend.text.position = "bottom", # Adjust font size
    legend.title = element_text(size = 12)  # Adjust title size
  )


ggsave("figures/TimeSpaceVar/4SEASONS/gam_carbon_subduction_stippled_4seasons.png",
       final_carb_plot, width = 7, height = 5)
