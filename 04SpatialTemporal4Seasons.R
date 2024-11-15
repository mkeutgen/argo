# Load necessary libraries
library(conflicted)
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
library(segmented)
library(dplyr)
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

# first investigate structure of data
df_complete_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
df_complete_clean %>% head()
df_argo_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_argo_loc.csv")



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
plot_argo_profiles(df_argo_djf_sf, "DJF Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/map_djf_argo_prof.jpg")
plot_argo_profiles(df_argo_mam_sf, "MAM Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/map_mam_argo_prof.jpg")
plot_argo_profiles(df_argo_jja_sf, "JJA Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/map_jja_argo_prof.jpg")
plot_argo_profiles(df_argo_son_sf, "SON Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/map_son_argo_prof.jpg")

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
plot_kde(kde_djf_df, "Density of DJF Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/kde_djf_argo_profiles.png")
plot_kde(kde_mam_df, "Density of MAM Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/kde_mam_argo_profiles.png")
plot_kde(kde_jja_df, "Density of JJA Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/kde_jja_argo_profiles.png")
plot_kde(kde_son_df, "Density of SON Argo Profiles", "/data/GLOBARGO/figures/TimeSpaceVar/kde_son_argo_profiles.png")


# Bin data
# Define bin size for longitude and latitude
bin_size <- 10
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

# Compute merged_counts for each season
merged_counts_djf <- compute_counts(df_argo_djf, df_complete_djf)
merged_counts_mam <- compute_counts(df_argo_mam, df_complete_mam)
merged_counts_jja <- compute_counts(df_argo_jja, df_complete_jja)
merged_counts_son <- compute_counts(df_argo_son, df_complete_son)



# Create proportion maps without saving them to files directly
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
ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/4SEASONS/combined_proportion_maps_res10.png", plot = combined_proportion_maps, width = 15, height = 12)


threshold <- 50  # Example threshold
filtered_counts_djf <- merged_counts_djf %>%
  filter(count_total >= threshold)  # Filter out low count cells

filtered_counts_mam <- merged_counts_mam %>%
  filter(count_total >= threshold)  # Filter out low count cells

filtered_counts_jja <- merged_counts_jja %>%
  filter(count_total >= threshold)  # Filter out low count cells

filtered_counts_son <- merged_counts_son %>%
  filter(count_total >= threshold)  # Filter out low count cells





# Define a binned color scale for the proportion maps
binned_color_scale <- scale_fill_viridis_b(
  name = "Proportion",
  breaks = seq(0, max(c(merged_counts_djf$proportion, merged_counts_mam$proportion, 
                        merged_counts_jja$proportion, merged_counts_son$proportion)), by = 0.05), 
  limits = c(0, max(c(merged_counts_djf$proportion, merged_counts_mam$proportion, 
                      merged_counts_jja$proportion, merged_counts_son$proportion))),
  oob = scales::squish # Squishes out-of-bound values into the closest bin
)

# Create proportion maps with the binned scale
prop_map_djf <- ggplot(merged_counts_djf, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  binned_color_scale +
  labs(title = "Proportion of Subduction Events (DJF)", x = "Longitude", y = "Latitude") +
  theme_minimal()

prop_map_mam <- ggplot(merged_counts_mam, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_mam$lon_bin), ylim = range(merged_counts_mam$lat_bin)) +
  binned_color_scale +
  labs(title = "Proportion of Subduction Events (MAM)", x = "Longitude", y = "Latitude") +
  theme_minimal()

prop_map_jja <- ggplot(merged_counts_jja, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_jja$lon_bin), ylim = range(merged_counts_jja$lat_bin)) +
  binned_color_scale +
  labs(title = "Proportion of Subduction Events (JJA)", x = "Longitude", y = "Latitude") +
  theme_minimal()

prop_map_son <- ggplot(merged_counts_son, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_son$lon_bin), ylim = range(merged_counts_son$lat_bin)) +
  binned_color_scale +
  labs(title = "Proportion of Subduction Events (SON)", x = "Longitude", y = "Latitude") +
  theme_minimal()

# Combine proportion maps with the discontinuous scale
combined_proportion_maps <- prop_map_djf + prop_map_mam + prop_map_jja + prop_map_son + 
  plot_layout(ncol = 2, nrow = 2, guides = "collect") + 
  plot_annotation(title = "Proportion of Subduction Events Across Seasons (Discontinuous Scale)")

# Save combined figure
ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/4SEASONS/combined_proportion_maps_discontinuous_k1000_res10.png", 
       plot = combined_proportion_maps, width = 15, height = 12)


# Fit the GAM models for each season
gam_djf <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lon_bin,lat_bin , bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = filtered_counts_djf,
  weights = log(count_total),
  method = "REML"
)

# Spatial autocorrelation investigation
# Extract the residuals from your fitted GAM model
residuals_djf <- residuals(gam_djf, type = "pearson")

# Add the residuals to your original data
filtered_counts_djf$residuals <- residuals_djf

# Quick summary of the residuals
summary(filtered_counts_djf$residuals)
# Create a spatial weights matrix using k-nearest neighbors
coords <- filtered_counts_djf %>% select(lon_bin, lat_bin) %>% as.matrix()
nb <- knearneigh(coords, k = 4)
listw <- nb2listw(knn2nb(nb), style = "W")

# Calculate Moran's I for the residuals
moran_test <- moran.test(filtered_counts_djf$residuals, listw)
# No spatial autocorellation

gam_djf %>% summary() # R-sq.(adj) =  0.837   Deviance explained = 96.5%

gam_mam <- gam(
  formula = cbind(count_anomaly, count_total - count_anomaly) ~ s(lon_bin, lat_bin, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = filtered_counts_mam,
  weights = log(count_total),
  method = "REML"
)

gam_mam %>% summary()


gam_jja <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = filtered_counts_jja,
  weights = log(count_total),
  method = "REML"
)

gam_son <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = filtered_counts_son,
  weights = log(count_total),
  method = "REML"
)

# Create prediction grids for each season
create_prediction_grid <- function(merged_counts) {
  lon_seq <- seq(min(merged_counts$lon_bin), max(merged_counts$lon_bin), by = 1)
  lat_seq <- seq(min(merged_counts$lat_bin), max(merged_counts$lat_bin), by = 1)
  expand.grid(lon_bin = lon_seq, lat_bin = lat_seq)
}

# Prediction grids
prediction_grid_djf <- create_prediction_grid(merged_counts_djf)
prediction_grid_mam <- create_prediction_grid(merged_counts_mam)
prediction_grid_jja <- create_prediction_grid(merged_counts_jja)
prediction_grid_son <- create_prediction_grid(merged_counts_son)

# Predict the smoothed proportions for each season
prediction_grid_djf$proportion <- predict(gam_djf, newdata = prediction_grid_djf, type = "response")
prediction_grid_mam$proportion <- predict(gam_mam, newdata = prediction_grid_mam, type = "response")
prediction_grid_jja$proportion <- predict(gam_jja, newdata = prediction_grid_jja, type = "response")
prediction_grid_son$proportion <- predict(gam_son, newdata = prediction_grid_son, type = "response")

# Remove NA values
prediction_grid_djf <- prediction_grid_djf %>% filter(!is.na(proportion))
prediction_grid_mam <- prediction_grid_mam %>% filter(!is.na(proportion))
prediction_grid_jja <- prediction_grid_jja %>% filter(!is.na(proportion))
prediction_grid_son <- prediction_grid_son %>% filter(!is.na(proportion))


# Define a binned color scale for the GAM maps
binned_color_scale_gam <- scale_fill_viridis_b(
  name = "Proportion",
  breaks = seq(0, max(c(prediction_grid_djf$proportion, prediction_grid_mam$proportion, 
                        prediction_grid_jja$proportion, prediction_grid_son$proportion)), by = 0.05),
  limits = c(0, max(c(prediction_grid_djf$proportion, prediction_grid_mam$proportion, 
                      prediction_grid_jja$proportion, prediction_grid_son$proportion))),
  oob = scales::squish
)

# Create GAM maps using the binned color scale
gam_map_djf <- ggplot() +
  geom_tile(data = prediction_grid_djf, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid_djf, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "gray80") +
  coord_sf(xlim = range(prediction_grid_djf$lon_bin), ylim = range(prediction_grid_djf$lat_bin), expand = FALSE) +
  binned_color_scale_gam +
  labs(title = "Estimated Proportion of Anomalous Argo Profiles (DJF)", x = "Longitude", y = "Latitude") +
  theme_minimal()

gam_map_mam <- ggplot() +
  geom_tile(data = prediction_grid_mam, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid_mam, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "gray80") +
  coord_sf(xlim = range(prediction_grid_mam$lon_bin), ylim = range(prediction_grid_mam$lat_bin), expand = FALSE) +
  binned_color_scale_gam +
  labs(title = "Estimated Proportion of Anomalous Argo Profiles (MAM)", x = "Longitude", y = "Latitude") +
  theme_minimal()

gam_map_jja <- ggplot() +
  geom_tile(data = prediction_grid_jja, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid_jja, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "gray80") +
  coord_sf(xlim = range(prediction_grid_jja$lon_bin), ylim = range(prediction_grid_jja$lat_bin), expand = FALSE) +
  binned_color_scale_gam +
  labs(title = "Estimated Proportion of Anomalous Argo Profiles (JJA)", x = "Longitude", y = "Latitude") +
  theme_minimal()

gam_map_son <- ggplot() +
  geom_tile(data = prediction_grid_son, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid_son, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "gray80") +
  coord_sf(xlim = range(prediction_grid_son$lon_bin), ylim = range(prediction_grid_son$lat_bin), expand = FALSE) +
  binned_color_scale_gam +
  labs(title = "Estimated Proportion of Anomalous Argo Profiles (SON)", x = "Longitude", y = "Latitude") +
  theme_minimal()

# Combine GAM maps with the discontinuous scale
combined_gam_maps <- gam_map_djf + gam_map_mam + gam_map_jja + gam_map_son + 
  plot_layout(ncol = 2, nrow = 2, guides = "collect") + 
  plot_annotation(title = "GAM Estimated Proportion of Anomalous Argo Profiles Across Seasons (dx = 10, k = 150, )")

# Save combined figure
ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/4SEASONS/combined_gam_maps_discontinuous_k1000_res10_weighted.png", 
       plot = combined_gam_maps, width = 15, height = 12)




