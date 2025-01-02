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

model <- readRDS("/data/GLOBARGO/src/data/gam_model_21nov_N_2_MLD_LAT_LON_Season.Rds")

summary(model)
# first investigate structure of data
df_complete_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
df_complete_clean %>% head()
df_argo_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_argo_loc.csv")
df_argo_clean %>% head()

# Define months for the two periods: DJFMAM and JJASON
winter_early_spring_months <- c(12, 1, 2, 3, 4)
late_spring_summer_months <- c(5, 6, 7, 8, 9, 10)

# Filter data for DJFMAM and JJASON
df_argo_winter_early_spring <- df_argo_clean %>%
  filter(month(TIME) %in% winter_early_spring_months)

df_argo_late_spring_summer <- df_argo_clean %>%
  filter(month(TIME) %in% late_spring_summer_months)

df_complete_winter_early_spring <- df_complete_clean %>%
  filter(month(TIME) %in% winter_early_spring_months)

df_complete_late_spring_summer <- df_complete_clean %>%
  filter(month(TIME) %in% late_spring_summer_months)

# Convert data frames to spatial objects and transform to Robinson projection
df_argo_winter_early_spring_sf <- df_argo_winter_early_spring %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_late_spring_summer_sf <- df_argo_late_spring_summer %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_complete_winter_early_spring_sf <- df_complete_winter_early_spring %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_complete_late_spring_summer_sf <- df_complete_late_spring_summer %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")


# Map Generations 
# Load world map data and transform to Robinson projection
world <- ne_countries(scale = "medium", returnclass = "sf")
meridian <- -180
wld.new <- st_break_antimeridian(world, lon_0 = meridian)

wld.rob.sf <- st_transform(wld.new, paste("+proj=robin +lon_0=", meridian, "+k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# Plot maps for DJFMAM Argo profiles
map_winter_early_spring_argo_prof <- ggplot() +
  geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
  geom_sf(data = df_argo_winter_early_spring_sf, aes(geometry = geometry), color = "blue", alpha = 0.1, size = 0.1) +
  coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
  labs(title = "DJFMAM Argo Profiles",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"))

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_winter_early_spring_argo_prof.jpg", map_winter_early_spring_argo_prof)

# Plot maps for JJASON Argo profiles
map_late_spring_summer_argo_prof <- ggplot() +
  geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
  geom_sf(data = df_argo_late_spring_summer_sf, aes(geometry = geometry), color = "blue", alpha = 0.1, size = 0.1) +
  coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
  labs(title = "JJASON Argo Profiles",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"))

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_late_spring_summer_argo_prof.jpg", map_late_spring_summer_argo_prof)

# Plot maps for DJFMAM detected subduction events
map_winter_early_spring_detected_events <- ggplot() +
  geom_sf(data = wld.rob.sf, fill = "gray80", color = "gray") +
  geom_sf(data = df_complete_winter_early_spring_sf, aes(geometry = geometry), color = "red", alpha = 0.5, size = 0.3) +
  coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
  labs(title = "DJFMAM Detected Subduction Events",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"))

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_winter_early_spring_detected_events.jpg", map_winter_early_spring_detected_events)

# Plot maps for JJASON detected subduction events
map_late_spring_summer_detected_events <- ggplot() +
  geom_sf(data = wld.rob.sf, fill = "gray80", color = "gray") +
  geom_sf(data = df_complete_late_spring_summer_sf, aes(geometry = geometry), color = "red", alpha = 0.5, size = 0.3) +
  coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
  labs(title = "JJASON Detected Subduction Events",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"))

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_late_spring_summer_detected_events.jpg", map_late_spring_summer_detected_events)

##############
## KDE #######
##############
# Run kde2d on the separated Argo profile data
kde_argo_winter <- MASS::kde2d(df_argo_winter_early_spring$LONGITUDE, df_argo_winter_early_spring$LATITUDE, n = 180)
kde_argo_summer <- MASS::kde2d(df_argo_late_spring_summer$LONGITUDE, df_argo_late_spring_summer$LATITUDE, n = 180)

# Run kde2d on the separated subduction events data
kde_events_winter <- MASS::kde2d(df_complete_winter_early_spring$LONGITUDE, df_complete_winter_early_spring$LATITUDE, n = 180)
kde_events_summer <- MASS::kde2d(df_complete_late_spring_summer$LONGITUDE, df_complete_late_spring_summer$LATITUDE, n = 180)

# Convert the KDE results to a data frame
kde_to_dataframe <- function(kde_result) {
  expand.grid(x = kde_result$x, y = kde_result$y) %>%
    mutate(z = as.vector(kde_result$z))
}

kde_argo_winter_df <- kde_to_dataframe(kde_argo_winter)
kde_argo_summer_df <- kde_to_dataframe(kde_argo_summer)
kde_events_winter_df <- kde_to_dataframe(kde_events_winter)
kde_events_summer_df <- kde_to_dataframe(kde_events_summer)

# Bin the longitude and latitude for each Argo profile and subduction event
bin_size <- 5
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE)), ceiling(max(df_argo_clean$LONGITUDE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE)), ceiling(max(df_argo_clean$LATITUDE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins to Argo profiles and anomalies for each season
df_argo_winter_early_spring <- df_argo_winter_early_spring %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
  )

df_argo_late_spring_summer <- df_argo_late_spring_summer %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
  )

df_complete_winter_early_spring <- df_complete_winter_early_spring %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
  )

df_complete_late_spring_summer <- df_complete_late_spring_summer %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
  )

# Convert bin labels to numeric for calculations
df_argo_winter_early_spring <- df_argo_winter_early_spring %>%
  mutate(
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

df_argo_late_spring_summer <- df_argo_late_spring_summer %>%
  mutate(
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

df_complete_winter_early_spring <- df_complete_winter_early_spring %>%
  mutate(
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

df_complete_late_spring_summer <- df_complete_late_spring_summer %>%
  mutate(
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

# Compute counts for each bin and time period
total_counts_winter <- df_argo_winter_early_spring %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count_total = n(), .groups = 'drop')

total_counts_summer <- df_argo_late_spring_summer %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count_total = n(), .groups = 'drop')

anomaly_counts_winter <- df_complete_winter_early_spring %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count_anomaly = n(), .groups = 'drop')

anomaly_counts_summer <- df_complete_late_spring_summer %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count_anomaly = n(), .groups = 'drop')

# Merge counts and compute proportions for each time period
merged_counts_winter <- full_join(total_counts_winter, anomaly_counts_winter, by = c("lon_bin", "lat_bin")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  ) %>%
  filter(count_total > 0) %>% 
  filter(proportion<0.8)  %>% filter(!is.na(lat_bin))


merged_counts_summer <- full_join(total_counts_summer, anomaly_counts_summer, by = c("lon_bin", "lat_bin")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  ) %>%
  filter(count_total > 0)%>% 
  filter(proportion<0.8) %>% filter(!is.na(lat_bin))


# Function to plot density map
plot_density_map <- function(kde_df, title, output_file) {
  ggplot() +
    geom_tile(data = kde_df, aes(x = x, y = y, fill = z)) +
    geom_contour(data = kde_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
    coord_fixed() +
    geom_sf(data = world, fill = "gray80", color = "black") +
    scale_fill_viridis_c() +
    labs(title = title,
         x = "Longitude",
         y = "Latitude",
         fill = "Density") +
    theme_minimal() -> p
  
  # Save the plot
  ggsave(filename = output_file, plot = p, width = 10, height = 6)
}

# Function to plot proportion map
plot_proportion_map <- function(merged_counts, title, output_file) {
  ggplot(merged_counts, aes(x = lon_bin, y = lat_bin)) +
    geom_tile(aes(fill = proportion)) +
    geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
    geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
    coord_sf(xlim = range(merged_counts$lon_bin), ylim = range(merged_counts$lat_bin)) +
    scale_fill_viridis_c() +
    labs(
      title = title,
      x = "Longitude",
      y = "Latitude",
      fill = "Proportion"
    ) +
    theme_minimal() -> p
  
  # Save the plot
  ggsave(output_file, plot = p, width = 10, height = 6)
}

# Plotting density maps for DJFMAM and JJASON

# Winter/Early Spring (DJFMAM) Argo profile density
plot_density_map(kde_argo_winter_df,
                 "Density of Winter/Early Spring Argo Profiles",
                 "/data/GLOBARGO/figures/TimeSpaceVar/density_winter_argo_profiles.png")

# Late Spring/Summer (JJASON) Argo profile density
plot_density_map(kde_argo_summer_df,
                 "Density of Late Spring/Summer Argo Profiles",
                 "/data/GLOBARGO/figures/TimeSpaceVar/density_summer_argo_profiles.png")

# Winter/Early Spring (DJFMAM) Subduction events density
plot_density_map(kde_events_winter_df,
                 "Density of Winter/Early Spring Subduction Events",
                 "/data/GLOBARGO/figures/TimeSpaceVar/density_winter_subduction_events.png")

# Late Spring/Summer (JJASON) Subduction events density
plot_density_map(kde_events_summer_df,
                 "Density of Late Spring/Summer Subduction Events",
                 "/data/GLOBARGO/figures/TimeSpaceVar/density_summer_subduction_events.png")

# Plotting proportion maps for DJFMAM and JJASON

# Winter/Early Spring (DJFMAM) Proportion of subduction events to Argo profiles
plot_proportion_map(merged_counts_winter,
                    "Proportion of Subduction Events (Winter/Early Spring)",
                    "/data/GLOBARGO/figures/TimeSpaceVar/prop_winter_subduction_events.png")

# Late Spring/Summer (JJASON) Proportion of subduction events to Argo profiles
plot_proportion_map(merged_counts_summer,
                    "Proportion of Subduction Events (Late Spring/Summer)",
                    "/data/GLOBARGO/figures/TimeSpaceVar/prop_summer_subduction_events.png")


# Fit the GAM for Winter/Early Spring (DJFMAM)
gam_winter <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos", k = 300),
  family = binomial(link = "logit"),
  data = merged_counts_winter,
  method = "REML"
)

# Fit the GAM for Late Spring/Summer (JJASON)
gam_summer <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos", k = 300),
  family = binomial(link = "logit"),
  data = merged_counts_summer,
  method = "REML"
)

# Create prediction grids for each time period
# Define grid resolution for winter/early spring
lon_seq_winter <- seq(min(merged_counts_winter$lon_bin), max(merged_counts_winter$lon_bin), by = 1)
lat_seq_winter <- seq(min(merged_counts_winter$lat_bin), max(merged_counts_winter$lat_bin), by = 1)

prediction_grid_winter <- expand.grid(lon_bin = lon_seq_winter, lat_bin = lat_seq_winter)

# Define grid resolution for late spring/summer
lon_seq_summer <- seq(min(merged_counts_summer$lon_bin), max(merged_counts_summer$lon_bin), by = 1)
lat_seq_summer <- seq(min(merged_counts_summer$lat_bin), max(merged_counts_summer$lat_bin), by = 1)
prediction_grid_summer <- expand.grid(lon_bin = lon_seq_summer, lat_bin = lat_seq_summer)

# Predict the smoothed probabilities for winter/early spring
prediction_grid_winter$proportion <- predict(
  gam_winter,
  newdata = prediction_grid_winter,
  type = "response"
)

# Predict the smoothed probabilities for late spring/summer
prediction_grid_summer$proportion <- predict(
  gam_summer,
  newdata = prediction_grid_summer,
  type = "response"
)

# Remove any NA values from both grids
prediction_grid_winter <- prediction_grid_winter %>% filter(!is.na(proportion))
prediction_grid_summer <- prediction_grid_summer %>% filter(!is.na(proportion))

# Determine x and y limits for winter
x_limits_winter <- c(min(prediction_grid_winter$lon_bin), max(prediction_grid_winter$lon_bin))
y_limits_winter <- c(min(prediction_grid_winter$lat_bin), max(prediction_grid_winter$lat_bin))

# Determine x and y limits for summer
x_limits_summer <- c(min(prediction_grid_summer$lon_bin), max(prediction_grid_summer$lon_bin))
y_limits_summer <- c(min(prediction_grid_summer$lat_bin), max(prediction_grid_summer$lat_bin))

# Plot the smoothed proportions for Winter/Early Spring (DJFMAM)
ggplot() +
  geom_tile(data = prediction_grid_winter, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid_winter, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3, inherit.aes = TRUE) +
  geom_sf(data = world, fill = "gray80", color = "gray80", inherit.aes = FALSE) +
  coord_sf(
    xlim = x_limits_winter,
    ylim = y_limits_winter,
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(
    title = "Estimated Proportion of Anomalous Argo Profiles (Winter/Early Spring)",
    x = "Longitude",
    y = "Latitude",
    fill = "Proportion"
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/gam_winter_subduction_events.png")

# Plot the smoothed proportions for Late Spring/Summer (JJASON)
ggplot() +
  geom_tile(data = prediction_grid_summer, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid_summer, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3, inherit.aes = TRUE) +
  geom_sf(data = world, fill = "gray80", color = "gray80", inherit.aes = FALSE) +
  coord_sf(
    xlim = x_limits_summer,
    ylim = y_limits_summer,
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(
    title = "Estimated Proportion of Anomalous Argo Profiles (Late Spring/Summer)",
    x = "Longitude",
    y = "Latitude",
    fill = "Proportion"
  ) +
  scale_fill_viridis_c() +
  theme_minimal()

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/gam_summer_subduction_events.png")



# Function to fit the GAM with different k values and return model evaluation metrics (AIC and GCV), biggest drop in AIC and GCV at around k = 300
fit_gam_and_evaluate <- function(k_value, data, season) {
  gam_model <- gam(
    cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos", k = k_value),
    family = binomial(link = "logit"),
    data = data,
    method = "REML"
  )
  
  # Return AIC and GCV scores
  list(
    k = k_value,
    AIC = AIC(gam_model),
    GCV = gam_model$gcv.ubre,
    model = gam_model,
    season = season
  )
}
# Find the number of unique lat/lon bins in the winter dataset
n_unique_winter <- nrow(unique(merged_counts_winter[, c("lat_bin", "lon_bin")]))
n_unique_summer <- nrow(unique(merged_counts_summer[, c("lat_bin", "lon_bin")]))

# Print the number of unique bins
print(paste("Number of unique lat/lon bins (Winter):", n_unique_winter))
print(paste("Number of unique lat/lon bins (Summer):", n_unique_summer))


# Define a range of k values to test
k_values <- seq(100, 1000, by = 50)  # You can modify this range as needed

# Fit models for Winter/Early Spring (DJFMAM)
results_winter <- lapply(k_values, function(k) fit_gam_and_evaluate(k, merged_counts_winter, "Winter"))

# Fit models for Late Spring/Summer (JJASON)
results_summer <- lapply(k_values, function(k) fit_gam_and_evaluate(k, merged_counts_summer, "Summer"))

# Extract AIC and GCV results for comparison
extract_results <- function(results) {
  do.call(rbind, lapply(results, function(res) {
    data.frame(k = res$k, AIC = res$AIC, GCV = res$GCV, Season = res$season)
  }))
}

# Combine results for winter and summer
results_df <- rbind(extract_results(results_winter), extract_results(results_summer))

# View the results
print(results_df)

# Plot AIC and GCV vs k to visualize the best k
library(ggplot2)

ggplot(results_df, aes(x = k, y = AIC, color = Season)) +
  geom_line() +
  geom_point() +
  labs(title = "AIC vs k for Winter and Summer", x = "k (Number of Basis Functions)", y = "AIC") +
  theme_minimal()

ggplot(results_df, aes(x = k, y = GCV, color = Season)) +
  geom_line() +
  geom_point() +
  labs(title = "GCV vs k for Winter and Summer", x = "k (Number of Basis Functions)", y = "GCV") +
  theme_minimal()

# Observe that AIC and GCV drop significantly around $k = 300$.
