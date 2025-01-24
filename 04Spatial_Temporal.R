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

# Two datasets, one dataset df_complete_clean which contains time, location of the 4122 anomalous Argo profiles detected
# df_argo_clean which contains time and location of all the 125,836 Argo profiles (anomalous and non-anomalous)
# Goal of this script is to investigate the time-space variability of the anomalies, showing hotspots of anomalies.
# Of course, more anomalies are expected in regions where there are more floats, which is why we need to find a sensible way
# to normalize the anomalies rate by the number of non-anomalous present in the vicinity. A gross way to do that might be to compute
# anomalies proportion per square of degree LON/LAT but there might be a more elegant way to do that by computing a constant area radius arround
# the anomalies...

# first investigate structure of data
df_complete_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
df_complete_clean %>% head()
df_argo_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_argo_loc.csv")
df_argo_clean %>% head()

df_complete_clean
df_argo_clean

# so now that you've this dataframe of detected events, 
# we are going to study the spatial and temporal variation of those 
# subduction events of filaments of AOU and salinity,
# most likely due to symmetric and MLI instabilities 
# in the mixed layer of the ocean. Since those 3545 events are detected out 
# of the much larger Argo Database, we begin by importing the Argo Database to
# compute proportions of events out of total number of argo floats present in each 5*5 degree gridcell.
# To do this, we have at our disposal df_argo with 126,591 rows, one row for each profile,
# uniquely identified by float number (WMO) and TIME, LAT, LON, and CYCLE NUMBER# 

# Load world map data and transform to Robinson projection
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define new meridian
meridian <- -180

# Split world at new meridian
wld.new <- st_break_antimeridian(world, lon_0 = meridian)


wld.rob.sf <-  st_transform(wld.new, 
                            paste("+proj=robin +lon_0=", meridian,
                                  "+k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") )




# Convert data frames to spatial objects and transform to Robinson projection
df_argo_sf <- df_argo_clean %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%  # Remove missing values
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_complete_sf <- df_complete_clean %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%  # Remove missing values
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

# Plot total Argo profiles
map_total_argo_prof <- ggplot() +
  geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
  geom_sf(data = df_argo_sf, aes(geometry = geometry), color = "blue", alpha = 0.1, size = 0.1) +
  coord_sf(crs = "+proj=robin", datum = NA) +
  labs(title = "125,836 Argo Profiles",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90")) +
  annotation_north_arrow(location = "tl")

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_total_argo_prof.jpg",map_total_argo_prof)

# Plot detected subduction events

map_detected_events <- ggplot() +
  geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
  geom_sf(data = df_complete_sf, aes(geometry = geometry), color = "red", alpha = 0.1, size = 0.1) +
  coord_sf(crs = "+proj=robin", datum = NA) +
  labs(title = "4,390 subduction anomalies",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90")) +
  annotation_north_arrow(location = "tl")


ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_detected_events.jpg",map_detected_events)


####################""
### KDE ##############""
#######################

# Run kde2d on cleaned data
kde_argo <- MASS::kde2d(df_argo_clean$LONGITUDE, df_argo_clean$LATITUDE, n = 180)


# Run kde2d on cleaned subduction events data
kde_events <- MASS::kde2d(df_complete_clean$LONGITUDE, df_complete_clean$LATITUDE, n = 180)

# Convert kde2d results to a data frame for ggplot2
kde_to_dataframe <- function(kde_result) {
  expand.grid(x = kde_result$x, y = kde_result$y) %>%
    mutate(z = as.vector(kde_result$z))
}


# Convert the KDE outputs
kde_argo_df <- kde_to_dataframe(kde_argo)
kde_events_df <- kde_to_dataframe(kde_events)



# Plot density of Argo profiles with continents using coord_sf()
ggplot() +
  geom_sf(data = kde_argo_df, aes(x = x, y = y, fill = z)) +
  geom_contour(data = kde_argo_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
  coord_fixed() +
  geom_sf(data = world, fill = "gray80", color = "black") +  # Add continents
  coord_sf(xlim = range(kde_argo_df$x), ylim = range(kde_argo_df$y)) +  # Use coord_sf() to properly render sf data
  scale_fill_viridis_c() +
  labs(title = "Density of Argo Floats Events",
       x = "Longitude",
       y = "Latitude",
       fill = "Density") +
  theme_minimal()

ggsave(a,filename = "/data/GLOBARGO/figures/TimeSpaceVar/kde_argo_n1000.png")

# Plot density of argo events with continents using coord_sf()
map_argo_float_density <- ggplot() +
  geom_tile(data = kde_argo_df, aes(x = x, y = y, fill = z)) +
  geom_contour(data = kde_argo_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
  coord_fixed() +
  geom_sf(data = world, fill = "gray80", color = "black") +  # Add continents
  coord_sf(xlim = range(kde_events_df$x), ylim = range(kde_events_df$y)) +  # Use coord_sf() to properly render sf data
  scale_fill_viridis_c() +
  labs(title = "Density of Argo Floats",
       x = "Longitude",
       y = "Latitude",
       fill = "Density") +
  theme_minimal()

ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/map_argo_float_density.jpg",map_argo_float_density)

# Shift lon : 
# Shift longitudes in kde_argo_df
# Shift longitudes in your data


# Plot density of subduction events with continents using coord_sf()
map_event_density <- ggplot() +
  geom_tile(data = kde_events_df, aes(x = x, y = y, fill = z)) +
  geom_contour(data = kde_events_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
  coord_fixed() +
  geom_sf(data = world, fill = "gray80", color = "black") +  # Add continents
  coord_sf(xlim = range(kde_events_df$x), ylim = range(kde_events_df$y)) +  # Use coord_sf() to properly render sf data
  scale_fill_viridis_c() +
  labs(title = "Density of Subduction Events",
       x = "Longitude",
       y = "Latitude",
       fill = "Density") +
  theme_minimal()



# Define bin size (e.g., 1-degree bins)
bin_size <- 5
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE)), ceiling(max(df_argo_clean$LONGITUDE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE)), ceiling(max(df_argo_clean$LATITUDE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins to total Argo profiles
df_argo_clean <- df_argo_clean %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
  )

# Assign bins to anomalies
df_complete_clean <- df_complete_clean %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
  )

# Convert bin labels to numeric
df_argo_clean <- df_argo_clean %>%
  mutate(
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

df_complete_clean <- df_complete_clean %>%
  mutate(
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

# Compute counts of total profiles per bin
total_counts <- df_argo_clean %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count_total = n(), .groups = 'drop')

# Compute counts of anomalies per bin
anomaly_counts <- df_complete_clean %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count_anomaly = n(), .groups = 'drop')

# Merge counts and compute proportions
merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  ) %>%
  filter(count_total > 0)  # Remove bins with zero total profiles



# Remove rows where count_anomaly > count_total
merged_counts <- merged_counts %>%
  filter(count_anomaly <= count_total) %>% filter(!is.na(lat_bin) )
merged_counts$proportion <- ifelse(merged_counts$proportion > 0.5, NA, merged_counts$proportion)


# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")


blr <- function(x) (log(x/(1-x)))

merged_counts$proportion %>% blr() %>% hist()

# Plot the proportion of anomalies in proportion space
prop_map <- ggplot(merged_counts, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black",inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts$lon_bin), ylim = range(merged_counts$lat_bin)) +
  scale_fill_viridis_c() +
  labs(
    title = "Proportion of Anomalous Argo Profiles",
    x = "Longitude",
    y = "Latitude",
    fill = "Proportion"
  ) +
  theme_minimal()


ggsave(prop_map,filename = "/data/GLOBARGO/figures/TimeSpaceVar/prop_map.png")



# Plot the proportion of anomalies in logit space 
ggplot(merged_counts, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = blr(proportion))) +
  geom_contour(aes(z = blr(proportion)), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black",inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts$lon_bin), ylim = range(merged_counts$lat_bin)) +
  scale_fill_viridis_c() +
  labs(
    title = "Proportion of Anomalous Argo Profiles",
    x = "Longitude",
    y = "Latitude",
    fill = "Proportion"
  ) +
  theme_minimal()
library(spatstat)




# Fit the GAM model



gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos",k = 100),
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
)


# Define grid resolution
lon_seq <- seq(min(merged_counts$lon_bin), max(merged_counts$lon_bin), by = 1)
lat_seq <- seq(min(merged_counts$lat_bin), max(merged_counts$lat_bin), by = 1)

# Create grid
prediction_grid <- expand.grid(lon_bin = lon_seq, lat_bin = lat_seq)

# Predict the smoothed probabilities
prediction_grid$proportion <- predict(
  gam_model,
  newdata = prediction_grid,
  type = "response"
)


# Remove any NA values
prediction_grid <- prediction_grid %>% filter(!is.na(proportion))


# Determine x and y limits
x_limits <- c(min(prediction_grid$lon_bin), max(prediction_grid$lon_bin))
y_limits <- c(min(prediction_grid$lat_bin), max(prediction_grid$lat_bin))


prediction_grid %>% head()

# Plot the smoothed proportions
ggplot() +
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = (proportion))) +
  geom_contour(data = prediction_grid,aes(z = proportion,x=lon_bin,y=lat_bin), color = "white", alpha = 0.3,inherit.aes = TRUE) +
  geom_sf(data = world, fill = "gray80", color = "gray80", inherit.aes = FALSE) +
  coord_sf(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(
    title = "Estimated Proportion of Anomalous Argo Profiles from Generalized Additive Model with Spherical Spline Smoothing",
    x = "Longitude",
    y = "Latitude",
    fill = "Proportion"
  )  + scale_fill_viridis_c()+
  theme_minimal()







# Comparative analysis
# Example: Assign regions based on latitude and longitude boundaries
df_complete <- df_complete %>%
  mutate(Region = case_when(
    LATITUDE < -50 & LATITUDE > -60 ~ "Southern Ocean",
    LONGITUDE > -80 & LONGITUDE < -60 & LATITUDE > 30 & LATITUDE < 45 ~ "Gulf Stream",
    LONGITUDE > 140 & LONGITUDE < 160 & LATITUDE > 30 & LATITUDE < 45 ~ "Kuroshio Extension",
    TRUE ~ "Other"
  ))


df_argo <- df_argo %>%
  mutate(Region = case_when(
    LATITUDE < -50 & LATITUDE > -60 ~ "Southern Ocean",
    LONGITUDE > -80 & LONGITUDE < -60 & LATITUDE > 30 & LATITUDE < 45 ~ "Gulf Stream",
    LONGITUDE > 140 & LONGITUDE < 160 & LATITUDE > 30 & LATITUDE < 45 ~ "Kuroshio Extension",
    TRUE ~ "Other"
  ))

# # Total profiles per region
total_profiles_region <- df_argo %>%
  group_by(Region) %>%
  summarise(TotalProfiles = n())

# Subduction events per region
events_region <- df_complete %>%
  group_by(Region) %>%
  summarise(SubductionEvents = n())

# Combine and calculate proportions
region_summary <- total_profiles_region %>%
  left_join(events_region, by = "Region") %>%
  mutate(
    SubductionEvents = coalesce(SubductionEvents, 0),
    Proportion = SubductionEvents / TotalProfiles
  )

# View results
print(region_summary)


# b Statistical Testing

# Create contingency table
contingency_table <- region_summary %>%
  select(Region, SubductionEvents, TotalProfiles) %>%
  mutate(NonEvents = TotalProfiles - SubductionEvents) %>%
  select(-TotalProfiles) %>%
  pivot_longer(cols = c(SubductionEvents, NonEvents), names_to = "Event", values_to = "Count") %>%
  pivot_wider(names_from = Region, values_from = Count)

# Perform chi-squared test
chisq_result <- chisq.test(as.matrix(contingency_table[,-1]))
