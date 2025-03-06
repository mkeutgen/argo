library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)          # for coord_sf, st_crs, etc.
library(lubridate)   # if you need date manipulation

#---------------------------------------------------------
# Example: Show 10 km bins that contain >= 1 Argo profile
#---------------------------------------------------------

# --- 1) Create WEEK and YEAR columns, then define bin edges ---
df_argo_clean <- df_argo_clean %>%
  mutate(
    WEEK = week(TIME),
    YEAR = year(TIME),
    MONTH = year(TIME)
  )


# 1. Define bin size (0.2° ~ 20 km at mid-latitudes)
resolution <- 0.2

# 2. Bin the data & aggregate
argo_bins <- df_argo_clean %>%
  mutate(
    # create a 'bin origin' by flooring longitude and latitude
    lon_bin = floor(LONGITUDE / resolution) * resolution,
    lat_bin = floor(LATITUDE / resolution) * resolution
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count = n(), .groups = "drop")

nb_argo_total <- argo_bins %>% nrow()

argo_mam <- df_argo_clean %>% filter(MONTH %in% c(12,1,2)) %>%
  mutate(
    # create a 'bin origin' by flooring longitude and latitude
    lon_bin = floor(LONGITUDE / resolution) * resolution,
    lat_bin = floor(LATITUDE / resolution) * resolution
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count = n(), .groups = "drop")


argo_mam_2023 <- df_argo_clean %>% filter(YEAR == 2023, MONTH %in% c(3,4,5)) %>%
  mutate(
    # create a 'bin origin' by flooring longitude and latitude
    lon_bin = floor(LONGITUDE / resolution) * resolution,
    lat_bin = floor(LATITUDE / resolution) * resolution
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count = n(), .groups = "drop")

nb_argo_mam_2023 <- argo_mam_2023 %>% nrow()

nb_argo_mam <- argo_mam %>% nrow()

argo_submeso <- df_argo_clean %>% filter(YEAR == 2024,WEEK == 5) %>%
  mutate(
    # create a 'bin origin' by flooring longitude and latitude
    lon_bin = floor(LONGITUDE / resolution) * resolution,
    lat_bin = floor(LATITUDE / resolution) * resolution
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(count = n(), .groups = "drop") 

nb_argo_submeso <- argo_submeso %>% nrow()

# Optional: If you want a full grid from -180 to 180, -90 to 90, do:
full_grid <- expand.grid(
   lon_bin = seq(-180, 180 - resolution, by = resolution),
   lat_bin = seq(-90, 90 - resolution, by = resolution)
 ) %>% 
   as_tibble() %>%
   left_join(argo_bins, by = c("lon_bin", "lat_bin")) %>%
   mutate(count = replace_na(count, 0)) 

full_grid_mam <- expand.grid(
  lon_bin = seq(-180, 180 - resolution, by = resolution),
  lat_bin = seq(-90, 90 - resolution, by = resolution)
) %>% 
  as_tibble() %>%
  left_join(argo_mam, by = c("lon_bin", "lat_bin")) %>%
  mutate(count = replace_na(count, 0)) 



# Optional: If you want a full grid from -180 to 180, -90 to 90, do:
full_grid_submeso <- expand.grid(
  lon_bin = seq(-180, 180 - resolution, by = resolution),
  lat_bin = seq(-90, 90 - resolution, by = resolution)
) %>% 
  as_tibble() %>%
  left_join(argo_submeso, by = c("lon_bin", "lat_bin")) %>%
  mutate(count = replace_na(count, 0)) 

# Compute proportions of ocean cells that don't have floats using rule of thumb that
# 71 \% of globe area is ocean
ocean_cells_without_argo <- ((full_grid) %>% nrow())*0.71

ocean_cells_with_argo <- (full_grid %>% filter(full_grid$count> 0) %>% nrow())
ocean_cells_with_argo_in_mam <- (full_grid_mam %>% filter(full_grid_mam$count> 0) %>% nrow())
ocean_cells_with_argo_submeso <- (full_grid_submeso %>% filter(full_grid_submeso$count> 0) %>% nrow())


nb_argo/ocean_cells_without_argo * 100


nb_argo_mam/ocean_cells_without_argo * 100

nb_argo_mam_2023/ocean_cells_without_argo * 100

nb_argo_submeso/ocean_cells_without_argo * 100



# 3. Filter to bins that DO have at least one float
sampled_bins <- argo_bins %>% filter(count > 0)

# 4. For each bin, generate its 4 corners so we can plot squares with geom_polygon
#    (each row -> 4 corners)
sampled_polygons <- sampled_bins %>%
  rowwise() %>%
  mutate(
    corners = list(
      data.frame(
        LON = c(lon_bin,
                lon_bin + resolution,
                lon_bin + resolution,
                lon_bin),
        LAT = c(lat_bin,
                lat_bin,
                lat_bin + resolution,
                lat_bin + resolution),
        # group identifier for geom_polygon
        bin_id = paste0(lon_bin, "_", lat_bin)
      )
    )
  ) %>%
  ungroup() %>%
  unnest(cols = corners)

#---------------------------------------------------------
# 5. Plot: show red squares where Argo floats exist
#---------------------------------------------------------

# Suppose you have an sf object "world_data" in EPSG:4326 
# (e.g., from rnaturalearth or a shapefile read with st_read)
# If you do not have an sf object, you can get quick data with:
#   library(rnaturalearth)
  world_data <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  # Draw land
  geom_sf(data = world_data, fill = "gray", color = "gray") +
  
  # Draw squares in red (group by bin_id)
  geom_polygon(
    data = sampled_polygons,
    aes(x = LON, y = LAT, group = bin_id),
    fill = "red",
    color = NA,
    alpha = 1
  ) +
  
  # Use coordinate system in EPSG:4326
  coord_sf(
    crs    = st_crs(4326),
    xlim   = c(-180, 180),
    ylim   = c(-90, 90),
    expand = FALSE
  ) +
  
  labs(
    title = "Argo Float Locations (~20 km resolution)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

# Show where Argo mam floats exist in mam 2023

# 3. Filter to bins that DO have at least one float
sampled_bins <- argo_mam_2023 %>% filter(count > 0)

# 4. For each bin, generate its 4 corners so we can plot squares with geom_polygon
#    (each row -> 4 corners)
sampled_polygons <- sampled_bins %>%
  rowwise() %>%
  mutate(
    corners = list(
      data.frame(
        LON = c(lon_bin,
                lon_bin + resolution,
                lon_bin + resolution,
                lon_bin),
        LAT = c(lat_bin,
                lat_bin,
                lat_bin + resolution,
                lat_bin + resolution),
        # group identifier for geom_polygon
        bin_id = paste0(lon_bin, "_", lat_bin)
      )
    )
  ) %>%
  ungroup() %>%
  unnest(cols = corners)


undersampled_data <- ggplot() +
  # Draw land
  geom_sf(data = world_data, fill = "gray", color = "gray") +
  geom_point(data = sampled_bins,aes(x=lon_bin,y=lat_bin),color="black",alpha=.5)+
  # Draw squares in red (group by bin_id), with thicker outlines
  geom_polygon(
    data = sampled_polygons,
    aes(
      x = LON, 
      y = LAT, 
      group = bin_id,
      fill = "0.2° gridboxes where at least one float was present between March and May 2023"
    ),
    color = NA,
    alpha = 1,
    size = 5  # multiplied the original size=5 by 10
  ) +
  
  # Use coordinate system in EPSG:4326
  coord_sf(
    crs    = st_crs(4326),
    xlim   = c(-180, 180),
    ylim   = c(-90, 90),
    expand = FALSE
  ) +
  labs(
    title = "At the submesoscale (0.2° resolution) Argo floats are blind to over 99 % of the ocean",
    subtitle = "Argo floats location between 1 March and 30 May 2023",
    x = "Longitude",
    y = "Latitude",
    fill = NULL  # you can adjust or remove the fill label as you wish
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "blue", size = 16),
    axis.text.y.right  = element_text(color = "blue", size = 12),
    legend.text = element_text(size = 12),
    title = element_text(size = 15),
    legend.position = "bottom"
  )
ggsave(undersampled_data,
       filename = "figures/map_undersampled.png",
       width = 20,height = 15)


# Since it can accumulate for 3 months, relevant frequency estimates are : 
# How many carbon subduction events happen in (0.1)^2 gridcell per season. Dimension ?
# 3 subduction per month
