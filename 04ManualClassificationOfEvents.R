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

# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load argo dataset
df_argo <- read_csv("/data/GLOBARGO/data/argo_profiles_df.csv")


# Load detected events :
df_abs_sal <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v4.csv")
df_spic <- read_csv("/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")



# Classified datasets
df_spic_class <- read_csv("/data/GLOBARGO/data/anom_in_spic_and_sal_cat1_and2.csv")
df_partial_salinity_class <- read_csv("/data/GLOBARGO/data/classification_results_salinity_fig_not_in_spic.csv")


df_partial_salinity_class <- df_partial_salinity_class %>%
  filter(Category %in% c(1, 2)) %>%
  mutate(WMO = gsub("_plot", "", WMO))  # Remove the 'plot' suffix

df_partial_salinity_class$CYCLE_NUMBER <- df_partial_salinity_class$Cycle

df_partial_salinity_class$WMO <- as.numeric(df_partial_salinity_class$WMO )
df_spic_class$WMO             <- as.numeric(df_spic_class$WMO             )    

# Perform the join, 
df_combined <- bind_rows(df_partial_salinity_class, df_spic_class) %>% select(WMO,CYCLE_NUMBER,Category)


# Next we wish to complete this dataframe by retriving
# LON/LAT/TIME info for each of those individual subduction events, df_abs_sal has the required info


# Check for duplicates in df_abs_sal based on WMO and CYCLE_NUMBER
df_abs_sal_unique <- df_abs_sal %>%
  distinct(WMO, CYCLE_NUMBER, .keep_all = TRUE)

df_abs_sal_unique$CYCLE_NUMBER

# Perform the join again with the unique rows
df_complete <- df_combined %>%
  left_join(df_abs_sal_unique %>%
              mutate(WMO = as.character(WMO)) %>% 
              select(WMO, CYCLE_NUMBER, LATITUDE, LONGITUDE, TIME), 
            by = c("WMO", "CYCLE_NUMBER"))




# so now that you've this dataframe of detected events, 
# we are going to study the spatial and temporal variation of those 
# subduction events of filaments of AOU and salinity,
# most likely due to symmetric and MLI instabilities 
# in the mixed layer of the ocean. Since those 3545 events are detected out 
# of the much larger Argo Database, we begin by importing the Argo Database to
# compute proportions of events out of total number of argo floats present in each 5*5 degree gridcell.
# To do this, we have at our disposal df_argo with 126,591 rows, one row for each profile,
# uniquely identified by float number (WMO) and TIME, LAT, LON, and CYCLE NUMBER# 

dx <- 10

# Step 1: Bin Argo profiles into 5x5 degree grid cells in df_argo
df_argo_binned <- df_argo %>%
  mutate(
    LAT_bin = floor(LATITUDE / dx) * dx,   # Bin latitude into 5-degree intervals
    LON_bin = floor(LONGITUDE / dx) * dx   # Bin longitude into 5-degree intervals
  ) %>%
  group_by(LAT_bin, LON_bin) %>%
  summarise(total_profiles = n(), .groups = 'drop')  # Count total profiles in each grid cell

# Step 2: Bin detected subduction events into 5x5 degree grid cells
df_events_binned <- df_complete %>%
  mutate(
    LAT_bin = floor(LATITUDE / dx) * dx,   # Bin latitude into 5-degree intervals
    LON_bin = floor(LONGITUDE / dx) * dx   # Bin longitude into 5-degree intervals
  ) %>%
  group_by(LAT_bin, LON_bin) %>%
  summarise(subduction_events = n(), .groups = 'drop')  # Count events in each grid cell

# Step 3: Join the two datasets and compute the proportion of subduction events
df_proportion <- df_argo_binned %>%
  left_join(df_events_binned, by = c("LAT_bin", "LON_bin")) %>%
  mutate(
    subduction_events = coalesce(subduction_events, 0),  # Replace NAs with 0 (no events)
    proportion = subduction_events / total_profiles      # Calculate proportion of events
  )

# Step 4: View the result
df_proportion %>% head()

# Load world map data and transform to Robinson projection
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define new meridian
meridian <- -180

# Split world at new meridian
wld.new <- st_break_antimeridian(world, lon_0 = meridian)

wld.rob.sf <-  st_transform(wld.new, 
                            paste("+proj=robin +lon_0=", meridian2 ,
                                  "+k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") )




# Convert data frames to spatial objects and transform to Robinson projection
df_argo_sf <- df_argo %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%  # Remove missing values
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_complete_sf <- df_complete %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%  # Remove missing values
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

# Plot total Argo profiles
ggplot() +
  geom_sf(data = wld.rob.sf, fill = "grey", color = "white") +
  geom_sf(data = df_argo_sf, aes(geometry = geometry), color = "blue", alpha = 0.1, size = 0.1) +
  coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
  labs(title = "Total Argo Profiles",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90")) +
  annotation_north_arrow(location = "tl")

# Plot detected subduction events
ggplot() +
  geom_sf(data = wld.rob.sf, fill = "gray80", color = "black") +
  geom_sf(data = df_complete_sf, aes(geometry = geometry), color = "red", alpha = 0.5, size = 0.3) +
  coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
  labs(title = "Detected Subduction Events",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90")) +
  annotation_north_arrow(location = "tl")



 ####################""
### KDE ##############""
#######################

# Clean Argo data
df_argo_clean <- df_argo %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))


# Run kde2d on cleaned data
kde_argo <- MASS::kde2d(df_argo_clean$LONGITUDE, df_argo_clean$LATITUDE, n = 100)

# Clean subduction events data
df_complete_clean <- df_complete %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))


# Run kde2d on cleaned subduction events data
kde_events <- MASS::kde2d(df_complete_clean$LONGITUDE, df_complete_clean$LATITUDE, n = 100)

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
  geom_tile(data = kde_argo_df, aes(x = x, y = y, fill = z)) +
  geom_contour(data = kde_argo_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
  coord_fixed() +
  geom_sf(data = world, fill = "gray80", color = "black") +  # Add continents
  coord_sf(xlim = range(kde_events_df$x), ylim = range(kde_events_df$y)) +  # Use coord_sf() to properly render sf data
  scale_fill_viridis_c() +
  labs(title = "Density of Argo Floats Events",
       x = "Longitude",
       y = "Latitude",
       fill = "Density") +
  theme_minimal()



# Plot density of subduction events with continents using coord_sf()
ggplot() +
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


# Ripley's K Function
library(spatstat)

# Define observation window using adjusted longitudes
win <- owin(
  xrange = range(df_complete_clean$LONGITUDE),
  yrange = range(df_complete_clean$LATITUDE)
)

# Create point pattern object
ppp_events <- ppp(
  x = df_complete_clean$LONGITUDE,
  y = df_complete_clean$LATITUDE,
  window = win
)

# Perform Ripley's K analysis
K <- Kest(ppp_events)
plot(K)



# d. Hotspot Analysis

#Objective: Identify regions where subduction events are significantly concentrated.

#Actions:
  
#  Getis-Ord Gi Statistic:* r ??


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
