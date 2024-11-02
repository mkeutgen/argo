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
library(multcomp)

logitInv <- function(x){1/(1+exp(-x))}
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
# Read data
df_complete_clean <- read_csv("/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
df_argo_clean <- read_csv("/data/GLOBARGO/src/data/df_argo_loc.csv")

# Define bin size
bin_size <- 5

# Define longitude and latitude bins
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE, na.rm = TRUE)), ceiling(max(df_argo_clean$LONGITUDE, na.rm = TRUE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE, na.rm = TRUE)), ceiling(max(df_argo_clean$LATITUDE, na.rm = TRUE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins to Argo profiles and anomalies
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

df_argo_clean <- assign_bins(df_argo_clean)
df_complete_clean <- assign_bins(df_complete_clean)


df_argo_clean <- df_argo_clean %>%
  mutate(
    hemisphere = ifelse(LATITUDE >= 0, "Northern", "Southern"),
    local_season = case_when(
      hemisphere == "Northern" & month(TIME) %in% c(12, 1, 2) ~ "Winter",
      hemisphere == "Northern" & month(TIME) %in% c(3, 4, 5) ~ "Spring",
      hemisphere == "Northern" & month(TIME) %in% c(6, 7, 8) ~ "Summer",
      hemisphere == "Northern" & month(TIME) %in% c(9, 10, 11) ~ "Autumn",
      hemisphere == "Southern" & month(TIME) %in% c(12, 1, 2) ~ "Summer",
      hemisphere == "Southern" & month(TIME) %in% c(3, 4, 5) ~ "Autumn",
      hemisphere == "Southern" & month(TIME) %in% c(6, 7, 8) ~ "Winter",
      hemisphere == "Southern" & month(TIME) %in% c(9, 10, 11) ~ "Spring"
    )
  )


df_complete_clean <- df_complete_clean %>%
  mutate(
    hemisphere = ifelse(LATITUDE >= 0, "Northern", "Southern"),
    local_season = case_when(
      hemisphere == "Northern" & month(TIME) %in% c(12, 1, 2) ~ "Winter",
      hemisphere == "Northern" & month(TIME) %in% c(3, 4, 5) ~ "Spring",
      hemisphere == "Northern" & month(TIME) %in% c(6, 7, 8) ~ "Summer",
      hemisphere == "Northern" & month(TIME) %in% c(9, 10, 11) ~ "Autumn",
      hemisphere == "Southern" & month(TIME) %in% c(12, 1, 2) ~ "Summer",
      hemisphere == "Southern" & month(TIME) %in% c(3, 4, 5) ~ "Autumn",
      hemisphere == "Southern" & month(TIME) %in% c(6, 7, 8) ~ "Winter",
      hemisphere == "Southern" & month(TIME) %in% c(9, 10, 11) ~ "Spring"
    )
  )



# Assign regions
assign_region <- function(lat, lon) {
  # Adjust longitude to 0-360 degrees
  if (lon < 0) {
    lon <- lon + 360
  }
  
  if (lat < -30) {
    return("Southern Ocean")
  } else if (lat > 30 & ((lon >= 260 & lon <= 360) | (lon >= 0 & lon <= 20))) {
    # North Atlantic: lat > 23.5°, lon between 260°-360° or 0°-20°
    return("North Atlantic")
  } else if (lat > 30 & lon >= 120 & lon < 260) {
    # North Pacific: lat > 23.5°, lon between 120°-260°
    return("North Pacific")
  } else {
    # All other regions
    return("Tropics")
  }
}

# Apply the function to Argo profiles
df_argo_clean <- df_argo_clean %>%
  mutate(region = mapply(assign_region, LATITUDE, LONGITUDE),
         year = year(TIME))

# Apply the function to subduction events
df_complete_clean <- df_complete_clean %>%
  mutate(region = mapply(assign_region, LATITUDE, LONGITUDE),
         year = year(TIME))

# Compute counts for each bin, season, and region
total_counts <- df_argo_clean %>%
  group_by(lon_bin, lat_bin, local_season, region,year) %>%
  summarize(count_total = n(), .groups = 'drop')

anomaly_counts <- df_complete_clean %>%
  group_by(lon_bin, lat_bin, local_season, region,year) %>%
  summarize(count_anomaly = n(), .groups = 'drop')

# Merge counts and compute proportions
merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin", "local_season", "region")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  )%>%
  filter(count_total > 0) %>%
  filter(proportion < 0.6) %>%  # Filter out high proportions (anomalies > 60%)
  filter(!is.na(lat_bin))


# Set 'Winter' as the reference season
merged_counts$local_season <- factor(merged_counts$local_season, levels = c("Winter", "Spring", "Summer", "Autumn"))

# Compute counts for each bin, season, and region
total_counts_yr <- df_argo_clean %>%
  group_by(local_season, region,year) %>%
  summarize(count_total = n(), .groups = 'drop')

anomaly_counts_yr <- df_complete_clean %>%
  group_by(local_season, region,year) %>%
  summarize(count_anomaly = n(), .groups = 'drop')



# Merge counts and compute proportions
merged_counts_yr <- full_join(total_counts_yr, anomaly_counts_yr, by = c("year", "local_season", "region")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  )%>%
  filter(count_total > 0)

# Ensure region and local_season are factors
merged_counts_yr <- merged_counts_yr %>%
  mutate(
    region = factor(region),
    local_season = factor(local_season)
  )


merged_counts_yr$local_season <- factor(merged_counts_yr$local_season, levels = c("Autumn", "Winter", "Spring", "Summer"))

merged_counts_yr$region <- factor(merged_counts_yr$region, levels = c("Tropics","North Atlantic", "North Pacific" , "Southern Ocean"))

# Create a GLM with proportion as the response variable, using binomial family (logit link)
glm_model <- glm(
  cbind(count_anomaly, count_total - count_anomaly) ~ region + local_season,
  data = merged_counts_yr,
  family = binomial(link = "logit")
)

summary(glm_model)

baseline_odds <- glm_model$coefficients[1]

baseline_proba <- exp(baseline_odds)

adjustedlogodds <- baseline_odds + glm_model$coefficients
# Obtain probabilities from the GLM
probability <- exp(adjustedlogodds)/(1+exp(adjustedlogodds))

# Conduct Tukey's post-hoc test for pairwise region comparisons within each season
posthoc_results <- glht(glm_model, linfct = mcp(region = "Tukey"))

# Summarize results
summary(posthoc_results)



# Fit the model
gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ 
    s(lat_bin, lon_bin, bs = "sos", k = 300) + 
    local_season + region,
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
)

summary(gam_model)

# Check model diagnostics
gam.check(gam_model)

# Fit the model (interaction)
gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ 
    s(lat_bin, lon_bin, bs = "sos", k = 5) + 
    local_season * region,
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
) # R^2 = 0.56

gam_model %>% summary()
# No interaction
gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ 
    s(lat_bin, lon_bin, bs = "sos", k = 5) + 
    local_season + region,
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
) # R^2 = 0.56



gam_model %>% summary()
# Convert 'region' to a factor
merged_counts$region <- factor(merged_counts$region)

gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ 
    s(lat_bin, lon_bin, by = region, bs = "sos", k = 5) + 
    local_season+region,
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
)

summary(gam_model)

gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ 
    s(lat_bin, lon_bin, bs = "sos", k = 900) + 
    local_season + region,
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
) 

summary(gam_model)
logitInv <- function(x){1/(1+exp(-x))}

# Set 'Winter' as the reference season
merged_counts$local_season <- factor(merged_counts$local_season, levels = c("Winter", "Spring", "Summer", "Autumn"))

# Compute counts for each bin, season, and region
total_counts_yr <- df_argo_clean %>%
  group_by(local_season, region,year) %>%
  summarize(count_total = n(), .groups = 'drop')

anomaly_counts_yr <- df_complete_clean %>%
  group_by(local_season, region,year) %>%
  summarize(count_anomaly = n(), .groups = 'drop')




# Merge counts and compute proportions
merged_counts_yr <- full_join(total_counts_yr, anomaly_counts_yr, by = c("year", "local_season", "region")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  )%>%
  filter(count_total > 0)



# Define logit and inverse logit functions
logit <- function(p) { log(p / (1 - p)) }
logitInv <- function(x) { 1 / (1 + exp(-x)) }


# Calculate yearly mean proportions and standard errors using the logit scale
region_season_summary <- merged_counts_yr %>%
  filter(proportion > 0 & proportion < 1) %>%  # Avoids logit issues at bounds
  group_by(region, local_season) %>%
  summarize(
    mean_proportion = mean(proportion),
    se_proportion = sd(logit(proportion)) / sqrt(n()),  # Standard error in logit scale
    .groups = 'drop'
  ) %>%
  mutate(
    # Transform back to proportion scale for plotting error bars
    lower_bound = logitInv(logit(mean_proportion) - se_proportion),
    upper_bound = logitInv(logit(mean_proportion) + se_proportion)
  )

# Set 'Winter' as the reference season
region_season_summary$local_season <- factor(region_season_summary$local_season, levels = c("Winter", "Spring", "Summer", "Autumn"))


# Plot with error bars
ggplot(region_season_summary, aes(x = region, y = mean_proportion, fill = local_season)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black",alpha=.8) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),
                position = position_dodge(width = 0.9), width = 0.25, color = "black") +
  scale_fill_viridis_d() +
  labs(
    title = "Mean Proportion of Anomalous Profiles by Region and Season (with Standard Error)",
    x = "Region",
    y = "Mean Proportion of Anomalies",
    fill = "Season"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"))

ggsave(filename="/data/GLOBARGO/figures/TimeSpaceVar/4SEASONS/prop_map.png",width = 10,height = 6)



# Plot
ggplot(region_season_summary, aes(x = region, y = proportion, fill = local_season)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  scale_fill_viridis_d() +
  labs(
    title = "Proportion of Anomalous Profiles by Region and Season",
    x = "Region",
    y = "Proportion of Anomalies",
    fill = "Season"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"))


# Fit a GLM model with polynomial terms
glm_model <- glm(
  cbind(count_anomaly, count_total - count_anomaly) ~ 
     local_season + region,
  family = binomial(link = "logit"),
  data = merged_counts
)
# Summary of the model
summary(glm_model)


# Visualization
prediction_grid <- expand.grid(
  lon_bin = seq(min(merged_counts$lon_bin), max(merged_counts$lon_bin), by = 1),
  lat_bin = seq(min(merged_counts$lat_bin), max(merged_counts$lat_bin), by = 1),
  local_season = unique(merged_counts$local_season),
  region = unique(merged_counts$region)
) %>%
  filter(!is.na(region))

prediction_grid$proportion <- predict(gam_model, newdata = prediction_grid, type = "response")
prediction_grid <- prediction_grid %>% filter(!is.na(proportion))
prediction_grid <- prediction_grid %>% filter(proportion <= 0.5)

# Plot predicted probabilities with continents
ggplot() +
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(
    data = prediction_grid,
    aes(x = lon_bin, y = lat_bin, z = proportion),
    color = "white",
    alpha = 0.3
  ) +
  geom_sf(data = world, fill = "gray80", color = "black") +
  coord_sf(
    xlim = c(min(prediction_grid$lon_bin), max(prediction_grid$lon_bin)),
    ylim = c(min(prediction_grid$lat_bin), max(prediction_grid$lat_bin)),
    expand = FALSE
  ) +
  facet_grid(region ~ local_season) +
  scale_fill_viridis_c() +
  labs(
    title = "Estimated Proportion of Subduction Events by Region and Season",
    x = "Longitude",
    y = "Latitude",
    fill = "Proportion"
  ) +
  theme_minimal()









##########################

# Define seasons based on months
df_argo_clean <- df_argo_clean %>%
  mutate(
    season = case_when(
      month(TIME) %in% djf_months ~ "DJF",
      month(TIME) %in% mam_months ~ "MAM",
      month(TIME) %in% jja_months ~ "JJA",
      month(TIME) %in% son_months ~ "SON"
    )
  )

df_complete_clean <- df_complete_clean %>%
  mutate(
    season = case_when(
      month(TIME) %in% djf_months ~ "DJF",
      month(TIME) %in% mam_months ~ "MAM",
      month(TIME) %in% jja_months ~ "JJA",
      month(TIME) %in% son_months ~ "SON"
    )
  )

# Compute counts including season
total_counts <- df_argo_clean %>%
  group_by(lon_bin, lat_bin, season) %>%
  summarize(count_total = n(), .groups = 'drop')

anomaly_counts <- df_complete_clean %>%
  group_by(lon_bin, lat_bin, season) %>%
  summarize(count_anomaly = n(), .groups = 'drop')

# Merge counts and compute proportions
merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin", "season")) %>%
  mutate(
    count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
    proportion = count_anomaly / count_total
  ) %>%
  filter(count_total > 0) %>%
  filter(proportion < 0.8) %>%  # Filter out high proportions (anomalies > 80%)
  filter(!is.na(lat_bin))


# Fit one GAM model including season
gam_model <- gam(
  cbind(count_anomaly, count_total - count_anomaly) ~ s(lat_bin, lon_bin, bs = "sos", k = 300) + season,
  family = binomial(link = "logit"),
  data = merged_counts,
  method = "REML"
)


# Create prediction grid including season
prediction_grid <- expand.grid(
  lon_bin = seq(min(merged_counts$lon_bin), max(merged_counts$lon_bin), by = 1),
  lat_bin = seq(min(merged_counts$lat_bin), max(merged_counts$lat_bin), by = 1),
  season = unique(merged_counts$season)
)


# Predict the smoothed proportions
prediction_grid$proportion <- predict(gam_model, newdata = prediction_grid, type = "response")

# Remove NA values
prediction_grid <- prediction_grid %>% filter(!is.na(proportion))



# Plot predicted proportions for each season
gam_map <- ggplot() +
  geom_tile(data = prediction_grid, aes(x = lon_bin, y = lat_bin, fill = proportion)) +
  geom_contour(data = prediction_grid, aes(x = lon_bin, y = lat_bin, z = proportion), color = "white", alpha = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "black") +
  coord_sf(xlim = range(prediction_grid$lon_bin), ylim = range(prediction_grid$lat_bin), expand = FALSE) +
  facet_wrap(~ season, ncol = 2) +
  labs(title = "Estimated Proportion of Subduction Events by Season", x = "Longitude", y = "Latitude", fill = "Proportion") +
  scale_fill_viridis_c() +
  theme_minimal()

# Save the plot
ggsave(filename = "/data/GLOBARGO/figures/TimeSpaceVar/gam_map_by_season.png", plot = gam_map, width = 15, height = 12)
