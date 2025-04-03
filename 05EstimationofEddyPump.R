library(tidyverse)
library(lubridate)
library(mgcv)
library(sf)
library(ggspatial)
library(viridis)
library(patchwork)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)

# Read the data
df_carbon_with_poc <- read_csv("data/df_carbon_subduction_anom_with_poc.csv") %>% filter(integrated_poc > 0.01) 

# Implementing Siegel's correction : 

fseq_df <- read_csv("co2_sequestration_50years.csv")
fseq_df$fseq %>% summary()

fseq_df <- fseq_df %>% na.omit()


library(RANN)  # For fast nearest neighbor search

# 1. Prepare your data frames
df_carbon <- df_carbon_with_poc %>%
  mutate(depth_meters = critical_pres) %>%  # Approximate pressure to depth
  select(LATITUDE, LONGITUDE, depth_meters, integrated_poc,TIME)

# 2. Create 3D coordinate matrices
coords_carbon <- df_carbon %>% 
  select(LATITUDE, LONGITUDE, depth_meters) %>% 
  as.matrix()

coords_fseq <- fseq_df %>% 
  select(latitude, longitude, depth_meters) %>% 
  as.matrix()

# 3. Find nearest neighbors (takes <1 second even for large datasets)
nn <- nn2(coords_fseq, coords_carbon, k = 1)

# 4. Add the matched fseq values
df_carbon$fseq_50yr <- fseq_df$fseq[nn$nn.idx]

# 5. Calculate remaining POC
df_carbon_with_poc <- df_carbon %>%
  mutate(poc_remaining_50yr = integrated_poc * fseq_50yr) %>%
  bind_cols(df_carbon_with_poc %>% select(WMO, CYCLE_NUMBER))




# --- 1. Split the Data by Season --------------------------------------------
df_carbon_poc_djf <- df_carbon_with_poc %>% filter(month(TIME) %in% c(12, 1, 2))
df_carbon_poc_mam <- df_carbon_with_poc %>% filter(month(TIME) %in% c(3, 4, 5))
df_carbon_poc_jja <- df_carbon_with_poc %>% filter(month(TIME) %in% c(6, 7, 8))
df_carbon_poc_son <- df_carbon_with_poc %>% filter(month(TIME) %in% c(9, 10, 11))





df_carbon_with_poc <- df_carbon_with_poc %>%
  mutate(
    month_val = month(TIME),
    season = case_when(
      month_val %in% c(12, 1, 2)  ~ "DJF",
      month_val %in% c(3, 4, 5)   ~ "MAM",
      month_val %in% c(6, 7, 8)   ~ "JJA",
      month_val %in% c(9, 10, 11) ~ "SON"
    ),
    # Make it a factor in a logical order
    season = factor(season, levels = c("DJF", "MAM", "JJA", "SON"))
  )

# s(LATITUDE, LONGITUDE, by=season) allows each season to have its own spatial surface
# s(season, bs="re") is an optional random intercept per season.
m_full_st <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "tp", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_tp200 <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "tp", k = 200) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_tp600 <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "tp", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_tp1200 <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "tp", k = 1200) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_sos200 <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 200) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_sos600 <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_sos600_poc_remaining <- gam(
  log(poc_remaining_50yr) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_sos600 %>% summary()
m_full_st_sos1200 <- gam(
  log(integrated_poc) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 1200) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)


# TP k = 200 R-sq.(adj) =  0.178   Deviance explained = 21.9%
# TP k = 600 
# SOS k = 200
# SOS k = 600

m_full_st <- m_full_st_sos600

summary(m_full_st)

# Create base lat-lon grid at 1° resolution
lat_seq <- seq(-90, 89.5, by = 1) + 0.5
lon_seq <- seq(-180, 179.5, by = 1) + 0.5

# Expand for all seasons
pred_grid <- expand.grid(
  LATITUDE  = lat_seq,
  LONGITUDE = lon_seq,
  season    = c("DJF","MAM","JJA","SON")
) %>%
  mutate(season = factor(season, levels = c("DJF","MAM","JJA","SON")))

pred_grid$log_poc <- predict(m_full_st_sos600, newdata = pred_grid, type = "response")
pred_grid$poc     <- exp(pred_grid$log_poc)
pred_grid$log_poc_50yrs <- predict(m_full_st_sos600_poc_remaining, newdata = pred_grid, type = "response")
pred_grid$poc_50yrs     <- exp(pred_grid$log_poc_50yrs)



saveRDS(pred_grid,"data/pred_grid_poc_season.Rds")

library(sf)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

plot_by_season <- function(season_code) {
  # Filter to one season
  df_plot <- pred_grid %>% filter(season == season_code)
  
  ggplot() +
    geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc)) +
    geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc),
                 color = "white", alpha = 0.3) +
    geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(xlim=c(-180,180), ylim=c(-90,90), expand=FALSE, crs=st_crs(4326)) +
    scale_fill_viridis(name="Average POC [mg C/m²]") +
    labs(title = paste("Season:", season_code),
         x="Longitude", y="Latitude") +
    theme_minimal(base_size=20) +
    theme(legend.position='right', panel.grid=element_blank())
}

p_djf <- plot_by_season("DJF")
p_mam <- plot_by_season("MAM")
p_jja <- plot_by_season("JJA")
p_son <- plot_by_season("SON")

# Combine with patchwork or ggarrange
library(patchwork)
combined <- p_djf + p_mam + p_jja + p_son + 
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(title = "Average POC content of a carbon subduction event per region and season (factor-smooth model)")

ggsave("figures/TimeSpaceVar/4SEASONS/gam_poc_per_area.png",
       combined, width = 18, height = 10)

print(combined)


# Siegel Correction





# Estimating POC export :

# Carb prob gives probability that Argo float captures a subduction event per 1 degree gridcell and season
# For instance, if p = 0.20, every 5 T, with T a typical timescale, a carbon subduction event happens.
# To get a flux in (mg C m^2 / day ) f
# How to get f from carb_prob_df ? If T is 1 week, then 0.20 * 1/7
carb_prob_df <- readRDS("data/pred_grid_carb_prob.Rds")
pred_grid <- readRDS("data/pred_grid_poc_season.Rds")

carb_prob_df$LONGITUDE <- carb_prob_df$lon_bin
carb_prob_df$LATITUDE <- carb_prob_df$lat_bin

pred_grid$LATITUDE 
pred_grid$LONGITUDE

merged_df <- merge(
  carb_prob_df,
  pred_grid,
  by = c("LONGITUDE", "LATITUDE", "season")
)
# POC flux assuming T = 3 week
merged_df <- merged_df %>% mutate(poc_flux1week = poc * proportion * 1/3,
                                  poc_flux1week50yrs = poc_50yrs * proportion *1/3)

make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.1) {
  # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
  breaks_vec <- seq(prob_min, prob_max, by = binwidth)
  
  scale_fill_viridis_b(
    name = "POC flux (mg C m⁻² day⁻¹)",
    breaks = breaks_vec,
    limits = c(prob_min, prob_max),
    oob = scales::squish
  )
}

poc_scale <- make_discrete_scale(0, 10, binwidth = 1)

pred_grid$poc_50yrs
pred_grid$poc

plot_by_season_pocflux <- function(pred_grid, world_data, season_label, common_scale, argo_months) {
  # 1. Aggregate Argo float locations into 5° bins for stippling
  argo_bins <- df_argo_clean %>%
    filter(month(TIME) %in% argo_months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  full_grid <- expand.grid(
    lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
    lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
  ) %>%
    as_tibble()
  
  # Fill in any missing bins with count = 0
  argo_bins_full <- full_grid %>%
    left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # 2. Identify undersampled areas (e.g., fewer than 1 profile)
  undersampled <- argo_bins_full %>% filter(count < 1)
  
  # 3. Generate corner points for each undersampled grid cell (4 corners per cell)
  undersampled_corners <- undersampled %>%
    rowwise() %>%
    mutate(corners = list(
      data.frame(
        LON = c(lon_bin_stipple,
                lon_bin_stipple + stipple_resolution,
                lon_bin_stipple,
                lon_bin_stipple + stipple_resolution),
        LAT = c(lat_bin_stipple,
                lat_bin_stipple,
                lat_bin_stipple + stipple_resolution,
                lat_bin_stipple + stipple_resolution)
      )
    )) %>%
    ungroup() %>%
    unnest(corners)
  
  # 4. Filter the prediction grid to the desired season
  df_plot <- pred_grid %>% filter(season == season_label)
  
  # 5. Build the plot using the provided discrete scale (common_scale)
  ggplot() +
    geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc_flux1week)) +
    geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc_flux1week),
                 color = "white", alpha = 0.3,binwidth = 1) +
    geom_sf(data = world_data, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE, crs = st_crs(4326)) +
    common_scale +   # Apply the discrete scale passed as argument
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT),
               color = "white", alpha = 1, size = 1, shape = 20) +
    labs(title = paste("Average POC flux (",season_label,")",sep = ""),
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.position = 'top', 
          legend.title = element_text(size = 25, face = "bold", margin = margin(r = 100)),  # Improve legend title
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(8,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) 
}

p_djf <- plot_by_season_pocflux(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "DJF", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(12, 1, 2)
)

p_mam <- plot_by_season_pocflux(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "MAM", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(3, 4, 5)
)

p_jja <- plot_by_season_pocflux(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "JJA", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(6, 7, 8)
)

p_son <- plot_by_season_pocflux(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "SON", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(9, 10, 11)
)


combined_poc <- ggarrange(p_djf,p_mam,p_jja,p_son,ncol = 2,nrow=2,
                           common.legend = T,legend="bottom")
ggsave("figures/TimeSpaceVar/4SEASONS/poc_map.png",
       combined_poc, width = 16, height = 12)

make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.1) {
  # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
  breaks_vec <- seq(prob_min, prob_max, by = binwidth)
  
  scale_fill_viridis_b(
    name = "POC flux sequestered for > 50 years (mg C m⁻² day⁻¹)",
    breaks = breaks_vec,
    limits = c(prob_min, prob_max),
    oob = scales::squish
  )
}
poc_scale <- make_discrete_scale(0, 10, binwidth = 1)

plot_by_season_pocflux_50years <- function(pred_grid, world_data, season_label, common_scale, argo_months) {
  # 1. Aggregate Argo float locations into 5° bins for stippling
  argo_bins <- df_argo_clean %>%
    filter(month(TIME) %in% argo_months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  full_grid <- expand.grid(
    lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
    lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
  ) %>%
    as_tibble()
  
  # Fill in any missing bins with count = 0
  argo_bins_full <- full_grid %>%
    left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # 2. Identify undersampled areas (e.g., fewer than 1 profile)
  undersampled <- argo_bins_full %>% filter(count < 1)
  
  # 3. Generate corner points for each undersampled grid cell (4 corners per cell)
  undersampled_corners <- undersampled %>%
    rowwise() %>%
    mutate(corners = list(
      data.frame(
        LON = c(lon_bin_stipple,
                lon_bin_stipple + stipple_resolution,
                lon_bin_stipple,
                lon_bin_stipple + stipple_resolution),
        LAT = c(lat_bin_stipple,
                lat_bin_stipple,
                lat_bin_stipple + stipple_resolution,
                lat_bin_stipple + stipple_resolution)
      )
    )) %>%
    ungroup() %>%
    unnest(corners)
  
  # 4. Filter the prediction grid to the desired season
  df_plot <- pred_grid %>% filter(season == season_label)
  
  # 5. Build the plot using the provided discrete scale (common_scale)
  ggplot() +
    geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc_flux1week50yrs)) +
    geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc_flux1week50yrs),
                 color = "white", alpha = 0.3,binwidth = 1) +
    geom_sf(data = world_data, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE, crs = st_crs(4326)) +
    common_scale +   # Apply the discrete scale passed as argument
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT),
               color = "white", alpha = 1, size = 1, shape = 20) +
    labs(title = paste("Average POC flux (",season_label,")",sep = ""),
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.position = 'top', 
          legend.title = element_text(size = 20, face = "bold", margin = margin(r = 150)),  # Improve legend title
          panel.grid = element_blank(),
          legend.text = element_text(size = 18),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(5,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) 
}

p_djf <- plot_by_season_pocflux_50years(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "DJF", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(12, 1, 2)
)

p_mam <- plot_by_season_pocflux_50years(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "MAM", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(3, 4, 5)
)

p_jja <- plot_by_season_pocflux_50years(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "JJA", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(6, 7, 8)
)

p_son <- plot_by_season_pocflux_50years(
  pred_grid = merged_df, 
  world_data = world, 
  season_label = "SON", 
  common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
  argo_months = c(9, 10, 11)
)


combined_poc_50yrs <- ggarrange(p_djf,p_mam,p_jja,p_son,ncol = 2,nrow=2,
                          common.legend = T,legend="bottom")
ggsave("figures/TimeSpaceVar/4SEASONS/poc_map_50yrs.png",
       combined_poc_50yrs, width = 16, height = 12)


w# Estimating POC export :
# Read in the carb probability data (if not already loaded)
carb_prob_df <- readRDS("data/pred_grid_carb_prob.Rds")
carb_prob_df$LONGITUDE <- carb_prob_df$lon_bin
carb_prob_df$LATITUDE <- carb_prob_df$lat_bin

# Assuming pred_grid is already loaded and has matching LATITUDE and LONGITUDE columns
merged_df <- merge(
  carb_prob_df,
  pred_grid,
  by = c("LONGITUDE", "LATITUDE", "season")
)

# Save a "base" merged data frame (only the necessary columns)
base_merged_df <- merged_df %>% 
  select(LONGITUDE, LATITUDE, season, proportion, poc,poc_50yrs)

# Earth radius in meters (for cell area calculation)
R <- 6371000

# Define season lengths (in days); adjust if you have more accurate numbers
season_days <- data.frame(
  season = c("DJF", "MAM", "JJA", "SON"),
  days = c(90, 92, 92, 91)
)

# -----------------------------------------
# 2. Sensitivity Analysis for Different T
# -----------------------------------------
# Define T values (in days) for sensitivity analysis
# T=1 (1 day), T=3 (3 days), T=7 (1 week), T=14 (2 weeks)
T_values <- c(1:14)



# Create a results data frame to store total export for each T
sensitivity_results <- data.frame(
  T_days = T_values,
  total_export_PgC = NA_real_,
  total_export_PgC_50yrs = NA_real_
)

# Loop over each T value
for (i in seq_along(T_values)) {
  T_val <- T_values[i]
  
  # For each T, compute the daily POC flux: poc_flux = poc * proportion * (1/T)
  # (i.e., if T = 7 days then 1/7 per day, if T = 1 day then 1/1, etc.)
  sensitivity_df <- base_merged_df %>%
    mutate(poc_flux = poc * proportion * (1 / T_val),
           poc_flux_50yrs = poc_50yrs * proportion * (1 / T_val)) %>%
    # Merge the season lengths (in days)
    left_join(season_days, by = "season") %>%
    # Calculate the cell area for a 1° x 1° grid cell (m²)
    mutate(
      lat_rad = LATITUDE * pi / 180,
      cell_area = (R^2) * (pi/180)^2 * cos(lat_rad),
      # Multiply the daily flux by the number of days in the season to get mg C/m² per season
      export_season_mg_m2 = poc_flux * days,
      export_season_mg_m2_50years = poc_flux_50yrs * days,
      # Multiply by the grid cell area to get mg C per grid cell per season
      export_season_mg = export_season_mg_m2 * cell_area,
      export_season_mg_50years = export_season_mg_m2_50years * cell_area
    )
  
  # Sum over seasons to obtain annual export per grid cell
  annual_export_df <- sensitivity_df %>%
    group_by(LONGITUDE, LATITUDE) %>%
    summarise(export_annual_mg = sum(export_season_mg, na.rm = TRUE),
              export_annual_mg_50years = sum(export_season_mg_50years, na.rm = TRUE)) %>%
    ungroup()
  
  # Total global annual export in mg C
  total_export_mg <- sum(annual_export_df$export_annual_mg, na.rm = TRUE)
  total_export_mg_50years <- sum(annual_export_df$export_annual_mg_50years, na.rm = TRUE)

  # Convert mg C to Pg C (1 mg = 1e-18 Pg)
  total_export_PgC <- total_export_mg * 1e-18
  total_export_PgC_50yrs <- total_export_mg_50years * 1e-18
  
  # Store the result for this T value
  sensitivity_results$total_export_PgC[i] <- total_export_PgC
  sensitivity_results$total_export_PgC_50yrs[i] <- total_export_PgC_50yrs
}

# Print the sensitivity results
print(sensitivity_results)
sensitivity_results$total_export_PgC %>% hist()
library(ggplot2)
library(viridis)

plot_esp <- sensitivity_results %>%
  ggplot(aes(x = T_days)) +
  # Original export line
  geom_line(aes(y = total_export_PgC, color = "Export rate"), 
            size = 1.5) +
  geom_point(aes(y = total_export_PgC, color = "Export rate"), 
             size = 3) +
  # 50-year sequestered line
  geom_line(aes(y = total_export_PgC_50yrs, color = "Export rate that remains sequestered for > 50 years"), 
            size = 1.5, linetype = "dashed") +
  geom_point(aes(y = total_export_PgC_50yrs, color = "Export rate that remains sequestered for > 50 years"), 
             size = 3) +
  # Styling
  scale_color_manual(name = "",
                     values = c("Export rate" = "darkgreen",
                                "Export rate that remains sequestered for > 50 years" = "darkblue")) +
  scale_x_continuous(breaks = seq(1, 14, by = 1)) +
  labs(
    title = "Sensitivity of ESP Export Rate to Timescale T",
    x = "Characteristic Timescale (T) [Days]",
    y = "Carbon Export Rate (Pg C/year)",
    color = "Metric"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(color = "black", size = 25),
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 25),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    legend.position = "bottom"
  )

plot_esp
ggsave(plot_esp,filename = "figures/sensitivity_ESP_magnitude.png",width=15,height=10)


