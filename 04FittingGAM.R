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
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(viridis)
library(sp)
library(spdep)
library(mgcv)
library(patchwork)
library(gratia)

# Read data to fit the gam
df_full <- read_csv("dataframe_full_mld_and_n2.csv") 

# Standardize the log-transformed predictors
df_full <- df_full %>%
  mutate(
    log_cleaned_mld_std = scale(log(cleaned_mld)),
    log_cleaned_N2_std = scale(log(cleaned_N2))
  )


# MLD is significant (positive)
m_mld <- glm(formula = Anomaly ~ log(cleaned_mld),
         family = binomial(link = "logit"),
         data = df_full) 

m_mld %>% summary() # 0.25844 log odds for log of MLD
# Stratification is significant (negative)
m_N2 <- glm(formula = Anomaly ~ log(cleaned_N2) ,
         family = binomial(link = "logit"),
         data = df_full) 

m_N2 %>% summary() 

# Compare the 2 models
# Compare AIC values
AIC_mld <- AIC(m_mld)
AIC_N2 <- AIC(m_N2)

print(paste("AIC for MLD model:", AIC_mld))
print(paste("AIC for N2 model:", AIC_N2))
# -0.30904 for log odds of N2

# Stratification is not significant if combined with mld
m_mld_N2 <- glm(formula = Anomaly ~ log(cleaned_N2) + log(cleaned_mld) ,
         family = binomial(link = "logit"),
         data = df_full) 

m_mld_N2 %>% summary() 
# This indicates multicollinearity
m_gam <- gam(formula = Anomaly ~ s(log(cleaned_N2)) + s(log(cleaned_mld)),
             family = binomial(link = "logit"),
             data = df_full)


# Extract smooth effects for log(cleaned_mld)
smooth_mld <- smooth_estimates(m_gam, smooth = "s(log(cleaned_mld))")
# Extract smooth effects for log(cleaned_N2)
smooth_n2 <- smooth_estimates(m_gam, smooth = "s(log(cleaned_N2))")

# # Calculate the range of the smooth effect
mld_range <- max(smooth_mld$.estimate) - min(smooth_mld$.estimate)
N2_range <- max(smooth_n2$.estimate) - min(smooth_n2$.estimate)

# Print the ranges
cat("Effect size for log(cleaned_mld):", mld_range, "\n")
cat("Effect size for log(cleaned_N2):", N2_range, "\n")

# Back transform
smooth_mld$mld <- exp(smooth_mld$`log(cleaned_mld)`)  # If cleaned_mld was log-transformed earlier
smooth_n2$N2 <- exp(smooth_n2$`log(cleaned_N2)`)     # If cleaned_N2 was log-transformed earlier

# Transform the smooth effects for cleaned_mld to probabilities

smooth_mld <- smooth_mld %>%
  mutate(
    prob = 1 / (1 + exp(-.estimate)),        # Transform to probability
    prob_lower = 1 / (1 + exp(-(.estimate - .se))),  # Lower bound of confidence interval
    prob_upper = 1 / (1 + exp(-(.estimate + .se)))   # Upper bound of confidence interval
  )

# Transform the smooth effects for N2 to probabilities
smooth_n2 <- smooth_n2 %>%
  mutate(
    prob = 1 / (1 + exp(-.estimate)),        # Transform to probability
    prob_lower = 1 / (1 + exp(-(.estimate - .se))),  # Lower bound of confidence interval
    prob_upper = 1 / (1 + exp(-(.estimate + .se)))   # Upper bound of confidence interval
  )

library(ggplot2)

mld_effect <- ggplot(smooth_mld, aes(x = mld, y = prob)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), alpha = 0.2, fill = "blue") +
  labs(
    title = "Effect of MLD on Subduction Probability",
    x = "Mixed Layer Depth (MLD)",
    y = "Predicted Probability"
  ) +
  theme_minimal()+xlim(0,500)


n2_effect <- ggplot(smooth_n2, aes(x = N2, y = prob)) +
  geom_line(color = "red", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), alpha = 0.2, fill = "red") +
  labs(
    title = "Effect of Stratification (N²) on Subduction Probability over Likely Range of N2 (0 to 0.01)",
    x = "Brunt-Väisälä Frequency (N²) log10 scale",
    y = "Predicted Probability"
  ) +
  theme_minimal()+scale_x_log10()+xlim(0,0.01)

# Combine the two plots side by side
combined_plot <- mld_effect + n2_effect + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save the combined plot
ggsave("/data/GLOBARGO/figures/GAM_combined_effects_plot.png", combined_plot, width = 12, height = 6, dpi = 300)

m_gam %>% summary
# R-sq.(adj) =  0.0103   Deviance explained = 2.53%, so clearly, not much if just using MLD and N2


# Add spatial effect
gam_model_spatial <- bam(
  Anomaly ~  s(log(cleaned_N2)) + s(log(cleaned_mld)) + s(LATITUDE, LONGITUDE, bs = "sos", k = 5),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)

gam_model_spatial %>% summary()
# Deviance explained = 19.4%

gam_model_spatial_temporal <- gam(
  Anomaly ~  s(log(cleaned_N2),k=15) + s(log(cleaned_mld),k=15) + s(LATITUDE, LONGITUDE, bs = "sos", k = 600)+
    s(AdjustedDayOfYear, bs = 'cc', k = 20),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)
gam_model_spatial_temporal %>% summary()

# Deviance explained = 19.5%

gam_model_spatial_temporal_int <- gam(
  Anomaly ~  s(log(cleaned_N2)) + s(log(cleaned_mld)) + s(LATITUDE, LONGITUDE, bs = "sos", k = 600)+
    s(AdjustedDayOfYear, bs = 'cc', k = 20) + ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(20,20)),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)


# Add seasonal effect :
# Create Season factor
df_full$Season <- cut(
  df_full$AdjustedDayOfYear,
  breaks = c(0, 90, 180, 270, 365),
  labels = c("Winter", "Spring", "Summer", "Fall"),
  include.lowest = TRUE
)

gam_model_spatial_temporal_int_hemi_season <- gam(
  Anomaly ~ Hemisphere + Season + s(log(cleaned_N2)) + s(log(cleaned_mld)) + s(LATITUDE, LONGITUDE, bs = "sos", k = 600)+
    s(AdjustedDayOfYear, bs = 'cc', k = 20) + ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(20,20)),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)




# Define grid size to approximately 25 km
# Define bin size
bin_size <- 0.25

# Define longitude and latitude bins
longitude_bins <- seq(floor(min(df_full$LONGITUDE, na.rm = TRUE)), ceiling(max(df_full$LONGITUDE, na.rm = TRUE)), by = bin_size)
latitude_bins <- seq(floor(min(df_full$LATITUDE, na.rm = TRUE)), ceiling(max(df_full$LATITUDE, na.rm = TRUE)), by = bin_size)

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

df_full_binned <- assign_bins(df_full)



# Optionally, average predictors within grid cells and months
df_agg <- df_full_binned %>%
  group_by(lat_bin, lon_bin, Month) %>%
  summarise(
    Anomaly = sum(Anomaly),
    Total =  n(),
    cleaned_mld = mean(cleaned_mld),
    cleaned_N2 = mean(cleaned_N2),
    LATITUDE = mean(LATITUDE),
    LONGITUDE = mean(LONGITUDE),
  ) %>%
  ungroup()

df_agg <- df_agg %>%
  mutate(
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "DJF",
      Month %in% c(3, 4, 5) ~ "MAM",
      Month %in% c(6, 7, 8) ~ "JJA",
      Month %in% c(9, 10, 11) ~ "SON"
    ),
    Season = factor(Season, levels = c("DJF", "MAM", "JJA", "SON"))
  )

gam_model <- gam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 500) +
    Season,
  family = binomial,
  data = df_agg
)
# R-sq.(adj) =  0.271   Deviance explained = 42.3%
gam_model %>% summary()

write_rds(gam_model,file = "models/MLD_LAT_LON_Season.Rds")

gam_model <- read_rds(file = "models/MLD_LAT_LON_Season.Rds")
gam_model %>% gam.check()

gam_model_weights_high_k <- gam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 200) +
    Season,
  family = binomial,
  weights = Total,
  data = df_agg
)
write_rds(gam_model_weights_high_k,file = "gam_model_weights_high_k.Rds")


gam_model_weights_high_k <- bam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 1000) +
    Season,
  family = binomial,
  weights = Total,
  data = df_agg
)


write_rds(gam_model_weights_high_k,file = "gam_model_weights_high_k.Rds")



gam_model_weights_high_k <- NULL 

gam_model_weights_high_k_int <- bam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(20,20))+
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 1000) +
    Season,
  family = binomial,
  weights = Total,
  data = df_agg
)

write_rds(gam_model_weights_high_k_int,file = "high_int_mod.Rds")



# Compare AIC values
AIC(gam_model_no_bin, gam_model_bin)

# Plot for non-binned model

plot(gam_model_no_bin, scheme = 2)

# Plot for binned model
plot(gam_model_bin, scheme = 2)



# 2 Spatial Autocorrelation for model without binning

library(spdep)

# Create spatial weights matrix
coordinates <- cbind(df_full_djf$LONGITUDE, df_full_djf$LATITUDE)
neighbors <- knearneigh(coordinates, k = 8,longlat = TRUE)
weights <- nb2listw(knn2nb(neighbors))

# Extract residuals
residuals_no_bin <- residuals(gam_model_no_bin, type = "pearson")

# Moran's I test
moran_test <- moran.test(residuals_no_bin, weights)
print(moran_test)

# Fit a GAMM with spatial correlation
gamm_model <- gamm(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = df_full_djf,
  correlation = corExp(form = ~ LONGITUDE + LATITUDE),
  method = "REML"
)
# 2 Spatial Autocorrelation for model WITH binning

library(spdep)

# Create spatial weights matrix
coordinates <- cbind(df_binned$LON_BIN, df_binned$LAT_BIN)
neighbors <- knearneigh(coordinates, k = 8,longlat = TRUE)
weights <- nb2listw(knn2nb(neighbors))

# Extract residuals
residuals_no_bin <- residuals(gam_model_no_bin, type = "pearson")

# Moran's I test
moran_test <- moran.test(residuals_no_bin, weights)
print(moran_test)

# Fit a GAMM with spatial correlation
gamm_model <- gamm(
  Anomaly ~ s(LAT_BIN, LON_BIN, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = df_binned,
  correlation = corExp(form = ~ LON_BIN + LAT_BIN),
  method = "REML"
)
