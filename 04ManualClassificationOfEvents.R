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

df_partial_salinity_class$WMO <- as.numeric(df_partial_salinity_class$WMO)
df_spic_class$WMO             <- as.numeric(df_spic_class$WMO)    

# Join two datasets, assuming that if a row is common to spic and partial_salinity in terms of their unique
# combination of WMO and Cycle, partial_salinity should be prefered :


# Perform the join, 
df_combined <- bind_rows(df_partial_salinity_class, df_spic_class) %>% select(WMO,CYCLE_NUMBER,Category)


# Next we wish to complete this dataframe by retriving
# LON/LAT/TIME info for each of those individual subduction events, df_abs_sal has the required info


# Check for duplicates in df_abs_sal based on WMO and CYCLE_NUMBER
df_abs_sal_unique <- df_abs_sal %>%
  distinct(WMO, CYCLE_NUMBER, .keep_all = TRUE)


# Perform the join again with the unique rows
df_complete <- df_combined %>%
  left_join(df_abs_sal_unique %>% 
              select(WMO, CYCLE_NUMBER, LATITUDE, LONGITUDE, TIME), 
            by = c("WMO", "CYCLE_NUMBER"))

# Clean subduction events data
df_complete_clean <- df_complete %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))

df_argo_clean <- df_argo %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))


write_csv(df_complete_clean,file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
write_csv(df_argo_clean, file = "/data/GLOBARGO/src/data/df_argo_loc.csv")
