# Load necessary libraries
library(dplyr)
library(gsw)
library(oce)
library(robustbase)
# Read the list of WMO IDs
wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")

# Initialize a data frame to store results
mld_results <- data.frame(WMO = integer(), Latitude = numeric(), Longitude = numeric(), Time = as.POSIXct(character()), MLD = numeric())

# Main processing loop over each WMO ID
for (j in seq_along(wmolist)) {
  try({
    # Get current WMO ID
    wmo <- wmolist[j]
    
    # Load data for the current float (assuming 'load_float_data' is an internal function)
    float_data <- load_float_data(
      float_ids = wmo,
      variables = c(
        "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
        "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
        "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
        "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
        "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
        "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
        "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
        "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
      ),
      format = "dataframe"
    )
    
    # Preprocess data: calculate additional oceanographic parameters
    float_data <- float_data %>%
      filter(!is.na(PRES_ADJUSTED) & !is.na(PSAL_ADJUSTED) & !is.na(TEMP_ADJUSTED)) %>%
      group_by(CYCLE_NUMBER) %>%
      mutate(
        ABS_SAL = gsw::gsw_SA_from_SP(
          SP = PSAL_ADJUSTED,
          p = PRES_ADJUSTED,
          longitude = first(LONGITUDE),
          latitude = first(LATITUDE)
        ),
        CONS_TEMP = gsw::gsw_CT_from_t(
          SA = ABS_SAL,
          t = TEMP_ADJUSTED,
          p = PRES_ADJUSTED
        ),
        SIGMA0 = gsw::gsw_sigma0(ABS_SAL, CONS_TEMP)  # Calculate density anomaly
      )
    
    # Define density threshold for MLD calculation
    density_threshold <- 0.03  # kg/m³ increase from the surface density
    
    # Loop over each profile (cycle) to compute MLD
    for (cycle_data in split(float_data, float_data$CYCLE_NUMBER)) {
      # Surface density
      surface_density <- cycle_data$SIGMA0[1]
      
      # Find the depth where density exceeds the surface density by the threshold
      mld_index <- which(cycle_data$SIGMA0 >= surface_density + density_threshold)[1]
      
      # Determine MLD based on the threshold, or set to NA if no MLD found
      mld <- if (!is.na(mld_index)) cycle_data$PRES_ADJUSTED[mld_index] else NA
      
      # Store WMO, Latitude, Longitude, Time, and MLD for the profile
      mld_results <- rbind(mld_results, data.frame(
        WMO = wmo,
        Latitude = first(cycle_data$LATITUDE),
        Longitude = first(cycle_data$LONGITUDE),
        Time = first(cycle_data$TIME),
        MLD = mld,
        CYCLE_NUMBER = first(cycle_data$CYCLE_NUMBER)
      ))
    }
  }, silent = TRUE)  # Continue to the next float if there's an error
}

# Display or save the results
print(head(mld_results))

mld_results$MLD %>% hist(bins=100)



# Calculate the lower and upper percentiles
lower_bound <- quantile(mld_results$MLD, 0.001, na.rm = TRUE)
upper_bound <- quantile(mld_results$MLD, 0.999, na.rm = TRUE)



# Create a new column 'cleaned_mld' that sets outliers to NA
mld_results <- mld_results %>%
  mutate(cleaned_mld = ifelse(MLD >= lower_bound & MLD <= upper_bound, MLD, NA))

# Define bins for 'binned_mld' column following  https://doi.org/10.1029/2004JC002378
breaks <- c(0,10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 400, 500, Inf)
labels <- c("0-10","10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", 
            "90-100", "100-125", "125-150", "150-200", "200-300", "300-400", "400-500", "500+")
# Create the 'binned_mld' column using cut
mld_results <- mld_results %>%
  mutate(binned_mld = cut(MLD, breaks = breaks, labels = labels, right = FALSE))

# Save to file
write.csv(mld_results, "mld_results.csv", row.names = FALSE)

mld_results <- read.csv(mld_results, "mld_results.csv", row.names = FALSE)


# What about change in MLD? Assuming Argo are Lagrangian, we can compute the difference 
# between the 3 bin centered median MLD  in 
# the Argo floats and the median of the last 5 preceding cycle numbers timesteps, to know if MLD is restatifying or not.
mld_results %>% filter(WMO == 1902303) %>% View()
