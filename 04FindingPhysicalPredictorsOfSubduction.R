# Load necessary libraries
library(dplyr)
library(gsw)
library(oce)

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
        MLD = mld
      ))
    }
  }, silent = TRUE)  # Continue to the next float if there's an error
}

# Display or save the results
print(head(mld_results))
# Optionally save to a file
write.csv(mld_results, "mld_results.csv", row.names = FALSE)
