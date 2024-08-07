library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
library(segmented)
# Load it last
library(dplyr)
# IDEAS : downscale the residuals, not the original data BUT 
# BUT filter out profiles that are too monotonic, use segmentate ?

wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")
detected.events.list <- list()


# cycle_number <- results %>% filter(WMO == wmo) %>% select(CYCLE_NUMBER) %>% unique() %>% as_vector()

# Necessary functions : 
calculate_derivative <- function(data) {
  # Calculate the difference in VALUE and PRES_ADJUSTED
  data %>%
    arrange(PRES_ADJUSTED) %>%
    mutate(
      dVALUE = c(NA, diff(VALUE) / diff(PRES_ADJUSTED))
    )
}
# Post-processing function that check if there's indeed a peak at the detected_level : 
check_sign_change <- function(derivatives, target_level, check_depth = 50) {
  # Find the index of the target level
  target_index <- which.min(abs(derivatives$PRES_ADJUSTED - target_level))
  
  # Get indices to check around the target level
  lower_index <- which.min(abs(derivatives$PRES_ADJUSTED - (target_level - check_depth)))
  upper_index <- which.min(abs(derivatives$PRES_ADJUSTED - (target_level + check_depth)))
  
  # Check for sign changes in the derivative around the target level
  sign_changes <- sign(derivatives$dVALUE[lower_index:upper_index])
  change_detected <- any(diff(sign_changes) != 0, na.rm = TRUE)
  
  return(change_detected)
}

downscale_data_fun <- function(df, b = 20, cutoff = 1.96) {
  data <- df %>%
    dplyr::select(PRES_ADJUSTED, SCALE.RES.ROB, VAR, CYCLE_NUMBER, LONGITUDE, LATITUDE,TIME) %>%
    pivot_wider(names_from = VAR, values_from = SCALE.RES.ROB)
  
  bin_width <- b
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(from = floor(pressure_range[1] / bin_width) * bin_width,
              to = ceiling(pressure_range[2] / bin_width) * bin_width,
              by = bin_width)
  
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(-c(PRES_ADJUSTED, LATITUDE, LONGITUDE, CYCLE_NUMBER,TIME), \(x) mean(x, na.rm = TRUE)))
  
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- downscaled_data %>%
    mutate(OUT.S = ifelse(abs(AOU) > cutoff & abs(SPIC) > cutoff & AOU < 0, 1, 0))
  
  return(downscaled_data)
}


# Functions for SPIC anomaly 
find_closest_to_surface <- function(data) {
  data %>%
    filter(!is.na(SPIC)) %>%    # Filter out rows where SPIC is NA
    arrange(PRES_ADJUSTED) %>%  # Sort by pressure
    slice(1) %>%
    pull(SPIC)
}

mean_spic_at_min_max_levels <- function(data, target_pressure) {
  # Filter for rows within 100 units of the target pressure
  filtered_data <- data %>%
    filter(!is.na(SPIC)) %>%  # Filter out rows where SPIC is NA
    filter(PRES_ADJUSTED >= (target_pressure - 100) & PRES_ADJUSTED <= (target_pressure + 100))  # Filter for the range
  
  # Find the minimum and maximum pressure levels within the range
  min_pressure <- min(filtered_data$PRES_ADJUSTED, na.rm = TRUE)
  max_pressure <- max(filtered_data$PRES_ADJUSTED, na.rm = TRUE)
  
  # Filter the data for the minimum and maximum pressure levels
  spic_values <- filtered_data %>%
    filter(PRES_ADJUSTED == min_pressure | PRES_ADJUSTED == max_pressure) %>%
    pull(SPIC)  # Extract SPIC values
  
  # Calculate the mean SPIC
  mean_spic <- mean(spic_values, na.rm = TRUE)
  
  return(mean_spic)
}


for (j in seq_along(wmolist)) {
  try({
    
    wmo <- wmolist[j]
    
    # wmo <- 5904677
    # wmo <- 1902455
    # wmo <- 1902593
    
    data_df = load_float_data(float_ids = wmo,
                              variables = c("DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
                                            "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
                                            "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
                                            "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
                                            "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
                                            "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
                                            "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
                                            "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"),
                              format = "dataframe")
    
    
    data_df <- data_df %>% filter(!is.na(DOXY)) %>% group_by(CYCLE_NUMBER) %>%
      mutate(SPIC = swSpice(salinity = PSAL_ADJUSTED, temperature = TEMP_ADJUSTED,
                            latitude = first(LATITUDE), longitude = first(LONGITUDE), eos = "unesco"),
             ABS_SAL = gsw::gsw_SA_from_SP(SP = PSAL_ADJUSTED, p = PRES_ADJUSTED,
                                           longitude = first(LONGITUDE), latitude = first(LATITUDE)),
             CONS_TEMP = gsw::gsw_CT_from_t(SA = ABS_SAL, t = TEMP_ADJUSTED, p = PRES_ADJUSTED),
             SAT_DOXY = gsw_O2sol(SA = ABS_SAL, CT = CONS_TEMP, p = PRES_ADJUSTED,
                                  longitude = first(LONGITUDE), latitude = first(LATITUDE)),
             AOU = SAT_DOXY - DOXY_ADJUSTED,
             SIGMA0 = gsw::gsw_sigma0(SA = ABS_SAL, CT = CONS_TEMP))
    
    downscale_data_fun_wo_out <- function(df, b = 20) {
      data <- df %>%
        dplyr::select(PRES_ADJUSTED, AOU, SPIC, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME)
      
      bin_width <- b
      pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
      bins <- seq(from = floor(pressure_range[1] / bin_width) * bin_width,
                  to = ceiling(pressure_range[2] / bin_width) * bin_width,
                  by = bin_width)
      
      data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
      
      downscaled_data <- data %>%
        group_by(bin) %>%
        mutate(across(-c(PRES_ADJUSTED, LATITUDE, LONGITUDE, CYCLE_NUMBER, TIME), mean, na.rm = TRUE))
      
      downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
      
      downscaled_data <- unique(downscaled_data)
      return(downscaled_data)
    }
    
    list_of_tibbles <- data_df %>%
      group_by(CYCLE_NUMBER) %>%
      group_split()
    
    
    mld.vec <- c()
    for (i in seq_along(list_of_tibbles)) {
      temperature <- list_of_tibbles[[i]]$TEMP_ADJUSTED
      pressure <- list_of_tibbles[[i]]$PRES_ADJUSTED
      sigma0 <- list_of_tibbles[[i]]$SIGMA0
      
      # Calculate N2 and find the index of the maximum value
      n2 <- swN2(pressure, sigmaTheta = sigma0)
      mid <- which.max(n2)
      
      # Check if mid is not empty and assign value to mld.vec, otherwise assign NA
      if (length(mid) > 0 && !is.na(mid)) {
        mld.vec[i] <- pressure[mid]
      } else {
        mld.vec[i] <- NA
      }
    }
    
    
    mld.df <- data.frame(CYCLE_NUMBER = unique(data_df$CYCLE_NUMBER), MLD_DEPTH = mld.vec)
    
    # In this implementation, we downscale before computing residuals :
    
    downscaled_ds_list <- lapply(list_of_tibbles, downscale_data_fun_wo_out)
    df <- downscaled_ds_list %>% bind_rows()
    
    A <- df %>%  #df/data_df should be data_df if you don't want to downscale residuals
      group_by(CYCLE_NUMBER) %>%
      group_modify(~ .x %>%
                     dplyr::select(PRES_ADJUSTED, AOU, SPIC, LATITUDE, LONGITUDE, TIME) %>%
                     pivot_longer(!c(PRES_ADJUSTED, LATITUDE, LONGITUDE, TIME), names_to = "VAR", values_to = "VALUE") %>%
                     group_by(VAR) %>%
                     mutate(
                       MA_3 = rollmean(VALUE, 3, fill = NA),
                       TM_11 = rollapply(VALUE, 11,
                                         function(x) {
                                           x_subset <- x[x >= quantile(x, 0.2, na.rm = TRUE) & x <= quantile(x, 0.8, na.rm = TRUE)]
                                           if (length(x_subset) > 0) {
                                             mean(x_subset, na.rm = TRUE)
                                           } else {
                                             NA
                                           }
                                         },
                                         fill = NA),
                       TM_9 = rollapply(VALUE, 9,
                                        function(x) {
                                          x_subset <- x[x >= quantile(x, 0.2, na.rm = TRUE) & x <= quantile(x, 0.8, na.rm = TRUE)]
                                          if (length(x_subset) > 0) {
                                            mean(x_subset, na.rm = TRUE)
                                          } else {
                                            NA
                                          }
                                        },
                                        fill = NA),
                       MM_11 = rollmedian(VALUE, 11, fill = NA),
                       ROB.RES = MA_3 - TM_9,
                       ROB.RES.RAW = VALUE - TM_9
                     ) %>%
                     mutate(
                       IQRN = IQR(ROB.RES.RAW, na.rm = TRUE) / 1.349,
                       MEDIAN_RES = median(ROB.RES.RAW[ROB.RES.RAW != 0], na.rm = TRUE)
                     ) %>%
                     mutate(
                       SCALE.RES.ROB = ifelse(ROB.RES.RAW == 0, 0,
                                              (ROB.RES.RAW - MEDIAN_RES) / IQRN)
                     )
      ) %>%
      ungroup()
    
    # Just comparing trimmed mean and median behavior : 
    
    # Plot cycle number 19, float id wmo5904677 : 
    # A %>% filter(CYCLE_NUMBER==19) %>%
    #   ggplot(aes(x=PRES_ADJUSTED,y=VALUE))+facet_grid(.~VAR,scales="free")+
    #   geom_line(aes(y=TM_11,color="trimmed mean k=11"))+
    #   geom_line(aes(y=TM_9,color="trimmed mean k=9"))+
    #   geom_line(aes(y=MM_11,color="moving median"))+
    #   geom_line(aes(y=VALUE,color="observations"))+geom_point(aes(y=VALUE,color="observations"))+
    #   coord_flip()+scale_x_reverse()
    
    
    
    
    
    
    
    
    
    
    
    
    list_of_tibbles <- A %>%
      group_by(CYCLE_NUMBER) %>%
      group_split()
    
    
    # How to remove all those floats where one is clearly monotonic ? 
    output <- lapply(list_of_tibbles, downscale_data_fun, b = 40)
    output <- output %>% bind_rows()
    
    B <- output %>% bind_rows()
    
    prop_zero <- B %>% dplyr::select(AOU, SPIC, CYCLE_NUMBER, PRES_ADJUSTED,LATITUDE,LONGITUDE,TIME) %>%
      unique() %>%
      pivot_longer(cols = c("AOU", "SPIC"), names_to = "VAR", values_to = "VALUE") %>%
      group_by(CYCLE_NUMBER, VAR) %>% summarize(zero_proportion = mean(VALUE == 0, na.rm = TRUE), .groups = 'drop')
    
    pivot_proportion <- prop_zero %>%
      pivot_wider(names_from = VAR, values_from = zero_proportion)
    
    selected_cycles <- pivot_proportion %>%
      filter(AOU <= 0.5, SPIC <= 0.5) %>%
      pull(CYCLE_NUMBER)
    
    
    
    carb_eddy.id <- B %>% filter(OUT.S == 1) %>% dplyr::select(CYCLE_NUMBER, PRES_ADJUSTED,LATITUDE,LONGITUDE,TIME) %>% unique()
    
    carb_eddy.id <- carb_eddy.id %>% filter(PRES_ADJUSTED <= 700) %>% filter(PRES_ADJUSTED >= 200)
    
    
    carb_eddy.id <- carb_eddy.id %>%
      inner_join(mld.df, by = "CYCLE_NUMBER") %>%
      filter(PRES_ADJUSTED > MLD_DEPTH)
    
    const.vec <- c()
    for (i in 1:nrow(carb_eddy.id)){
      spic <- data_df %>% filter(CYCLE_NUMBER==carb_eddy.id$CYCLE_NUMBER[i]) %>% dplyr::select(SPIC,PRES_ADJUSTED) %>% ungroup()
      
      # The target pressure level
      target_pressure <- carb_eddy.id$PRES_ADJUSTED[[i]]
      

      
      
      # Find the value of spiciness at the pressure level closest to target pressure
      closest_spic <- spic %>%
        slice_min(abs(PRES_ADJUSTED - target_pressure), n = 1) %>%
        pull(SPIC) 
      
      # If spic at outlying level is closer to surface value than the mean spic is of the surface value, anomaly is consistent
      const.vec[i] <- ifelse( 
        abs(closest_spic- find_closest_to_surface(spic) )  < 
                               abs(mean_spic_at_min_max_levels(spic,target_pressure = target_pressure)-
                                     find_closest_to_surface(spic) ) ,
        1,0)
    }
    
    carb_eddy.id$CONSISTENT_ANOM <- const.vec
    
    # Split the detection dataframe in a list of dataframe with a separate dataframe for each cycle number :
    list.carb_eddy <- carb_eddy.id %>% group_by(CYCLE_NUMBER) %>%
      group_split()
    
    
    # Create an empty list to store plots
    prof.plot <- list()
    res.plot <- list()
    list.plots <- list()
    current_eddy.l <- list()
    
    
    
    
    # Iterate over each cycle number in carb_eddy.id
    for (i in seq_along(carb_eddy.id$CYCLE_NUMBER) ) {
      current_cycle <- carb_eddy.id$CYCLE_NUMBER[i]
      current_data <- A %>% filter(CYCLE_NUMBER == current_cycle)
      current_eddy <- carb_eddy.id %>% filter(CYCLE_NUMBER == current_cycle)
      pres_level <-  carb_eddy.id$PRES_ADJUSTED[i]
      
      # Filter data for AOU and SPIC
      data_aou <- current_data %>% filter(VAR == "AOU")
      data_spic <- current_data %>% filter(VAR == "SPIC")
      
      # Calculate derivatives for AOU and SPIC
      data_aou_deriv <- calculate_derivative(data_aou)
      data_spic_deriv <- calculate_derivative(data_spic)
      
      # Check for sign change around the target level (100 meters around)
      sign_change_detected_aou <- check_sign_change(data_aou_deriv, pres_level, check_depth = 100)
      sign_change_detected_spic <- check_sign_change(data_spic_deriv, pres_level, check_depth = 100)
      
      # Add sign change detection results to current_eddy
      current_eddy.l[[i]] <- current_eddy %>%
        mutate(SIGN_AOU = sign_change_detected_aou,
               SIGN_SPIC = sign_change_detected_spic)
    }
    
    current_eddy.l <- current_eddy.l %>% bind_rows()  %>% ungroup() %>% dplyr::select(-bin) %>% unique()
    # ONLY pick both SPIC and AOU change sign
    eddy_dataframe <- current_eddy.l %>%
      filter(SIGN_AOU == TRUE & SIGN_SPIC == TRUE)
    
    
    
    # Plotting profile
    for (i in seq_along(eddy_dataframe$CYCLE_NUMBER)) {
      current_cycle <- eddy_dataframe$CYCLE_NUMBER[i]
      current_data <- A %>% filter(CYCLE_NUMBER == current_cycle)
      current_eddy <- eddy_dataframe %>% filter(CYCLE_NUMBER == current_cycle)
      pres_level <-  eddy_dataframe$PRES_ADJUSTED[i]
      
      
      prof.plot[[i]] <- current_data %>%
        ggplot(aes(x = PRES_ADJUSTED, y = VALUE)) +
        facet_grid(. ~ VAR, scales = "free") +
        coord_flip() +
        scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40)) +
        geom_line(aes(y = VALUE, color = "observed values")) +
        geom_point(aes(y = VALUE, color = "observed values"), size = .3) +
        geom_point(aes(y = TM_11, color = "Trimmed Mean (k=9)"), size = .3) +
        geom_line(aes(y = TM_11, color = "Trimmed Mean (k=9)")) +
        theme_bw() +
        labs(x = "Adjusted pressure (dbar)", y = "") +
        theme(legend.position = "bottom") +
        geom_vline(xintercept = current_eddy$MLD_DEPTH[1], color = "red", alpha = .3, size = 1) +
        # This line currently displays all PRES_ADJUSTED levels where an anomaly is detected but this needs to be
        # modified so that the green lines are only shown if the SPIKE test is passed successfully. 
        geom_vline(xintercept = current_eddy$PRES_ADJUSTED, color = "darkgreen", alpha = .3, size = 1)
      
      df <- current_data
      df.ds <- downscale_data_fun(df, b = 40)
      df.ds <- df.ds %>% ungroup() %>% dplyr::select(AOU, SPIC, PRES_ADJUSTED) %>%
        pivot_longer(cols = !PRES_ADJUSTED, names_to = "VAR", values_to = "VALUE")
      
      hline_data <- data.frame(
        VAR = c("AOU", "AOU", "SPIC", "SPIC"),
        hline = c(-2, 2, -2, 2),
        label = c("-2 sigma", "+2 sigma", "-2 sigma", "+2 sigma")
      )
      
      # Plotting residuals
      res.plot[[i]] <- df %>%
        ggplot(aes(x = PRES_ADJUSTED, y = SCALE.RES.ROB)) +
        scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40)) +
        facet_grid(. ~ VAR, scales = "free") +
        coord_flip() +
        theme_bw() +
        geom_point(data = df.ds, aes(x = PRES_ADJUSTED, y = VALUE, color = "downscaled residuals to 40 meters")) +
        labs(x = "Adjusted pressure (dbar)", y = "") +
        geom_hline(data = hline_data, aes(yintercept = hline, color = label), inherit.aes = TRUE) +
        scale_color_manual(name = "Threshold", values = c("-3 sigma" = "red", "+3 sigma" = "red", "-2 sigma" = "blue", "+2 sigma" = "blue")) +
        theme(legend.position = "bottom") +
        geom_vline(xintercept = current_eddy$PRES_ADJUSTED, color = "darkgreen", alpha = .3, size = 1) +
        geom_vline(xintercept = current_eddy$MLD_DEPTH[i], color = "red", alpha = .3, size = 1)
      
      annotation_text <- paste("Cycle Number: ", current_cycle,
                               "\nFloat ID (WMO): ", wmo,
                               "\nLongitude: ", current_eddy$LONGITUDE,
                               "\nLatitude: ", current_eddy$LATITUDE,
                               "\nTime: ", format(as.POSIXct(current_eddy$TIME, origin = "1970-01-01"), "%Y-%m-%d"))
      
      # Combine plots with annotation
      combined_plot <- ggarrange(prof.plot[[i]], res.plot[[i]], common.legend = FALSE, legend = "bottom", nrow = 2)
      combined_plot <- annotate_figure(combined_plot, top = text_grob(annotation_text, face = "bold", size = 10))
      
      list.plots[[i]] <- combined_plot
    }
    
    if (length(list.plots) > 0) {
      dir <- paste0("/data/GLOBARGO/figures/EddySubductionFiguresSensSpecIncr/", wmo)
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
      
      for (k in seq_along(list.plots) ) {
        cycle_number <- eddy_dataframe$CYCLE_NUMBER[k]
        file_name <- paste0(dir, "/", wmo, "_plot_cycle_", cycle_number, ".png")
        ggsave(file_name, list.plots[[k]], width = 10, height = 10)
      }
    }
    
    eddy_dataframe$WMO <- wmo
    detected.events.list[[j]] <- eddy_dataframe
    
  }, silent = TRUE)
}


detected.events.list <- detected.events.list[!sapply(detected.events.list, is.null)]
detected.events.list <- lapply(detected.events.list, function(x) {
  x$WMO <- as.character(x$WMO)
  return(x)
})

detected.events.df <- detected.events.list %>% bind_rows()

detected.events.df %>% dplyr::select(CYCLE_NUMBER,WMO) %>% unique()



write_csv(detected.events.df, "/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")


# Note on the sensitivity increase
# If we don't downscale prior to the computation of residuals : 
#For wmo5904677 we have that carb_eddy.id$CYCLE_NUMBER %>% unique()
# [1]  4  5 13 16 18 19 20 21 23 25 26 28 32 34 35 38 40 46 47 49 50 51 53 55 57 59


# We should add a subdetection test to check that the first derivative changes sign : 

#CYCLE_NUMBER = 77

