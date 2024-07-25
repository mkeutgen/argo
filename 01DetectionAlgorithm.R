library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")
detected.events.list <- list()
for (j in seq_along(wmolist)) {
  try({
    
    wmo <- wmolist[j]
    #wmo <- 5904677
    # cycle_number <- results %>% filter(WMO == wmo) %>% select(CYCLE_NUMBER) %>% unique() %>% as_vector()
    
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
        select(PRES_ADJUSTED, AOU, SPIC, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME)
      
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
    
    downscaled_ds_list <- lapply(list_of_tibbles, downscale_data_fun_wo_out)
    df <- downscaled_ds_list %>% bind_rows()
    
    A <- df %>%
      group_by(CYCLE_NUMBER) %>%
      group_modify(~ .x %>%
                     select(PRES_ADJUSTED, AOU, SPIC, LATITUDE, LONGITUDE, TIME) %>%
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
                       ROB.RES = MA_3 - TM_11,
                       ROB.RES.RAW = VALUE - TM_11
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
    
    list_of_tibbles <- A %>%
      group_by(CYCLE_NUMBER) %>%
      group_split()
    
    
    downscale_data_fun <- function(df, b = 20, cutoff = 3.3) {
      data <- df %>%
        select(PRES_ADJUSTED, SCALE.RES.ROB, VAR, CYCLE_NUMBER, LONGITUDE, LATITUDE,TIME) %>%
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
    
    output <- lapply(list_of_tibbles, downscale_data_fun, b = 40)
    output <- output %>% bind_rows()
    
    B <- output %>% bind_rows()
    
    prop_zero <- B %>% select(AOU, SPIC, CYCLE_NUMBER, PRES_ADJUSTED,LATITUDE,LONGITUDE,TIME) %>%
      unique() %>%
      pivot_longer(cols = c("AOU", "SPIC"), names_to = "VAR", values_to = "VALUE") %>%
      group_by(CYCLE_NUMBER, VAR) %>% summarize(zero_proportion = mean(VALUE == 0, na.rm = TRUE), .groups = 'drop')
    
    pivot_proportion <- prop_zero %>%
      pivot_wider(names_from = VAR, values_from = zero_proportion)
    
    selected_cycles <- pivot_proportion %>%
      filter(AOU <= 0.5, SPIC <= 0.5) %>%
      pull(CYCLE_NUMBER)
    
    
    
    carb_eddy.id <- B %>% filter(OUT.S == 1) %>% select(CYCLE_NUMBER, PRES_ADJUSTED,LATITUDE,LONGITUDE,TIME) %>% unique()
    
    carb_eddy.id <- carb_eddy.id %>% filter(PRES_ADJUSTED <= 700) %>% filter(PRES_ADJUSTED >= 200)
    
    
    carb_eddy.id <- carb_eddy.id %>%
      inner_join(mld.df, by = "CYCLE_NUMBER") %>%
      filter(PRES_ADJUSTED > MLD_DEPTH)
    
    const.vec <- c()
    for (i in 1:nrow(carb_eddy.id)){
      spic <- data_df %>% filter(CYCLE_NUMBER==carb_eddy.id$CYCLE_NUMBER[i]) %>% select(SPIC,PRES_ADJUSTED) %>% ungroup()
      
      # The target pressure level
      target_pressure <- carb_eddy.id$PRES_ADJUSTED[[i]]
      
      # Find the value of spiciness at the pressure level closest to 700
      closest_spic <- spic %>%
        slice_min(abs(PRES_ADJUSTED - target_pressure), n = 1) %>%
        select(SPIC) %>% ungroup()
      
      # If spic at outlying level is closer to surface value than the mean spic is of the surface value, anomaly is consistent
      const.vec[i] <- ifelse(abs(closest_spic-spic$SPIC[1]) < abs(mean(spic$SPIC)-spic$SPIC[1]),1,0)
    }
    carb_eddy.id$CONSISTENT_ANOM <- const.vec
    
    # in this case it's not consistent
    prof.plot <- list()
    res.plot <- list()
    list.plots <- list()
    for (i in seq_along(carb_eddy.id$CYCLE_NUMBER)) {
      current_cycle <- carb_eddy.id$CYCLE_NUMBER[i]
      current_data <- A %>% filter(CYCLE_NUMBER == current_cycle)
      current_eddy <- carb_eddy.id %>% filter(CYCLE_NUMBER == current_cycle)
      
      prof.plot[[i]] <- current_data %>%
        ggplot(aes(x = PRES_ADJUSTED, y = VALUE)) +
        facet_grid(. ~ VAR, scales = "free") +
        coord_flip() +
        scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40)) +
        geom_line(aes(y = VALUE, color = "observed values")) +
        geom_point(aes(y = VALUE, color = "observed values"), size = .3) +
        geom_point(aes(y = TM_11, color = "Trimmed Mean (k=11)"), size = .3) +
        geom_line(aes(y = TM_11, color = "Trimmed Mean (k=11)")) +
        theme_bw() +
        labs(x = "Adjusted pressure (dbar)", y = "") +
        theme(legend.position = "bottom") +
        geom_vline(xintercept = current_eddy$MLD_DEPTH[1], color = "red", alpha = .3, size = 1) +
        geom_vline(xintercept = current_eddy$PRES_ADJUSTED, color = "darkgreen", alpha = .3, size = 1)
      
      df <- current_data
      df.ds <- downscale_data_fun(df, b = 40)
      df.ds <- df.ds %>% ungroup() %>% select(AOU, SPIC, PRES_ADJUSTED) %>%
        pivot_longer(cols = !PRES_ADJUSTED, names_to = "VAR", values_to = "VALUE")
      
      hline_data <- data.frame(
        VAR = c("AOU", "AOU", "SPIC", "SPIC"),
        hline = c(-3, 3, -3, 3),
        label = c("-3 sigma", "+3 sigma", "-3 sigma", "+3 sigma")
      )
      
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
      
      combined_plot <- ggarrange(prof.plot[[i]], res.plot[[i]], common.legend = F, legend = "bottom", nrow = 2)
      combined_plot <- annotate_figure(combined_plot, top = text_grob(annotation_text, face = "bold", size = 10))
      
      list.plots[[i]] <- combined_plot
    }
    
    if (length(list.plots) > 0) {
      dir <- paste0("~/Documents/ARGO/EddySubductionFigures/", wmo)
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
      
      for (k in seq_along(list.plots)) {
        cycle_number <- carb_eddy.id$CYCLE_NUMBER[k]
        file_name <- paste0(dir, "/", wmo, "_plot_cycle_", cycle_number, ".png")
        ggsave(file_name, list.plots[[k]], width = 10, height = 10)
      }
    }
    
    carb_eddy.id$WMO <- wmo
    detected.events.list[[j]] <- carb_eddy.id
    
  }, silent = TRUE)
}


detected.events.list <- detected.events.list[!sapply(detected.events.list, is.null)]
detected.events.list <- lapply(detected.events.list, function(x) {
  x$WMO <- as.character(x$WMO)
  return(x)
})

detected.events.df <- detected.events.list %>% bind_rows()

detected.events.df %>% select(CYCLE_NUMBER,WMO) %>% unique()



write_csv(detected.events.df, "~/Documents/GLOBARGO/data/detected_events.csv")

# Finding out which profiles were already classified
class.df.part <- read_csv("~/Documents/GLOBARGO/data/classification_results_portion_data.csv")
class.df.part <- class.df.part %>%  mutate(WMO = str_replace(WMO, "_plot$", ""))
class.df.part$CYCLE_NUMBER <- class.df.part$Cycle

class.df.part %>% select(CYCLE_NUMBER,WMO)
detected.events.df %>% select(CYCLE_NUMBER,WMO)

matched_rows <- class.df.part %>%
  semi_join(detected.events.df, by = c("CYCLE_NUMBER", "WMO"))

write_csv(matched_rows,"~/Documents/GLOBARGO/data/classification_results_portion_data.csv")
       
# Manually classify data with Python script 05-ClassificationV2.ipynb :

detected.events.df <- read_csv("~/Documents/GLOBARGO/data/detected_events.csv")

manual.class.df <- read_csv("~/Documents/GLOBARGO/data/classification_results_manually_classified.csv")
manual.class.df <- manual.class.df %>%  mutate(WMO = str_replace(WMO, "_plot$", ""))
manual.class.df$CYCLE_NUMBER <- manual.class.df$Cycle

# Find if there are outlying profiles 
# Function to check if a group contains only category 0 or 4
only_category_0_3_or_4 <- function(categories) {
  all(categories %in% c(0, 4,3))
}

# Find WMOs with only category 0 or 4
wm_only_0_3_or_4 <- manual.class.df %>%
  group_by(WMO) %>%
  filter(only_category_0_3_or_4(Category)) %>%
  summarize(unique_categories = unique(Category))

# Display the result
print(wm_only_0_or_4)









# Function to prioritize the category values
prioritize_category <- function(categories) {
  if (4 %in% categories) {
    return(4)
  } else if (1 %in% categories) {
    return(1)
  } else if (2 %in% categories) {
    return(2)
  } else if (3 %in% categories) {
    return(3)
  } else {
    return(0)
  }
}

# Find WMOs with only category 0 or 4
anomalous_wmos <- manual.class.df %>%
  group_by(WMO) %>%
  filter(only_category_0_3_or_4(Category)) %>%
  distinct(WMO)

# Remove anomalous floats and apply the function to the data frame
manual.class.df_filtered <- manual.class.df %>%
  filter(!WMO %in% anomalous_wmos$WMO) %>%
  group_by(WMO, Cycle) %>%
  summarize(Category = prioritize_category(Category), .groups = 'drop') %>%
  filter(Category != 4)

# Calculate the proportion of each category
category_proportions <- manual.class.df_filtered %>%
  group_by(Category) %>%
  summarize(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

# Display the result
print(category_proportions)

# Remove anomalous floats and apply the function to the data frame
manual.class.df_filtered <- manual.class.df %>%
  filter(!WMO %in% anomalous_wmos$WMO) %>%
  group_by(WMO, Cycle) %>%
  summarize(Category = prioritize_category(Category), .groups = 'drop') %>%
  filter(Category != 4)

# Calculate the proportion of each category
category_proportions <- manual.class.df_filtered %>%
  group_by(Category) %>%
  summarize(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))



# Combine both dataframes
manual.class.df_filtered <- manual.class.df_filtered %>%
  rename(CYCLE_NUMBER = Cycle)

# Convert the WMO column in detected.events.df to character type
detected.events.df <- detected.events.df %>%
  mutate(WMO = as.character(WMO))


# Perform the full join
combined_df <- full_join(manual.class.df_filtered, detected.events.df, by = c("WMO", "CYCLE_NUMBER"))

# Display the result
combined_df_filtered <- combined_df %>%
  filter(!is.na(CONSISTENT_ANOM))


write_csv(combined_df_filtered,file = "~/Documents/GLOBARGO/data/classification_results_manually_classified_merged.csv")

# Subduction  events df :
df <- combined_df_filtered %>% filter(Category %in% c(1,2)) 
write_csv(df,"~/Documents/GLOBARGO/data/subduction_events.csv")
combined_df_filtered <- combined_df_filtered %>% filter(Category %in% c(1,2,3,0))
category_proportions_final <- combined_df_filtered %>%
  group_by(Category) %>%
  summarize(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))
# A tibble: 4 × 3
### Category Count Proportion
###
### 0   1977    0.497 
### 1   852     0.214 
### 2   958     0.241 
### 3   194     0.0487

