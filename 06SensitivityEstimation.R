##############################
### Sensitivity Estimation ###
##############################


library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")

detected.events.list <- list()
for (j in c(1:200)) {
  try({
    
    wmo <- wmolist[j]


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
    list_of_tibbles <- sample(list_of_tibbles, 20)
    data_df <- list_of_tibbles %>% bind_rows()
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
    
    # filter(OUT.S == 1) %>% is REMOVED so all images are considered 
    
    carb_eddy.id <- B %>%  select(CYCLE_NUMBER, PRES_ADJUSTED,LATITUDE,LONGITUDE,TIME,OUT.S) %>% unique()
    

    
    carb_eddy.id <- carb_eddy.id %>%
      inner_join(mld.df, by = "CYCLE_NUMBER") %>%
      filter(PRES_ADJUSTED > MLD_DEPTH)
    # Assuming you have the required data and functions loaded
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
        geom_vline(xintercept = current_eddy$MLD_DEPTH[1], color = "red", alpha = .3, size = 1)
      
      if (  sum((current_eddy$OUT.S),na.rm = TRUE )  >= 1 ) {
        prof.plot[[i]] <- prof.plot[[i]] +
          geom_vline(xintercept = current_eddy$PRES_ADJUSTED, color = "darkgreen", alpha = .3, size = 1)
      }
      
      df <- current_data
      df.ds <- downscale_data_fun(df, b = 40)
      df.ds <- df.ds %>% ungroup() %>% select(AOU, SPIC, PRES_ADJUSTED) %>%
        pivot_longer(cols = !PRES_ADJUSTED, names_to = "VAR", values_to = "VALUE")
      
      hline_data <- data.frame(
        VAR = c("AOU", "AOU", "SPIC", "SPIC"),
        hline = c(-2, 3, -3, 3),
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
        geom_vline(xintercept = current_eddy$MLD_DEPTH[i], color = "red", alpha = .3, size = 1)
      
      if (!is.na(current_eddy$OUT.S[1]) && current_eddy$OUT.S[1] == 1) {
        res.plot[[i]] <- res.plot[[i]] +
          geom_vline(xintercept = current_eddy$PRES_ADJUSTED, color = "darkgreen", alpha = .3, size = 1)
      }
      
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
      dir <- paste0("~/Documents/GLOBARGO/figures/SensitivityEstimation/", wmo)
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

# DO Sensitivity, run algorithm 
