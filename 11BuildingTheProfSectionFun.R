# Testing with Llort Float
# TWO FUNCTIONS, ONE TO BUILD THE DATAFRAME WITH SPIC AOU AND BBP
# ONE TO PLOT

WMO <- 5904677




unique_sub_loc <- read_csv("/data/GLOBARGO/data/unique_sub_loc.csv")
det_event_cat1 <- read_csv("/data/GLOBARGO/data/detected_events_unique_with_carbon_cat1.csv")

det_event_cat1 %>% filter(WMO == 5904677)

anomalies_det <- det_event_cat1 %>% filter(WMO == 5904677) %>% filter(CYCLE_NUMBER %in% c(25,27)) %>% select(PRES_ADJUSTED,TIME)
anomalies_det$TIME <- as.POSIXct(anomalies_det$TIME)

variables <- c("DOXY","BBP_700","PSAL","NITRATE")
float_profs <- NULL
plot_isopyc=1
plot_mld=1
raw="no"
obs="off"
qc_flags=0:9

# Set the max depth to 1000 meters

max_depth <- 1000

# Downscaling function
downscale_data_fun_wo_out <- function(df,b=20) {
  # Select and pivot data
  data <- df %>%
    select(PRES_ADJUSTED, AOU,BBP700_ADJUSTED, SPIC,CYCLE_NUMBER,LONGITUDE,LATITUDE,TIME,BBP700_ADJUSTED)
  
  # Determine bin edges
  bin_width <- b
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(from = floor(pressure_range[1]/bin_width)*bin_width,
              to = ceiling(pressure_range[2]/bin_width)*bin_width,
              by = bin_width)
  
  # Assign data points to bins
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  # Calculate mean for each bin
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(-c(PRES_ADJUSTED,LATITUDE,LONGITUDE,CYCLE_NUMBER,TIME), mean, na.rm = TRUE))
  
  # Add bin's central value
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- unique(downscaled_data)
  # Return downscaled data
  return(downscaled_data)
}


interpol_df_fun <- function(WMO){
    
    
    # make sure Setting is initialized
    if (exists("Setting")==F) {
      initialize_argo()
    }
    
    # download Sprof files if necessary
    good_float_ids = download_multi_floats(float_ids)
    
    loaded = load_float_data(good_float_ids,variables = NULL, float_profs)
    Data = loaded$Data
    Mdata = loaded$Mdata

    
    
    
    # vertical interpolation to depths with regular intervals
    
      prs_res = 20
      
      Datai = depth_interp(Data[[floats[f]]], qc_flags, calc_dens=calc_dens, 
                           calc_mld_temp=(plot_mld==1), calc_mld_dens=(plot_mld==2),
                           prs_res=prs_res)
      
    # computation of spiciness (SPIC) and apparent oxygen utilization (AOU)
    # based on adjusted variables if available 
      
      SPIC <- matrix(NA, nrow = nrow(Datai$PRES), ncol = ncol(Datai$PRES))
      
      # Loop over each profile (column in the matrix)
      for (i in 1:ncol(Datai$PRES_ADJUSTED)) {
        
        # Extract salinity, temperature, latitude, and longitude for the current profile
        salinity <- if (all(is.na(Datai$PSAL_ADJUSTED[, i]) | is.nan(Datai$PSAL_ADJUSTED[, i])) ) {
          Datai$PSAL[, i]
        } else {
          Datai$PSAL_ADJUSTED[, i]
        }
        
        temperature <- if (all(is.na(Datai$TEMP_ADJUSTED[, i]) | is.nan(Datai$TEMP_ADJUSTED[, i]) ) ) {
          Datai$TEMP[, i]
        } else {
          Datai$TEMP_ADJUSTED[, i]
        }
    
        latitude <- Datai$LATITUDE[1, i]  # Latitude is constant for a given profile
        longitude <- Datai$LONGITUDE[1, i]  # Longitude is constant for a given profile
        
        
        # Compute spiciness for each depth level in the profile
        for (j in 1:nrow(Datai$PRES_ADJUSTED)) {
          SPIC[j, i] <- swSpice(salinity = salinity[j], temperature = temperature[j],
                                latitude = latitude, longitude = longitude, eos = "unesco")
        }
      }
      
      # Add SPIC to Datai list
      Datai$SPIC <- SPIC
      
      
      # Initialize an empty matrix to store AOU values
      AOU <- matrix(NA, nrow = nrow(Datai$PRES_ADJUSTED), ncol = ncol(Datai$PRES_ADJUSTED))
      
      # Loop over each profile (column in the matrix)
      for (i in 1:ncol(Datai$PRES_ADJUSTED)) {
        
        # Conditional handling of adjusted vs non-adjusted variables
        salinity <- if (all(is.na(Datai$PSAL_ADJUSTED[, i]) | is.nan(Datai$PSAL_ADJUSTED[, i])) ) {
          Datai$PSAL[, i]
        } else {
          Datai$PSAL_ADJUSTED[, i]
        }
        
        temperature <- if (all(is.na(Datai$TEMP_ADJUSTED[, i]) | is.nan(Datai$TEMP_ADJUSTED[, i]))) {
          Datai$TEMP[, i]
        } else {
          Datai$TEMP_ADJUSTED[, i]
        }
        
        pressure <- if (all(is.na(Datai$PRES_ADJUSTED[, i]) | is.nan(Datai$PRES_ADJUSTED[, i]))) {
          Datai$PRES[, i]
        } else {
          Datai$PRES_ADJUSTED[, i]
        }
        
        doxy <- if (all(is.na(Datai$DOXY_ADJUSTED[, i]) | is.nan(Datai$DOXY_ADJUSTED[, i]))) {
          Datai$DOXY[, i]
        } else {
          Datai$DOXY_ADJUSTED[, i]
        }
        
        latitude <- Datai$LATITUDE[1, i]  # Latitude is constant for a given profile
        longitude <- Datai$LONGITUDE[1, i]  # Longitude is constant for a given profile
      
        
        # Compute Absolute Salinity (ABS_SAL) and Conservative Temperature (CONS_TEMP) for each depth level
        ABS_SAL <- gsw_SA_from_SP(SP = salinity, p = pressure, longitude = longitude, latitude = latitude)
        CONS_TEMP <- gsw_CT_from_t(SA = ABS_SAL, t = temperature, p = pressure)
        
        # Compute Saturation Dissolved Oxygen (SAT_DOXY)
        SAT_DOXY <- gsw_O2sol(SA = ABS_SAL, CT = CONS_TEMP, p = pressure, longitude = longitude, latitude = latitude)
        
        # Compute AOU for each depth level
        AOU[, i] <- SAT_DOXY - doxy
      }
      
      # Add AOU to Datai list
      Datai$AOU <- AOU
      

      # Convert time to universal time coordinates UTC
      
      df = NULL
      for (name in names(Datai)) {
        if (name == "TIME") {
          df[[name]] = as.POSIXct(Datai[[name]], tz="UTC")
        } else {
          df[[name]] = as.vector(Datai[[name]])
        }
      }
      
      df = data.frame(df)
      
      # Downscaling AOU SPIC and BBP700 to uniform 20 meters vertical res 
      list_of_tibbles <- df  %>%
        group_by(CYCLE_NUMBER) %>%
        group_split()
      
      list_of_tibbles <- lapply(list_of_tibbles, function(df) {
        df$BBP700_ADJUSTEDORIG <- df$BBP700_ADJUSTED
        df$BBP700_ADJUSTED <- rollmedian(df$BBP700_ADJUSTED, k = 3, fill = NA)
        return(df)
      })    
      # Downscale data to a uniform resolution of 20 meters then compute residuals at 20 meters
      
      downscaled_ds_list <- lapply(list_of_tibbles,downscale_data_fun_wo_out)
      
      df <- downscaled_ds_list %>% bind_rows()
      
      
      

      
      # Create min/max parameters to use the geom_rect() function that can plot
      # rectangles with variable width/height
      
      if ("PRES_ADJUSTED" %in% names(Datai) ) {
        df$ymin = df$PRES_ADJUSTED - prs_res/2
        df$ymax = df$PRES_ADJUSTED + prs_res/2
      } else {
        df$ymin = df$PRES - prs_res/2
        df$ymax = df$PRES + prs_res/2
      }
      # Not NA time
      nna = which(!is.na(Datai$TIME[1,]))
      
      # Number of not na
      nc = length(nna)
      
      xvec = as.POSIXct(Datai$TIME[1,nna], tz="UTC")
      xmin = as.POSIXct(rep(NA, nc), tz="UTC")
      xmax = as.POSIXct(rep(NA, nc), tz="UTC")
      
      xmin[2:(nc-1)] = xvec[2:(nc-1)] - ( xvec[2:(nc-1)] - xvec[1:(nc-2)] ) / 2
      xmax[2:(nc-1)] = xvec[2:(nc-1)] + ( xvec[3:nc] - xvec[2:(nc-1)] ) / 2
      
      xmin[1] = xvec[1] - ( xvec[2] - xvec[1] ) / 2
      xmax[nc] = xvec[nc] + ( xvec[nc] - xvec[nc-1] ) / 2
      xmin[nc] = xmax[nc-1]
      xmax[1] = xmin[2]
      
      full_xmin = as.POSIXct(rep(NA, ncol(Datai$TIME)), tz="UTC")
      full_xmax = as.POSIXct(rep(NA, ncol(Datai$TIME)), tz="UTC")
      full_xmin[nna] = xmin
      full_xmax[nna] = xmax
      # Debugging : 198283 = nrow(Datai$TIME is 851) and full_xmin is 233 where does those number come from ?
      #data_df <- load_float_data(float_ids = 5904183,format="dataframe")
      # data_df$TIME %>% unique() %>% length()
      
      # So 233 is unique number of dates in that float
      # Where does 851 comes from ? format="dataframe"
      list_of_tibbles <- df  %>%
        group_by(TIME) %>%
        group_split()
      list.df <- list()
      
      for (i in seq_along(list_of_tibbles)) {
        tibble_i <- list_of_tibbles[[i]]
        tibble_i$BBP700_ADJUSTED <- rollmedian(tibble_i$BBP700_ADJUSTED, k = 3, fill = NA)
        # Assign xmin and xmax to each tibble
        tibble_i$xmin <- rep(full_xmin[i], nrow(tibble_i))
        tibble_i$xmax <- rep(full_xmax[i], nrow(tibble_i))
        
        # Replace the tibble in the list with the modified one
        list.df[[i]] <- tibble_i
      }
      
      

      df <- list.df %>% bind_rows()
      
      ### END OF DATAFRAME GEN
      
      return(df)
      
}


  
df.25 <- df %>% filter(CYCLE_NUMBER == 25)
df.25 %>% ggplot(aes(x=PRES_ADJUSTED,y=BBP700_ADJUSTED))+geom_point()+geom_line()+
  scale_x_reverse()+coord_flip()

  variables <- c("SPIC")
  #  if(is.null(Data[[floats[f]]][[variables[v] ]])){ # Check if the float has variable
  #    print(paste("No",variables[v],"available for float",floats[f]))
  #    next
  #  }

  
   # For 5904677 : 
  # 5 September to 25 November 2016
  

df$TIME %>% str()  
  
  
plotting_fun <- function(df)  
  
start_date <- as.POSIXct("2016-01-01")
end_date <- as.POSIXct("2016-12-31")

# Filter the dataframe for the specified date range
df_filtered <- df %>%
  filter(TIME >= start_date & TIME <= end_date)


  if ("PRES_ADJUSTED" %in%   colnames(df)
 ) {
    g1 = ggplot(df.wf, aes(x=TIME, y=PRES_ADJUSTED)) + facet_wrap(. ~ name)
  } else {
    g1 = ggplot(df.wf, aes(x=TIME, y=PRES)) + facet_wrap(. ~ name,scales = "free")
  }
  


  # Create individual plots for each variable
  plot_AOU <- ggplot(df_filtered, aes(x=TIME, y=PRES_ADJUSTED)) +
    geom_rect(aes(fill=AOU, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    scale_fill_viridis_c() + theme_bw() + ggtitle("AOU")+ 
    scale_y_reverse(limits = c(max_depth, 0))+
    geom_point(data = anomalies_det, aes(x=TIME, y=PRES_ADJUSTED), color="red", size=3)
  
  plot_SPIC <- ggplot(df_filtered, aes(x=TIME, y=PRES_ADJUSTED)) +
    geom_rect(aes(fill=SPIC, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    scale_fill_viridis_c() + theme_bw() + ggtitle("SPIC") +
    scale_y_reverse(limits = c(max_depth, 0))+
    geom_point(data = anomalies_det, aes(x=TIME, y=PRES_ADJUSTED), color="red", size=3)
  
  plot_BBP700_ADJUSTED <- ggplot(df_filtered, aes(x=TIME, y=PRES_ADJUSTED)) +
    geom_rect(aes(fill=BBP700_ADJUSTED, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    scale_fill_viridis_c() + theme_bw() + ggtitle("BBP700_ADJUSTED") +
    scale_y_reverse(limits = c(max_depth, 0))+
    geom_point(data = anomalies_det, aes(x=TIME, y=PRES_ADJUSTED), color="red", size=3)
  
  
  plot_BBP700 <- ggplot(df, aes(x=TIME, y=PRES_ADJUSTED)) +
    geom_rect(aes(fill=BBP700, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    scale_fill_viridis_c() + theme_bw() + ggtitle("BBP700")+
    geom_point(data = anomalies_det, aes(x=TIME, y=PRES_ADJUSTED), color="red", size=3)
  
  
  # Use ggarrange to combine the plots
  combined_plots <- ggarrange(plot_AOU, plot_SPIC, plot_BBP700_ADJUSTED, 
                              ncol = 2, nrow = 2, common.legend = FALSE, legend = "right")


  
  name_units = get_var_name_units(variables[v])
  long_name = name_units$long_name
  units = name_units$units
  
  if (obs == "on") {
    index = which(!is.na(Data[[floats[f]]][[variables[v]]]))
    g1 = g1 + geom_point(aes(x = x, y = y),
                         data = data.frame(x = as.POSIXct(Data[[floats[f]]]$TIME[index], tz="UTC"),
                                           y = as.vector(Data[[floats[f]]]$PRES[index])),
                         size=0.1, alpha=0.2
    )
  }
  # Plot the mixed layer
  if (plot_mld == 1) {
    g1 = g1 + geom_line(aes(x = x, y = y),
                        data = data.frame(x = as.POSIXct(Datai$TIME[1,], tz="UTC"),
                                          y = as.vector(Datai$MLD_TEMP)),
                        size=0.4,color="red",alpha=.5
    )
  } else if (plot_mld == 2 & "MLD_DENS" %in% names(Datai)) {
    # some old core floats don't have PSAL and therefore
    # density cannot be computed
    g1 = g1 + geom_line(aes(x = x, y = y),
                        data = data.frame(x = as.POSIXct(Datai$TIME[1,], tz="UTC"),
                                          y = as.vector(Datai$MLD_DENS)),
                        size=1
    )
  }
  # If you want to plot isopycnals
  if ( plot_isopyc ) {
    if ("DENS_ADJUSTED" %in% names(Datai)) {
      g1 = g1 + geom_contour(aes(z = DENS_ADJUSTED), color="red") +
        geom_label_contour(aes(z = DENS_ADJUSTED))
    } else {
      g1 = g1 + geom_contour(aes(z = DENS), color="red") +
        geom_label_contour(aes(z = DENS))
    }
  }
  # Reverse y axis so that top of the ocean is on top of the plot
  if ( !is.null(max_depth) ) {
    g1 = g1 + scale_y_reverse(limits = c(max_depth, 0))
  } else {
    g1 = g1 + scale_y_reverse()
  }
  
  g1 = g1 +
    labs(title = paste0("Float ", Mdata[[float_ids[f]]]$WMO_NUMBER,": ",
                        long_name, title_add),
         x = "Time",
         y = "Pressure (dbar)",face = 'bold',family = "serif",
         fill = units)
  g1= g1+theme (axis.title.y = element_text(size=16,colour = "black",face = "bold",family = "serif") ) 
  g1= g1+theme (axis.title.x = element_text(size=16,colour = "black",face = "bold",family = "serif") ) 
  g1= g1+theme (axis.text.y = element_text(size=16,colour = "black",face = "bold",family = "serif") ) 
  g1= g1+theme (axis.text.x = element_text(size=16,colour = "black",face = "bold",family = "serif") )
  g1=g1+theme(legend.text = element_text(size = 16,face = 'bold',family = "serif"),
              legend.title  = element_text(size = 16,face = 'bold',family = "serif"),
              legend.key.width=unit(1,'cm'),
              legend.key.height=unit(1,'cm'))#
  g1=g1+theme(plot.title = element_text(size = 16, face = "bold",family = "serif"))
  
  
  x11()
  plot(g1)


