library(tidyverse)
library(gsw)
library(zoo)
library(oce)
library(robustbase)
library(parallel)

setwd("~/Documents/ARGO/")
save_dir <- "~/Documents/ARGO/Partie1/"

data.list <- readRDS("ToTheGlobe/results_865.rds")
 


outlier.detection.fun.row <- function(WMO){
  
  data_df = load_float_data( float_ids= WMO, # specify WMO number
                             variables=c(
                               "DATA_TYPE" 
                               ,"PLATFORM_NUMBER"
                               ,"BBP700"
                               ,"BBP700_dPRES"
                               ,"BBP700_ADJUSTED_QC"
                               ,"LATITUDE"
                               ,"LONGITUDE"
                               ,"PROFILE_TEMP_QC"                              
                               ,"PROFILE_DOXY_QC"                                 
                               ,"PROFILE_BBP700_QC"                                                
                               ,"PRES_QC"            
                               ,"PRES"  
                               ,"PRES_ADJUSTED"       
                               ,"PROFILE_PSAL_QC"
                               ,"CHLA_QC"  
                               ,"CHLA_ADJUSTED"
                               ,"CHLA_ADJUSTED_ERROR"
                               ,"DOXY"
                               ,"DOXY_QC"
                               ,"DOXY_ADJUSTED"
                               ,"DOXY_ADJUSTED_QC"
                               ,"DOXY_ADJUSTED_ERROR"
                               ,"PSAL"
                               ,"PSAL_dPRES"
                               ,"PSAL_ADJUSTED"
                               ,"PSAL_ADJUSTED_QC"
                               ,"TEMP"                            
                               ,"TEMP_QC"                        
                               ,"TEMP_dPRES"
                               ,"TEMP_ADJUSTED"                  
                               ,"TEMP_ADJUSTED_QC"
                               ,"TEMP_ADJUSTED_ERROR"), # specify variables,
                             format="dataframe" # specify format;  
  )
  
  
  # Check if ADJUSTED variables are available
  # Calculate percentage of missing values for each column
  missing_percent <- data_df %>% select(PRES_ADJUSTED,TEMP_ADJUSTED) %>% 
    summarize_all(~ mean(is.na(.))) %>% 
    gather(column, missing_percentage)
  
  # Check if any column has more than 50% missing values
  has_high_missing <- any(missing_percent$missing_percentage > 0.5)
  
  # Logical check and perform action
  if (has_high_missing) {
    output <- NULL
  } else {
    # Use adjusted variables to compute AOU & SPIC
    
    data_df <- data_df %>% group_by(CYCLE_NUMBER) %>% 
      mutate(SPIC = swSpice(salinity = PSAL_ADJUSTED,
                            temperature = TEMP_ADJUSTED,
                            latitude = first(LATITUDE),
                            longitude = first(LONGITUDE),
                            eos = "unesco"),
             ABS_SAL = gsw::gsw_SA_from_SP(SP = PSAL_ADJUSTED, p = PRES_ADJUSTED, 
                                           longitude = first(LONGITUDE),
                                           latitude = first(LATITUDE)),
             CONS_TEMP =  gsw::gsw_CT_from_t(SA = ABS_SAL, t = TEMP_ADJUSTED,
                                             p = PRES_ADJUSTED),
             SAT_DOXY = gsw_O2sol(SA= ABS_SAL,CT= CONS_TEMP,p= PRES_ADJUSTED,
                                  longitude= first(LONGITUDE),latitude= first(LATITUDE)),
             AOU =  SAT_DOXY - DOXY_ADJUSTED
      )
    
    
    
    result_df <- data_df %>%
      # Filter out missing DOXY values
      filter(!is.na(DOXY)) %>%
      # Group by CYCLE_NUMBER
      group_by(CYCLE_NUMBER) %>%
      # Apply the computations for each unique CYCLE_NUMBER
      group_modify(~ .x %>%
                     select(TIME,PRES_ADJUSTED, AOU, SPIC, BBP700_ADJUSTED,LATITUDE,LONGITUDE) %>%
                     pivot_longer(!c(TIME,PRES_ADJUSTED,LATITUDE,LONGITUDE), names_to = "VAR", values_to = "VALUE") %>%
                     group_by(VAR) %>%
                     mutate(
                       # moving average k = 3
                       MA_3 = rollmean(VALUE, 3, fill = NA),
                       MM_11 = rollmedian(VALUE, 11, fill = NA),
                       ROB.RES = MA_3 - MM_11,
                       SCALE.RES.ROB = ifelse(abs(VALUE-MM_11) == 0, 0, 
                                              (ROB.RES - median(ROB.RES, na.rm = TRUE)) / mad(ROB.RES, na.rm = TRUE))
                     )) %>% 
      ungroup()  # To remove the grouping
    
    # Check the result
    
    
    result.df.wf <- result_df %>% select(TIME,PRES_ADJUSTED,VAR,SCALE.RES.ROB,LONGITUDE,LATITUDE,VALUE,CYCLE_NUMBER) %>%
      pivot_wider(values_from = c(SCALE.RES.ROB,VALUE),names_from=VAR) 
    
    result.df.wf <- result.df.wf %>%
      mutate(OUT.S = ifelse(abs(SCALE.RES.ROB_AOU) > 3 & abs(SCALE.RES.ROB_SPIC) > 3 & SCALE.RES.ROB_AOU < 0,1,0),
             OUT.T = ifelse(OUT.S == 1 & abs(SCALE.RES.ROB_BBP700_ADJUSTED) > 1.96 &
                              SCALE.RES.ROB_BBP700_ADJUSTED > 0,1,0))
    
    
    output <- result.df.wf %>% filter(OUT.S == 1)
    # To avoid ties in the dataset pick only largest BBP700
    
    output$WMO <- WMO
    output <- output %>%
      select(WMO,CYCLE_NUMBER,TIME,LONGITUDE,LATITUDE,PRES_ADJUSTED,OUT.T) %>%
      unique() 
  }
  return(output)
  
  
}

full.var.sel <- c(
  "BBP700"
  ,"PRES"  
  ,"CHLA"
  ,"DOXY"
  ,"PSAL"
  ,"TEMP")

contains_all <- function(x) {
  all(full.var.sel %in% x)
}

indices <- sapply(Sprof$split_sens,contains_all)

wmolist <- Sprof$wmo[indices] %>% unique()

# Define a directory to save the results
save_dir <- "~/Documents/ARGO/ToTheGlobev3/"






# Example of saving intermediate results
result_list <- list()
for (i in seq_along(wmolist)) {
  
  # Using tryCatch to handle potential errors
  processed_result <- tryCatch({
    outlier.detection.fun.row(wmolist[i])
  }, 
  error = function(e) {
    message("Error processing result ", i, ": ", e$message)
    return(NULL)  # or some other default value if needed
  })
  
  result_list[[i]] <- processed_result
  
  # Save every 50th result as before
  if (i %% 100 == 0) {
    saveRDS(result_list, paste0(save_dir,"results_", i, ".rds"))
  }
}


saveRDS(r)












####################### OTHER CODE


result_list <- list()

for (i in seq_along(wmolist)) {
  
  # Using tryCatch to handle potential errors
  processed_result <- tryCatch({
    fun_rowwise_det(wmolist[i])
  }, 
  error = function(e) {
    message("Error processing result ", i, ": ", e$message)
    return(NULL)  # or some other default value if needed
  })
  
  result_list[[i]] <- processed_result
  
  # Save and clear the list every 100th result
  if (i %% 10 == 0) {
    save_file_name <- paste0(save_dir,"results_", i, ".rds")
    saveRDS(result_list, save_file_name)
    
    # Reinitialize result_list after saving
    result_list <- list()
    
    # Optional: message to indicate saving progress
    message("Results up to ", i, " saved to ", save_file_name)
  }
}

# Don't forget to save the final batch if it doesn't end exactly on a multiple of 100
if (length(result_list) > 0) {
  save_file_name <- paste0(save_dir,"results_final.rds")
  saveRDS(result_list, save_file_name)
  message("Final batch of results saved to ", save_file_name)
}

