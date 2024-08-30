# Begin by reverse engineering "show_sections.R"


setwd("~/Documents/GLOBARGO")
library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
wmo <- 5904677
data_df = load_float_data( float_ids= wmo, # specify WMO number
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

