library(readr)
library(tidyverse)

data_class <- read_csv("~/Documents/GLOBARGO/data/classification_results_manually_classified_merged.csv")

data_class %>% select(WMO,CYCLE_NUMBER) %>% unique() 
