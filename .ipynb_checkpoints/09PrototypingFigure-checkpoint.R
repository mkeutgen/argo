# 9 PROTOTYPING FIGURE 
library(tidyverse)


# Let WMO : 
WMO = 5904105


# Find right times : 
# Cycle 20 : 2013-06-09 
# Cycle 23 : "2013-06-24"  
# Cycle 57 : "2013-12-16"
# Cycle 83 : "2014-06-15"
# Cycle 99 : "2014-10-08"
# Cycle 108 : "2014-12-12"
# Cycle 142 : "2015-08-13"


show_sections(float_ids=5904105, 
              variables= c('DOXY'),
              plot_mld=1,       # tells the function to plot mixed layer depth
              raw="yes") # tells the function to plot raw data

