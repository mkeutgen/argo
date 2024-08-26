#CYCLE_NUMBER = 77
library(readr)
classification_results_v3 <- read_csv("/data/GLOBARGO/data/classification_results_v3.csv")
classification_results_v3 <- classification_results_v3 %>%  mutate(WMO = str_replace(WMO, "_plot$", ""))
classification_results_v3$CYCLE_NUMBER <- classification_results_v3$Cycle
classification_results_v3$WMO <- as.numeric(classification_results_v3$WMO)

detected.events.df <- read_csv("/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")
detected.events.df %>% select(LATITUDE,LONGITUDE,WMO,CYCLE_NUMBER) %>% unique()

merged_df <- detected.events.df %>%
  left_join(classification_results_v3, by = c("WMO", "CYCLE_NUMBER"))

# Remove spurrious Cycle column
merged_df <- merged_df %>% dplyr::select(-Cycle)

# Create a dataframe for the mapping of events, pick lowest value for PRES_ADJUSTED :
unique_sub_loc_df <- merged_df  %>%
  group_by(WMO, CYCLE_NUMBER) %>%
  slice_min(PRES_ADJUSTED) %>%
  ungroup() %>%
  select(CYCLE_NUMBER, LATITUDE, LONGITUDE, WMO, Category, TIME, PRES_ADJUSTED,CONSISTENT_ANOM)


unique_sub_loc_df %>% filter(Category %in% c(1,2))

unique_sub_loc_df %>% filter(Category %in% c(1,2)) %>% filter(CONSISTENT_ANOM == 1)

unique_sub_loc_df %>% filter(Category %in% c(1,2)) %>% filter(CONSISTENT_ANOM == 0)


write_csv(x = unique_sub_loc_df,file="/data/GLOBARGO/src/data/unique_sub_loc.csv")

unique_sub_loc_df %>% filter(WMO==5904677)

merged_df %>% filter(Category == 0) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()
merged_df %>% filter(Category == 1) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()
merged_df %>% filter(Category == 2) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()
merged_df %>% filter(Category == 3) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()


merged_df %>% filter(Category == 0) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
merged_df %>% filter(Category == 1) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
merged_df %>% filter(Category == 2) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
merged_df %>% filter(Category == 3) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()


unique_sub_loc_df %>% filter(Category == 0) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
unique_sub_loc_df %>% filter(Category == 1) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
unique_sub_loc_df %>% filter(Category == 2) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
unique_sub_loc_df %>% filter(Category == 3) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()



category_proportions_final <- classification_results_v3 %>%
  group_by(Category) %>%
  summarize(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

