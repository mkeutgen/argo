# Load necessary libraries
library(tidyverse)
library(fs)
library(conflicted)

# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load datasets
df_abs_sal <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v4.csv")
df_spic <- read_csv("/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")
df_spic_class <- read_csv("/data/GLOBARGO/data/classification_results_v4.csv")

# Preprocess df_spic_class
df_spic_class <- df_spic_class %>%
  mutate(WMO = gsub("_plot", "", WMO),
         CYCLE_NUMBER = Cycle)

# Count unique profiles
df_abs_sal %>% distinct(WMO, CYCLE_NUMBER) %>% count()  # ABS_SAL unique profiles
df_spic %>% distinct(WMO, CYCLE_NUMBER) %>% count()     # SPIC unique profiles

# Find profiles in ABS_SAL but not in SPIC
df_not_in_spic <- df_abs_sal %>%
  anti_join(df_spic, by = c("WMO", "CYCLE_NUMBER"))

# Find profiles in SPIC but not in ABS_SAL
df_not_in_sal <- df_spic %>%
  anti_join(df_abs_sal, by = c("WMO", "CYCLE_NUMBER"))

# Create destination folder if it doesn't exist
destination_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresInSalinityButNotSpic"
dir_create(destination_folder)

# Copy files that are in ABS_SAL but not in SPIC, preserving folder structure
df_not_in_spic %>%
  rowwise() %>%
  mutate(
    source_subfolder = file.path("/data/GLOBARGO/figures/EddySubductionFiguresSalinityVarV4", as.character(WMO)),
    source_file = file.path(source_subfolder, paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")),
    destination_subfolder = file.path(destination_folder, as.character(WMO)),
    destination_file = file.path(destination_subfolder, paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png"))
  ) %>%
  pwalk(function(source_file, destination_subfolder, destination_file) {
    if (file_exists(source_file)) {
      dir_create(destination_subfolder, recurse = TRUE)
      file_copy(source_file, destination_file, overwrite = TRUE)
    }
  })

# Filter anomalies classified as Category 1 or 2
df_spic_anomalies <- df_spic_class %>%
  filter(Category %in% c(1, 2))

# Check for corresponding files in ABS_SAL and save those
df_spic_anomalies_with_files <- df_spic_anomalies %>%
  mutate(
    source_subfolder = file.path("/data/GLOBARGO/figures/EddySubductionFiguresSalinityVarV4", as.character(WMO)),
    source_file = file.path(source_subfolder, paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")),
    file_exists = file_exists(source_file)
  ) %>%
  filter(file_exists == TRUE) %>%
  select(WMO, CYCLE_NUMBER, Category)

# Save the anomalies that exist in both SPIC and ABS_SAL
write_csv(df_spic_anomalies_with_files, "/data/GLOBARGO/data/anom_in_spic_and_sal_cat1_and2.csv")

# Create destination folder for classified anomalies
destination_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresSpicClass1and2"
dir_create(destination_folder)

# Copy files of classified anomalies, preserving folder structure
df_spic_anomalies %>%
  rowwise() %>%
  mutate(
    source_subfolder = file.path("/data/GLOBARGO/figures/EddySubductionFiguresSalinityVarV4", as.character(WMO)),
    source_file = file.path(source_subfolder, paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")),
    destination_subfolder = file.path(destination_folder, as.character(WMO)),
    destination_file = file.path(destination_subfolder, paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png"))
  ) %>%
  pwalk(function(source_file, destination_subfolder, destination_file) {
    if (file_exists(source_file)) {
      dir_create(destination_subfolder, recurse = TRUE)
      file_copy(source_file, destination_file, overwrite = TRUE)
    }
  })

cat("Process complete.")
