# directory paths for file architecture

# home directory for project
dir_home <- here::here() %>% 
  dirname() %>% 
  dirname() %>% 
  dirname()

# Baseline PFAS - NAFLD_final project directory
dir_project <- fs::path(
  here::here() %>% dirname(),
  "8_Finalizing results",
  "Baseline PFAS - NAFLD_final project")

# data folder 
# dir_cleaned_data <- fs::path(dir_home,"2_Cleaned Data") # original data (data is: tl_covariates_outcomes_environmentalchemicals_w.rds)
dir_data <- fs::path(dir_home,"2_Cleaned Data")

#Merged data folder path
dir_project_data <- fs::path(dir_project,"0_data")

# report folder
dir_report <- fs::path(dir_project,"2_reports")

# figure folder
dir_figure <- fs::path(dir_project, "3_figures")

dir_temp <- fs::path(dir_home,
                     "4_Projects",
                     "TL_PFAS NAFLD",
                     "1_2_Plasma Metabolomics Scripts",
                     "0_data",
                     "temp")

dir_final <- fs::path(dir_home,
                      "4_Projects",
                      "TL_PFAS NAFLD",
                      "1_2_Plasma Metabolomics Scripts",
                      "0_data")

dir_results <- fs::path(dir_final %>% dirname(),
                        "2_results")
