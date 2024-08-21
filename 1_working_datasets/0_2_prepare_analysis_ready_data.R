# The purpose of this script is to restrict the variables in the dataframe to just those used in this project
covars <- c("bmi_0", 
            "race_binary", 
            "smoke_0", 
            "age_0", 
            "sex", 
            
            "parents_income_0")

sensitivity_variables <- c("alb_0", 
                           "e_gfr_0",
                           "site")

cat_outcomes <- c("nafld_type_0", "nafld_nash_mul_0", "nafld_di_0",
                  "steato_0", "steato_mul_0", "steatgrd_0","steatgrd_mul_0", "bhepa_0",
                  "fibrostg_0", "fibrostg_di_0","lob_0","lob_di_0" ,"nash_0","nash_mul_0",
                  "htn_0", "dyslipid_0", "dyslipid_di_0")

cont_outcomes <- c("alt_0", "ast_0", "ggt_0",
                   "dbp_0", "sbp_0", "bmi_0", 
                   "hdl_0", "ldl_0", "tc_0", 
                   "trig_0", "iwaist_0", "uwaist_0")

data <- data_w_complete %>%
  mutate(age_0 = agemos_0/12) %>%
  dplyr::select(key,
                all_of(covars),
                all_of(sensitivity_variables),
                all_of(cat_outcomes),
                all_of(cont_outcomes),
                contains("_targeted_plasma"),
                contains("_liver"),
                contains("_imputed"),
                -contains("batch_number"),
                -contains("niddk"),
                -contains("run_order"))


# Saving analysis ready data

write_rds(data, fs::path(dir_project_data, "tl_analysis_ready_data.rds"))

rm(data_w_complete, data)               
