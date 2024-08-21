## This script is to get imputed pfas and add them back to the original wide format data(data_w)

## Load LOD data file------
LOD <- read_csv(fs::path(dir_home, 
                         "4_Projects", 
                         "TL_PFAS NAFLD", 
                         "0_Data",
                         "LOD.csv"))
## Create a variable that contains all pfas names (targeted and untargeted)------
plasma_pfas_names<-
  data_w %>% dplyr::select(contains("_plasma"),
                         -contains("run_order"), -contains("batch_number"), -contains("niddk")) %>% 
  colnames()
liver_pfas_names<- 
  data_w %>% dplyr::select(matches("_liver$"),
                         -contains("identifier_targeted_liver"), -contains("weightofbiopsy_targeted_liver")) %>% 
  colnames()

all_pfas_names <- c(plasma_pfas_names,liver_pfas_names)

## Imputed data preparing------- 
###(1) Converting wide format data to long format data based on all exposures-------
###(2) Adding pfas names, measurement type, matrices and visit number -------
data_impute <- data_w %>%
  dplyr::mutate(row = row_number()) %>%
  labelled::remove_labels() %>% 
  pivot_longer(names_to="pfas", 
               values_to="pfas_value",
               cols =  all_pfas_names) %>%
  tidyr::separate(col = "pfas" , 
                  into = c("pfas_names", "visit_new"), 
                  sep = "_plasma_|_liver|_plasma", 
                  remove = FALSE) %>%
  mutate(measurement_type = ifelse(grepl("_targeted",pfas),
                            'targeted',
                            'untargeted'), 
        matrices = ifelse(grepl("plasma", pfas),
                          "plasma",
                          "liver")) 
### (3) merge with LOD----
data_impute2 <- data_impute %>% 
  tidylog::left_join(LOD, by = c("pfas_names"="pfas_names_id")) %>%
  mutate(is_less_lod = ifelse((pfas_value - lod)<0, 1, 0)) 

### (4) adding two variables to the data_impute set----
#### is_less_lod: 1, plasma value is less than LOD; 0, plasma value is equal and more than LOD
#### percent_less_than: percentage of each pfas variable in data_impute less than LOD
data_impute3<- data_impute2 %>%
  group_by(pfas) %>% 
  dplyr::mutate(percent_less_lod =
                sum(is_less_lod,na.rm=TRUE)/n())


### (5) Adding one variable: pfas_value_imputed that contains the imputed value of pfas-----
#### getting the name of the continuous pfas
continuous_pfas_names <- data_impute3 %>% 
  filter(percent_less_lod <= 0.3) %>% 
  dplyr::select(pfas) %>%
  unique()
  
continuous_pfas_names <- continuous_pfas_names$pfas

#### getting the name of the categorical pfas
categorical_pfas_names<-
  setdiff(all_pfas_names,
          continuous_pfas_names)
  
continuous_untargeted_pfas_names <- continuous_pfas_names[grepl("_untargeted",continuous_pfas_names)]
#### notes: all the liver pfas and the target plasma pfas are continuous variable.
data_impute4 <- data_impute3 %>%
  mutate(pfas_value_imputed =
           case_when(
             (pfas %in%
                continuous_untargeted_pfas_names) &
               is_less_lod == 1 &
               matrices == "plasma" ~ lod/2, 
             (pfas %in%
                categorical_pfas_names) &
               is_less_lod == 1 ~ -1,
             (pfas %in%
                categorical_pfas_names) & 
               is_less_lod == 0 ~ 9999,
             TRUE ~ pfas_value
             ))

data_impute5 <- data_impute4%>%
  mutate(visit_new = ifelse(visit_new == "", "0", visit_new))

data_descriptive_stats <- data_impute5

###(6) adding either "imputed" or "dichotomous" to plasma names-----
data_impute6 <- data_impute5 %>% mutate(
  pfas_value_imputed = pfas_value_imputed %>% as.character(),
  suffix = case_when(
    pfas %in% categorical_pfas_names ~ "_dichotomous",
    pfas %in% c("pf_hp_a_untargeted_plasma_6", "pf_hp_a_untargeted_plasma_12") ~ "_dichotomous",
    pfas %in% continuous_untargeted_pfas_names ~ "_imputed",
    TRUE ~ ""
  ),
  pfas = paste(pfas, suffix, sep = ""),
  pfas_value_imputed = ifelse(pfas_value_imputed == 9999, "Detected", 
                              ifelse(pfas_value_imputed == -1, "Not Detected", pfas_value_imputed)),
  pfas_value_imputed = case_when(
    pfas_names == "pf_hp_a_untargeted"& visit_new == "6" & is_less_lod == 1 ~ "Not Detected",
    pfas_names == "pf_hp_a_untargeted"& visit_new == "12" & is_less_lod == 1 ~ "Not Detected",
    pfas_names == "pf_hp_a_untargeted"& visit_new == "6" & is_less_lod == 0 ~ "Detected",
    pfas_names == "pf_hp_a_untargeted"& visit_new == "12" & is_less_lod == 0 ~ "Detected",
    TRUE ~ pfas_value_imputed
  )) %>%
  dplyr::select(-c(suffix,is_less_lod, pfas_value, 
                   percent_less_lod, lod, visit_new, 
                   pfas_names, matrices, measurement_type))

###(7) converting longer format back to wide format------
data_w_imputed <- data_impute6 %>% 
  tidyr::pivot_wider(names_from = pfas, values_from = pfas_value_imputed) %>%
  mutate_at(.vars = vars(contains("imputed")),
            .funs = as.numeric)

###(8) Merging imputed pfas with original wide format data(data_w)
data_w_all_pfas <- data_w %>% tidylog::left_join(data_w_imputed %>% 
                                                   dplyr::select(key,
                                                                 contains("_imputed"),
                                                                 contains("_dichotomous")), 
                                                 by = "key")

rm(LOD, plasma_pfas_names, liver_pfas_names, all_pfas_names,
   data_impute, data_impute2, data_impute3, data_impute4, 
   data_impute5, data_impute6, data_w_imputed, data_w)




