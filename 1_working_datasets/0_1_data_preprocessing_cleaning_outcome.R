## The aims of this script is to make sure all the outcomes are in the desired format

# Check categorical outcome variable type (Ordinal)
prep_data_fun <- function(df){
  # steato
  df$steato_0 <- df$steato_0 %>%
    factor(levels = c("No steatohepatitis",
                      "Possible/borderline steatohepatitis (Type 1, typical zone 3 pattern)",
                      "Possible/borderline steatohepatitis (Type 2, zone 1 pattern)",
                      "Definite steatohepatitis"), ordered = TRUE)
  #steatgrd
  df$steatgrd_0 <- df$steatgrd_0%>% 
    factor(levels = c("0%", "<5%", "5-33%", "33-67%", ">67%" ), ordered = TRUE)
  
  # bhepa
  df <- df %>% mutate(bhepa_0 = case_when(bhepa_0 == 0 ~ "None",
                                          bhepa_0 == 1 ~ "Fedata_w, less characteristic",
                                          bhepa_0 == 2 ~ "Many, prominent"
                                          ))
  df$bhepa_0 <- df$bhepa_0 %>%
    factor(levels = c("None", 
                      "Fedata_w, less characteristic",
                      "Many, prominent"
                      ), ordered = TRUE)
  
  #Fibrostg
  df$fibrostg_0 <- df$fibrostg_0 %>%
    factor(levels = c("None", "Periportal OR Perisinusoidal", "Periportal AND Perisinusoidal","Bridging fibrosis","Cirrhosis"), ordered = TRUE)
  
  #lob
  df$lob_0 <- df$lob_0%>% 
    factor(levels = c("None", "Mild", "Moderate", "Marked"), ordered = TRUE)
  
  # NASH
  df$nash_0 <- df$nash_0 %>%
    factor(levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8), ordered = TRUE)
  
  # htn
  df$htn_0 <- df$htn_0 %>%
    factor(levels = c("No", "Yes"), ordered = TRUE)
  
  # dyslipid
  df$dyslipid_0 <- df$dyslipid_0 %>%
    factor(levels = c("Not Present",
                      "No pharmacologic treatment of dyslipidemia",
                      "Treatment data_with single medication",
                      "Treatment data_with tdata_wo or more medications"), ordered = TRUE)
  # hyptn
  df$hyptn_0 <- df$hyptn_0 %>%
    factor(levels = c("No BP elevation diagnosed",
                      "Hypertension, no pharmacologic treatment",
                      "Hypertension, data_with single pharmacologic medication",
                      "Hypertension, data_with tdata_wo or more pharmacologic medications"),
           ordered = TRUE)
  return(df)
}

data_w_complete <- prep_data_fun(data_w_all_pfas)

rm(data_w_all_pfas)