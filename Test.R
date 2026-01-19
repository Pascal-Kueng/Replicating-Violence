library(tidyverse)
library(janitor)
library(survival)
library(sandwich)
library(lmtest)

# 1. Ensure Data is Loaded and Cleaned
df <- readxl::read_xlsx("S1_file_combined.xlsx")
df_model <- df %>%
  clean_names() %>%
  mutate(
    # Reconstruct Period (Handle inconsistent underscores from raw data)
    period_raw = case_when(
      t0n_ == 1 ~ "Neolithic",       # Has underscore
      t1c_ == 1 ~ "Chalcolithic",    # Has underscore
      t2be == 1 ~ "Early Bronze",    # No underscore
      t3bm == 1 ~ "Middle Bronze",   # No underscore
      t4bl == 1 ~ "Late Bronze",     # No underscore
      t5i_ == 1 ~ "Iron Age",        # Has underscore
      TRUE ~ NA_character_
    ),
    period_factor = factor(period_raw, levels = c(
      "Chalcolithic", "Neolithic", "Early Bronze", 
      "Middle Bronze", "Late Bronze", "Iron Age"
    )),
    
    # Reconstruct Region
    region = case_when(
      co2 == "l_" ~ "Levant",
      co2 == "ms" ~ "Mesopotamia",
      co2 == "tr" ~ "Turkey",
      co2 == "ir" ~ "Iran",
      TRUE ~ "Unknown"
    ),
    
    # Prepare Interval Data (Log-transformed)
    # tr_nop is the variable used in Stata Model 1
    y_lower = ifelse(tr_nop > 0, log(tr_nop), -Inf),
    y_upper = ifelse(tr_nop > 0, log(tr_nop), -3) # -3 is approx log(0.05)
  )

# 2. FORCE creation of status_var (Fixes your error)
df_model$status_var <- 3  # Code 3 tells Surv this is "Interval" data

# 3. Create Survival Object
surv_obj <- Surv(time  = df_model$y_lower, 
                 time2 = df_model$y_upper, 
                 event = df_model$status_var, 
                 type  = "interval")

# 4. Run Model 1 (Replicating Stata's intreg)
m1_intreg <- survreg(
  surv_obj ~ period_factor,
  data = df_model,
  dist = "gaussian", # Gaussian interval regression = Tobit
  weights = ncases   # Weight by number of skeletons
)

# 5. Clustered Standard Errors (by 'co2lf')
# Note: if co2lf is missing/character, ensure it's a factor/ID
coeftest(m1_intreg, vcov = vcovCL(m1_intreg, cluster = df_model$co2lf))