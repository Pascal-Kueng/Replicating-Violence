library(readxl)
library(tidyverse)

viol_data <- read_excel("S1_file_combined.xlsx")

# Explanation: For Figure 2, the authors excluded all cases with less than 30 observations except for Late Bronze Turkey (n = 9). 
# In the supplementary materials, they justify keeping Turkey-LBA because it aligns with interregional trends of increasing violence during Late Bronze Age,
# while excluding Mesopotamia-Chalcolithic (n = 12) because it contradicts these trends (shows 0% violence when they claim violence "peaked" in Chalcolithic).
# They acknowledge this exclusion could reflect either real regional variation or selectivity bias from small sample size.
# To assess robustness of their choice, we recreate their plot and compare it to versions with all exclusions (n < 30) and no exclusions.


##############################################################################
##### Table 1 (Number of cases for violence shares by region and period) #####
##############################################################################

viol_data %>%
  group_by(co2, LF) %>%
  summarise(ncases = sum(ncases, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = LF, values_from = ncases) %>%
  bind_rows(summarise(., across(where(is.numeric), sum), co2 = "ME")) %>%
  mutate(Total = rowSums(select(., where(is.numeric)), na.rm = TRUE),
         co2 = recode(co2, ir = "Iran", l_ = "Levant", ms = "Mesopotamia", 
                      tr = "Turkey", ME = "Middle East")) %>%
  select(Region = co2, Neolithic = `0n_`, Chalcolithic = `1c_`, 
         EBA = `2be`, MBA = `3bm`, LBA = `4bl`, IA = `5i_`, Total)


#################################################
##### Violence shares by region and period) #####
#################################################

viol_data %>%
  group_by(co2, LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = LF, values_from = tr_nop) %>%
  bind_rows(
    viol_data %>%
      group_by(LF) %>%
      summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = LF, values_from = tr_nop) %>%
      mutate(co2 = "ME")
  ) %>%
  mutate(co2 = recode(co2, ir = "Iran", l_ = "Levant", ms = "Mesopotamia", 
                      tr = "Turkey", ME = "Middle East")) %>%
  select(Region = co2, Neolithic = `0n_`, Chalcolithic = `1c_`, 
         EBA = `2be`, MBA = `3bm`, LBA = `4bl`, IA = `5i_`)


###################################################
##### Figure 2 as displayed in the manuscript #####
###################################################

# Inclusion of Turkey-Late Bronze (n = 9) but exclusion of all other cases with less than 30 observations 

# Collapse by region and time period
collapsed_data_paper <- viol_data %>%
  group_by(co2, LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop')

# Drop observations with low case numbers except Turkey LBA
collapsed_data_paper <- collapsed_data_paper %>%
  filter(!(co2 == "ir" & LF == "2be"),   # Iran EBA (n = 2)
         !(co2 == "l_" & LF == "4bl"),   # Levant LBA (n = 5)
         !(co2 == "ms" & LF == "1c_"))   # Mesopotamia Chalcolithic (n = 12) 

# Calculate Middle East average (all regions combined)
middle_east_paper <- viol_data %>%
  filter(!(co2 == "ir" & LF == "2be"),
         !(co2 == "l_" & LF == "4bl"),
         !(co2 == "ms" & LF == "1c_")) %>%
  group_by(LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop') %>%
  mutate(co2 = "ME")

# Combine individual regions with Middle East summary
figure3_data_paper <- bind_rows(collapsed_data_paper, middle_east_paper)

# Create complete data with NAs for missing periods
all_combinations_paper <- expand.grid(
  co2 = unique(figure3_data_paper$co2),
  LF = c("0n_", "1c_", "2be", "3bm", "4bl", "5i_")
)

figure3_data_complete_paper <- all_combinations_paper %>%
  left_join(figure3_data_paper, by = c("co2", "LF")) %>%
  mutate(
    period = factor(LF, 
                    levels = c("0n_", "1c_", "2be", "3bm", "4bl", "5i_"),
                    labels = c("Neolithic/\nMesolithic", "Chalcolithic", "Early\nBronze", 
                               "Middle\nBronze", "Late\nBronze", "Iron Age")),
    region = factor(co2,
                    levels = c("ir", "l_", "ms", "tr", "ME"),
                    labels = c("Iran", "Levant", "Mesopotamia", "Turkey", "Middle East"))
  )

# Create line plot
ggplot(figure3_data_complete_paper, aes(x = period, y = tr_nop, color = region, group = region)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3.5, alpha = 0.9, data = figure3_data_complete_paper %>% filter(!is.na(tr_nop))) +
  scale_color_brewer(palette = "Set1") +
  labs(x = NULL, 
       y = "Violence share (%)", 
       color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 11)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.30, 0.05),
    labels = function(x) paste0(x, "%")
  ) +
  ggtitle("Figure 2 as displayed in the manuscript") +
  labs(subtitle = "Inclusion of Turkey-Late Bronze (n = 9) but exclusion of all other cases with less than 30 observations")


###############################################################################
##### Figure 2 with exclusion of all cases with less than 30 observations #####
###############################################################################

# As displayed in supplementary figure 2

# Collapse by region and time period
collapsed_data_exclusion <- viol_data %>%
  group_by(co2, LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop')

# Drop observations with low case numbers
collapsed_data_exclusion <- collapsed_data_exclusion %>%
  filter(!(co2 == "ir" & LF == "2be"),   # Iran EBA (n = 2)
         !(co2 == "l_" & LF == "4bl"),   # Levant LBA (n = 5)
         !(co2 == "ms" & LF == "1c_"),   # Mesopotamia Chalcolithic (n = 12)
         !(co2 == "tr" & LF == "4bl"))   # Turkey LBA (n = 9)

# Calculate Middle East average (all regions combined)
middle_east_exclusion <- viol_data %>%
  filter(!(co2 == "ir" & LF == "2be"),
         !(co2 == "l_" & LF == "4bl"),
         !(co2 == "ms" & LF == "1c_"),
         !(co2 == "tr" & LF == "4bl")) %>%
  group_by(LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop') %>%
  mutate(co2 = "ME")

# Combine individual regions with Middle East summary
figure3_data_exclusion <- bind_rows(collapsed_data_exclusion, middle_east_exclusion)

# Create complete data with NAs for missing periods
all_combinations_exclusion <- expand.grid(
  co2 = unique(figure3_data_exclusion$co2),
  LF = c("0n_", "1c_", "2be", "3bm", "4bl", "5i_")
)

figure3_data_complete_exclusion <- all_combinations_exclusion %>%
  left_join(figure3_data_exclusion, by = c("co2", "LF")) %>%
  mutate(
    period = factor(LF, 
                    levels = c("0n_", "1c_", "2be", "3bm", "4bl", "5i_"),
                    labels = c("Neolithic/\nMesolithic", "Chalcolithic", "Early\nBronze", 
                               "Middle\nBronze", "Late\nBronze", "Iron Age")),
    region = factor(co2,
                    levels = c("ir", "l_", "ms", "tr", "ME"),
                    labels = c("Iran", "Levant", "Mesopotamia", "Turkey", "Middle East"))
  )

# Create line plot
ggplot(figure3_data_complete_exclusion, aes(x = period, y = tr_nop, color = region, group = region)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3.5, alpha = 0.9, data = figure3_data_complete_exclusion %>% filter(!is.na(tr_nop))) +
  scale_color_brewer(palette = "Set1") +
  labs(x = NULL, 
       y = "Violence share (%)", 
       color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 11)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.30, 0.05),
    labels = function(x) paste0(x, "%")
  ) +
  ggtitle("Figure 2") +
  labs(subtitle = "Exclusion of all cases with less than 30 observations")

# Comparison to paper: 
# The rise in violence share during Late Bronze in the Middle East becomes slightly less pronounced due to the exclusion of Turkey-LBA
# Although the absolute value is only slightly affected, the visual perception of the graph changes noticeably 


#######################################
##### Figure 2 without exclusions #####
#######################################

# Collapse by region and time period
collapsed_data_complete <- viol_data %>%
  group_by(co2, LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop')

# Calculate Middle East average (all regions combined)
middle_east_complete <- viol_data %>%
  group_by(LF) %>%
  summarise(tr_nop = weighted.mean(tr_nop, w = ncases, na.rm = TRUE), .groups = 'drop') %>%
  mutate(co2 = "ME")

# Combine individual regions with Middle East summary
figure3_data_complete <- bind_rows(collapsed_data_complete, middle_east_complete)

# Create complete data with NAs for missing periods
all_combinations_complete <- expand.grid(
  co2 = unique(figure3_data_complete$co2),
  LF = c("0n_", "1c_", "2be", "3bm", "4bl", "5i_")
)

figure3_data_complete_with_na <- all_combinations_complete %>%
  left_join(figure3_data_complete %>% select(co2, LF, tr_nop), by = c("co2", "LF")) %>%
  mutate(
    period = factor(LF, 
                    levels = c("0n_", "1c_", "2be", "3bm", "4bl", "5i_"),
                    labels = c("Neolithic/\nMesolithic", "Chalcolithic", "Early\nBronze", 
                               "Middle\nBronze", "Late\nBronze", "Iron Age")),
    region = factor(co2,
                    levels = c("ir", "l_", "ms", "tr", "ME"),
                    labels = c("Iran", "Levant", "Mesopotamia", "Turkey", "Middle East"))
  )

# Create line plot
ggplot(figure3_data_complete_with_na, aes(x = period, y = tr_nop, color = region, group = region)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3.5, alpha = 0.9, data = figure3_data_complete_with_na %>% filter(!is.na(tr_nop))) +
  scale_color_brewer(palette = "Set1") +
  labs(x = NULL, 
       y = "Violence share (%)", 
       color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 11)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.30, 0.05),
    labels = function(x) paste0(x, "%")
  ) +  
  ggtitle("Figure 2") +
  labs(subtitle = "No exclusions")

# Comparison to paper: 
# Including all observations reveals the contradictions that motivated the authors' exclusion decisions
# Violence share is 0% in Mesopotamia-Chalcolithic (n = 12), which directly contradicts the authors' claim that "interpersonal violence peaked during the Chalcolithic period"
# Violence share is also 0% in Levant-LBA (n = 5), inconsistent with the stated increase during the Bronze Age
# It remains unclear whether these low-case observations reflect true (regional) variation or are artifacts of small sample sizes and selectivity bias
# The authors' exclusion criteria (keeping Turkey-LBA but dropping Mesopotamia-Chalcolithic) create a "cleaner" narrative that aligns with their main claims
# But also effectively removes contradictory evidence while retaining supporting evidence, even when both have similarly small sample sizes (n = 9 vs n = 12)
# Overall values are only slightly affected, but the visual perception of the figure and narrative coherence change noticeably
