### SPRUCE Litter Decomposition Study 2015 - 2018
## Supplementary Table 2 - Index prediction of decomposition rate 

## Libraries ------
library(tidyverse)
library(readxl)
library(kableExtra)

# Settings ------
decomp_order <- c("MAG", "SPR", "LTR", "ANG", "SPL", "LTL")  # slowest → fastest

# Input Data -----
df <- read_xlsx("Input/SPRUCE_decomposition_chem_data.xlsx")
df$Litter <- factor(df$Litter, levels = decomp_order)

fticr_raw <- read.csv("Input/Metabodirect_All_Sum_Time/1_preprocessing_output/Report_processed_MolecFormulas.csv")

metadata <- read.csv("Input/metadata.csv") %>%
  select(SampleID, Pickup_t, Litter)

# Decomposition Rate Calculation ------
data_k <- df %>%
  select(Litter, Temp, CO2, Pickup_t, percent_mass) %>%
  filter(Pickup_t != "T_0") %>%
  pivot_wider(names_from = Pickup_t, values_from = percent_mass) %>%
  mutate(
    M0 = 100,
    k2 = -log(T_2 / M0) / 740    # 0–2 yr
  )

mean_k2 <- data_k %>%
  group_by(Litter) %>%
  summarise(
    mean_K2 = mean(k2, na.rm = TRUE),
    .groups = "drop"
  )

# C:N Ratio ------
cn_mean_litter <- df %>%
  filter(Pickup_t == "T_0") %>%
  mutate(C_N_ratio = percent_C / percent_N) %>%
  group_by(Litter) %>%
  summarise(
    C_N_ratio = mean(C_N_ratio, na.rm = TRUE),
    .groups = "drop"
  )

# Recalculate AI with boundary conditions ------
fticr_raw <- fticr_raw %>%
  mutate(
    AI_numerator   = 1 + C - O - S - 0.5 * (N + P + H),
    AI_denominator = C - O - N - S - P,
    AI_recalculated = ifelse(AI_numerator <= 0, 0,
                             ifelse(AI_denominator <= 0, 0,
                                    AI_numerator / AI_denominator))
  )

# Define indices from MetaboDirect output + AI_recalculated ------
indices <- c("NOSC", "GFE", "DBE", "DBE_O", "AI_mod", "AI_recalculated")

# Calculate sample-level mean indices ------
sample_cols <- names(fticr_raw)[str_starts(names(fticr_raw), "Kelly_")]

sample_indices <- fticr_raw %>%
  select(all_of(c(indices, sample_cols))) %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "SampleID",
    values_to = "intensity"
  ) %>%
  filter(intensity > 0) %>%
  group_by(SampleID) %>%
  summarise(
    across(
      all_of(indices),
      list(mean = ~ mean(.x, na.rm = TRUE)),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  )

sample_indices_meta <- sample_indices %>%
  left_join(metadata, by = "SampleID")

# Mean indices per litter at T_0 ------
indices_mean_litter <- sample_indices_meta %>%
  filter(Pickup_t == "T_0") %>%
  group_by(Litter) %>%
  summarise(
    across(
      ends_with("_mean"),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Combine all predictors with mean K2 ------
indices_k2 <- indices_mean_litter %>%
  left_join(cn_mean_litter, by = "Litter") %>%
  left_join(mean_k2, by = "Litter") %>%
  mutate(Litter = factor(Litter, levels = decomp_order))

# Regression function ------
calc_fit_stats <- function(data, index_col) {
  formula <- as.formula(paste("mean_K2 ~", index_col))
  model   <- lm(formula, data = data)
  
  data.frame(
    Index    = index_col,
    R2       = summary(model)$r.squared,
    Adj_R2   = summary(model)$adj.r.squared,
    Slope    = summary(model)$coefficients[2, 1],
    P_value  = summary(model)$coefficients[2, 4]
  )
}

# Run regressions for all indices ------
predictors <- c("C_N_ratio", "NOSC_mean", "GFE_mean", 
                "DBE_mean", "DBE_O_mean", "AI_mod_mean", "AI_recalculated_mean")

regression_table <- map_dfr(predictors, ~ calc_fit_stats(indices_k2, .x)) %>%
  mutate(
    Index_name = case_when(
      Index == "C_N_ratio"           ~ "C:N ratio",
      Index == "NOSC_mean"           ~ "NOSC",
      Index == "GFE_mean"            ~ "GFE",
      Index == "DBE_mean"            ~ "DBE",
      Index == "DBE_O_mean"          ~ "DBE-O",
      Index == "AI_mod_mean"         ~ "AI (Modified)",
      Index == "AI_recalculated_mean"~ "AI"
    ),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  arrange(desc(R2)) %>%
  select(Index_name, R2, Adj_R2, Slope, P_value, Significance)

# Display and save ------
print(regression_table)

write_csv(regression_table, "Index/indices_K2_regression_table.csv")


