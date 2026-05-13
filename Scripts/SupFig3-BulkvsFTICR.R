### Supplementary Figure 3 - Carbohydrate and Aromatic/Lignin Percent Change
### Bulk vs FT-ICR MS Comparison

## Libraries ------
library(tidyverse)
library(readxl)
library(ggpubr)

# Settings ------
decomp_order <- c("MAG", "LTR", "SPR", "ANG", "LTL", "SPL")

my_colors <- c(
  "LTL" = "#e6ab02",
  "LTR" = "#d95f02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291",
  "SPR" = "#1b9e77"
)

# Input Data ------
df <- read_xlsx("Input/SPRUCE_decomposition_chem_data.xlsx")
df$Litter <- factor(df$Litter, levels = decomp_order)

fticr_raw <- read.csv("Input/Metabodirect_All_Sum_Time/1_preprocessing_output/Report_processed_MolecFormulas.csv")
metadata  <- read.csv("Input/metadata.csv") %>%
  select(SampleID, Pickup_t, Litter)

# PANEL A: Carbohydrate Percent Change (Bulk vs FT-ICR)  ------

# 1. Bulk carbohydrate percent change
bulk_initial <- df %>%
  filter(Pickup_t == "T_0") %>%
  group_by(Litter) %>%
  summarise(initial_carbs = mean(percent_carbs, na.rm = TRUE), .groups = "drop")

bulk_carbs <- df %>%
  left_join(bulk_initial, by = "Litter") %>%
  mutate(
    carbs_change = ((percent_carbs - initial_carbs) / initial_carbs) * 100,
    Pickup_t     = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2")),
    Litter       = factor(Litter, levels = decomp_order),
    Source       = "Bulk"
  ) %>%
  select(Litter, Pickup_t, carbs_change, Source)

# 2. FT-ICR MS carbohydrate-like relative abundance percent change
sample_cols <- names(fticr_raw)[str_starts(names(fticr_raw), "Kelly_")]

fticr_carbs <- fticr_raw %>%
  filter(Class == "Carbohydrate") %>%
  select(all_of(sample_cols)) %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "SampleID",
    values_to = "intensity"
  ) %>%
  filter(intensity > 0) %>%
  group_by(SampleID) %>%
  summarise(carbs_intensity = sum(intensity, na.rm = TRUE), .groups = "drop")

# Total intensity per sample for relative abundance
fticr_total <- fticr_raw %>%
  select(all_of(sample_cols)) %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "SampleID",
    values_to = "intensity"
  ) %>%
  filter(intensity > 0) %>%
  group_by(SampleID) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE), .groups = "drop")

fticr_carbs_rel <- fticr_carbs %>%
  left_join(fticr_total, by = "SampleID") %>%
  mutate(carbs_rel = carbs_intensity / total_intensity * 100) %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2")),
    Litter   = factor(Litter, levels = decomp_order)
  )

# Percent change from T_0
fticr_carbs_initial <- fticr_carbs_rel %>%
  filter(Pickup_t == "T_0") %>%
  group_by(Litter) %>%
  summarise(initial_carbs_rel = mean(carbs_rel, na.rm = TRUE), .groups = "drop")

fticr_carbs_change <- fticr_carbs_rel %>%
  left_join(fticr_carbs_initial, by = "Litter") %>%
  mutate(
    carbs_change = ((carbs_rel - initial_carbs_rel) / initial_carbs_rel) * 100,
    Source       = "FT-ICR MS"
  ) %>%
  select(Litter, Pickup_t, carbs_change, Source)

# 3. Combine bulk and FT-ICR carbohydrate data
carbs_combined <- bind_rows(bulk_carbs, fticr_carbs_change) %>%
  mutate(Source = factor(Source, levels = c("Bulk", "FT-ICR MS")))
# Time point labels
time_labels <- c(
  "T_0.5" = "0.5 year",
  "T_1"   = "1 year",
  "T_2"   = "2 years"
)

# 4. Panel A plot
# Panel A
panel_A <- ggplot(
  carbs_combined %>% filter(Pickup_t != "T_0"),
  aes(x = Litter, y = carbs_change, fill = Litter)
) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  facet_grid(Source ~ Pickup_t, scales = "free_y",
             labeller = labeller(Pickup_t = time_labels)) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Carbohydrate Percent Change: FT-ICR vs Bulk",
    x     = "Litter",
    y     = "Percent Change (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    panel.grid       = element_blank()
  )

panel_A

# PANEL B: Aromatics (Bulk) vs Lignin (FT-ICR) Percent Change -----

# 1. Bulk aromatics percent change
bulk_arom_initial <- df %>%
  filter(Pickup_t == "T_0") %>%
  group_by(Litter) %>%
  summarise(initial_arom = mean(percent_aromatics, na.rm = TRUE), .groups = "drop")

bulk_arom <- df %>%
  left_join(bulk_arom_initial, by = "Litter") %>%
  mutate(
    arom_change = ((percent_aromatics - initial_arom) / initial_arom) * 100,
    Pickup_t    = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2")),
    Litter      = factor(Litter, levels = decomp_order),
    Source      = "Aromatics (Bulk)"
  ) %>%
  select(Litter, Pickup_t, arom_change, Source)

# 2. FT-ICR MS lignin-like relative abundance percent change
fticr_lignin <- fticr_raw %>%
  filter(Class == "Lignin") %>%
  select(all_of(sample_cols)) %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "SampleID",
    values_to = "intensity"
  ) %>%
  filter(intensity > 0) %>%
  group_by(SampleID) %>%
  summarise(lignin_intensity = sum(intensity, na.rm = TRUE), .groups = "drop")

fticr_lignin_rel <- fticr_lignin %>%
  left_join(fticr_total, by = "SampleID") %>%
  mutate(lignin_rel = lignin_intensity / total_intensity * 100) %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2")),
    Litter   = factor(Litter, levels = decomp_order)
  )

fticr_lignin_initial <- fticr_lignin_rel %>%
  filter(Pickup_t == "T_0") %>%
  group_by(Litter) %>%
  summarise(initial_lignin_rel = mean(lignin_rel, na.rm = TRUE), .groups = "drop")

fticr_lignin_change <- fticr_lignin_rel %>%
  left_join(fticr_lignin_initial, by = "Litter") %>%
  mutate(
    arom_change = ((lignin_rel - initial_lignin_rel) / initial_lignin_rel) * 100,
    Source      = "Lignin (FT-ICR MS)"
  ) %>%
  select(Litter, Pickup_t, arom_change, Source)

# 3. Combine bulk aromatics and FT-ICR lignin
arom_combined <- bind_rows(bulk_arom, fticr_lignin_change) %>%
  mutate(Source = factor(Source, levels = c("Aromatics (Bulk)", "Lignin (FT-ICR MS)")))

# 4. Panel B plot
panel_B <- ggplot(
  arom_combined %>% filter(Pickup_t != "T_0"),
  aes(x = Litter, y = arom_change, fill = Litter)
) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  facet_grid(Source ~ Pickup_t, scales = "free_y",
             labeller = labeller(Pickup_t = time_labels)) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Aromatics (Bulk) vs Lignin (FT-ICR) Percent Change",
    x     = "Litter",
    y     = "Percent Change (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    panel.grid       = element_blank()
  )

panel_B

# Combine panels and save ------
supp_fig <- ggpubr::ggarrange(
  panel_A, panel_B,
  labels  = c("A", "B"),
  ncol    = 1,
  nrow    = 2
)

supp_fig

ggsave(
  "Plots/SuppFig3_Carbs_Arom_Bulk_vs_FTICR.png",
  supp_fig,
  width  = 7,
  height = 7,
  dpi    = 300,
  bg     = "white"
)
