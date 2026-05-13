### SPRUCE Litter Decomposition Study 2015 - 2018
## Bulk Characterization

library(tidyverse)
library(ggpubr)
library(rstatix)
library(forcats)
library(readxl)

# Read data ---------------------------------------------------------------
df <- read_xlsx("Input/SPRUCE_decomposition_chem_data.xlsx")

# Litter orderings --------------------------------------------------------
decomp_order  <- c("MAG", "SPR", "LTR", "ANG", "SPL", "LTL")  # by k
litter_order  <- c("MAG", "ANG", "SPR", "SPL", "LTR", "LTL")  # taxonomic-ish
df$Litter     <- factor(df$Litter, levels = decomp_order)

# Color palette (consistent across figures)
my_colors <- c(
  "LTL" = "#e6ab02",
  "LTR" = "#d95f02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291",
  "SPR" = "#1b9e77"
)

custom_labels <- c(
  "MAG" = "S. magellanicum",
  "ANG" = "S. angustifolium",
  "SPR" = "Spruce roots",
  "LTR" = "Labrador tea roots",
  "LTL" = "Labrador tea leaves",
  "SPL" = "Spruce needles"
)

# Decomposition Rate ------------------------------------------------------

data_k <- df %>%
  select(Litter, Temp, CO2, Pickup_t, percent_mass) %>%
  filter(Pickup_t != "T_0") %>%
  pivot_wider(names_from = Pickup_t, values_from = percent_mass)

data_k <- data_k %>%
  mutate(
    M0 = 100,
    #k0.5 = -log(T_0.5 / M0) / 252,   # 0–0.5 yr (not used here)
    #k1   = -log(T_1   / M0) / 384,   # 0–1 yr (not used here)
    k2   = -log(T_2   / M0) / 740    # 0–2 yr
  )

# Long format for rate comparisons ----------------------------------------
data_long <- data_k %>%
  pivot_longer(cols = starts_with("k"),
               names_to = "Rate",
               values_to = "Decomposition") %>%
  mutate(
    Litter = factor(Litter, levels = decomp_order),
    Rate   = factor(Rate, levels = c("k0.5", "k1", "k2"))
  )

# C:N ratio at T0 ---------------------------------------------------------
meta_with_cn <- df %>%
  filter(Pickup_t == "T_0") %>%
  select(SampleID, Litter, Plot, Temp, CO2, percent_C, percent_N) %>%
  mutate(
    C_N_ratio = percent_C / percent_N,
    Litter    = factor(Litter, levels = decomp_order)
  )


# PANEL A (for Sup Fig 1): Temp × CO2 × Litter line plot of k2 ------

fig1A <- ggplot(
  data_k,
  aes(x = Temp, y = k2,
      group = interaction(Litter, CO2),
      color = as.factor(CO2))
) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ Litter,
             labeller = labeller(Litter = as_labeller(custom_labels))) +
  labs(
    x     = "Temperature (°C)",
    y     = "Decomposition rate (k/day)",
    color = "CO₂ treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )


# EFFECT SIZE: export as TABLE ------

# Prepare data for effect-size model
es_dat <- data_k %>%
  filter(!is.na(k2)) %>%
  mutate(
    Temp_f = factor(Temp),
    CO2_f  = factor(CO2)
  )

# Linear model with k2
mod_k2 <- lm(k2 ~ Litter + Temp_f + CO2_f, data = es_dat)

# Partial eta^2 with 95% CIs
eta_tbl <- effectsize::eta_squared(
  mod_k2,
  partial = TRUE,
  ci = 0.95
) %>%
  as_tibble() %>%
  filter(Parameter %in% c("Litter", "Temp_f", "CO2_f")) %>%
  mutate(
    Parameter_clean = recode(Parameter,
                             "Litter" = "Litter",
                             "Temp_f" = "Temp",
                             "CO2_f"  = "CO₂"),
    Parameter_clean = factor(Parameter_clean,
                             levels = c("Litter", "Temp", "CO₂"))
  )

# Create a clean table for export (Supplementary Table)
eta_table_export <- eta_tbl %>%
  select(Parameter_clean,
         Eta2_partial,
         CI_low,
         CI_high) %>%
  arrange(desc(Eta2_partial))

# Export as CSV 
readr::write_csv(
  eta_table_export,
  "Tables/SuppTable_effect_sizes_k_partial_eta2.csv"
)


# PANEL 1A: Decomposition rate boxplot (by litter) -----

fig1C <- ggplot(
  data_long %>% filter(Rate == "k2"),
  aes(x = fct_reorder(Litter, Decomposition, .fun = median, na.rm = TRUE),
      y = Decomposition,
      fill = Litter)
) +
  geom_boxplot(alpha = 0.85) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = my_colors) +
  labs(
    x = NULL,
    y = "Decomposition rate (k,/day)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 25, vjust = 1, hjust = 1),
    legend.position = "none"
  )


# PANEL 1B C:N ratio boxplot---------

fig1D <- ggplot(meta_with_cn,
                aes(x = Litter, y = C_N_ratio, fill = Litter)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  # geom_jitter(width = 0.2, alpha = 0.5,
  #             color = "black", size = 1.8) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = custom_labels) +
  labs(
    x = NULL,
    y = "Initial C:N ratio"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 25, vjust = 1, hjust = 1),
    legend.position = "none"
  )


# Main Fig 1 ------

fig1_main <- ggpubr::ggarrange(
  fig1C, fig1D,
  labels = c("A", "B"),
  ncol   = 2,
  nrow   = 1,
  widths = c(1, 1)
)

fig1_main

ggsave(
  "Plots/Figure1_main_k_CNR.png",
  plot   = fig1_main,
  dpi    = 300,
  width  = 7.5,
  height = 3.8
)

# SUPPLEMENTARY FIGURE 1 ------
#   A: Temp × CO₂ × Litter, k2
#   B: initial %C
#   C: initial %N
#   D: initial %P


# Initial averages per litter ---------------------------------------------
df_initial <- df %>%
  filter(Pickup_t == "T_0")

average_values <- df_initial %>%
  group_by(Litter) %>%
  summarise(
    across(c("percent_mass", "percent_C", "percent_N",
             "percent_P", "percent_carbs", "percent_aromatics"),
           mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename_with(~ paste0("Initial_", .), -Litter)

df_merged <- df %>%
  left_join(average_values, by = "Litter")

# Litter ordering by median k2 --------------------------------------------
k2_order <- data_long %>%
  filter(Rate == "k2") %>%
  group_by(Litter) %>%
  summarise(med_k2 = median(Decomposition, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(med_k2)) %>%
  pull(Litter)

df_T0 <- df %>%
  filter(Pickup_t == "T_0") %>%
  mutate(Litter = factor(Litter, levels = rev(k2_order)))

summary_stats <- function(df, column) {
  df %>%
    group_by(Litter) %>%
    summarise(
      Mean = mean(.data[[column]], na.rm = TRUE),
      SE   = sd(.data[[column]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

stats_C <- summary_stats(df_T0, "percent_C")
stats_N <- summary_stats(df_T0, "percent_N")
stats_P <- summary_stats(df_T0, "percent_P")

p_init_C <- ggplot(stats_C,
                   aes(x = Litter, y = Mean, fill = Litter)) +
  geom_col(width = 0.7, alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0.2) +
  scale_fill_manual(values = my_colors) +
  labs(x = NULL, y = "%C") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1),
        legend.position = "none")

p_init_N <- ggplot(stats_N,
                   aes(x = Litter, y = Mean, fill = Litter)) +
  geom_col(width = 0.7, alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0.2) +
  scale_fill_manual(values = my_colors) +
  labs(x = NULL, y = "%N") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1),
        legend.position = "none")

p_init_P <- ggplot(stats_P,
                   aes(x = Litter, y = Mean, fill = Litter)) +
  geom_col(width = 0.7, alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0.2) +
  scale_fill_manual(values = my_colors) +
  labs(x = NULL, y = "%P") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1),
        legend.position = "none")

# --- Supplementary Figure 1 layout update ---
# Row 1: panel A (fig1A)
# Row 2: panels B–D (p_init_C, p_init_N, p_init_P)

row2_supp1 <- ggpubr::ggarrange(
  p_init_C, p_init_N, p_init_P,
  ncol   = 3,
  nrow   = 1,
  labels = c("B", "C", "D")  # bottom row labels
)

supp_fig1 <- ggpubr::ggarrange(
  fig1A,          # Panel A on its own row
  row2_supp1,     # Panels B–D on second row
  ncol    = 1,
  nrow    = 2,
  heights = c(1.5, 1)       
)

supp_fig1

ggsave("Plots/SupFig1_tempCO2_k_plus_initial_CNP.png",
       plot   = supp_fig1,
       dpi    = 300,
       width  = 7.5,
       height = 6.5)



# SUPPLEMENTARY FIGURE 2 – Percent Change  ------------

# Percent Change ----------------------------------------------------------
df_initial <- df %>% filter(Pickup_t == "T_0")

average_values <- df_initial %>%
  group_by(Litter) %>%
  summarise(across(c("percent_mass", "percent_C", "percent_N",
                     "percent_P", "percent_carbs", "percent_aromatics"),
                   mean, na.rm = TRUE),
            .groups = 'drop')

average_values <- average_values %>%
  rename_with(~ paste0("Initial_", .), -Litter)

df_merged <- df %>%
  left_join(average_values, by = "Litter")

for (col in c("percent_mass", "percent_C", "percent_N",
              "percent_P", "percent_carbs", "percent_aromatics")) {
  initial_col <- paste0("Initial_", col)
  change_col  <- paste0(col, "_change")
  
  df_merged[[change_col]] <-
    ((df_merged[[col]] - df_merged[[initial_col]]) /
       df_merged[[initial_col]]) * 100
}

# Long format
df_long <- df_merged %>%
  select(Litter, Pickup_t, ends_with("_change")) %>%
  pivot_longer(cols = ends_with("_change"),
               names_to = "Metric",
               values_to = "Change") %>%
  filter(Pickup_t != "T_0") %>%
  mutate(
    Litter   = factor(Litter, levels = decomp_order),
    Pickup_t = factor(Pickup_t, levels = c("T_0.5", "T_1", "T_2"))
  )

### C/N/P
df_long_CNP <- df_long %>%
  filter(Metric %in% c("percent_C_change",
                       "percent_N_change",
                       "percent_P_change")) %>%
  mutate(
    Metric = recode(Metric,
                    "percent_C_change" = "C",
                    "percent_N_change" = "N",
                    "percent_P_change" = "P")
  )

supp2_CNP <- ggplot(df_long_CNP,
                    aes(x = Litter, y = Change, fill = Litter)) +
  geom_boxplot(alpha = 0.9, outlier.size = 0.6) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(
    x = "Litter type",
    y = "Percent change from initial (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 25, vjust = 1, hjust = 1),
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )

# Bulk carbohydrates and aromatics line plot --------
litter_order_bulk  <- c("ANG", "MAG", "SPL", "LTL")

my_colors_bulk <- c(
  "LTL" = "#e6ab02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291"
)

litter_labels_bulk <- c(
  "ANG" = "ANG (S. angustifolium)",
  "MAG" = "MAG (S. magellanicum)",
  "SPL" = "SPL (Spruce needles)",
  "LTL" = "LTL (Labrador tea leaves)"
)

parse_time <- function(x) as.numeric(gsub("^T_", "", x))

Bulk_df <- df %>%
  select(SampleID, Name, Litter, Pickup_t,
         percent_carbs, percent_aromatics) %>%
  mutate(
    Pickup_t = factor(
      Pickup_t,
      levels  = c("T_0", "T_0.5", "T_1", "T_2"),
      ordered = TRUE
    ),
    t_num = parse_time(as.character(Pickup_t))
  ) %>%
  filter(percent_carbs != 0,
         percent_aromatics != 0) %>%
  pivot_longer(
    c(percent_carbs, percent_aromatics),
    names_to  = "Bulk_measure",
    values_to = "percent"
  ) %>%
  mutate(
    Litter       = factor(Litter, levels = litter_order_bulk),
    Bulk_measure = recode(
      Bulk_measure,
      "percent_carbs"     = "Carbohydrates",
      "percent_aromatics" = "Aromatics"
    )
  ) %>%
  filter(!Litter %in% c("SPR", "LTR"))

sum_bulk <- Bulk_df %>%
  group_by(Bulk_measure, Litter, Pickup_t, t_num) %>%
  summarise(
    n            = dplyr::n(),
    mean_measure = mean(percent, na.rm = TRUE),
    sd_measure   = sd(percent,   na.rm = TRUE),
    se_measure   = sd_measure / sqrt(n),
    .groups      = "drop"
  )%>%
  filter(!is.na(Litter))

facet_labels_bulk <- c(
  "Carbohydrates" = "Bulk carbohydrates (%)",
  "Aromatics"     = "Bulk aromatics (%)"
)

supp2_bulk <- ggplot(
  sum_bulk,
  aes(x = t_num, y = mean_measure,
      color = Litter, fill = Litter, group = Litter)
) +
  geom_smooth(
    method = "loess",
    span   = 1.2,
    se     = FALSE,
    linewidth = 1.5
  ) +
  geom_point(size = 3) +
  facet_wrap(
    ~ Bulk_measure,
    scales   = "free_y",
    labeller = as_labeller(facet_labels_bulk)
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 2),
    labels = c("0", "0.5", "1", "2")
  ) +
  labs(
    x = "Collection time (years)",
    y = "Bulk percentage (%)"
  ) +
  scale_color_manual(
    values = my_colors_bulk,
    limits = litter_order_bulk,
    breaks = litter_order_bulk,
    labels = litter_labels_bulk,
    drop   = FALSE,
    name   = "Litter type"
  ) +
  scale_fill_manual(
    values = my_colors_bulk,
    limits = litter_order_bulk,
    breaks = litter_order_bulk,
    labels = litter_labels_bulk,
    drop   = FALSE,
    name   = "Litter type"
  ) +
  guides(color = guide_legend(order = 1),
         fill  = guide_legend(order = 1)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

# Combined Supplementary Figure 2 ----------------------------------------
supp_fig2 <- ggpubr::ggarrange(
  supp2_CNP, supp2_bulk,
  ncol    = 1,
  nrow    = 2,
  labels  = c("A", "B"),
  heights = c(1, 1)
)
supp_fig2


ggsave("Plots/SuppFig2_percent_change_CNP_and_bulk_carbs_aromatics.png",
       plot = supp_fig2, dpi = 300, width = 7, height = 5.5)
