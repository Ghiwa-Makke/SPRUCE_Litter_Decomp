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









# Main figure: initial molecular indices vs. decomposition rate ------
library(tidyverse)

dir.create("Index", showWarnings = FALSE, recursive = TRUE)

my_colors <- c(
  "LTL" = "#e6ab02",
  "LTR" = "#d95f02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291",
  "SPR" = "#1b9e77"
)

custom_labels_expr <- list(
  MAG = bquote(italic("S. magellanicum")),
  ANG = bquote(italic("S. angustifolium")),
  SPR = "Spruce roots",
  LTR = "Labrador tea roots",
  LTL = "Labrador tea leaves",
  SPL = "Spruce needles"
)

# Preserve the litter order used in the original analysis
decomp_order <- c("MAG", "SPR", "LTR", "ANG", "SPL", "LTL")

legend_labels <- as.expression(
  unname(custom_labels_expr[decomp_order])
)

panel_names <- c(
  "C_N_ratio"            = "C:N ratio",
  "NOSC_mean"            = "NOSC",
  "DBE_mean"             = "DBE",
  "DBE_O_mean"           = "DBE-O",
  "AI_recalculated_mean" = "AI",
  "AI_mod_mean"          = "AI-mod"
)

# Prepare litter means; convert k from day⁻¹ to year⁻¹
plot_data <- indices_k2 %>%
  mutate(
    Litter = factor(Litter, levels = decomp_order),
    k_year = mean_K2 * 365
  ) %>%
  select(Litter, k_year, all_of(names(panel_names))) %>%
  pivot_longer(
    cols = all_of(names(panel_names)),
    names_to = "Index",
    values_to = "Initial_value"
  ) %>%
  filter(
    is.finite(Initial_value),
    is.finite(k_year)
  ) %>%
  mutate(
    Panel = factor(
      unname(panel_names[Index]),
      levels = unname(panel_names)
    )
  )

# Fit regressions and prepare annotations
plot_stats <- plot_data %>%
  group_by(Panel) %>%
  group_modify(~ {
    fit <- lm(k_year ~ Initial_value, data = .x)
    
    tibble(
      R2 = summary(fit)$r.squared,
      Slope = unname(coef(fit)["Initial_value"]),
      p = coef(summary(fit))["Initial_value", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_text = if_else(
      p < 0.001,
      "p < 0.001",
      paste0("p = ", formatC(p, format = "f", digits = 3))
    ),
    slope_text = formatC(Slope, format = "g", digits = 3),
    Label = paste0(
      "R\u00B2 = ", formatC(R2, format = "f", digits = 3),
      "\nSlope = ", slope_text,
      "\n", p_text
    )
  )

# Plot
p_indices <- ggplot(
  plot_data,
  aes(x = Initial_value, y = k_year)
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    level = 0.95,
    color = "grey25",
    fill = "grey75",
    linewidth = 0.7,
    alpha = 0.35
  ) +
  geom_point(
    aes(fill = Litter),
    shape = 21,
    size = 3.5,
    color = "white",
    stroke = 0.5
  ) +
  geom_text(
    data = plot_stats,
    aes(x = -Inf, y = Inf, label = Label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.12,
    size = 3.3,
    lineheight = 1.05
  ) +
  facet_wrap(
    ~ Panel,
    ncol = 3,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = my_colors,
    breaks = decomp_order,
    labels = legend_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.30))
  ) +
  labs(
    x = "Initial litter-mean index value",
    y = expression(
      "Mean decomposition constant, " * k ~ "(" * year^{-1} * ")"
    ),
    fill = "Litter type"
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing = grid::unit(1.2, "lines")
  )

print(p_indices)

# Save figure
ggsave(
  "Index/Initial_indices_decomposition_main_figure.png",
  plot = p_indices,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "Plots_frontier/tiff/Figure2_index.tiff",
  p_indices,
  width       = 180,
  height      = 130,
  units       = "mm",
  dpi         = 300,
  device      = agg_tiff,
  compression = "lzw"
)



# Save plotted values and regression statistics
write_csv(
  plot_data,
  "Index/Initial_indices_decomposition_plot_data.csv"
)

write_csv(
  plot_stats,
  "Index/Initial_indices_decomposition_plot_statistics.csv"
)




# Quantify the DBE-O relationship across the observed litter range ------

# Use the same complete observations as the DBE-O regression
dbeo_data <- indices_k2 %>%
  select(Litter, DBE_O_mean, mean_K2) %>%
  filter(
    is.finite(DBE_O_mean),
    is.finite(mean_K2)
  )

# Show the underlying means for each litter type
print(dbeo_data, n = Inf)

# Fit the original regression: k in day⁻¹
dbeo_fit <- lm(mean_K2 ~ DBE_O_mean, data = dbeo_data)

slope_day <- unname(coef(dbeo_fit)["DBE_O_mean"])
slope_year <- slope_day * 365

# Observed ranges
dbeo_min <- min(dbeo_data$DBE_O_mean)
dbeo_max <- max(dbeo_data$DBE_O_mean)
dbeo_span <- dbeo_max - dbeo_min

k_min <- min(dbeo_data$mean_K2)
k_max <- max(dbeo_data$mean_K2)
k_span <- k_max - k_min

stopifnot(k_span > 0)

# Magnitude of fitted change across the observed DBE-O range
predicted_k_difference <- abs(slope_day) * dbeo_span

percent_k_range <- 100 * predicted_k_difference / k_span

# Also calculate using the rounded slope requested
difference_using_rounded_slope <- 0.00055 * dbeo_span
percent_using_rounded_slope <- (
  100 * difference_using_rounded_slope / k_span
)

range_summary <- tibble(
  n_litter_types = nrow(dbeo_data),
  DBE_O_min = dbeo_min,
  DBE_O_max = dbeo_max,
  DBE_O_span = dbeo_span,
  k_min_day = k_min,
  k_max_day = k_max,
  k_span_day = k_span,
  fitted_slope_day = slope_day,
  fitted_slope_year = slope_year,
  predicted_k_difference_day = predicted_k_difference,
  percent_of_observed_k_range = percent_k_range,
  difference_using_rounded_slope = difference_using_rounded_slope,
  percent_using_rounded_slope = percent_using_rounded_slope
)

print(range_summary, width = Inf)


