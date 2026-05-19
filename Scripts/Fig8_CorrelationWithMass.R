### SPRUCE Litter decomposition 2015-2017
### Correlations of classes with mass loss
### Figure 8

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
  library(GGally)
  library(ragg)
})

# Data --------------------------------------------------------------------
meta <- read.csv("Input/metadata.csv") %>%
  transmute(
    SampleID,
    Litter,
    t_num = as.numeric(gsub("^T_", "", Pickup_t)),
    Temp = case_when(
      Temp == 0   ~ "+0°C",
      Temp == 4.5 ~ "+4.5°C",
      Temp == 9   ~ "+9°C",
      TRUE ~ as.character(Temp)
    ),
    CO2 = case_when(
      CO2 == "aCO2" ~ "Ambient",
      CO2 == "eCO2" ~ "Elevated",
      TRUE ~ as.character(CO2)
    ),
    percent_mass
  )

div <- read.csv("Tables/Class_Shannon_Diversity_PerSample.csv") %>%
  left_join(meta %>% select(SampleID, Temp, CO2), by = "SampleID") %>%
  transmute(
    SampleID, Litter, Temp, CO2, t_num,
    Class = recode(
      Class,
      "Amino sugar"             = "AminoSugar",
      "Condensed hydrocarbon"   = "CondensedHC",
      "Unsaturated hydrocarbon" = "Unsaturated_HC"
    ),
    Shannon
  )

temp_levels <- c("+0°C", "+4.5°C", "+9°C")
co2_levels  <- c("Ambient", "Elevated")

duplicate_T0_none <- function(df) {
  df %>%
    uncount(if_else(t_num == 0 & Temp == "None", length(temp_levels), 1), .id = "iT") %>%
    mutate(Temp = if_else(t_num == 0 & Temp == "None", temp_levels[iT], Temp)) %>%
    select(-iT) %>%
    uncount(if_else(t_num == 0 & CO2 == "None", length(co2_levels), 1), .id = "iC") %>%
    mutate(CO2 = if_else(t_num == 0 & CO2 == "None", co2_levels[iC], CO2)) %>%
    select(-iC) %>%
    mutate(
      Temp = factor(Temp, levels = temp_levels),
      CO2  = factor(CO2,  levels = co2_levels)
    )
}

div <- duplicate_T0_none(div)

mydata <- div %>%
  pivot_wider(names_from = Class, values_from = Shannon) %>%
  left_join(meta, by = c("SampleID", "Litter", "Temp", "CO2", "t_num")) %>%
  mutate(percent_mass = if_else(t_num == 0, 100, percent_mass))

# Slopes (weighted by n per timepoint) ------------------------------------
slope_table_all <- data.frame(
  Litter = character(),
  Temp = character(),
  CO2 = character(),
  Slope_percent_mass = numeric(),
  Slope_AminoSugar = numeric(),
  Slope_Carbohydrate = numeric(),
  Slope_CondensedHC = numeric(),
  Slope_Lignin = numeric(),
  Slope_Lipid = numeric(),
  Slope_Other = numeric(),
  Slope_Peptide = numeric(),
  Slope_Tannin = numeric(),
  Slope_Unsaturated_HC = numeric(),
  stringsAsFactors = FALSE
)

for (litter_type in unique(mydata$Litter)) {
  for (temp in unique(mydata$Temp)) {
    for (co2 in unique(mydata$CO2)) {
      
      subset_data <- mydata %>%
        filter(Litter == litter_type, Temp == temp, CO2 == co2) %>%
        filter(!is.na(t_num))
      
      if (nrow(subset_data) < 2 || dplyr::n_distinct(subset_data$t_num) < 2) next
      
      subset_data <- subset_data %>%
        group_by(t_num) %>%
        mutate(w = 1 / n()) %>%
        ungroup()
      
      slopes_list <- list()
      
      for (col_name in colnames(subset_data)[6:15]) {
        
        df_xy <- subset_data %>%
          select(t_num, y = all_of(col_name), w) %>%
          drop_na()
        
        if (nrow(df_xy) < 2 || dplyr::n_distinct(df_xy$t_num) < 2) {
          slopes_list[[paste0("Slope_", col_name)]] <- NA_real_
        } else {
          slopes_list[[paste0("Slope_", col_name)]] <-
            unname(coef(lm(y ~ t_num, data = df_xy, weights = w))[2])
        }
      }
      
      slope_table_all <- rbind(
        slope_table_all,
        data.frame(
          Litter = litter_type,
          Temp = temp,
          CO2 = co2,
          do.call(data.frame, slopes_list),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

rownames(slope_table_all) <- paste(
  slope_table_all$Litter,
  slope_table_all$Temp,
  slope_table_all$CO2,
  sep = "_"
)

write.csv(
  slope_table_all,
  "Tables/Slope_table_all_WLS.csv",
  row.names = TRUE
)

# Correlations -------------------------------------------------------------
data_matrix <- slope_table_all[, 4:13]

spearman_corr_matrix <- cor(
  data_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

p_values <- sapply(names(data_matrix), function(v) {
  xy <- data_matrix[, c("Slope_percent_mass", v)]
  xy <- xy[complete.cases(xy), , drop = FALSE]
  if (nrow(xy) < 3) return(NA_real_)
  suppressWarnings(
    cor.test(xy[[1]], xy[[2]], method = "spearman", exact = FALSE)$p.value
  )
})

mass_corr  <- spearman_corr_matrix["Slope_percent_mass", ]
mass_pvals <- p_values[names(mass_corr)]

mass_corr_df <- data.frame(
  Compound = names(mass_corr),
  Correlation = as.numeric(mass_corr),
  P_value = as.numeric(mass_pvals),
  stringsAsFactors = FALSE
) %>%
  filter(Compound != "Slope_percent_mass")

# Labels + colors ----------------------------------------------------------
my_colors <- c(
  "LTL" = "#e6ab02",
  "LTR" = "#d95f02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291",
  "SPR" = "#1b9e77"
)

clean_class_names <- function(x) {
  x %>%
    gsub("^Slope_", "", .) %>%
    gsub("Unsaturated_HC", "Unsaturated HC", .) %>%
    gsub("AminoSugar", "Amino Sugar", .) %>%
    gsub("CondensedHC", "Condensed HC", .)
}

mass_corr_df <- mass_corr_df %>%
  mutate(
    Compound = clean_class_names(Compound),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

compound_order <- c(
  "Tannin", "Lignin", "Carbohydrate", "Condensed HC", "Lipid",
  "Peptide", "Amino Sugar", "Other", "Unsaturated HC"
)

mass_corr_df <- mass_corr_df %>%
  mutate(Compound = factor(Compound, levels = compound_order))

# Panel A ------------------------------------------------------------------
bar_plot <- ggplot(
  mass_corr_df,
  aes(x = reorder(Compound, Correlation), y = Correlation, fill = Correlation)
) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(
    aes(label = Significance),
    nudge_y = ifelse(mass_corr_df$Correlation >= 0, 0.08, -0.08),
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#40B0A6",
    mid = "white",
    high = "yellow4",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_flip() +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.title  = element_text(size = 10),
    legend.position = "none",
    plot.title  = element_text(size = 10, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Universal",
    x = "Compound class",
    y = "Spearman correlation (ρ)"
  )

# Panel B ------------------------------------------------------------------
long_df <- slope_table_all %>%
  as_tibble() %>%
  mutate(Mass = Slope_percent_mass) %>%
  select(Litter, Temp, CO2, Mass, starts_with("Slope_")) %>%
  pivot_longer(
    cols = starts_with("Slope_"),
    names_to = "Class",
    values_to = "ClassSlope"
  ) %>%
  filter(Class != "Slope_percent_mass") %>%
  mutate(Class = clean_class_names(Class))

heat_df <- long_df %>%
  group_by(Litter, Class) %>%
  summarise(
    n = sum(complete.cases(Mass, ClassSlope)),
    rho = ifelse(
      n < 3,
      NA_real_,
      cor(Mass, ClassSlope, method = "spearman", use = "complete.obs")
    ),
    p = ifelse(n < 3, NA_real_, {
      xy <- cur_data() %>% select(Mass, ClassSlope) %>% drop_na()
      suppressWarnings(
        cor.test(xy$Mass, xy$ClassSlope, method = "spearman", exact = FALSE)$p.value
      )
    }),
    .groups = "drop"
  ) %>%
  group_by(Litter) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    Stars = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup() %>%
  mutate(
    Class  = factor(Class, levels = compound_order),
    Litter = factor(Litter, levels = names(my_colors))
  )

heat_plot <- ggplot(heat_df, aes(x = Litter, y = Class, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Stars), size = 4.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#40B0A6",
    mid = "white",
    high = "yellow4",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation\n(ρ)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.title  = element_text(size = 10),
    plot.title  = element_text(size = 10, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "By litter",
    x = "Litter",
    y = "Compound class"
  )

# Combine ------------------------------------------------------------------
figure_8 <- ggarrange(
  bar_plot,
  heat_plot,
  nrow = 1,
  labels = c("A", "B"),
  widths = c(1, 1.2)
)

figure_8


ggsave(
  "Plots/Updated/Figure8_MassLoss_vs_ClassCorrelations.tiff",
  figure_8,
  width       = 190,
  height      = 85,
  units       = "mm",
  dpi         = 300,
  device      = agg_tiff,
  compression = "lzw",
  bg          = "white"
)

# ggpairs ------------------------------------------------------------------
litter_order <- c("LTL", "LTR", "ANG", "MAG", "SPL", "SPR")

litter_labels_expr <- c(
  "LTL" = "LTL (Labrador tea leaves)",
  "LTR" = "LTR (Labrador tea fine roots)",
  "ANG" = expression("ANG (" * italic("S. angustifolium") * ")"),
  "MAG" = expression("MAG (" * italic("S. magellanicum") * ")"),
  "SPL" = "SPL (Spruce needles)",
  "SPR" = "SPR (Spruce fine roots)"
)

plot_pairs <- ggpairs(
  slope_table_all,
  columns = 4:13,
  ggplot2::aes(color = Litter),
  upper = list(
    continuous = wrap(
      "cor",
      method = "spearman",
      use = "pairwise.complete.obs",
      size = 3
    )
  )
) +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = my_colors,
    breaks = litter_order,
    labels = litter_labels_expr,
    name   = "Litter type"
  ) +
  scale_fill_manual(
    values = my_colors,
    breaks = litter_order,
    labels = litter_labels_expr,
    name   = "Litter type"
  )

ggsave(
  "Plots/ggpairs_plot.png",
  plot_pairs,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Plots/ggpairs_plot.tiff",
  plot_pairs,
  width       = 190,
  height      = 150,
  units       = "mm",
  dpi         = 300,
  device      = agg_tiff,
  compression = "lzw",
  bg          = "white"
)
