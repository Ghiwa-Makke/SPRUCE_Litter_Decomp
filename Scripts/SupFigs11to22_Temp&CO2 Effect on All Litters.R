### SPRUCE Litter decomposition 2015–2017
### Shannon diversity of compound classes (from All.csv intensity matrix)

### Sup Figs 11 to 22

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(janitor)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------- Helpers ---------------------------------
parse_time <- function(x) as.numeric(gsub("^T_", "", x))

time_colors <- c(
  "T_0"   = "#4DBBD5",
  "T_0.5" = "#3C5488",
  "T_1"   = "#B09C85",
  "T_2"   = "#503A2C"
)

tp_levels <- c("T_0", "T_0.5", "T_1", "T_2")
tp_labels <- c("T_0" = "0", "T_0.5" = "0.5", "T_1" = "1", "T_2" = "2")

temp_levels <- c("+0°C", "+4.5°C", "+9°C")
co2_levels  <- c("Ambient", "Elevated")

# Panel order 
class_levels <- c(
  "AminoSugar",
  "Carbohydrate",
  "CondensedHC",
  "Lignin",
  "Lipid",
  "Other",
  "Peptide",
  "Tannin",
  "Unsaturated_HC"
)

# Titles shown on each panel
class_titles <- c(
  "AminoSugar"     = "AminoSugar",
  "Carbohydrate"   = "Carbohydrate",
  "CondensedHC"    = "Condensed HC",
  "Lignin"         = "Lignin",
  "Lipid"          = "Lipid",
  "Other"          = "Other",
  "Peptide"        = "Peptide",
  "Tannin"         = "Tannin",
  "Unsaturated_HC" = "Unsaturated HC"
)

theme_shannon <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 10),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text   = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

make_shannon_plot <- function(df, facet_var, title) {
  ggplot(df, aes(x = Pickup_t, y = Shannon, fill = Pickup_t)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = Pickup_t), width = 0.2, size = 2, alpha = 0.9) +
    labs(title = title, x = "", y = "Shannon Diversity") +
    scale_x_discrete(labels = tp_labels, drop = FALSE) +
    scale_fill_manual(
      values = time_colors,
      breaks = tp_levels,
      drop = FALSE,
      name = "Collection time\n(years)"
    ) +
    scale_color_manual(
      values = time_colors,
      breaks = tp_levels,
      drop = FALSE,
      name = "Collection time\n(years)"
    )+
    facet_wrap(vars({{ facet_var }})) +
    theme_shannon
}

# Duplicate "None" rows into each target level (Temp or CO2)
duplicate_none <- function(df, var, levels_vec) {
  var_chr <- rlang::as_name(rlang::ensym(var))
  df %>%
    mutate("{var_chr}" := as.character(.data[[var_chr]])) %>%
    filter(.data[[var_chr]] %in% c("None", levels_vec)) %>%
    tidyr::uncount(
      weights = if_else(.data[[var_chr]] == "None", length(levels_vec), 1),
      .id = "rep_id"
    ) %>%
    mutate(
      "{var_chr}" := if_else(.data[[var_chr]] == "None", levels_vec[rep_id], .data[[var_chr]]),
      "{var_chr}" := factor(.data[[var_chr]], levels = levels_vec)
    ) %>%
    select(-rep_id)
}

# Build 3x3 patchwork panel set (per litter, per effect)
make_panel_grid <- function(df_litter_effect, facet_var) {
  
  # enforce class ordering (layout order)
  df_litter_effect <- df_litter_effect %>%
    mutate(Class = factor(Class, levels = class_levels))
  
  p_list <- vector("list", length(class_levels))
  names(p_list) <- class_levels
  
  for (cl in class_levels) {
    df_cl <- df_litter_effect %>% filter(Class == cl)
    
    # Keep layout stable even if a class is missing
    if (nrow(df_cl) == 0) {
      p_list[[cl]] <- ggplot() + theme_void() + labs(title = class_titles[[cl]])
    } else {
      p_list[[cl]] <- make_shannon_plot(df_cl, {{ facet_var }}, title = class_titles[[cl]])
    }
  }
  
  (p_list[[1]] + p_list[[2]] + p_list[[3]] +
      p_list[[4]] + p_list[[5]] + p_list[[6]] +
      p_list[[7]] + p_list[[8]] + p_list[[9]]) +
    plot_layout(nrow = 3, byrow = TRUE)
}

# Inputs: All.csv pipeline -----------------------
all <- read.csv("Input/Shannon/All.csv", check.names = FALSE)

classes <- all %>%
  column_to_rownames("Mass") %>%
  t() %>%
  as.data.frame() %>%
  mutate(Class = recode(Class, "Protein" = "Peptide"))

df <- classes %>%
  rownames_to_column("Mass") %>%
  select(-Class)

class_meta <- classes %>%
  rownames_to_column("Mass") %>%
  select(Mass, Class)

rm(all, classes)

metadata <- read.csv("Input/metadata.csv") %>%
  select(SampleID, Litter, Pickup_t, Temp, CO2) %>%
  mutate(
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
    Pickup_t = factor(Pickup_t, levels = tp_levels, ordered = TRUE)
  )

# Build intensity matrix (sample x Mass)
intensity_matrix <- df %>%
  pivot_longer(-Mass, names_to = "SampleID", values_to = "intensity") %>%
  mutate(intensity = suppressWarnings(as.numeric(intensity))) %>%
  pivot_wider(names_from = "SampleID", values_from = "intensity", values_fill = 0) %>%
  select(Mass, all_of(metadata$SampleID)) %>%
  column_to_rownames("Mass") %>%
  t()

# Total-sum normalize
norm_intensity_matrix <- decostand(intensity_matrix, method = "total")

# Long format
df_long <- as.data.frame(t(norm_intensity_matrix)) %>%
  rownames_to_column("Mass") %>%
  pivot_longer(-Mass, names_to = "SampleID", values_to = "intensity") %>%
  left_join(class_meta, by = "Mass") %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    t_num = parse_time(as.character(Pickup_t))
  )

# Shannon per sample × class
div_by_sample <- df_long %>%
  group_by(SampleID, Litter, Pickup_t, t_num, Temp, CO2, Class) %>%
  summarise(Shannon = diversity(intensity, index = "shannon"), .groups = "drop") %>%
mutate(
  Class = recode(Class,
                 "Amino sugar"            = "AminoSugar",
                 "Condensed hydrocarbon"  = "CondensedHC",
                 "Unsaturated hydrocarbon"= "Unsaturated_HC"
  )
)


# Output directory ----------------------------
if (!dir.exists("Plots")) dir.create("Plots", recursive = TRUE)

# Generate plots for ALL litters -----------------
litters <- div_by_sample %>%
  filter(!is.na(Litter)) %>%
  distinct(Litter) %>%
  pull(Litter) %>%
  as.character()

for (lit in litters) {
  
  #  TEMP effect 
  div_temp <- div_by_sample %>%
    filter(Litter == lit) %>%
    duplicate_none(Temp, temp_levels)
  
  temp_grid <- make_panel_grid(div_temp, Temp) +
    plot_layout(guides = "collect") +   
    plot_annotation(
      title = paste0(lit, " — Temperature effect"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
  
  ggsave(
    filename = paste0("Plots/Shannon/SupFig_", lit, "_shannon_temp_fromAll.png"),
    plot = temp_grid,
    width = 12, height = 8, dpi = 300
  )
  
  # CO2 effect 
  div_co2 <- div_by_sample %>%
    filter(Litter == lit) %>%
    duplicate_none(CO2, co2_levels)
  
  co2_grid <- make_panel_grid(div_co2, CO2) +
    plot_layout(guides = "collect") +  
    plot_annotation(
      title = paste0(lit, " — CO₂ effect"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
  
  ggsave(
    filename = paste0("Plots/Shannon/SupFig_", lit, "_shannon_CO2_fromAll.png"),
    plot = co2_grid,
    width = 10, height = 6.5, dpi = 300
  )
}

message("Done. Saved figures to: Plots/")
