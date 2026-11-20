### SPRUCE Litter decomposition 2015–2017
### Supplementary figures — Shannon diversity by litter and treatment

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(patchwork)
  library(ragg)
})

time_colors <- c(
  "T_0" = "#4DBBD5",
  "T_0.5" = "#3C5488",
  "T_1" = "#B09C85",
  "T_2" = "#503A2C"
)
tp_levels <- names(time_colors)
tp_labels <- setNames(c("0", "0.5", "1", "2"), tp_levels)
temp_levels <- c("+0°C", "+4.5°C", "+9°C")
co2_levels <- c("Ambient", "Elevated")

class_titles <- c(
  AminoSugar = "Amino sugar",
  Carbohydrate = "Carbohydrate",
  CondensedHC = "Condensed HC",
  Lignin = "Lignin",
  Lipid = "Lipid",
  Other = "Other",
  Peptide = "Peptide",
  Tannin = "Tannin",
  Unsaturated_HC = "Unsaturated HC"
)

theme_shannon <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(margin = margin(t = 2)),
    axis.ticks.length = unit(1, "mm"),
    panel.spacing.x = unit(1, "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12, margin = margin(2, 1, 2, 1)),
    legend.position = "none",
    plot.margin = margin(3, 3, 3, 3)
  )

make_shannon_plot <- function(data, facet_var, title, y_title) {
  ggplot(data, aes(Pickup_t, Shannon, fill = Pickup_t)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(
      aes(colour = Pickup_t),
      position = position_jitter(width = 0.16, height = 0, seed = 1),
      size = 1.5,
      alpha = 0.9
    ) +
    labs(title = title, x = NULL, y = y_title) +
    scale_x_discrete(labels = tp_labels, drop = FALSE) +
    scale_y_continuous(n.breaks = 4) +
    scale_fill_manual(values = time_colors, guide = "none") +
    scale_colour_manual(values = time_colors, guide = "none") +
    facet_wrap(vars({{ facet_var }}), nrow = 1, drop = FALSE) +
    theme_shannon
}

# Include untreated baseline rows in every treatment facet.
duplicate_none <- function(data, var, levels_vec) {
  var_name <- rlang::as_name(rlang::ensym(var))

  data %>%
    mutate("{var_name}" := as.character(.data[[var_name]])) %>%
    filter(.data[[var_name]] %in% c("None", levels_vec)) %>%
    tidyr::uncount(
      weights = if_else(.data[[var_name]] == "None", length(levels_vec), 1L),
      .id = "rep_id"
    ) %>%
    mutate(
      "{var_name}" := if_else(
        .data[[var_name]] == "None", levels_vec[rep_id], .data[[var_name]]
      ),
      "{var_name}" := factor(.data[[var_name]], levels = levels_vec)
    ) %>%
    select(-rep_id)
}

make_panel_grid <- function(data, facet_var, ncol) {
  panels <- lapply(seq_along(class_titles), function(i) {
    class_name <- names(class_titles)[i]
    class_data <- filter(data, Class == class_name)
    title <- class_titles[[i]]
    y_title <- if ((i - 1L) %% ncol == 0L) "Shannon diversity" else NULL

    if (nrow(class_data) == 0L) {
      return(ggplot() + theme_void() + labs(title = title) + theme_shannon)
    }

    make_shannon_plot(class_data, {{ facet_var }}, title, y_title)
  })

  wrap_plots(panels, ncol = ncol, byrow = TRUE)
}

# Read mass intensities and sample metadata.
all_intensities <- read.csv("Input/Shannon/All.csv", check.names = FALSE)

classes <- all_intensities %>%
  column_to_rownames("Mass") %>%
  t() %>%
  as.data.frame() %>%
  mutate(Class = recode(Class, "Protein" = "Peptide"))

intensities <- classes %>%
  rownames_to_column("Mass") %>%
  select(-Class)

class_meta <- classes %>%
  rownames_to_column("Mass") %>%
  select(Mass, Class)

rm(all_intensities, classes)

metadata <- read.csv("Input/metadata.csv") %>%
  select(SampleID, Litter, Pickup_t, Temp, CO2) %>%
  mutate(
    Temp = case_when(
      Temp == 0 ~ "+0°C",
      Temp == 4.5 ~ "+4.5°C",
      Temp == 9 ~ "+9°C",
      TRUE ~ as.character(Temp)
    ),
    CO2 = case_when(
      CO2 == "aCO2" ~ "Ambient",
      CO2 == "eCO2" ~ "Elevated",
      TRUE ~ as.character(CO2)
    ),
    Pickup_t = factor(Pickup_t, levels = tp_levels, ordered = TRUE)
  )

intensity_matrix <- intensities %>%
  pivot_longer(-Mass, names_to = "SampleID", values_to = "intensity") %>%
  mutate(intensity = suppressWarnings(as.numeric(intensity))) %>%
  pivot_wider(names_from = "SampleID", values_from = "intensity", values_fill = 0) %>%
  select(Mass, all_of(metadata$SampleID)) %>%
  column_to_rownames("Mass") %>%
  t()

norm_intensity_matrix <- decostand(intensity_matrix, method = "total")

df_long <- as.data.frame(t(norm_intensity_matrix)) %>%
  rownames_to_column("Mass") %>%
  pivot_longer(-Mass, names_to = "SampleID", values_to = "intensity") %>%
  left_join(class_meta, by = "Mass") %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(t_num = as.numeric(sub("^T_", "", as.character(Pickup_t))))

div_by_sample <- df_long %>%
  group_by(SampleID, Litter, Pickup_t, t_num, Temp, CO2, Class) %>%
  summarise(Shannon = diversity(intensity, index = "shannon"), .groups = "drop") %>%
  mutate(
    Class = recode(
      Class,
      "Amino sugar" = "AminoSugar",
      "Condensed hydrocarbon" = "CondensedHC",
      "Unsaturated hydrocarbon" = "Unsaturated_HC"
    )
  )

output_dir <- "Plots_frontier/Shannon"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

litters <- div_by_sample %>%
  filter(!is.na(Litter)) %>%
  distinct(Litter) %>%
  pull(Litter) %>%
  as.character()

annotation_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
  plot.caption = element_text(hjust = 0.5, size = 12)
)

save_shannon_figure <- function(plot, litter, effect, height_mm) {
  ggsave(
    filename = file.path(
      output_dir, paste0("SupFig_", litter, "_shannon_", effect, "_fromAll.tiff")
    ),
    plot = plot,
    width = 180,
    height = height_mm,
    units = "mm",
    dpi = 300,
    device = ragg::agg_tiff,
    compression = "lzw",
    bg = "white"
  )
}

for (lit in litters) {
  div_litter <- filter(div_by_sample, Litter == lit)
  div_temp <- duplicate_none(div_litter, Temp, temp_levels)
  div_co2 <- duplicate_none(div_litter, CO2, co2_levels)

  temp_grid <- make_panel_grid(div_temp, Temp, ncol = 2) +
    plot_annotation(
      title = paste0(lit, " - Temperature effect"),
      caption = "Collection time (years)",
      theme = annotation_theme
    )

  co2_grid <- make_panel_grid(div_co2, CO2, ncol = 3) +
    plot_annotation(
      title = bquote(.(lit) ~ "- " * CO[2] ~ "effect"),
      caption = "Collection time (years)",
      theme = annotation_theme
    )

  save_shannon_figure(temp_grid, lit, "temp", height_mm = 230)
  save_shannon_figure(co2_grid, lit, "CO2", height_mm = 190)
}

message("Done. Saved figures to: ", output_dir)
