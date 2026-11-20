### SPRUCE litter decomposition 2015–2017
## Supplementary Figures 8–10 — Transformations

library(tidyverse)
library(patchwork)
library(ragg)

## Input data
transformations_summary_df <- read.csv(
  "../Input/Metabodirect_All_Sum_Time/6_transformations/Transformations_summary_counts.csv"
)
metadata <- read.csv("../Input/metadata.csv")

litter_order <- c("ANG", "MAG", "SPL", "SPR", "LTL", "LTR")
custom_labels_facet <- c(
  ANG = "italic('S. angustifolium')",
  MAG = "italic('S. magellanicum')",
  SPL = "'Spruce needles'",
  SPR = "'Spruce roots'",
  LTL = "'Labrador tea leaves'",
  LTR = "'Labrador tea roots'"
)
time_breaks <- c(0, 0.5, 1, 2)

metadata <- metadata %>%
  mutate(
    Time = recode(
      Pickup_t,
      "T_0" = 0,
      "T_0.5" = 0.5,
      "T_1" = 1,
      "T_2" = 2
    )
  )

trans_df <- transformations_summary_df %>%
  left_join(
    metadata %>% select(SampleID, Litter, Time),
    by = "SampleID"
  )

group_per_sample <- trans_df %>%
  group_by(SampleID, Litter, Time, Group) %>%
  summarise(Perc = sum(Perc_Counts), .groups = "drop")

## Figure-wide typography
figure_theme <- theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
    #plot.margin = margin(4, 5, 4, 5)
  )

## Shared axis and colors for the stacked bar charts
group_levels <- sort(unique(na.omit(group_per_sample$Group)))
group_fill <- scale_fill_brewer(
  palette = "Set2",
  limits = group_levels,
  name = "Transformation group"
)
abundance_axis <- scale_y_continuous(
  breaks = c(0, 0.5, 1),
  expand = expansion(mult = 0)
)
abundance_view <- coord_cartesian(ylim = c(0, 1))

## Panel A — Groups over time
group_time <- group_per_sample %>%
  group_by(Time, Group) %>%
  summarise(MeanPerc = mean(Perc), .groups = "drop") %>%
  mutate(Time = factor(Time, levels = time_breaks))

pA_overall <- ggplot(group_time, aes(Time, MeanPerc, fill = Group)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.9) +
  abundance_axis +
  abundance_view +
  group_fill +
  labs(
    title = "A) Transformation groups over time",
    x = "Time (years)",
    y = "Relative abundance"
  ) +
  figure_theme +
  theme(
    legend.position = "right",
    legend.key.size = unit(4.5, "mm")
  )

## Panel B — Groups by litter and time
group_litter_time <- group_per_sample %>%
  group_by(Litter, Time, Group) %>%
  summarise(MeanPerc = mean(Perc), .groups = "drop") %>%
  mutate(
    Litter = factor(Litter, levels = litter_order),
    Time = factor(Time, levels = time_breaks)
  )

pB_litter <- ggplot(group_litter_time, aes(Time, MeanPerc, fill = Group)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.9) +
  facet_wrap(
    ~ Litter,
    ncol = 3,
    labeller = as_labeller(custom_labels_facet, label_parsed)
  ) +
  abundance_axis +
  abundance_view +
  group_fill +
  labs(
    title = "B) Transformation groups by litter",
    x = "Time (years)",
    y = "Relative abundance"
  ) +
  figure_theme +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.major = element_blank(),
    panel.spacing = unit(2, "mm"),
    legend.position = "none"
  )

## Transformation-level means for C and D
major_groups <- trans_df %>%
  filter(Group %in% c("Abiotic", "Amino Acid")) %>%
  group_by(SampleID, Litter, Time, Group, Transformation) %>%
  summarise(Perc = sum(Perc_Counts), .groups = "drop")

major_trans <- major_groups %>%
  group_by(Group, Transformation, Time) %>%
  summarise(MeanPerc = mean(Perc), .groups = "drop")

abiotic_df <- major_trans %>% filter(Group == "Abiotic")
amino_df <- major_trans %>% filter(Group == "Amino Acid")

## Supplementary Figures 9 and 10 — Heatmaps
make_heatmap <- function(data, title, palette) {
  transformation_order <- data %>%
    group_by(Transformation) %>%
    summarise(total = mean(MeanPerc), .groups = "drop") %>%
    arrange(total) %>%
    pull(Transformation)

  data <- data %>%
    mutate(
      Time = factor(Time, levels = time_breaks),
      Transformation = factor(Transformation, levels = transformation_order)
    )

  ggplot(data, aes(Time, Transformation, fill = MeanPerc)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    scale_fill_viridis_c(
      option = palette,
      name = "Relative abundance",
      n.breaks = 3
    ) +
    labs(title = title, x = "Time (years)", y = NULL) +
    figure_theme +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title.position = "top",
      legend.key.width = unit(5, "mm"),
      legend.key.height = unit(10, "mm")
    )
}

pC_abiotic <- make_heatmap(abiotic_df, "Abiotic transformations", "C")
pD_amino <- make_heatmap(amino_df, "Amino acid-related transformations", "D")

## Export Supplementary Figures 8–10
# Each heatmap gets its own page to accommodate all transformation labels.
Supplementary_Figure_8 <- wrap_plots(
  wrap_elements(full = pA_overall),
  wrap_elements(full = pB_litter),
  ncol = 1,
  heights = c(1, 1.65)
)

output_dir <- "../Plots_frontier/tiff"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_supplementary_figure <- function(plot, filename, height_mm) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = 180,
    height = height_mm,
    units = "mm",
    dpi = 300,
    device = agg_tiff,
    bg = "white",
    compression = "lzw"
  )
}

save_supplementary_figure(Supplementary_Figure_8, "SupFig8.tiff", height_mm = 150)
save_supplementary_figure(pC_abiotic, "SupFig9.tiff", height_mm = 190)
save_supplementary_figure(pD_amino, "SupFig10.tiff", height_mm = 200)

