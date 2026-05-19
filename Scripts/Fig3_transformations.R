### SPRUCE Litter decomposition 2015-2017
## Fig 3 - Transformations
## Revised layout for clearer GCA-style presentation

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ragg)
library(scales)

## input data -------------------------------------------------------------
transformations_summary_df <- read.csv(
  "Input/Metabodirect_All_Sum_Time/6_transformations/Transformations_summary_counts.csv"
)
metadata <- read.csv("Input/metadata.csv")

## custom labels for facets ----------------------------------------------
custom_labels_facet <- c(
  ANG = "italic('S. angustifolium')",
  MAG = "italic('S. magellanicum')",
  SPL = "'Spruce needles'",
  SPR = "'Spruce roots'",
  LTL = "'Labrador tea leaves'",
  LTR = "'Labrador tea roots'"
)

## time variable ----------------------------------------------------------
metadata <- metadata %>%
  mutate(
    Time = recode(
      Pickup_t,
      "T_0"   = 0,
      "T_0.5" = 0.5,
      "T_1"   = 1,
      "T_2"   = 2
    )
  )

time_breaks <- c(0, 0.5, 1, 2)

## merge ------------------------------------------------------------------
trans_df <- transformations_summary_df %>%
  left_join(
    metadata %>% select(SampleID, Litter, Time),
    by = "SampleID"
  )

## summarize groups per sample -------------------------------------------
group_per_sample <- trans_df %>%
  group_by(SampleID, Litter, Time, Group) %>%
  summarise(
    Perc = sum(Perc_Counts),
    .groups = "drop"
  )

## Panel A: Groups x Time -------------------------------------------------
group_time <- group_per_sample %>%
  group_by(Time, Group) %>%
  summarise(
    MeanPerc = mean(Perc),
    .groups  = "drop"
  ) %>%
  mutate(Time = factor(Time, levels = time_breaks))

pA_overall <- ggplot(
  group_time,
  aes(x = Time, y = MeanPerc, fill = Group)
) +
  geom_col(color = "black", linewidth = 0.2, width = 0.9) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "A) Transformation groups over time",
    x     = "Time (years)",
    y     = "Transformation relative abundance"
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(hjust = 0, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 8, margin = margin(b = 2)),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(4, "mm")
  ) +
  scale_fill_brewer(palette = "Set2")

## Panel B: Groups x Litter x Time ---------------------------------------
litter_order <- c("ANG", "MAG", "SPL", "SPR", "LTL", "LTR")

group_litter_time <- group_per_sample %>%
  group_by(Litter, Time, Group) %>%
  summarise(
    MeanPerc = mean(Perc),
    .groups  = "drop"
  ) %>%
  mutate(
    Litter = factor(Litter, levels = litter_order),
    Time   = factor(Time, levels = time_breaks)
  )

pB_litter <- ggplot(
  group_litter_time,
  aes(x = Time, y = MeanPerc, fill = Group)
) +
  geom_col(color = "black", linewidth = 0.2, width = 0.9) +
  facet_wrap(
    ~ Litter,
    ncol = 3,
    labeller = as_labeller(custom_labels_facet, label_parsed)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "B) Transformation groups by litter",
    x     = "Time (years)",
    y     = "Transformation relative abundance"
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(hjust = 0, face = "bold"),
    strip.text       = element_text(size = 7),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "none", 
    plot.margin    = margin(t = 5, r = 0, b = 5, l = 5)) +
  scale_fill_brewer(palette = "Set2")

## Major groups: Abiotic and Amino Acid ----------------------------------
major_groups <- trans_df %>%
  filter(Group %in% c("Abiotic", "Amino Acid")) %>%
  group_by(SampleID, Litter, Time, Group, Transformation) %>%
  summarise(
    Perc = sum(Perc_Counts),
    .groups = "drop"
  )

major_trans <- major_groups %>%
  group_by(Group, Transformation, Time) %>%
  summarise(
    MeanPerc = mean(Perc),
    .groups  = "drop"
  ) %>%
  ungroup()

abiotic_df <- major_trans %>%
  filter(Group == "Abiotic")

amino_df <- major_trans %>%
  filter(Group == "Amino Acid")

## reorder transformations by overall abundance --------------------------
abiotic_order <- abiotic_df %>%
  group_by(Transformation) %>%
  summarise(total = mean(MeanPerc), .groups = "drop") %>%
  arrange(total) %>%
  pull(Transformation)

amino_order <- amino_df %>%
  group_by(Transformation) %>%
  summarise(total = mean(MeanPerc), .groups = "drop") %>%
  arrange(total) %>%
  pull(Transformation)

abiotic_df <- abiotic_df %>%
  mutate(
    Time = factor(Time, levels = time_breaks),
    Transformation = factor(Transformation, levels = abiotic_order)
  )

amino_df <- amino_df %>%
  mutate(
    Time = factor(Time, levels = time_breaks),
    Transformation = factor(Transformation, levels = amino_order)
  )

## common heatmap theme ---------------------------------------------------
heatmap_theme <- theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(hjust = 1, face = "bold"),
    panel.grid      = element_blank(),
    axis.text.y     = element_text(size = 7),
    axis.text.x     = element_text(size = 8),
    legend.position = "bottom",
    legend.title    = element_text(size = 7, margin = margin(b = 2)),
    legend.text     = element_text(size = 7),
    legend.key.width = unit(12, "mm"),
    legend.key.height = unit(3, "mm"),
    plot.margin     = margin(t = 5, r = 2, b = 2, l = 5)
  )

## Panel C: Abiotic heatmap ----------------------------------------------
pC_abiotic <- ggplot(
  abiotic_df,
  aes(x = Time, y = Transformation, fill = MeanPerc)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "C",
    name   = "Relative abundance"
  ) +
  labs(
    title = "C) Abiotic transformations",
    x     = "Time (years)",
    y     = NULL
  ) +
  heatmap_theme +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = unit(30, "mm"),
      barheight = unit(3, "mm")
    )
  )

## Panel D: Amino acid-related heatmap -----------------------------------
pD_amino <- ggplot(
  amino_df,
  aes(x = Time, y = Transformation, fill = MeanPerc)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Relative abundance"
  ) +
  labs(
    title = "D) Amino acid-related transformations",
    x     = "Time (years)",
    y     = NULL
  ) +
  heatmap_theme +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = unit(30, "mm"),
      barheight = unit(3, "mm")
    )
  )

## layout -----------------------------------------------------------------
library(cowplot) 
## Top row with A and B (no patchwork)
empty_plot <- ggplot() + theme_void()
top_row <- cowplot::plot_grid(
  empty_plot, pA_overall, pB_litter, 
  nrow       = 1,
  rel_widths = c(0.2, 1, 1.5),  # last slot has almost no width
  align      = "v",
  axis       = "lr"
)
bottom_row <- cowplot::plot_grid(
  pC_abiotic, pD_amino,
  nrow       = 1,
  rel_widths = c(1, 1),
  align      = "v",
  axis       = "lr"
)

fig3_cow <- cowplot::plot_grid(
  top_row,
  bottom_row,
  ncol        = 1,
  rel_heights = c(1, 1.5),
  align       = "v"
)
fig3_cow

ggsave(
  "Plots/Updated/Fig3_transformations_cowplot.tiff",
  fig3_cow,
  width       = 190,
  height      = 200,
  units       = "mm",
  dpi         = 300,
  device      = agg_tiff,bg = "white",
  compression = "lzw"
)

