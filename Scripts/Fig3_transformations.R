### SPRUCE Litter decomposition 2015-2017 
## Fig 3 - Transformations

library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridisLite)

## input data: 
transformations_summary_df <- read.csv("Input/Metabodirect_All_Sum_Time/6_transformations/Transformations_summary_counts.csv")
metadata <- read.csv("Input/metadata.csv")

# Make sure Time is numeric
metadata <- metadata %>%
  mutate(Time = recode(Pickup_t,
                       "T_0"   = 0,
                       "T_0.5" = 0.5,
                       "T_1"   = 1,
                       "T_2"   = 2))

# Merge transformation df with metadata
trans_df <- transformations_summary_df %>%
  left_join(metadata %>% select(SampleID, Litter, Time),
            by = "SampleID")


## Summarize Groups per sample (using Perc_Counts)
# Collapse to Group level within each sample
group_per_sample <- trans_df %>%
  group_by(SampleID, Litter, Time, Group) %>%
  summarise(Perc = sum(Perc_Counts),
            .groups = "drop")


## Panel A: Groups × Time (all litters combined) -----
group_time <- group_per_sample %>%
  group_by(Time, Group) %>%
  summarise(MeanPerc = mean(Perc),
            SDPerc   = sd(Perc),
            .groups  = "drop")

pA_overall <- ggplot(group_time,
                     aes(x = factor(Time),
                         y = MeanPerc,
                         fill = Group)) +
  geom_col(color = "black", linewidth = 0.2) +
  labs(
    title = "A) Transformation groups over time",
    x     = "Time (years)",
    y     = "Tranformation Relative Abundance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set2")
pA_overall

## Panel B: Groups × Litter × Time ------
# Manually specify litter order
litter_order <- c("ANG", "MAG", "SPL", "SPR", "LTL", "LTR")

group_litter_time <- group_per_sample %>%
  group_by(Litter, Time, Group) %>%
  summarise(MeanPerc = mean(Perc),
            SDPerc   = sd(Perc),
            .groups  = "drop") %>%
  mutate(
    Litter = factor(Litter, levels = litter_order)
  )

pB_litter <- ggplot(group_litter_time,
                    aes(x = factor(Time),
                        y = MeanPerc,
                        fill = Group)) +
  geom_col(color = "black", linewidth = 0.2) +
  facet_wrap(~ Litter, ncol = 3) +
  labs(
    title = "B) Transformation groups by litter type",
    x     = "Time (years)",
    y     = "Tranformation Relative Abundance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0, face = "bold"),
    strip.text      = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2")
pB_litter

## Major groups: Abiotic & Sugar -----
major_groups <- trans_df %>%
  filter(Group %in% c("Abiotic","Amino Acid")) %>%
  group_by(SampleID, Litter, Time, Group, Transformation) %>%
  summarise(Perc = sum(Perc_Counts),
            .groups = "drop")

major_trans <- major_groups %>%
  group_by(Group, Transformation, Time) %>%
  summarise(
    MeanPerc = mean(Perc),
    SDPerc   = sd(Perc),
    .groups  = "drop"
  ) %>%
  group_by(Group, Time) %>%
  arrange(desc(MeanPerc), .by_group = TRUE) %>%
  mutate(Transformation_sorted = factor(Transformation,
                                        levels = unique(Transformation))) %>%
  ungroup()

abiotic_df <- major_trans %>% filter(Group == "Abiotic")
sugar_df   <- major_trans %>% filter(Group == "Sugar")
amino_df <- major_trans %>% filter(Group == "Amino Acid")

## Panel C: Abiotic transformations over time ---------
pC_abiotic <- ggplot(abiotic_df,
                     aes(x = factor(Time),
                         y = MeanPerc,
                         fill = Transformation_sorted)) +
  geom_col(color = "black", linewidth = 0.2) +
  labs(
    title = "C) Abiotic transformations",
    x     = "Time (years)",
    y     = "Tranformation Relative Abundance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_viridis_d(option = "C",
                       name = "Transformation")+
  guides(fill = guide_legend(ncol = 4))
pC_abiotic

## Panel D: Sugar-related transformations over time
pD_amino <- ggplot(amino_df,
                   aes(x = factor(Time),
                       y = MeanPerc,
                       fill = Transformation_sorted)) +
  geom_col(color = "black", linewidth = 0.2) +
  labs(
    title = "D) Amino Acid-related transformations",
    x     = "Time (years)",
    y     = "Tranformation Relative Abundance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_viridis_d(option = "D",
                       name = "Transformation")+
  guides(fill = guide_legend(ncol = 3))
pD_amino

## Fig 3 -------

pC_abiotic_for_legend <- pC_abiotic +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.3, "lines"),
    legend.margin     = margin(t = 1, r = 1, b = 1, l = 1)
  ) +
  guides(fill = guide_legend(ncol = 3))

# Extract the legend as a grob
legend_C <- cowplot::get_legend(pC_abiotic_for_legend)

# Create a version of the plot WITHOUT legend
pC_abiotic_noleg <- pC_abiotic +
  theme(legend.position = "none")
legend_C_wrap <- patchwork::wrap_elements(legend_C)

pC_amino_for_legend <- pD_amino +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.3, "lines"),
    legend.margin     = margin(t = 1, r = 1, b = 1, l = 1)
  ) +
  guides(fill = guide_legend(ncol = 3))
# 
# Extract the legend as a grob
legend_D <- cowplot::get_legend(pC_amino_for_legend)

# Create a version of the plot WITHOUT legend
pD_amino_noleg <- pD_amino +
  theme(legend.position = "none")

legend_D_wrap <- patchwork::wrap_elements(legend_D)




top <- (pA_overall + plot_spacer() + pB_litter) +
  plot_layout(widths = c(1, 0.3, 1))

bottom <- (pD_amino_noleg + plot_spacer() + legend_D_wrap) +
  plot_layout(widths = c(0.5, 0.1, 0.4))   

middle <- (pC_abiotic_noleg + plot_spacer() + legend_C_wrap) +
  plot_layout(widths = c(0.5, 0.1, 0.4))   

full <- (top / middle / bottom) +
  plot_layout(heights = c(1, 0.8, 0.8))  

full

ggsave("Plots/Fig3_transformations.png",
       full,
       width = 350, height = 250, units = "mm", dpi = 600)
