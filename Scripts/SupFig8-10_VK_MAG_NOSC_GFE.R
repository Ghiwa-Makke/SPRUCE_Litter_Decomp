### SPRUCE Litter decomposition 2015–2017
### Van krevelen Pairwise Consective timepoints – MAG
### SupFig8 - SupFig10

## Libraries: 
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(grid) 
library(ggrepel)
library(patchwork)
library(cowplot)
library(rstatix)
library(ggpubr)

tp_colors <- c(
  "T_0"   = "#4DBBD5", 
  "T_0.5" = "#3C5488",
  "T_1"   = "#B09C85", 
  "T_2"   = "#503A2C"
)

## Data: 
intensity_df <- read_csv('Input/Report_processed_noNorm_MAG.csv') # metabodirect output
metadata <- read_csv('Input/metadata.csv')%>%
  filter(Litter == "MAG")

# Data long format
df_longer <- intensity_df %>%
  pivot_longer(cols = metadata$SampleID, names_to = 'SampleID', values_to = 'Intensity') %>%
  left_join(metadata, by = 'SampleID')


# VK Pairwise Time --------------------------------------------------------

### 0-0.5 -------------------------------------------------------------------
sub_df_0_0.5 <- df_longer %>%
  filter(Pickup_t %in% c("T_0", "T_0.5")) %>%
  mutate(is_present = Intensity > 0)

# Masses present at each timepoint
masses_0 <- sub_df_0_0.5 %>%
  filter(Pickup_t == "T_0", is_present) %>%
  pull(Mass) %>% unique()

masses_0.5 <- sub_df_0_0.5 %>%
  filter(Pickup_t == "T_0.5", is_present) %>%
  pull(Mass) %>% unique()

# Unique sets only (no shared)
unique_0 <- setdiff(masses_0, masses_0.5)
unique_0.5 <- setdiff(masses_0.5, masses_0)

plot_df_0_0.5 <- sub_df_0_0.5 %>%
  filter(
    (Pickup_t == "T_0" & Mass %in% unique_0) |
      (Pickup_t == "T_0.5" & Mass %in% unique_0.5)
  ) %>%
  group_by(Mass, HC, OC, Pickup_t) %>%
  summarise(Intensity = mean(Intensity, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5")))

# VK plot
p_vk_0_0.5 <- ggplot(plot_df_0_0.5, aes(x = OC, y = HC)) +
  geom_point(aes(color = Pickup_t), size = 1.8, alpha = 0.9) +
  scale_color_manual(
    values = tp_colors,
    drop   = FALSE,
    name   = "Collection time (years)",
    labels = c(
      "T_0"   = "0 years",
      "T_0.5" = "0.5 years"
    )
  ) +
  # Vertical arrow
  annotate(
    "segment",
    x = 1.3, xend = 1.3,
    y = 2,  yend = 0.7,
    color = "grey40",
    linewidth = 2.8,
    arrow = arrow(length = unit(0.4, "cm"), type = "open")
  ) +
  # Text at top: "0 years"
  annotate("text",
           x = 1.3, y = 2.5,
           label = "0 years",
           fontface = "bold",
           size = 4,
           color = "grey40") +
  # Text at bottom: "0.5 years"
  annotate("text",
           x = 1.3, y = 0.3,
           label = "0.5 years",
           fontface = "bold",
           color = "grey40",
           size = 4) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0,1.4))+
  theme_bw(base_size = 14) +
  labs(
    title = "Unique compounds: 0 vs 0.5 years",
    x     = "O/C",
    y     = "H/C"
  ) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

p_vk_0_0.5


### 0.5-1 -------------------------------------------------------------------
sub_df_0.5_1 <- df_longer %>%
  filter(Pickup_t %in% c("T_0.5", "T_1")) %>%
  mutate(is_present = Intensity > 0)

# Masses present at each timepoint
masses_0.5 <- sub_df_0.5_1 %>%
  filter(Pickup_t == "T_0.5", is_present) %>%
  pull(Mass) %>% unique()

masses_1 <- sub_df_0.5_1 %>%
  filter(Pickup_t == "T_1", is_present) %>%
  pull(Mass) %>% unique()

# Unique sets only (no shared)
unique_0.5 <- setdiff(masses_0.5, masses_1)
unique_1 <- setdiff(masses_1, masses_0.5)

plot_df_0.5_1 <- sub_df_0.5_1 %>%
  filter(
    (Pickup_t == "T_0.5" & Mass %in% unique_0.5) |
      (Pickup_t == "T_1" & Mass %in% unique_1)
  ) %>%
  group_by(Mass, HC, OC, Pickup_t) %>%
  summarise(Intensity = mean(Intensity, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Pickup_t = factor(Pickup_t, levels = c("T_0.5", "T_1")))

# VK plot
p_vk_0.5_1 <- ggplot(plot_df_0.5_1, aes(x = OC, y = HC)) +
  geom_point(aes(color = Pickup_t), size = 1.8, alpha = 0.9) +
  scale_color_manual(
    values = tp_colors,
    drop   = FALSE,
    name   = "Collection time (years)",
    labels = c(
      "T_0.5"   = "0.5 years",
      "T_1" = "1 years"
    )
  ) +
  annotate(
    "segment",
    x = 0.5, xend = 1.1,
    y = 0.1,  yend = 0.1,
    color = "grey40",
    linewidth = 2.8,
    arrow = arrow(length = unit(0.4, "cm"), type = "open")
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           label = "0.5 years",
           fontface = "bold",
           size = 4,
           color = "grey40") +
  annotate("text",
           x = 1.3, y = 0.1,
           label = "1 year",
           fontface = "bold",
           color = "grey40",
           size = 4) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0,1.4))+
  theme_bw(base_size = 14) +
  labs(
    title = "Unique compounds: 0.5 vs 1 years",
    x     = "O/C",
    y     = "H/C"
  ) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom", 
    legend.title = element_blank()
  )

p_vk_0.5_1


### 1-2 ---------------------------------------------------------------------
sub_df_1_2 <- df_longer %>%
  filter(Pickup_t %in% c("T_1", "T_2")) %>%
  mutate(is_present = Intensity > 0)

# Masses present at each timepoint
masses_1 <- sub_df_1_2 %>%
  filter(Pickup_t == "T_1", is_present) %>%
  pull(Mass) %>% unique()

masses_2 <- sub_df_1_2 %>%
  filter(Pickup_t == "T_2", is_present) %>%
  pull(Mass) %>% unique()

# Unique sets only (no shared)
unique_1 <- setdiff(masses_1, masses_2)
unique_2 <- setdiff(masses_2, masses_1)

plot_df_1_2 <- sub_df_1_2 %>%
  filter(
    (Pickup_t == "T_1" & Mass %in% unique_1) |
      (Pickup_t == "T_2" & Mass %in% unique_2)
  ) %>%
  group_by(Mass, HC, OC, Pickup_t) %>%
  summarise(Intensity = mean(Intensity, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Pickup_t = factor(Pickup_t, levels = c("T_1", "T_2")))

# VK plot
p_vk_1_2 <- ggplot(plot_df_1_2, aes(x = OC, y = HC)) +
  geom_point(aes(color = Pickup_t), size = 1.8, alpha = 0.9) +
  scale_color_manual(
    values = tp_colors,
    drop   = FALSE,
    name   = "Collection time (years)",
    labels = c(
      "T_1"   = "1 years",
      "T_2" = "2 years"
    )
  ) +
  annotate(
    "segment",
    x = 0.4, xend = 1.1,
    y = 0.1,  yend = 0.1,
    color = "grey40",
    linewidth = 2.8,
    arrow = arrow(length = unit(0.4, "cm"), type = "open")
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           label = "1 year",
           fontface = "bold",
           size = 4,
           color = "grey40") +
  annotate("text",
           x = 1.35, y = 0.1,
           label = "2 years",
           fontface = "bold",
           color = "grey40",
           size = 4) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0,1.4))+
  theme_bw(base_size = 14) +
  labs(
    title = "Unique compounds: 1 vs 2 years",
    x     = "O/C",
    y     = "H/C"
  ) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom", 
    legend.title = element_blank()
  )

p_vk_1_2



#### All timepoints ---------------------------------------------------------
combined_vk_unique <- p_vk_0_0.5 +
  p_vk_0.5_1 +
  p_vk_1_2 +
  plot_layout(ncol = 3, 
              #guides = "collect"
              ) &
  theme(legend.position = "bottom")

combined_vk_unique

ggsave("Plots/Combined_VK_unique_0_0.5_1_2_MAG.png",
       combined_vk_unique, width = 12, height = 5, dpi = 300)


#### Summary VK -----------------------------------------------------------------
presence_df <- df_longer %>%
  filter(Intensity > 0)  

summary_df <- presence_df %>%
  group_by(Pickup_t) %>%
  summarize(
    mean_OC = mean(OC, na.rm = TRUE),
    sd_OC   = sd(OC, na.rm = TRUE),
    mean_HC = mean(HC, na.rm = TRUE),
    sd_HC   = sd(HC, na.rm = TRUE),
    n       = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2")),
    Pickup_label = recode(
      Pickup_t,
      "T_0"   = "0",
      "T_0.5" = "0.5",
      "T_1"   = "1",
      "T_2"   = "2"
    )
  )

vk_summary_plot2 <- ggplot(summary_df, aes(x = mean_OC, y = mean_HC, color = Pickup_t)) +
  geom_point(size = 4) +
  coord_cartesian(xlim = c(0.35, 0.8), ylim = c(0.8, 1.75)) +
  theme_bw(base_size = 12) +
  scale_color_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  labs(
    x = "O/C",
    y = "H/C"
  ) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  ) +
  geom_text_repel(
    aes(label = Pickup_label),
    size        = 6,
    fontface    = "bold",
    show.legend = FALSE
  )

vk_summary_plot2

VK_summary <- vk_summary_plot2 +
  geom_curve(
    aes(x = 0.58, y = 1.5, xend = 0.65, yend = 1.),
    curvature = 1.3,
    angle = 80,
    ncp = 20,                    # higher = smoother curve
    arrow = arrow(type = "closed", length = unit(0.5, "cm")),
    linewidth = 2,
    color = "grey40"
  )

VK_summary
ggsave("Plots/sum_VK_MAG.png",
       VK_summary, width = 6, height = 5, dpi = 300)


# GFE_NOSC ------------------------------------------------------------

presence_df <- df_longer %>%
  filter(Intensity > 0) %>%
  mutate(
    Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2"))
  )

gfe_stat <- presence_df %>%
  tukey_hsd(GFE ~ Pickup_t) %>%   
  add_y_position(step.increase = 1)  

gfe_violin <- ggplot(presence_df, aes(x = Pickup_t, y = GFE, fill = Pickup_t)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  scale_fill_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  scale_x_discrete(labels = c("T_0" = "0", "T_0.5" = "0.5", "T_1" = "1", "T_2" = "2")) +
  stat_pvalue_manual(
    gfe_stat,
    label      = 'p.adj.signif',
    inherit.aes = FALSE,
    hide.ns    = TRUE
  ) +
  theme_bw() +
  labs(
    x = "Collection time (years)",
    y = "Gibbs free energy (GFE)"
  ) +
  theme(
    plot.title      = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

gfe_violin


# NOSC
nosc_stat <- presence_df %>%
  tukey_hsd(NOSC ~ Pickup_t) %>%   
  add_y_position(step.increase = 1)  

nosc_violin <- ggplot(presence_df, aes(x = Pickup_t, y = NOSC, fill = Pickup_t)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  scale_fill_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  scale_x_discrete(labels = c("T_0" = "0", "T_0.5" = "0.5", "T_1" = "1", "T_2" = "2")) +
  stat_pvalue_manual(
    nosc_stat,
    label      = 'p.adj.signif',
    inherit.aes = FALSE,
    hide.ns    = TRUE
  ) +
  theme_bw() +
  labs(
    x = "Collection time (years)",
    y = "NOSC"
  ) +
  theme(
    plot.title      = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

nosc_violin

### GFE per class

presence_df <- presence_df %>%
  filter(!is.na(Class), !is.na(GFE)) %>%
  mutate(Class = recode(Class, "Protein" = "Peptide"))

gfe_stat <- presence_df %>%
  filter(!is.na(Class), !is.na(GFE)) %>%
  group_by(Class) %>%
  tukey_hsd(GFE ~ Pickup_t) %>%
  add_y_position(step.increase = 0.05)

gfe_violin_by_class <- presence_df %>%
  filter(!is.na(Class), !is.na(GFE)) %>%
  mutate(Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2"))) %>%
  ggplot(aes(x = Pickup_t, y = GFE, fill = Pickup_t)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, show.legend = FALSE) +
  scale_fill_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  scale_x_discrete(labels = c("T_0" = "0", "T_0.5" = "0.5", "T_1" = "1", "T_2" = "2")) +
  stat_pvalue_manual(
    gfe_stat,
    label      = "p.adj.signif",
    hide.ns    = TRUE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Class, scales = "free_y") +
  theme_bw() +
  labs(
    title = "GFE by compound class and collection time",
    x = "Collection time (years)",
    y = "Gibbs free energy (GFE)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

gfe_violin_by_class

# Mass distribution with time -------------------------------------------

selected_classes <- c("Lignin", "Tannin", "Lipid", "Condensed hydrocarbon")

mass_class_time <- presence_df %>%
  filter(Class %in% selected_classes, !is.na(Mass), !is.na(Pickup_t)) %>%
  mutate(
    Pickup_t = factor(Pickup_t, levels = c("T_0", "T_0.5", "T_1", "T_2"))
  )

rotated_density2 <- ggplot(mass_class_time, aes(x = Mass, fill = Pickup_t, color = Pickup_t)) +
  geom_density(alpha = 0.3, adjust = 1) +
  scale_fill_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  scale_color_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  facet_grid(
    Pickup_t ~ Class,
    scales   = "free_y",
    labeller = labeller(
      Pickup_t = c("T_0" = "0", "T_0.5" = "0.5", "T_1" = "1", "T_2" = "2")
    )
  ) +
  theme_bw() +
  labs(
    x    = "Mass (Da)",
    y    = "Density",
    fill = "Collection time",
    color = "Collection time"
  ) +
  theme(
    plot.title = element_blank(),
    strip.text = element_text(face = "bold")
  )

rotated_density2


library(scales)

rotated_density2 <- ggplot(mass_class_time, 
                           aes(x = Mass, fill = Pickup_t, color = Pickup_t)) +
  geom_density(alpha = 0.3, adjust = 1) +
  
scale_fill_manual(
  values = tp_colors,
  breaks = c("T_0", "T_0.5", "T_1", "T_2"),
  labels = c("0 years", "0.5 years", "1 year", "2 years"),
  name   = "Collection time"
) +
  scale_color_manual(
    values = tp_colors,
    breaks = c("T_0", "T_0.5", "T_1", "T_2"),
    labels = c("0 years", "0.5 years", "1 year", "2 years"),
    name   = "Collection time"
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.004, by = 0.002),
    labels = scales::label_number(accuracy = 0.001)
  )+
  
facet_grid(
  Pickup_t ~ Class,
  labeller = labeller(
    Pickup_t = c("T_0"="0", "T_0.5"="0.5", "T_1"="1", "T_2"="2"),
    Class = c("Lignin" = "Lignin", "Tannin" = "Tannin", "Lipid" = "Lipid", "Condensed hydrocarbon" = "Condensed HC")
  )
) +
  
theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
    strip.background = element_rect(fill = "grey", colour = "black"),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Mass (Da)",
    y = "Density",
    fill = "Collection time",
    color = "Collection time"
  )

rotated_density2
# SupFig8 --------------------------------------------------------------------

legend_time <- cowplot::get_legend(
  VK_summary +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 13)
    )
)

## Remove legends from all panels
p_vk_0_0.5_noleg <- p_vk_0_0.5  + theme(legend.position = "none", plot.title = element_blank())
p_vk_0.5_1_noleg <- p_vk_0.5_1  + theme(legend.position = "none", plot.title = element_blank())
p_vk_1_2_noleg   <- p_vk_1_2    + theme(legend.position = "none", plot.title = element_blank())
VK_summary_noleg <- VK_summary  + theme(legend.position = "none", plot.title = element_blank())
rotated_density2_noleg <- rotated_density2 +
  theme(legend.position = "none", plot.title = element_blank())

top_row <- p_vk_0_0.5_noleg + p_vk_0.5_1_noleg + p_vk_1_2_noleg +
  plot_layout(ncol = 3)

bottom_row <- (VK_summary_noleg | rotated_density2_noleg) +
  plot_layout(ncol = 2, widths = c(0.6, 1.4))

legend_row <- patchwork::wrap_elements(legend_time)

supfig8_full <- top_row / bottom_row / legend_row +
  plot_layout(heights = c(1.3, 2, 0.25)) +
  plot_annotation(tag_levels = "A")


main_panels <- top_row / bottom_row +
  plot_annotation(tag_levels = "A")

supfig8_full <- main_panels / legend_row +
  plot_layout(heights = c(1.8, 2, 0.25))

supfig8_full


ggsave(
  filename = "Plots/SupFig8_VK_mass_MAG_5panel.png",
  plot     = supfig8_full,
  width    = 11,
  height   = 10,
  dpi      = 300
)


# Sup Fig 10 ---------------------------------------------------------------
top_row <- nosc_violin + gfe_violin +
  plot_layout(ncol = 2)

bottom_row <- gfe_violin_by_class + theme(legend.position = "none", plot.title = element_blank())

Sup_fig10 <- top_row / bottom_row +
  plot_layout(heights = c(1, 3)) +
  plot_annotation(tag_levels = "A")
Sup_fig10


ggsave(
  filename = "Plots/SupFig10_GFE_NOSC_MAG.png",
  plot     = Sup_fig10,
  width    = 8,
  height   = 8,
  dpi      = 300
)
