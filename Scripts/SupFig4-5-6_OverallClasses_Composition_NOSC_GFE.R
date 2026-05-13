### SPRUCE Litter decomposition 2015-2017 
## Sup Fig 4-5-6 - overall shift in compound classes: - Metabodirect outputs 
##   Sup Fig 4: Significant molecular classes over time
##   Sup Fig 5: Significant elemental composition over time
##   Sup Fig 6: NOSC (A) and GFE (B) with time (weighted)


# Load Libraries
library(tidyverse)
library(ggpubr)
library(rstatix)
library(rlang)
source("Input/Metabodirect_All_Sum_Time/custom_functions.R")  # for calculate_weighted()

## Load tables
df <- readr::read_csv("Input/Metabodirect_All_Sum_Time/1_preprocessing_output/Report_processed_MolecFormulas.csv")
metadata <- readr::read_csv("Input/metadata.csv")
el_comp   <- readr::read_csv("Input/Metabodirect_All_Sum_Time/1_preprocessing_output/elemental_composition.csv")
class_comp<- readr::read_csv("Input/Metabodirect_All_Sum_Time/1_preprocessing_output/class_composition.csv")

group1      <- "Pickup_t"  # grouping variable

# Reformat data
metadata <- metadata %>%
  mutate(across(!SampleID, as.factor))

metadata <- metadata %>%
  mutate(!!group1 := forcats::fct_inorder(.data[[group1]]))

# Long intensity table, with metadata
df_longer <- df %>%
  pivot_longer(cols = metadata$SampleID,
               names_to  = "SampleID",
               values_to = "NormIntensity") %>%
  filter(NormIntensity != 0) %>%
  left_join(metadata, by = "SampleID")

# Class composition long, with metadata
class_comp <- class_comp %>%
  pivot_longer(cols = -SampleID,
               names_to  = "Class",
               values_to = "Count") %>%
  left_join(metadata, by = "SampleID")

# Elemental composition long, with metadata
el_comp <- el_comp %>%
  pivot_longer(cols = -SampleID,
               names_to  = "Element",
               values_to = "Count") %>%
  left_join(metadata, by = "SampleID")


my_colors <- c("#4DBBD5", "#3C5488", "#B09C85", "#503A2C")

# General formulas for stats
form_kw   <- reformulate(group1, response = "Rel")
form_dunn <- reformulate(group1, response = "Rel")


## Sup Fig 4 – Significant molecular classes over time ------


# Relative abundance of each class within sample
class_box_df <- class_comp %>%
  group_by(SampleID) %>%
  mutate(
    Total = sum(Count, na.rm = TRUE),
    Rel   = 100 * Count / ifelse(Total == 0, NA, Total)
  ) %>%
  ungroup() %>%
  mutate(
    Class = dplyr::recode(Class,
                          "Protein" = "Peptide")
  )

# Kruskal–Wallis per Class + BH across classes
kw_classes <- class_box_df %>%
  group_by(Class) %>%
  kruskal_test(form_kw) %>%
  adjust_pvalue(method = "BH") %>%
  arrange(p.adj)

sig_classes <- kw_classes %>%
  filter(p.adj < 0.05) %>%
  pull(Class) %>%
  unique()

# Fallback: if none pass FDR, show top few classes
if (length(sig_classes) == 0) {
  message("No classes significant at BH q<0.05; plotting top 6 by p.adj.")
  sig_classes <- kw_classes %>%
    slice_min(p.adj, n = min(6, n())) %>%
    pull(Class)
}

# Dunn post-hoc only for significant classes
dunn_classes <- class_box_df %>%
  filter(Class %in% sig_classes) %>%
  group_by(Class) %>%
  dunn_test(form_dunn, p.adjust.method = "BH") %>%
  add_xy_position(x = group1,
    step.increase = 0.05)   

# Plot – significant molecular classes over time
p_box_class_sig <- ggplot(
  class_box_df %>% filter(Class %in% sig_classes),
  aes(x = .data[[group1]], y = Rel, fill = .data[[group1]])
) +
  geom_boxplot(outlier.shape = 16, alpha = 0.9) +
  facet_wrap(~ Class, scales = "free_y") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = c(
    "T_0"   = "0",
    "T_0.5" = "0.5",
    "T_1"   = "1",
    "T_2"   = "2"
  )) +
  labs(
    x = "Collection time (years)",
    y = "Relative abundance (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "none",
    strip.text        = element_text(face = "bold"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  ) +
  ggpubr::stat_pvalue_manual(
    dunn_classes,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  )
p_box_class_sig
ggsave("Plots/SupFig4_Significant_molecular_classes_over_time.png",
  p_box_class_sig, dpi = 300, width = 10, height = 6
)


## Sup Fig 5 – Significant elemental composition over time ------

# Relative abundance of each element within sample
el_box_df <- el_comp %>%
  group_by(SampleID) %>%
  mutate(
    Total = sum(Count, na.rm = TRUE),
    Rel   = 100 * Count / ifelse(Total == 0, NA, Total)
  ) %>%
  ungroup()

# Kruskal–Wallis per Element + BH across elements
kw_elements <- el_box_df %>%
  group_by(Element) %>%
  kruskal_test(form_kw) %>%
  adjust_pvalue(method = "BH") %>%
  arrange(p.adj)

sig_elements <- kw_elements %>%
  filter(p.adj < 0.05) %>%
  pull(Element) %>%
  unique()

# Fallback: if none pass FDR, show top few elements
if (length(sig_elements) == 0) {
  message("No elements significant at BH q<0.05; plotting top 6 by p.adj.")
  sig_elements <- kw_elements %>%
    slice_min(p.adj, n = min(6, n())) %>%
    pull(Element)
}

# Dunn post-hoc only for significant elements
dunn_elements <- el_box_df %>%
  filter(Element %in% sig_elements) %>%
  group_by(Element) %>%
  dunn_test(form_dunn, p.adjust.method = "BH") %>%
  add_xy_position(
    x             = group1,
    step.increase = 0.05)

# Plot – significant elemental composition over time
p_box_el_sig <- ggplot(
  el_box_df %>% filter(Element %in% sig_elements),
  aes(x = .data[[group1]], y = Rel, fill = .data[[group1]])
) +
  geom_boxplot(outlier.shape = 16, alpha = 0.9) +
  facet_wrap(~ Element, scales = "free_y") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = c(
    "T_0"   = "0",
    "T_0.5" = "0.5",
    "T_1"   = "1",
    "T_2"   = "2"
  )) +
  labs(
    x = "Collection time (years)",
    y = "Relative abundance (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "none",
    strip.text        = element_text(face = "bold"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  ) +
  ggpubr::stat_pvalue_manual(
    dunn_elements,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  )
p_box_el_sig
ggsave(
  file.path("Plots/SupFig5_Significant_elemental_composition_over_time.png"),
  p_box_el_sig, dpi = 300, width = 10, height = 6
)


## Sup Fig 6 – NOSC (A) and GFE (B) with time ----
library(rstatix)
library(rlang)
library(ggpubr)

# 1) long df 
ng_df <- df_longer %>%
  filter(!is.na(NOSC), !is.na(GFE)) %>%
  select(SampleID, !!sym(group1), NOSC, GFE) %>%
  pivot_longer(cols = c(NOSC, GFE),
               names_to = "Index",
               values_to = "value") %>%
  mutate(
    Index = factor(Index, levels = c("NOSC", "GFE"))
  )

# Helper formulas
form_kw_ng   <- reformulate(group1, response = "value")
form_dunn_ng <- form_kw_ng


## A) NOSC plot: 

ng_NOSC <- ng_df %>% filter(Index == "NOSC")

# Stats for NOSC
kw_NOSC <- ng_NOSC %>%
  group_by(Index) %>%
  kruskal_test(form_kw_ng)

range_NOSC <- range(ng_NOSC$value, na.rm = TRUE)
step_NOSC  <- diff(range_NOSC) * 0.05  

dunn_NOSC <- ng_NOSC %>%
  group_by(Index) %>%
  dunn_test(form_dunn_ng, p.adjust.method = "BH") %>%
  add_xy_position(
    x             = group1,
    step.increase = 0   
  ) %>%
  arrange(y.position, .by_group = TRUE) %>%   
  mutate(
    y.position = y.position + (row_number() - 1) * step_NOSC
  )

p_NOSC <- ggplot(
  ng_NOSC,
  aes(x = .data[[group1]], y = value, fill = .data[[group1]])
) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8, color = "black") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = c(
    "T_0"   = "0",
    "T_0.5" = "0.5",
    "T_1"   = "1",
    "T_2"   = "2"
  )) +
  labs(
    x = "Collection time (years)",
    y = "NOSC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  ) +
  ggpubr::stat_pvalue_manual(
    dunn_NOSC,
    label        = "p.adj.signif",
    hide.ns      = TRUE,
    tip.length   = 0.01,
    bracket.size = 0.3
  )
p_NOSC

## B GFE plot -----------

ng_GFE <- ng_df %>% filter(Index == "GFE")

# Stats for GFE
kw_GFE <- ng_GFE %>%
  group_by(Index) %>%
  kruskal_test(form_kw_ng)

dunn_GFE <- ng_GFE %>%
  group_by(Index) %>%
  dunn_test(form_dunn_ng, p.adjust.method = "BH") %>%
  add_xy_position(
    x             = group1,
    step.increase = 1
  )

p_GFE <- ggplot(
  ng_GFE,
  aes(x = .data[[group1]], y = value, fill = .data[[group1]])
) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8, color = "black") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = c(
    "T_0"   = "0",
    "T_0.5" = "0.5",
    "T_1"   = "1",
    "T_2"   = "2"
  )) +
  labs(
    x = "Collection time (years)",
    y = "GFE"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  ) +
  ggpubr::stat_pvalue_manual(
    dunn_GFE,
    label        = "p.adj.signif",
    hide.ns      = TRUE,
    tip.length   = 0.01,
    bracket.size = 0.3
  )
p_GFE

## Sup Fig 6 ------
sup_fig6 <- ggpubr::ggarrange(
  p_NOSC, p_GFE,
  ncol   = 2, nrow = 1,
  labels = c("A", "B"),
  widths = c(1, 1)
)

sup_fig6

ggsave(
  file.path("Plots/SupFig6_NOSC_GFE_violin_boxplots_over_time.png"),
  sup_fig6, dpi = 300, width = 8.5, height = 4
)

