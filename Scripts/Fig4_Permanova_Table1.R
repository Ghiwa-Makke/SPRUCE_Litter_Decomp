### SPRUCE Litter decomposition 2015–2017
### Table 1 – Overall PERMANOVA (Adonis)
### Figure 4 – PERMANOVA per litter (Adonis bubble plot)

library(tidyverse)
library(vegan)
library(ggplot2)
library(viridis)
library(scales)
library(purrr)
library(ragg)

custom_labels <- c(
  "MAG" = "S. magellanicum",
  "ANG" = "S. angustifolium",
  "SPR" = "Spruce roots",
  "LTR" = "Labrador tea roots",
  "LTL" = "Labrador tea leaves",
  "SPL" = "Spruce needles"
)

## 1. Overall ------
# Metadata (explanatory variables)
explanatory_data <- read.csv("Input/metadata.csv") %>%
  select(SampleID, Litter, Pickup_t, Plot, Temp, CO2) %>%
  mutate(
    Pickup_t = as.factor(Pickup_t),
    Litter   = as.factor(Litter)
  )

# Response data = class composition matrix (rows = SampleID)
response_data <- read.csv(
  "Input/Metabodirect_All_Sum_Time/1_preprocessing_output/class_composition.csv",
  row.names = 1
)

# Align rows between explanatory and response
common_ids <- intersect(explanatory_data$SampleID, rownames(response_data))

explanatory_data_aligned <- explanatory_data %>%
  filter(SampleID %in% common_ids) %>%
  arrange(match(SampleID, common_ids))

response_data_aligned <- response_data[common_ids, , drop = FALSE]

stopifnot(identical(explanatory_data_aligned$SampleID, rownames(response_data_aligned)))

 
### PERMANOVA (Table 1) -----
adonis_overall <- adonis2(
  response_data_aligned ~ Litter + Pickup_t + Temp + CO2,
  data        = explanatory_data_aligned,
  permutations = 999,
  by          = "terms"
)

# Convert to tidy table
adonis_overall_tbl <- adonis_overall %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )

# Save Table 1
write.csv(
  adonis_overall_tbl,
  "Tables/Table1_Adonis_overall.csv",
  row.names = FALSE
)
 
 
## 2. PERMANOVA per litter -------------------
## Upload data tables:

### ANG ---------------------------------------------------------------------
ANG <- read.csv("Input/Permanova/class_composition_ANG.csv", row.names = 1)
explan_ANG <- explanatory_data %>%
  filter(Litter == "ANG")
  
  
ANG_ids <- intersect(explan_ANG$SampleID, rownames(ANG))

explanatory_ANG_aligned <- explan_ANG %>%
  filter(SampleID %in% ANG_ids) %>%
  arrange(match(SampleID, ANG_ids))

response_ANG_aligned <- ANG[ANG_ids, , drop = F]

adonis_ANG <- adonis2(
  response_ANG_aligned ~ Pickup_t * Temp * CO2,
  data        = explanatory_ANG_aligned,
  by          = "terms"
)

# Convert to tidy table
adonis_ANG_tbl <- adonis_ANG %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )

# Save Table 1
write.csv(
  adonis_ANG_tbl,
  "Tables/adonis_ANG_tbl.csv",
  row.names = FALSE
)
### MAG ---------------------------------------------------------------------
MAG <- read.csv("Input/Permanova/class_composition_MAG.csv", row.names = 1)
explan_MAG <- explanatory_data %>%
  filter(Litter == "MAG") %>%
  filter(SampleID != "Kelly_41")


MAG_ids <- intersect(explan_MAG$SampleID, rownames(MAG))

explanatory_MAG_aligned <- explan_MAG %>%
  filter(SampleID %in% MAG_ids) %>%
  arrange(match(SampleID, MAG_ids))

response_MAG_aligned <- MAG[MAG_ids, , drop = F]

adonis_MAG <- adonis2(
  response_MAG_aligned ~ Pickup_t * Temp * CO2,
  data        = explanatory_MAG_aligned,
  by          = "terms"
)


# Convert to tidy table
adonis_MAG_tbl <- adonis_MAG %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )
# Save Table 1
write.csv(
  adonis_MAG_tbl,
  "Tables/adonis_MAG_tbl.csv",
  row.names = FALSE
)

### LTL ---------------------------------------------------------------------
LTL <- read.csv("Input/Permanova/class_composition_LTL.csv", row.names = 1)
explan_LTL <- explanatory_data %>%
  filter(Litter == "LTL")


LTL_ids <- intersect(explan_LTL$SampleID, rownames(LTL))

explanatory_LTL_aligned <- explan_LTL %>%
  filter(SampleID %in% LTL_ids) %>%
  arrange(match(SampleID, LTL_ids))

response_LTL_aligned <- LTL[LTL_ids, , drop = F]

adonis_LTL <- adonis2(
  response_LTL_aligned ~ Pickup_t * Temp * CO2,
  data        = explanatory_LTL_aligned,
  by          = "terms"
)

# Convert to tidy table
adonis_LTL_tbl <- adonis_LTL %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )

# Save Table 1
write.csv(
  adonis_LTL_tbl,
  "Tables/adonis_LTL_tbl.csv",
  row.names = FALSE
)






### LTR ---------------------------------------------------------------------
LTR <- read.csv("Input/Permanova/class_composition_LTR.csv", row.names = 1)
explan_LTR <- explanatory_data %>%
  filter(Litter == "LTR")


LTR_ids <- intersect(explan_LTR$SampleID, rownames(LTR))

explanatory_LTR_aligned <- explan_LTR %>%
  filter(SampleID %in% LTR_ids) %>%
  arrange(match(SampleID, LTR_ids))

response_LTR_aligned <- LTR[LTR_ids, , drop = F]

adonis_LTR <- adonis2(
  response_LTR_aligned ~ Pickup_t * Temp * CO2,
  data        = explanatory_LTR_aligned,
  by          = "terms"
)

# Convert to tidy table
adonis_LTR_tbl <- adonis_LTR %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )

# Save Table 1
write.csv(
  adonis_LTR_tbl,
  "Tables/adonis_LTR_tbl.csv",
  row.names = FALSE
)


### SPL ---------------------------------------------------------------------
SPL <- read.csv("Input/Permanova/class_composition_SPL.csv", row.names = 1)
explan_SPL <- explanatory_data %>%
  filter(Litter == "SPL")


SPL_ids <- intersect(explan_SPL$SampleID, rownames(SPL))

explanatory_SPL_aligned <- explan_SPL %>%
  filter(SampleID %in% SPL_ids) %>%
  arrange(match(SampleID, SPL_ids))

response_SPL_aligned <- SPL[SPL_ids, , drop = F]

adonis_SPL <- adonis2(
  response_SPL_aligned ~ Pickup_t * Temp * CO2,
  data        = explanatory_SPL_aligned,
  by          = "terms"
)

# Convert to tidy table
adonis_SPL_tbl <- adonis_SPL %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )

# Save Table 1
write.csv(
  adonis_SPL_tbl,
  "Tables/adonis_SPL_tbl.csv",
  row.names = FALSE
)

### SPR ---------------------------------------------------------------------
SPR <- read.csv("Input/Permanova/class_composition_SPR.csv", row.names = 1)
explan_SPR <- explanatory_data %>%
  filter(Litter == "SPR")


SPR_ids <- intersect(explan_SPR$SampleID, rownames(SPR))

explanatory_SPR_aligned <- explan_SPR %>%
  filter(SampleID %in% SPR_ids) %>%
  arrange(match(SampleID, SPR_ids))

response_SPR_aligned <- SPR[SPR_ids, , drop = F]

adonis_SPR <- adonis2(
  response_SPR_aligned ~ Pickup_t * Temp * CO2,
  data        = explanatory_SPR_aligned,
  by          = "terms"
)

# Convert to tidy table
adonis_SPR_tbl <- adonis_SPR %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(Term != "Total") %>%
  mutate(
    R2                 = R2,
    Percent_explained  = 100 * R2
  ) %>%
  select(
    Term,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Percent_explained
  )

# Save Table 1
write.csv(
  adonis_SPR_tbl,
  "Tables/adonis_SPR_tbl.csv",
  row.names = FALSE
)



# Combine all per-litter tables and add a Litter column
adonis_per_litter_all <- list(
  ANG = adonis_ANG_tbl,
  MAG = adonis_MAG_tbl,
  LTL = adonis_LTL_tbl,
  LTR = adonis_LTR_tbl,
  SPL = adonis_SPL_tbl,
  SPR = adonis_SPR_tbl
) %>%
  imap_dfr(~ mutate(.x, Litter = .y, .before = 1))

adonis_per_litter_all

# Save combined PERMANOVA results (all litter types)
write.csv(
  adonis_per_litter_all,
  "Tables/Adonis_per_litter_all.csv",
  row.names = FALSE
)


# 3. Plotting Per litter Effect -------------------------------------------
plot_df <- adonis_per_litter_all %>%
  filter(Term %in% c(
    "Pickup_t",
    "Temp",
    "CO2",
    "Pickup_t:Temp",
    "Pickup_t:CO2",
    "Temp:CO2"
  )) %>%
  mutate(
    Type   = Litter,
    Factors = recode(
      Term,
      "Pickup_t"          = "Collection time",
      "Temp"              = "Temp",
      "CO2"               = "CO2",
      "Pickup_t:Temp"     = "Collection time & Temp",
      "Pickup_t:CO2"      = "Collection time & CO2",
      "Temp:CO2"          = "Temp & CO2"
    ),
    Percentage  = Percent_explained,
    P.value     = `Pr(>F)`,
    Significance = case_when(
      P.value < 0.001 ~ "***",
      P.value < 0.01  ~ "**",
      P.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  select(Type, Factors, Percentage, P.value, Significance)

sample_order <- c("ANG", "MAG", "LTL", "LTR", "SPL", "SPR")
factor_order <- c(
  "Collection time",
  "Temp",
  "CO2",
  "Collection time & Temp",
  "Collection time & CO2",
  "Temp & CO2"
)

plot_df <- plot_df %>%
  mutate(
    Type    = factor(Type,    levels = sample_order),
    Factors = factor(Factors, levels = factor_order)
  )

## axis labels with italics and CO2 ---------------------------------------

custom_labels_x_expr <- c(
  "ANG" = expression(italic("S. angustifolium")),
  "MAG" = expression(italic("S. magellanicum")),
  "LTL" = expression("Labrador tea leaves"),
  "LTR" = expression("Labrador tea roots"),
  "SPL" = expression("Spruce needles"),
  "SPR" = expression("Spruce roots")
)

custom_labels_y <- c(
  "Collection time"        = "Collection time",
  "Temp"                   = "Temp",
  "CO2"                    = expression(CO[2]),
  "Collection time & Temp" = "Collection time & Temp",
  "Collection time & CO2"  = expression("Collection time & " ~ CO[2]),
  "Temp & CO2"             = expression("Temp & " ~ CO[2])
)

plot_df <- plot_df %>%
  mutate(
    Significance = if_else(P.value < 0.05, "*", "")
  )

## Fig 4 bubble plot ------------------------------------------------------

adonis_plot <- ggplot(plot_df, aes(x = Type, y = Factors)) +
  geom_point(
    aes(fill = P.value, size = Percentage),
    shape  = 21,
    color  = "black",
    stroke = 0.7
  ) +
  geom_text(
    aes(label = Significance),
    color    = "black",
    size     = 4,
    vjust    = 0.7,
    fontface = "bold",
    na.rm    = TRUE
  ) +
  scale_size_continuous(
    name  = "Variance explained (%)",
    range = c(3, 10)
  ) +
  scale_fill_gradient(
    name   = "PERMANOVA p-value",
    low    = "#f5f5dc",   # beige
    high   = "#006d2c"
  ) +
  scale_y_discrete(labels = custom_labels_y) +
  scale_x_discrete(labels = custom_labels_x_expr) +
  labs(
    x = "Litter type",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(hjust = 0, face = "bold", size = 11),
    axis.text.x      = element_text(face = "bold", angle = 30, hjust = 1, vjust = 1),
    axis.text.y      = element_text(face = "bold", size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9),
    plot.margin      = margin(t = 5, r = 2, b = 5, l = 5)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth       = 0.4,
      barheight      = 4,
      title.position = "top",
      title.hjust    = 0.5
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust    = 0.5,
      override.aes   = list(shape = 21, fill = "grey90", color = "black")
    )
  )

adonis_plot

ggsave(
  filename = "Plots/Updated/Fig4_adonis_drivers_beige_green.tiff",
  plot     = adonis_plot,
  width    = 150,
  height   = 120,
  units    = "mm",
  dpi      = 300,
  device   = agg_tiff,
  compression = "lzw",
  bg       = "white"
)
