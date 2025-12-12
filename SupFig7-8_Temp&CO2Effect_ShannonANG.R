### SPRUCE Litter decomposition 2015–2017
### Diversity of compound classes – ANG
### Sup Figure -- Effect of CO2 and Temperature

# Libraries ---------------------------------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)

# Data ------
ANG_Shannon_Diversity <- read.csv("Tables/ANG_Shannon_Diversity.csv")
metadata <- read.csv("Input/metadata.csv")


# Temperature Effect ------------------------------------------------------
metadata_Temp <- metadata %>%
  filter(Litter == "ANG")%>%
  select(SampleID, Temp) %>%
  mutate(Temp = ifelse(Temp == 0, "+0°C", 
                       ifelse(Temp == 4.5, "+4.5°C",
                              ifelse(Temp == 9, "+9°C",Temp))))

ANG_Shannon_Diversity <- merge(ANG_Shannon_Diversity, metadata_Temp, by = "SampleID")

Diversity_Temp0 <- ANG_Shannon_Diversity %>%
  filter(Temp == "None" | Temp == "+0°C") %>%
  mutate(Temp = ifelse(Temp == "None", "+0°C", Temp))

Diversity_Temp4.5 <- ANG_Shannon_Diversity %>%
  filter(Temp == "None" | Temp == "+4.5°C") %>%
  mutate(Temp = ifelse(Temp == "None", "+4.5°C", Temp))

Diversity_Temp9 <- ANG_Shannon_Diversity %>%
  filter(Temp == "None" | Temp == "+9°C") %>%
  mutate(Temp = ifelse(Temp == "None", "+9°C", Temp))

Diversity_Temp_all <- bind_rows(Diversity_Temp0, Diversity_Temp4.5, Diversity_Temp9)


time_colors <- c(
  "T_0"   = "#4DBBD5",
  "T_0.5" = "#3C5488",
  "T_1"   = "#B09C85",
  "T_2"   = "#503A2C"
)


theme_shannon <- theme_classic(base_size = 10) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 10),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text   = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

make_shannon_plot <- function(df, class_var, title) {
  ggplot(df, aes(x = Pickup_t, y = {{ class_var }}, fill = Pickup_t)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = Pickup_t), width = 0.2, size = 2, alpha = 0.9) +
    labs(
      title = title,
      x = "",
      y = "Shannon Diversity"
    ) +
    scale_x_discrete(labels = c(
      "T_0"   = "0",
      "T_0.5" = "0.5",
      "T_1"   = "1",
      "T_2"   = "2"
    )) +
    scale_fill_manual(values = time_colors) +
    scale_color_manual(values = time_colors) +
    facet_wrap(~ Temp) +
    theme_shannon
}

ANG_plot_AminoSugar   <- make_shannon_plot(Diversity_Temp_all, AminoSugar,   "AminoSugar")
ANG_plot_Carbohydrate <- make_shannon_plot(Diversity_Temp_all, Carbohydrate, "Carbohydrate")
ANG_plot_CondensedHC <- make_shannon_plot(Diversity_Temp_all, CondensedHC , "Condensed HC")
ANG_plot_Lignin <- make_shannon_plot(Diversity_Temp_all, Lignin , "Lignin")
ANG_plot_Lipid <- make_shannon_plot(Diversity_Temp_all, Lipid , "Lipid")
ANG_plot_Other <- make_shannon_plot(Diversity_Temp_all, Other , "Other")
ANG_plot_Peptide <- make_shannon_plot(Diversity_Temp_all, Peptide , "Peptide")
ANG_plot_Tannin <- make_shannon_plot(Diversity_Temp_all, Tannin , "Tannin")
ANG_plot_UnsaturatedHC <- make_shannon_plot(Diversity_Temp_all, Unsaturated_HC , "Unsaturated HC")


library(patchwork)
ANG_Temp_Plot <- (ANG_plot_AminoSugar + ANG_plot_Carbohydrate + ANG_plot_CondensedHC + ANG_plot_Lignin + ANG_plot_Lipid + ANG_plot_Other + ANG_plot_Peptide + ANG_plot_Tannin + ANG_plot_UnsaturatedHC)+
  plot_layout(nrow = 3, byrow = TRUE) +
  plot_annotation(title = "Sphagnum ANG - Temp Effect",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
ANG_Temp_Plot

ggsave("Plots/SupFig_ANG_shannon_temp.png",
       plot   = ANG_Temp_Plot,
       width  = 12,
       height = 8,
       dpi    = 300)

# CO2 effect --------------------------------------------------------------
metadata_CO2 <- metadata %>%
  filter(Litter == "ANG")%>%
  select(SampleID, CO2) %>%
  mutate(CO2 = ifelse(CO2 == "aCO2", "Ambient", 
                       ifelse(CO2 == "eCO2", "Elevated", CO2)))

ANG_Shannon_Diversity <- merge(ANG_Shannon_Diversity, metadata_CO2, by = "SampleID")

Diversity_Ambient <- ANG_Shannon_Diversity %>%
  filter(CO2 == "None" | CO2 == "Ambient") %>%
  mutate(CO2 = ifelse(CO2 == "None", "Ambient", CO2))

Diversity_Elevated <- ANG_Shannon_Diversity %>%
  filter(CO2 == "None" | CO2 == "Elevated") %>%
  mutate(CO2 = ifelse(CO2 == "None", "Elevated", CO2))

Diversity_CO2_all <- bind_rows(Diversity_Ambient, Diversity_Elevated)


time_colors <- c(
  "T_0"   = "#4DBBD5",
  "T_0.5" = "#3C5488",
  "T_1"   = "#B09C85",
  "T_2"   = "#503A2C"
)


theme_shannon <- theme_classic(base_size = 10) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 10),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text   = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

make_shannon_plot <- function(df, class_var, title) {
  ggplot(df, aes(x = Pickup_t, y = {{ class_var }}, fill = Pickup_t)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = Pickup_t), width = 0.2, size = 2, alpha = 0.9) +
    labs(
      title = title,
      x = "",
      y = "Shannon Diversity"
    ) +
    scale_x_discrete(labels = c(
      "T_0"   = "0",
      "T_0.5" = "0.5",
      "T_1"   = "1",
      "T_2"   = "2"
    )) +
    scale_fill_manual(values = time_colors) +
    scale_color_manual(values = time_colors) +
    facet_wrap(~ CO2) +
    theme_shannon
}

ANG_plot_AminoSugar   <- make_shannon_plot(Diversity_CO2_all, AminoSugar,   "AminoSugar")
ANG_plot_Carbohydrate <- make_shannon_plot(Diversity_CO2_all, Carbohydrate, "Carbohydrate")
ANG_plot_CondensedHC <- make_shannon_plot(Diversity_CO2_all, CondensedHC , "Condensed HC")
ANG_plot_Lignin <- make_shannon_plot(Diversity_CO2_all, Lignin , "Lignin")
ANG_plot_Lipid <- make_shannon_plot(Diversity_CO2_all, Lipid , "Lipid")
ANG_plot_Other <- make_shannon_plot(Diversity_CO2_all, Other , "Other")
ANG_plot_Peptide <- make_shannon_plot(Diversity_CO2_all, Peptide , "Peptide")
ANG_plot_Tannin <- make_shannon_plot(Diversity_CO2_all, Tannin , "Tannin")
ANG_plot_UnsaturatedHC <- make_shannon_plot(Diversity_CO2_all, Unsaturated_HC , "Unsaturated HC")


library(patchwork)
ANG_CO2_Plot <- (ANG_plot_AminoSugar + ANG_plot_Carbohydrate + ANG_plot_CondensedHC + ANG_plot_Lignin + ANG_plot_Lipid + ANG_plot_Other + ANG_plot_Peptide + ANG_plot_Tannin + ANG_plot_UnsaturatedHC)+
  plot_layout(nrow = 3, byrow = TRUE) +
  plot_annotation(title = "Sphagnum ANG - Warming Effect",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
ANG_CO2_Plot

ggsave("Plots/SupFig8_ANG_shannon_CO2.png",
       plot   = ANG_CO2_Plot,
       width  = 10,
       height = 6.5,
       dpi    = 300)

