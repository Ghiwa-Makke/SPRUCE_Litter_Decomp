### SPRUCE Litter decomposition 2015–2017
### Diversity of compound classes – MAG
### Sup Fig 7

# Libraries ---------------------------------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)

# 1. Compute Shannon diversity per class ----------------------------------
## Amino Sugar
AminoSugar_MAG <- read.csv("Input/Shannon/AminoSugar_MAG.csv",
                           header = TRUE, row.names = 1)

mass_intensity_data <- AminoSugar_MAG[, 3:ncol(AminoSugar_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

AminoSugar_MAG$AminoSugar_shannon_Diversity <- shannon_diversity
AminoSugar_MAG <- AminoSugar_MAG %>%
  select(-contains("Mass_"))

## Carbohydrate
Carbohydrate_MAG <- read.csv("Input/Shannon/Carbohydrate_MAG.csv",
                             header = TRUE, row.names = 1)

mass_intensity_data <- Carbohydrate_MAG[, 3:ncol(Carbohydrate_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Carbohydrate_MAG$Carbohydrate_shannon_Diversity <- shannon_diversity
Carbohydrate_MAG <- Carbohydrate_MAG %>%
  select(-contains("Mass_"))

## Condensed HC
CondensedHC_MAG <- read.csv("Input/Shannon/CondensedHC_MAG.csv",
                            header = TRUE, row.names = 1)

mass_intensity_data <- CondensedHC_MAG[, 3:ncol(CondensedHC_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

CondensedHC_MAG$CondensedHC_shannon_Diversity <- shannon_diversity
CondensedHC_MAG <- CondensedHC_MAG %>%
  select(-contains("Mass_"))

## Lignin
Lignin_MAG <- read.csv("Input/Shannon/Lignin_MAG.csv",
                       header = TRUE, row.names = 1)

mass_intensity_data <- Lignin_MAG[, 3:ncol(Lignin_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Lignin_MAG$Lignin_shannon_Diversity <- shannon_diversity
Lignin_MAG <- Lignin_MAG %>%
  select(-contains("Mass_"))

## Lipid
Lipid_MAG <- read.csv("Input/Shannon/Lipid_MAG.csv",
                      header = TRUE, row.names = 1)

mass_intensity_data <- Lipid_MAG[, 3:ncol(Lipid_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Lipid_MAG$Lipid_shannon_Diversity <- shannon_diversity
Lipid_MAG <- Lipid_MAG %>%
  select(-contains("Mass_"))

## Other
Other_MAG <- read.csv("Input/Shannon/Other_MAG.csv",
                      header = TRUE, row.names = 1)

mass_intensity_data <- Other_MAG[, 3:ncol(Other_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Other_MAG$Other_shannon_Diversity <- shannon_diversity
Other_MAG <- Other_MAG %>%
  select(-contains("Mass_"))

## Peptide (Protein class)
Peptide_MAG <- read.csv("Input/Shannon/Peptides_MAG.csv",
                        header = TRUE, row.names = 1)

mass_intensity_data <- Peptide_MAG[, 3:ncol(Peptide_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Peptide_MAG$Peptide_shannon_Diversity <- shannon_diversity
Peptide_MAG <- Peptide_MAG %>%
  select(-contains("Mass_"))

## Tannin
Tannin_MAG <- read.csv("Input/Shannon/Tannin_MAG.csv",
                       header = TRUE, row.names = 1)

mass_intensity_data <- Tannin_MAG[, 3:ncol(Tannin_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Tannin_MAG$Tannin_shannon_Diversity <- shannon_diversity
Tannin_MAG <- Tannin_MAG %>%
  select(-contains("Mass_"))

## Unsaturated HC
Unsaturated_HC_MAG <- read.csv("Input/Shannon/UnsaturatedHC_MAG.csv",
                               header = TRUE, row.names = 1)

mass_intensity_data <- Unsaturated_HC_MAG[, 3:ncol(Unsaturated_HC_MAG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Unsaturated_HC_MAG$Unsaturated_HC_shannon_Diversity <- shannon_diversity
Unsaturated_HC_MAG <- Unsaturated_HC_MAG %>%
  select(-contains("Mass_"))

# 2. Combine into one Shannon table ---------------------------------------

AminoSugar_MAG    <- rownames_to_column(AminoSugar_MAG,    var = "SampleID")
Carbohydrate_MAG  <- rownames_to_column(Carbohydrate_MAG,  var = "SampleID")
CondensedHC_MAG   <- rownames_to_column(CondensedHC_MAG,   var = "SampleID")
Lignin_MAG        <- rownames_to_column(Lignin_MAG,        var = "SampleID")
Lipid_MAG         <- rownames_to_column(Lipid_MAG,         var = "SampleID")
Other_MAG         <- rownames_to_column(Other_MAG,         var = "SampleID")
Peptide_MAG       <- rownames_to_column(Peptide_MAG,       var = "SampleID")
Tannin_MAG        <- rownames_to_column(Tannin_MAG,        var = "SampleID")
Unsaturated_HC_MAG <- rownames_to_column(Unsaturated_HC_MAG, var = "SampleID")

all_shannon <- AminoSugar_MAG %>%
  rename(AminoSugar = AminoSugar_shannon_Diversity) %>%
  left_join(
    Carbohydrate_MAG %>%
      select(SampleID, Carbohydrate_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Carbohydrate = Carbohydrate_shannon_Diversity) %>%
  left_join(
    CondensedHC_MAG %>%
      select(SampleID, CondensedHC_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(CondensedHC = CondensedHC_shannon_Diversity) %>%
  left_join(
    Lignin_MAG %>%
      select(SampleID, Lignin_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Lignin = Lignin_shannon_Diversity) %>%
  left_join(
    Lipid_MAG %>%
      select(SampleID, Lipid_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Lipid = Lipid_shannon_Diversity) %>%
  left_join(
    Other_MAG %>%
      select(SampleID, Other_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Other = Other_shannon_Diversity) %>%
  left_join(
    Peptide_MAG %>%
      select(SampleID, Peptide_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Peptide = Peptide_shannon_Diversity) %>%
  left_join(
    Tannin_MAG %>%
      select(SampleID, Tannin_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Tannin = Tannin_shannon_Diversity) %>%
  left_join(
    Unsaturated_HC_MAG %>%
      select(SampleID, Unsaturated_HC_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Unsaturated_HC = Unsaturated_HC_shannon_Diversity)

write.csv(all_shannon, file = "Tables/MAG_Shannon_Diversity.csv", row.names = FALSE)

# 3. Long format + plot ----------------------------------------------------

all_pivoted_shannon <- all_shannon %>%
  pivot_longer(
    cols      = c(AminoSugar, Carbohydrate, CondensedHC, Lignin,
                  Lipid, Other, Peptide, Tannin, Unsaturated_HC),
    names_to  = "Type",
    values_to = "Shannon_Diversity"
  )

time_colors <- c(
  "T_0"   = "#4DBBD5",
  "T_0.5" = "#3C5488",
  "T_1"   = "#B09C85",
  "T_2"   = "#503A2C"
)

shannon_plot <- ggplot(all_pivoted_shannon,
                       aes(x = Pickup_t, y = Shannon_Diversity)) +
  geom_boxplot(aes(fill = Pickup_t)) +
  geom_jitter(aes(color = Pickup_t),
              width = 0.2, size = 2, alpha = 0.8) +
  labs(
    x = "Collection time (years)",
    y = "Shannon diversity"
  ) +
  theme_classic() +
  scale_x_discrete(labels = c(
    "T_0"   = "0",
    "T_0.5" = "0.5",
    "T_1"   = "1",
    "T_2"   = "2"
  )) +
  scale_fill_manual(
    values = time_colors,
    labels = c(
      "T_0"   = "0 years",
      "T_0.5" = "0.5 years",
      "T_1"   = "1 year",
      "T_2"   = "2 years"
    ),
    name = "Collection time"
  ) +
  scale_color_manual(
    values = time_colors,
    guide  = "none"
  ) +
  facet_wrap(~ Type, scales = "free")

shannon_plot

ggsave("Plots/SupFig7_MAG_shannon.png",
       plot   = shannon_plot,
       width  = 8,
       height = 4.5,
       dpi    = 300)
