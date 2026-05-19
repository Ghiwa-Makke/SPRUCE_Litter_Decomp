### SPRUCE Litter decomposition 2015–2017
### Diversity of compound classes – ANG
### Figure 5

# Libraries ---------------------------------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)

# 1. Compute Shannon diversity per class ----------------------------------
## Amino Sugar
AminoSugar_ANG <- read.csv("Input/Shannon/AminoSugar_ANG.csv",
                           header = TRUE, row.names = 1)

mass_intensity_data <- AminoSugar_ANG[, 3:ncol(AminoSugar_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

AminoSugar_ANG$AminoSugar_shannon_Diversity <- shannon_diversity
AminoSugar_ANG <- AminoSugar_ANG %>%
  select(-contains("Mass_"))

## Carbohydrate
Carbohydrate_ANG <- read.csv("Input/Shannon/Carbohydrate_ANG.csv",
                             header = TRUE, row.names = 1)

mass_intensity_data <- Carbohydrate_ANG[, 3:ncol(Carbohydrate_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Carbohydrate_ANG$Carbohydrate_shannon_Diversity <- shannon_diversity
Carbohydrate_ANG <- Carbohydrate_ANG %>%
  select(-contains("Mass_"))

## Condensed HC
CondensedHC_ANG <- read.csv("Input/Shannon/CondensedHC_ANG.csv",
                            header = TRUE, row.names = 1)

mass_intensity_data <- CondensedHC_ANG[, 3:ncol(CondensedHC_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

CondensedHC_ANG$CondensedHC_shannon_Diversity <- shannon_diversity
CondensedHC_ANG <- CondensedHC_ANG %>%
  select(-contains("Mass_"))

## Lignin
Lignin_ANG <- read.csv("Input/Shannon/Lignin_ANG.csv",
                       header = TRUE, row.names = 1)

mass_intensity_data <- Lignin_ANG[, 3:ncol(Lignin_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Lignin_ANG$Lignin_shannon_Diversity <- shannon_diversity
Lignin_ANG <- Lignin_ANG %>%
  select(-contains("Mass_"))

## Lipid
Lipid_ANG <- read.csv("Input/Shannon/Lipid_ANG.csv",
                      header = TRUE, row.names = 1)

mass_intensity_data <- Lipid_ANG[, 3:ncol(Lipid_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Lipid_ANG$Lipid_shannon_Diversity <- shannon_diversity
Lipid_ANG <- Lipid_ANG %>%
  select(-contains("Mass_"))

## Other
Other_ANG <- read.csv("Input/Shannon/Other_ANG.csv",
                      header = TRUE, row.names = 1)

mass_intensity_data <- Other_ANG[, 3:ncol(Other_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Other_ANG$Other_shannon_Diversity <- shannon_diversity
Other_ANG <- Other_ANG %>%
  select(-contains("Mass_"))

## Peptide (Protein class)
Peptide_ANG <- read.csv("Input/Shannon/Protein_ANG.csv",
                        header = TRUE, row.names = 1)

mass_intensity_data <- Peptide_ANG[, 3:ncol(Peptide_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Peptide_ANG$Peptide_shannon_Diversity <- shannon_diversity
Peptide_ANG <- Peptide_ANG %>%
  select(-contains("Mass_"))

## Tannin
Tannin_ANG <- read.csv("Input/Shannon/Tannin_ANG.csv",
                       header = TRUE, row.names = 1)

mass_intensity_data <- Tannin_ANG[, 3:ncol(Tannin_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Tannin_ANG$Tannin_shannon_Diversity <- shannon_diversity
Tannin_ANG <- Tannin_ANG %>%
  select(-contains("Mass_"))

## Unsaturated HC
Unsaturated_HC_ANG <- read.csv("Input/Shannon/UnsaturatedHC_ANG.csv",
                               header = TRUE, row.names = 1)

mass_intensity_data <- Unsaturated_HC_ANG[, 3:ncol(Unsaturated_HC_ANG)]
norm_mass_intensity_data <- decostand(mass_intensity_data, method = "total")
shannon_diversity <- diversity(norm_mass_intensity_data, index = "shannon")

Unsaturated_HC_ANG$Unsaturated_HC_shannon_Diversity <- shannon_diversity
Unsaturated_HC_ANG <- Unsaturated_HC_ANG %>%
  select(-contains("Mass_"))

# 2. Combine into one Shannon table ---------------------------------------

AminoSugar_ANG    <- rownames_to_column(AminoSugar_ANG,    var = "SampleID")
Carbohydrate_ANG  <- rownames_to_column(Carbohydrate_ANG,  var = "SampleID")
CondensedHC_ANG   <- rownames_to_column(CondensedHC_ANG,   var = "SampleID")
Lignin_ANG        <- rownames_to_column(Lignin_ANG,        var = "SampleID")
Lipid_ANG         <- rownames_to_column(Lipid_ANG,         var = "SampleID")
Other_ANG         <- rownames_to_column(Other_ANG,         var = "SampleID")
Peptide_ANG       <- rownames_to_column(Peptide_ANG,       var = "SampleID")
Tannin_ANG        <- rownames_to_column(Tannin_ANG,        var = "SampleID")
Unsaturated_HC_ANG <- rownames_to_column(Unsaturated_HC_ANG, var = "SampleID")

all_shannon <- AminoSugar_ANG %>%
  rename(AminoSugar = AminoSugar_shannon_Diversity) %>%
  left_join(
    Carbohydrate_ANG %>%
      select(SampleID, Carbohydrate_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Carbohydrate = Carbohydrate_shannon_Diversity) %>%
  left_join(
    CondensedHC_ANG %>%
      select(SampleID, CondensedHC_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(CondensedHC = CondensedHC_shannon_Diversity) %>%
  left_join(
    Lignin_ANG %>%
      select(SampleID, Lignin_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Lignin = Lignin_shannon_Diversity) %>%
  left_join(
    Lipid_ANG %>%
      select(SampleID, Lipid_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Lipid = Lipid_shannon_Diversity) %>%
  left_join(
    Other_ANG %>%
      select(SampleID, Other_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Other = Other_shannon_Diversity) %>%
  left_join(
    Peptide_ANG %>%
      select(SampleID, Peptide_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Peptide = Peptide_shannon_Diversity) %>%
  left_join(
    Tannin_ANG %>%
      select(SampleID, Tannin_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Tannin = Tannin_shannon_Diversity) %>%
  left_join(
    Unsaturated_HC_ANG %>%
      select(SampleID, Unsaturated_HC_shannon_Diversity),
    by = "SampleID"
  ) %>%
  rename(Unsaturated_HC = Unsaturated_HC_shannon_Diversity)

write.csv(all_shannon, file = "Tables/ANG_Shannon_Diversity.csv", row.names = FALSE)

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

strip_labels <- c(
  AminoSugar     = "Amino Sugar",
  Carbohydrate   = "Carbohydrate",
  CondensedHC    = "Condensed HC",
  Lignin         = "Lignin",
  Lipid          = "Lipid",
  Other          = "Other",
  Peptide        = "Peptide",
  Tannin         = "Tannin",
  Unsaturated_HC = "Unsaturated HC"
)

#my_colors <- c("#4DBBD5", "#3C5488", "#B09C85", "#7E6148")

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
  facet_wrap(
    ~ Type,
    scales  = "free_y",
    ncol    = 3,
    labeller = labeller(Type = strip_labels)
  )

shannon_plot

ggsave("Plots/Fig5_ANG_shannon.png",
       plot   = shannon_plot,
       width  = 8,
       height = 4.5,
       dpi    = 300)


# Export TIFF in GCA style
ggsave(
  "Plots/Updated/Figure5.tiff",
  plot   = shannon_plot,
  width  = 190,
  height = 120,
  units  = "mm",
  dpi    = 300,
  device = agg_tiff,
  compression = "lzw"
)
