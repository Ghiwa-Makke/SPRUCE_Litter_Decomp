### SPRUCE Litter decomposition 2015–2017
### Diversity of compound classes – All Litter comparison 
### Figure 7

## libraries:
library(tidyverse)
library(vegan)        
library(janitor)      
parse_time <- function(x) as.numeric(gsub("^T_", "", x))  

## data: 
all <- read.csv("Input/Shannon/All.csv", check.names = F)
classes <- all %>%
  column_to_rownames(var = "Mass")%>%
  t()%>%
  as.data.frame()

df <- classes%>%
  rownames_to_column(var = "Mass") %>%
  select(-Class)

class_meta <- classes%>%
  rownames_to_column(var = "Mass") %>%
  select(Mass, Class)

rm(all)
rm(classes)

sample_metadata <- read.csv("Input/metadata.csv")
sample_metadata <- sample_metadata%>%
  select(SampleID, Litter, Pickup_t)


# color palette
my_colors <- c(
  "LTL" = "#e6ab02",
  "LTR" =  "#d95f02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291",
  "SPR" = "#1b9e77"
)


## Intensity matrix
intensity_matrix <- df %>%
  pivot_longer(-Mass, names_to = "SampleID", values_to = "intensity") %>%
  mutate(intensity = suppressWarnings(as.numeric(intensity)))%>%
  pivot_wider(names_from = "SampleID", values_from = "intensity", values_fill = 0) %>%
  select(Mass, all_of(sample_metadata$SampleID)) %>%
  column_to_rownames(var = 'Mass') %>% 
  t()
## Sum normalize intensities
norm_intensity_matrix <- decostand(intensity_matrix, method = 'total')

df_long <- as.data.frame(t(norm_intensity_matrix)) %>%
  rownames_to_column(var = 'Mass')%>%
  pivot_longer(-Mass, names_to = "SampleID", values_to = "intensity") %>%
  left_join(class_meta, by = "Mass") %>%
  left_join(sample_metadata %>% select(SampleID, Litter, Pickup_t),
            by = "SampleID") %>%
  mutate(
    Pickup_t = factor(Pickup_t, levels = c("T_0","T_0.5","T_1","T_2"), ordered = TRUE),
    t_num    = parse_time(as.character(Pickup_t))
  )

div_by_sample <- df_long %>%
  group_by(SampleID, Litter, Pickup_t, t_num, Class) %>%
  summarise(Shannon = diversity(intensity, index = "shannon"),.groups = "drop")

sum_shannon <- div_by_sample %>%
  group_by(Class, Litter, Pickup_t, t_num) %>%
  summarise(
    n            = dplyr::n(),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon   = sd(Shannon,   na.rm = TRUE),
    se_shannon   = sd_shannon / sqrt(n),
    .groups = "drop"
  )


# Plotting ----------------------------------------------------------------
litter_order  <- c("ANG", "MAG", "SPL", "SPR", "LTL", "LTR")
litter_labels <- c(
  "ANG" = "ANG (S. angustifolium)",
  "MAG" = "MAG (S. magellanicum)",
  "SPL" = "SPL (Spruce needles)",
  "SPR" = "SPR (Spruce fine roots)",
  "LTL" = "LTL (Labrador tea leaves)",
  "LTR" = "LTR (Labrador tea fine roots)"
)
sum_shannon <- sum_shannon %>%
  mutate(Litter = factor(Litter, levels = litter_order))

p_shannon <- ggplot(sum_shannon,
                    aes(x = t_num, y = mean_shannon,
                        color = Litter, fill = Litter, group = Litter)) +
  geom_ribbon(aes(ymin = mean_shannon - se_shannon,
                  ymax = mean_shannon + se_shannon),
              alpha = 0.50, color = NA) +
  facet_wrap(~ Class, scales = "free") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2), labels = c("0", "0.5", "1", "2")) +
  labs(x = "Collection time (years)", y = "Shannon diversity",
       #title = "Shannon +/- st. error"
       ) +
  scale_color_manual(
    values = my_colors,
    limits = litter_order, breaks = litter_order, labels = litter_labels,
    drop = FALSE, name = "Litter type"
  ) +
  scale_fill_manual(
    values = my_colors,
    limits = litter_order, breaks = litter_order, labels = litter_labels,
    drop = FALSE, name = "Litter type"
  ) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 1)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

p_shannon

ggsave("Plots/Fig7_Shannon_All_Litter.png", p_shannon, width = 10, height = 5, dpi = 300)
