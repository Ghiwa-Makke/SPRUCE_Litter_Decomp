### SPRUCE Litter decomposition 2015 - 2017
### Figure 2 - Overall metabolic varaition and trajectory analysis
## Author: Ghiwa Makke

### SPRUCE Litter decomposition 2015 - 2017
### Figure 2 - Overall metabolic variation and trajectory analysis

## Load packages & settings
library(tidyverse)
library(vegan)
library(ecotraj)
library(patchwork)
library(ragg)

## Custom colors 
my_colors <- c(
  "LTL" = "#e6ab02",
  "LTR" = "#d95f02",
  "ANG" = "#7570b3",
  "MAG" = "#e7298a",
  "SPL" = "#96C291",
  "SPR" = "#1b9e77"
)

## Facet / axis labels
custom_labels <- c(
  "MAG" = "italic('S. magellanicum')",
  "ANG" = "italic('S. angustifolium')",
  "SPR" = "'Spruce roots'",
  "LTR" = "'Labrador tea roots'",
  "LTL" = "'Labrador tea leaves'",
  "SPL" = "'Spruce needles'"
)

axis_labels_expr <- c(
  "MAG" = expression(italic("S. magellanicum")),
  "ANG" = expression(italic("S. angustifolium")),
  "SPR" = expression("Spruce roots"),
  "LTR" = expression("Labrador tea roots"),
  "LTL" = expression("Labrador tea leaves"),
  "SPL" = expression("Spruce needles")
)

## Legend labels for litter (Figure 2 legends)
litter_labels_legend <- c(
  ANG = expression(italic("S. angustifolium")),
  MAG = expression(italic("S. magellanicum")),
  SPR = expression("Spruce roots"),
  LTR = expression("Labrador tea roots"),
  LTL = expression("Labrador tea leaves"),
  SPL = expression("Spruce needles")
)

## 1. Read & align data -----
mat_input <- readr::read_csv(
  "Input/Metabodirect_All_Sum_Time/1_preprocessing_output/matrix_features.csv",
  col_types = cols()
)
meta <- readr::read_csv("Input/metadata.csv", col_types = cols()) %>%
  mutate(SampleID = as.character(SampleID))

# Move Mass to rownames, then transpose so rows = samples, cols = features
mat <- mat_input %>%
  column_to_rownames("Mass") %>%
  t()

## Time variable & factors 
time_map <- c("T_0" = 0, "T_0.5" = 0.5, "T_1" = 1, "T_2" = 2)

meta <- meta %>%
  mutate(
    Pickup_t  = factor(Pickup_t, levels = names(time_map)),
    Time_yrs  = unname(time_map[as.character(Pickup_t)]),
    Litter    = factor(Litter)
  )

## 2. Distance + PCoA on samples -----
bc_dist <- vegdist(mat, method = "bray")

pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE, add = TRUE)

pcoa_scores <- as.data.frame(pcoa_res$points[, 1:2])
colnames(pcoa_scores) <- c("Axis1", "Axis2")

var_explained <- pcoa_res$eig[1:2] / sum(pcoa_res$eig) * 100

pcoa_scores$SampleID <- rownames(mat)
pcoa_scores <- pcoa_scores %>%
  left_join(meta, by = "SampleID")

## 3. Litter x time centroids -------
pcoa_centroids <- pcoa_scores %>%
  group_by(Litter, Time_yrs) %>%
  summarise(
    Axis1 = mean(Axis1, na.rm = TRUE),
    Axis2 = mean(Axis2, na.rm = TRUE),
    n     = n(),
    .groups = "drop"
  )

# Common time grid
time_grid <- c(0, 0.5, 1, 2)

## 5. Interpolate centroids (per litter) ------
interpolate_litter_centroids <- function(df_litter, time_grid) {
  times_obs <- df_litter$Time_yrs
  x_obs     <- df_litter$Axis1
  y_obs     <- df_litter$Axis2
  
  ord <- order(times_obs)
  times_obs <- times_obs[ord]
  x_obs     <- x_obs[ord]
  y_obs     <- y_obs[ord]
  
  x_interp <- approx(x = times_obs, y = x_obs, xout = time_grid, rule = 2)$y
  y_interp <- approx(x = times_obs, y = y_obs, xout = time_grid, rule = 2)$y
  
  tibble(
    Litter   = df_litter$Litter[1],
    Time_yrs = time_grid,
    Axis1    = x_interp,
    Axis2    = y_interp,
    observed = time_grid %in% times_obs
  )
}

pcoa_centroids_interp <- pcoa_centroids %>%
  group_by(Litter) %>%
  group_modify(~ interpolate_litter_centroids(.x, time_grid)) %>%
  ungroup()

## 6. Build ecotraj inputs from centroids ---------
xy_cent     <- as.matrix(pcoa_centroids_interp %>% select(Axis1, Axis2))
sites_cent  <- as.character(pcoa_centroids_interp$Litter)
survey_lvls <- sort(unique(time_grid))

surveys_cent <- as.integer(factor(pcoa_centroids_interp$Time_yrs,
                                  levels = survey_lvls))
times_cent   <- pcoa_centroids_interp$Time_yrs

## 7. ecotraj trajectory metrics ----------
lengths_cent <- trajectoryLengths2D(
  xy      = xy_cent,
  sites   = sites_cent,
  surveys = surveys_cent,
  relativeToInitial = FALSE,
  all     = FALSE
)

lengths_cent_df <- as.data.frame(lengths_cent) %>%
  rownames_to_column("Litter") %>%
  rename(total_length = Path)

speeds_cent <- trajectorySpeeds2D(
  xy      = xy_cent,
  sites   = sites_cent,
  surveys = surveys_cent,
  times   = times_cent
)

speeds_cent_df <- as.data.frame(speeds_cent) %>%
  rownames_to_column("Litter") %>%
  rename(path_speed = Path)

dist_cent <- dist(xy_cent)

traj_cent <- defineTrajectories(
  d     = dist_cent,
  sites = sites_cent,
  times = times_cent
)

stopifnot(is.synchronous(traj_cent))

traj_metrics_summary <- lengths_cent_df %>%
  left_join(speeds_cent_df, by = "Litter")

## 8. Convergence diagnostics 
pairwise_cent <- pcoa_centroids_interp %>%
  inner_join(pcoa_centroids_interp,
             by = "Time_yrs",
             suffix = c("_1", "_2")) %>%
  filter(as.character(Litter_1) < as.character(Litter_2)) %>%
  mutate(
    distance = sqrt((Axis1_1 - Axis1_2)^2 + (Axis2_1 - Axis2_2)^2)
  )

mean_pairwise_cent <- pairwise_cent %>%
  group_by(Time_yrs) %>%
  summarise(
    mean_distance = mean(distance),
    sd_distance   = sd(distance),
    .groups       = "drop"
  )

mean_pairwise_cent_no0 <- mean_pairwise_cent %>%
  filter(Time_yrs > 0)

grand_centroid_interp <- pcoa_centroids_interp %>%
  group_by(Time_yrs) %>%
  summarise(
    GC_Axis1 = mean(Axis1),
    GC_Axis2 = mean(Axis2),
    .groups  = "drop"
  )

centroids_with_gc <- pcoa_centroids_interp %>%
  left_join(grand_centroid_interp, by = "Time_yrs") %>%
  mutate(
    dist_from_gc = sqrt((Axis1 - GC_Axis1)^2 + (Axis2 - GC_Axis2)^2)
  )

segments_conv <- centroids_with_gc %>%
  group_by(Litter) %>%
  arrange(Time_yrs, .by_group = TRUE) %>%
  mutate(
    Axis1_next = lead(Axis1),
    Axis2_next = lead(Axis2),
    Time_next  = lead(Time_yrs),
    dist_next  = lead(dist_from_gc),
    delta_dist = dist_next - dist_from_gc,
    pattern    = case_when(
      is.na(delta_dist)          ~ NA_character_,
      delta_dist < -1e-6         ~ "Converging",
      delta_dist >  1e-6         ~ "Diverging",
      TRUE                       ~ "Stable"
    ),
    magnitude  = abs(delta_dist)
  ) %>%
  filter(!is.na(Axis1_next)) %>%
  ungroup()

max_mag <- max(segments_conv$magnitude, na.rm = TRUE)
segments_conv <- segments_conv %>%
  mutate(
    lw_scaled = ifelse(
      max_mag > 0,
      0.4 + 2.0 * (magnitude / max_mag),
      0.8
    )
  )

# NMDS --------------------------------------------------------------------

set.seed(123)
nmds_res <- metaMDS(mat, distance = "bray", k = 2, trymax = 100, trace = FALSE)

nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
colnames(nmds_scores) <- c("NMDS1", "NMDS2")
nmds_scores$SampleID <- rownames(mat)

nmds_scores <- nmds_scores %>%
  left_join(meta, by = "SampleID")

time_labels <- c(
  "T_0"   = "0 year",
  "T_0.5" = "0.5 year",
  "T_1"   = "1 year",
  "T_2"   = "2 years"
)

## Panel A: NMDS (color = Litter, shape = Pickup_t)
fig2A_NMDS <- ggplot(
  nmds_scores,
  aes(x = NMDS1, y = NMDS2,
      color = Litter,
      shape = Pickup_t)
) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(
    values = my_colors,
    name   = "Litter type",
    labels = litter_labels_legend
  ) +
  scale_shape_discrete(name = "Collection time",
                       labels = time_labels) +
  labs(
    title = "A) NMDS of litter metabolomes",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title      = element_text(hjust = 0, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.title        = element_text(margin = margin(b = 4)),
    legend.key.spacing.y = unit(1, "mm"),
    legend.key.height    = unit(4, "mm"),
    panel.grid.minor = element_blank()
  )
fig2A_NMDS
## Panel B: PCoA centroid trajectories w/ convergence coloring
fig2B_PCoA_conv <- ggplot() +
  geom_point(data = pcoa_scores,
             aes(x = Axis1, y = Axis2),
             color = "gray90", size = 0.8, alpha = 0.4) +
  geom_point(data = centroids_with_gc,
             aes(x = Axis1, y = Axis2, fill = Litter),
             shape = 21, size = 3, color = "black", stroke = 0.4) +
  geom_text(data = centroids_with_gc,
            aes(x = Axis1, y = Axis2, label = Time_yrs),
            size = 2.3, fontface = "bold", vjust = -0.9) +
  geom_segment(data = segments_conv,
               aes(x = Axis1, y = Axis2,
                   xend = Axis1_next, yend = Axis2_next,
                   color = pattern,
                   linewidth = lw_scaled),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               alpha = 0.9) +
  scale_fill_manual(
    values = my_colors,
    name   = "Litter type",
    labels = litter_labels_legend
  ) +
  scale_color_manual(
    values = c(
      "Converging" = "#1f78b4",
      "Diverging"  = "#e31a1c",
      "Stable"     = "grey40"
    ),
    name = "Segment pattern"
  ) +
  scale_linewidth(range = c(0.4, 1), guide = "none") +
  labs(
    title = "B) Litter molecular trajectories in PCoA space",
    x     = sprintf("PCoA Axis 1 [%.1f%%]", var_explained[1]),
    y     = sprintf("PCoA Axis 2 [%.1f%%]", var_explained[2])
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title      = element_text(hjust = 0, face = "bold"),
    legend.position = "right",
    legend.text     = element_text(size = 7),
    panel.grid      = element_blank(),
    legend.title        = element_text(margin = margin(b = 4)),
    legend.key.spacing.y = unit(1, "mm"),
    legend.key.height    = unit(4, "mm"),
  )

## Panel C: Mean pairwise distance among centroids
fig2C_mean_pair <- ggplot(mean_pairwise_cent_no0,
                          aes(x = Time_yrs, y = mean_distance)) +
  geom_ribbon(aes(ymin = mean_distance - sd_distance,
                  ymax = mean_distance + sd_distance),
              alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(0.5, 1, 2)) +
  labs(
    title    = "C) Mean pairwise distance among litter types",
    subtitle = "Declining distance indicates global convergence",
    x        = "Time (years)",
    y        = "Mean pairwise distance (PCoA space)"
  ) +
  theme_bw(base_size = 7) +
  theme(
    plot.title    = element_text(hjust = 0.3, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    panel.grid.minor = element_blank()
  )

## Panel D: Total trajectory length by litter
fig2D_length <- ggplot(traj_metrics_summary,
                       aes(x = reorder(Litter, total_length),
                           y = total_length,
                           fill = Litter)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  scale_fill_manual(values = my_colors, guide = "none") +
  coord_flip() +
  labs(
    title = "D) Total trajectory length by litter type",
    x     = "Litter type",
    y     = "Cumulative distance in PCoA space"
  ) +
  theme_bw(base_size = 7) +
  theme(
    plot.title         = element_text(hjust = 0, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 7)
  ) +
  scale_x_discrete(
    labels = axis_labels_expr[levels(reorder(traj_metrics_summary$Litter,
                                             traj_metrics_summary$total_length))]
  )

## Panel E: Mean trajectory speed by litter
fig2E_speed <- ggplot(traj_metrics_summary,
                      aes(x = reorder(Litter, path_speed),
                          y = path_speed,
                          fill = Litter)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  scale_fill_manual(values = my_colors, guide = "none") +
  coord_flip() +
  labs(
    title = "E) Mean trajectory speed by litter type",
    x     = "Litter type",
    y     = "Speed (distance per year)"
  ) +
  theme_bw(base_size = 7) +
  theme(
    plot.title         = element_text(hjust = 0, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 7),
    legend.position    = "none"
  ) +
  scale_x_discrete(
    labels = axis_labels_expr[levels(reorder(traj_metrics_summary$Litter,
                                             traj_metrics_summary$path_speed))]
  )

## final multi-panel figure -----------
row1 <- fig2A_NMDS | fig2B_PCoA_conv
row2 <- fig2C_mean_pair + fig2D_length + fig2E_speed

Figure2 <- row1 / row2 + plot_layout(heights = c(1.2, 1))

Figure2_ <- fig2A_NMDS / fig2B_PCoA_conv / row2 + plot_layout(heights = c(1.2, 1.2, 1))

ggsave(
  "Plots/Updated/Figure2.tiff",
  Figure2,
  dpi         = 300,
  width       = 190,
  height      = 145,
  units       = "mm",
  device      = agg_tiff,
  compression = "lzw"
)

