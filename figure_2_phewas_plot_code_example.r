
#example PheWAS figure plot adapted from https://www.nature.com/articles/s44220-024-00281-7
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(stringr)
#upload organized CSV
df <- read.csv("~/Downloads/PheWAS_plot_sleep_WHIIRS.csv")
df$pval <- as.numeric(df$pval)

timepoint_order <- c("Average", "Trimester1", "Trimester2", "Trimester3")

trimester_colors <- c(
  "Average"    = "grey55",
  "Trimester1" = "#1f78b4",
  "Trimester2" = "#d9bf6b",
  "Trimester3" = "#d67c1c"
)

timepoint_labels <- c(
  "Average"    = "Average",
  "Trimester1" = "T1",
  "Trimester2" = "T2",
  "Trimester3" = "T3"
)

# ensure that all imaging modalities are captured
modality_levels <- c(
  "GLOBAL",
  "SUB CORTICAL VOLUME",
  "CORTICAL VOLUME",
  "SURFACE AREA",
  "GYRIFICATION",
  "MEAN DIFFUSIVITY",
  "FRACTIONAL ANISOTROPY",
  "RADIAL DIFFUSIVITY",
  "L1",
  "RSFC ADULT ATLAS",
  "RSFC CHILD ATLAS"
)

# ordering and labels for Y axis 
modality_order <- c(
  "Global",
  "Subcortical volume",
  "Cortical volume",
  "Surface area",
  "Gyrification",
  "Mean diffusivity",
  "Fractional anisotropy",
  "Radial diffusivity",
  "Axial diffusivity",
  "Rsfc adult atlas",
  "Rsfc child atlas"
)

df <- df %>%
  mutate(
    neg_log10_p = -log10(pval),
    timepoint = str_squish(timepoint),
    timepoint = case_when(
      timepoint %in% c("average", "avg", "Average", "AVG") ~ "Average",
      TRUE ~ timepoint
    ),
    timepoint = factor(timepoint, levels = timepoint_order),
    
    modality_clean = str_to_upper(str_squish(modality)),
    modality_clean = case_when(
      modality_clean == "SUBCORTICAL VOLUME" ~ "SUB CORTICAL VOLUME",
      TRUE ~ modality_clean
    ),
    modality_clean = factor(modality_clean, levels = modality_levels),
    
    hypothesis = as.integer(hypothesis),
    exploratory = as.integer(exploratory),
    label_point = ifelse(hypothesis == 1 | exploratory == 1, descriptive, NA)
  )

modality_gap <- 1.55
timepoint_gap <- 0.18
n_mod <- length(modality_levels)

df <- df %>%
  mutate(
    modality_num = as.numeric(modality_clean),
    timepoint_num = as.numeric(timepoint),
    y_base = (n_mod - modality_num + 1) * modality_gap,
    y = y_base + (timepoint_num - 2.5) * timepoint_gap
  )

separator_df <- data.frame(
  yintercept = ((n_mod - 1):1 + 0.5) * modality_gap
)

#if needed for points that are significant either hypothesis or exploratory then label points manually  
sig_df <- df %>%
  filter(hypothesis == 1 | exploratory == 1) %>%
  mutate(
    label_nudge_y = case_when(
      str_detect(descriptive, "total cortical gray matter") ~ -0.28,
      str_detect(descriptive, "subcortical gray matter") ~ 0.18,
      str_detect(descriptive, "total-medial orbitofrontal") ~ 0.20,
      
      modality_clean == "SURFACE AREA" & str_detect(descriptive, "medial") ~ 0.30,
      modality_clean == "SURFACE AREA" & str_detect(descriptive, "pars triangularis") ~ 0.22,
      modality_clean == "SURFACE AREA" & str_detect(descriptive, "insula") ~ 0.18,
      modality_clean == "SURFACE AREA" & str_detect(descriptive, "lateral") ~ -0.38,
      modality_clean == "RSFC CHILD ATLAS" ~ 0.20,
      TRUE ~ 0
    )
  )

exploratory_df <- df %>%
  filter(exploratory == 1 & hypothesis != 1)

hypothesis_df <- df %>%
  filter(hypothesis == 1)

plot_example <- ggplot(df, aes(x = neg_log10_p, y = y)) +
  geom_hline(
    data = separator_df,
    aes(yintercept = yintercept),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "grey80"
  ) +
  geom_point(
    aes(fill = timepoint),
    shape = 21,
    size = 3,
    alpha = 0.65,
    color = "white",
    stroke = 0.25
  ) +
  geom_segment(
    data = exploratory_df,
    aes(
      x = neg_log10_p - 0.015,
      xend = neg_log10_p + 0.015,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.6,
    lineend = "round"
  ) +
  geom_point(
    data = hypothesis_df,
    aes(fill = timepoint),
    shape = 21,
    size = 3.5,
    alpha = 1,
    color = "black",
    stroke = 1.1
  ) +
  geom_text_repel(
    data = sig_df,
    aes(
      label = label_point,
      y = y + label_nudge_y
    ),
    size = 3.8,
    force = 6,
    force_pull = 0.08,
    box.padding = 0.45,
    point.padding = 0.7,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.size = 0.25,
    segment.alpha = 0.6,
    direction = "both",
    xlim = c(-0.30, 5.30),
    ylim = c(min(df$y, na.rm = TRUE), max(df$y, na.rm = TRUE))
  ) +
  scale_fill_manual(
    values = trimester_colors,
    breaks = timepoint_order,
    labels = timepoint_labels,
    name = "Timepoint"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        colour = "black",
        stroke = 0.6,
        alpha = 1
      )
    )
  ) +
  scale_y_continuous(
    breaks = ((n_mod:1) * modality_gap),
    labels = modality_order,
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  scale_x_continuous(
    breaks = 0:5,
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(
    xlim = c(-0.35, 5.35),
    clip = "off"
  ) +
  labs(
    title = "Associations of maternal insomnia symptoms with brain phenotypes",
    x = expression(-log[10](italic(p))),
    y = "Imaging modality"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 40, 10, 10)
  )

plot_example


