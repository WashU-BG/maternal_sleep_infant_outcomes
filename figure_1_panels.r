
library(dplyr)
library(ggplot2)
library(tibble)
library(grid)


data_path <- "~/Downloads/main_panel_brain.csv"

df <- read.csv(
  data_path,
  stringsAsFactors = FALSE
)


# clean up 
df <- df %>%
  mutate(
    outcome = trimws(as.character(outcome)),
    sleep_measure = trimws(as.character(sleep_measure)),
    timepoint = trimws(as.character(timepoint)),
    
    p_value = suppressWarnings(
      as.numeric(p_value)
    ),
    
    fdr_pvalue = suppressWarnings(
      as.numeric(fdr_pvalue)
    ),
    
    estimate = suppressWarnings(
      as.numeric(estimate)
    ),
    
    conf_low = suppressWarnings(
      as.numeric(conf_low)
    ),
    
    conf_high = suppressWarnings(
      as.numeric(conf_high)
    )
  ) %>%
  
  filter(
    sleep_measure != "RLS"
  ) %>%
  
  mutate(
    #ensure that facet style is connected between sleep quality and insommnia 
    sleep_key = case_when(
      tolower(sleep_measure) == "whiirs" ~
        "insomnia",
      
      tolower(sleep_measure) %in% c(
        "sleep_qual",
        "sleep qual",
        "sleepqual"
      ) ~
        "sleep_quality",
      
      TRUE ~ NA_character_
    ),
    
    sleep_measure_label = case_when(
      sleep_key == "insomnia" ~
        "Maternal insomnia symptoms",
      
      sleep_key == "sleep_quality" ~
        "Maternal sleep quality",
      
      TRUE ~ sleep_measure
    ),
    #all sig neural pheotypes 
    outcome_label = case_when(
      outcome == "ctx_lh_lateralorbitofrontal_gi" ~
        "L lateral orbitofrontal\n(gi)",
      
      outcome == "ctx_lateralorbitofrontal_surface_area_TOTAL" ~
        "Lateral orbitofrontal\n(sa)",
      
      outcome == "lt_cbif_fa_fsl" ~
        "L cingulum bundle-inferior\n(fa)",
      
      outcome == "subcortical_gray_matter_volume" ~
        "Subcortical gray matter\n(volume)",
      
      outcome == "total_cerebral_white_matter_volume" ~
        "Total cerebral white matter\n(volume)",
      
      outcome == "total_cortical_gray_matter_volume" ~
        "Total cortical gray matter\n(volume)",
      
      outcome == "rt_fx_fa_fsl" ~
        "R fornix\n(fa)",
      
      outcome == "ctx_lh_parstriangularis_surface_area" ~
        "L pars triangularis\n(sa)",
      
      outcome == "ctx_medialorbitofrontal_surface_area_TOTAL" ~
        "Total medial orbitofrontal\n(sa)",
      
      outcome == "lt_fx_fa_fsl" ~
        "L fornix\n(fa)",
      
      outcome == "CO_vs_aFPN_bcp" ~
        "Cingulo-opercular-anterior frontoparietal\n(fc)",
      
      outcome == "DMN_vs_pFPN_bcp" ~
        "Default mode-posterior frontoparietal\n(fc)",
      
      outcome == "RST_vs_AUD_adult" ~
        "Retrosplenial-auditory\n(fc)",
      
      outcome == "Thal_vs_AUD_adult" ~
        "Thalamus-auditory\n(fc)",
      
      outcome == "ctx_insula_surface_area_TOTAL" ~
        "Total insula\n(sa)",
      
      outcome == "intracranial_volume_volume" ~
        "Intracranial volume",
      
      outcome == "supratentorial_not_ventricle_volume" ~
        "Supratentorial not ventricle volume",
      
      outcome == "total_brain_volume_volume" ~
        "Total brain volume",
      
      TRUE ~ outcome
    )
  )


#order axis based on findings and imaging modality 
order_whiirs <- c(
  "ctx_lateralorbitofrontal_surface_area_TOTAL",
  "ctx_lh_lateralorbitofrontal_gi",
  "ctx_medialorbitofrontal_surface_area_TOTAL",
  "ctx_lh_parstriangularis_surface_area",
  "ctx_insula_surface_area_TOTAL",
  "lt_fx_fa_fsl",
  "rt_fx_fa_fsl",
  "lt_cbif_fa_fsl",
  "subcortical_gray_matter_volume",
  "total_cerebral_white_matter_volume",
  "total_cortical_gray_matter_volume",
  "intracranial_volume_volume",
  "supratentorial_not_ventricle_volume",
  "total_brain_volume_volume",
  "CO_vs_aFPN_bcp",
  "RST_vs_AUD_adult",
  "Thal_vs_AUD_adult"
)

order_sleepq <- c(
  "rt_fx_fa_fsl",
  "lt_fx_fa_fsl",
  "DMN_vs_pFPN_bcp"
)


# color and limits 
trimester_colors <- c(
  "avg" = "grey55",
  "T1"  = "#1f78b4",
  "T2"  = "#d9bf6b",
  "T3"  = "#d67c1c"
)

x_limits <- c(
  -1.0,
  0.8
)

x_breaks <- seq(
  from = -1.0,
  to = 0.5,
  by = 0.5
)

offset_a <- c(
  "avg" =  0.27,
  "T1"  =  0.09,
  "T2"  = -0.09,
  "T3"  = -0.27
)

# manually set size and y axis since both parts are connected
section_heading_size <- 12 / ggplot2::.pt


make_layout <- function(
    data_sub,
    desired_order,
    spacing = 1
) {
  
  present_outcomes <- data_sub %>%
    distinct(
      outcome,
      outcome_label
    )
  
  ordered <- tibble(
    outcome = desired_order
  ) %>%
    inner_join(
      present_outcomes,
      by = "outcome"
    )
  
  extras <- present_outcomes %>%
    anti_join(
      ordered,
      by = "outcome"
    ) %>%
    arrange(
      outcome_label
    )
  
  layout <- bind_rows(
    ordered,
    extras
  )
  
  if (nrow(layout) == 0) {
    return(
      layout %>%
        mutate(
          y = numeric(0)
        )
    )
  }
  
  layout %>%
    mutate(
      y = rev(
        seq_len(n())
      ) * spacing
    )
}


#organize plotting data
df_a <- df %>%
  filter(
    timepoint %in% c(
      "avg",
      "T1",
      "T2",
      "T3"
    )
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c(
        "avg",
        "T1",
        "T2",
        "T3"
      )
    ),
    
    sig_a = case_when(
      fdr_pvalue == 1 ~
        "fdr_sig",
      
      fdr_pvalue == 2 ~
        "nominal_only",
      
      TRUE ~
        "not_sig"
    )
  )


df_a_whiirs <- df_a %>%
  filter(
    sleep_key == "insomnia"
  )


df_a_sleepq <- df_a %>%
  filter(
    sleep_key == "sleep_quality"
  )


if (nrow(df_a_whiirs) == 0) {
  stop(
    paste0(
      "No maternal insomnia rows remained after filtering. ",
      "Observed sleep_measure values: ",
      paste(
        sort(unique(df$sleep_measure)),
        collapse = ", "
      )
    )
  )
}

#handle layout for spacing- depending on number of findings this will need to be changed
make_sectioned_layout <- function(
    insomnia_data,
    sleepq_data,
    insomnia_spacing = 1.9,
    sleepq_spacing = 1.9,
    section_gap = 3.4,
    insomnia_header_gap = 1.45,
    sleepq_header_gap = 1.00
) {
  
  insomnia_layout <- make_layout(
    data_sub = insomnia_data,
    desired_order = order_whiirs,
    spacing = insomnia_spacing
  ) %>%
    mutate(
      sleep_key = "insomnia"
    )
  
  sleepq_layout <- make_layout(
    data_sub = sleepq_data,
    desired_order = order_sleepq,
    spacing = sleepq_spacing
  ) %>%
    mutate(
      sleep_key = "sleep_quality"
    )
  
  if (nrow(insomnia_layout) == 0) {
    stop(
      "No WHIIRS"
    )
  }
  
  if (nrow(sleepq_layout) == 0) {
    stop(
      "No SQ"
    )
  }
  
  sleepq_top <- max(
    sleepq_layout$y,
    na.rm = TRUE
  )
  
 #ensure that the insommnia figure is on top 
  insomnia_layout <- insomnia_layout %>%
    mutate(
      y = y +
        sleepq_top +
        section_gap
    )
  
  separator_y <- sleepq_top +
    section_gap / 2
  
  # add section headers
  section_headers <- tibble(
    label = c(
      "Maternal insomnia symptoms",
      "Maternal sleep quality"
    ),
    
    x = rep(
      mean(x_limits),
      2
    ),
    
    y = c(
      max(
        insomnia_layout$y,
        na.rm = TRUE
      ) +
        insomnia_header_gap,
      
      sleepq_top +
        sleepq_header_gap
    )
  )
  
  list(
    layout = bind_rows(
      insomnia_layout,
      sleepq_layout
    ),
    
    separator_y = separator_y,
    
    headers = section_headers
  )
}


section_info <- make_sectioned_layout(
  insomnia_data = df_a_whiirs,
  sleepq_data = df_a_sleepq,
  
  insomnia_spacing = 1.9,
  
  sleepq_spacing = 1.9,
  
  section_gap = 3.4,
  
  insomnia_header_gap = 1.45,
  
  sleepq_header_gap = 1.00
)


#connect y axis
layout_all <- section_info$layout


plot_data_before_complete_filter <- bind_rows(
  df_a_whiirs,
  df_a_sleepq
) %>%
  left_join(
    layout_all,
    by = c(
      "sleep_key",
      "outcome",
      "outcome_label"
    )
  ) %>%
  mutate(
    y_plot = y +
      unname(
        offset_a[
          as.character(timepoint)
        ]
      )
  )



#y limit
y_limits <- c(
  min(
    plot_data$y_plot,
    na.rm = TRUE
  ) -
    0.70,
  
  max(
    section_info$headers$y,
    na.rm = TRUE
  ) +
    1.15
)


#keep theme consistent with other panels 
theme_brain <- function() {
  
  theme_classic(
    base_size = 15
  ) +
    theme(
      plot.title.position = "plot",
      
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        lineheight = 1.08,
        margin = margin(
          b = 10
        )
      ),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      axis.text = element_text(
        color = "black",
        size = 13
      ),
      
      axis.text.y = element_text(
        color = "black",
        size = 12,
        hjust = 1,
        lineheight = 0.92,
        margin = margin(
          r = 7
        )
      ),
      
      axis.text.x = element_text(
        color = "black",
        size = 13
      ),
      
      axis.title.x = element_blank(),
      
      axis.title.y = element_text(
        color = "black",
        size = 15,
        margin = margin(
          r = 20
        )
      ),
      
      axis.ticks.y = element_blank(),
      
      legend.position = "bottom",
      legend.justification = "center",
      legend.box = "horizontal",
      
      legend.title = element_text(
        size = 12,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.margin = margin(
        t = 8,
        r = 4,
        b = 4,
        l = 4
      ),
      
      legend.spacing.x = unit(
        0.45,
        "cm"
      ),
      
      plot.margin = margin(
        t = 12,
        r = 12,
        b = 8,
        l = 12
      )
    )
}


#generate brain panel 
main_title <- paste0(
  "Associations of maternal sleep problems during pregnancy\n",
  "with infant neural phenotypes"
)


brain_plot_panel <- ggplot(
  data = plot_data,
  mapping = aes(
    x = estimate,
    y = y_plot,
    color = timepoint
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black",
    show.legend = FALSE
  ) +
  
  geom_hline(
    yintercept = section_info$separator_y,
    linewidth = 0.8,
    color = "grey50",
    show.legend = FALSE
  ) +
  
  geom_segment(
    aes(
      x = conf_low,
      xend = conf_high,
      yend = y_plot
    ),
    linewidth = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  

  geom_point(
    shape = 21,
    fill = "white",
    size = 3.1,
    stroke = 1.0,
    show.legend = TRUE
  ) +
  
  geom_point(
    data = function(x) {
      dplyr::filter(
        x,
        sig_a == "fdr_sig"
      )
    },
    
    aes(
      fill = timepoint
    ),
    
    shape = 21,
    size = 3.1,
    stroke = 1.0,
    show.legend = FALSE
  ) +
  
  geom_point(
    data = function(x) {
      dplyr::filter(
        x,
        sig_a == "nominal_only"
      )
    },
    
    shape = 16,
    size = 1.2,
    stroke = 0,
    show.legend = FALSE
  ) +
  

  geom_label(
    data = section_info$headers,
    
    mapping = aes(
      x = x,
      y = y,
      label = label
    ),
    
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = section_heading_size,
    color = "black",
    fill = "white",
    linewidth = 0,
    
    label.padding = unit(
      0.16,
      "lines"
    ),
    
    label.r = unit(
      0,
      "lines"
    ),
    
    show.legend = FALSE
  ) +
  
  scale_y_continuous(
    breaks = axis_labels$y,
    labels = axis_labels$outcome_label,
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  scale_x_continuous(
    breaks = x_breaks,
    labels = sprintf(
      "%.1f",
      x_breaks
    )
  ) +
  
  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    clip = "off"
  ) +
  
  scale_color_manual(
    values = trimester_colors,
    
    limits = c(
      "avg",
      "T1",
      "T2",
      "T3"
    ),
    
    breaks = c(
      "avg",
      "T1",
      "T2",
      "T3"
    ),
    
    labels = c(
      "Average",
      "T1",
      "T2",
      "T3"
    ),
    
    name = "Trimester",
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = trimester_colors,
    
    limits = c(
      "avg",
      "T1",
      "T2",
      "T3"
    ),
    
    drop = FALSE,
    guide = "none"
  ) +
  
  labs(
    title = main_title,
    
    x = NULL,
    
    y = "Infant imaging"
  ) +
  
  theme_brain() +
  
  guides(
    color = guide_legend(
      title = "Trimester",
      order = 1,
      nrow = 1,
      byrow = TRUE,
      
      override.aes = list(
        shape = 21,
        size = 3.2,
        stroke = 1.0,
        
        fill = unname(
          trimester_colors[
            c(
              "avg",
              "T1",
              "T2",
              "T3"
            )
          ]
        )
      )
    )
  ) # make sure this is not removed 

brain_plot_panel

#----------------------------------------------------------

# behavior - sanity check and re-load packages- of note, need to manually set font unlike brain figure to match panels
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)

font_multiplier <- 1.5
figure_width <- 17.5
figure_height <- 6.5

df <- read.csv(
  "~/Downloads/main_panel_results_itsea.csv",
  stringsAsFactors = FALSE
)

df <- df %>%
  mutate(
    outcome = trimws(outcome),
    sleep_measure = trimws(sleep_measure),
    timepoint = trimws(timepoint),
    p_value = as.numeric(p_value),
    fdr_pvalue = as.numeric(fdr_pvalue)
  ) %>%
  mutate(
    outcome_label = case_when(
      outcome == "itsea_external_tscore" ~ "Externalizing",
      outcome == "itsea_internal_tscore" ~ "Internalizing",
      outcome == "itsea_dysreg_tscore" ~ "Dysregulation",
      outcome == "itsea_competence_tscore" ~ "Competence",
      TRUE ~ outcome
    ),
    
    sleep_measure = factor(
      sleep_measure,
      levels = c("WHIIRS", "RLS", "Sleep_qual"),
      labels = c(
        "Maternal\ninsomnia symptoms",
        "Maternal\nrestless legs\nsymptoms",
        "Maternal sleep\nquality"
      )
    ),
    
    sig = if_else(p_value < 0.05, "sig", "ns")
  )

outcome_order <- c(
  "Externalizing",
  "Internalizing",
  "Dysregulation",
  "Competence"
)

# keep colors identical to brain figure 
trimester_colors <- c(
  "avg" = "grey55",
  "T1"  = "#1f78b4",
  "T2"  = "#d9bf6b",
  "T3"  = "#d67c1c"
)

#x axis limits and breaks
x_limits <- c(-1.0, 0.8)
x_breaks <- seq(-1.0, 0.5, by = 0.5)

# fdr_pvalue:
# 1 = filled circle
# 2 = open circle + center dot-> nominal
# 0/other = open circle

df_a <- df %>%
  filter(timepoint %in% c("avg", "T1", "T2", "T3")) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("avg", "T1", "T2", "T3")
    ),
    
    outcome_label = factor(
      outcome_label,
      levels = outcome_order
    ),
    
    sig_a = case_when(
      fdr_pvalue == 1 ~ "fdr_sig",
      fdr_pvalue == 2 ~ "nominal_only",
      TRUE ~ "not_sig"
    )
  )

# facet wrap the sleep outcomes
facet_key_a <- df_a %>%
  distinct(sleep_measure) %>%
  arrange(sleep_measure) %>%
  mutate(facet_id = row_number())

y_key_a <- df_a %>%
  distinct(sleep_measure, outcome_label) %>%
  left_join(facet_key_a, by = "sleep_measure") %>%
  arrange(sleep_measure, outcome_label) %>%
  group_by(sleep_measure) %>%
  mutate(
    row_in_panel = dense_rank(outcome_label),
    row_rev = max(row_in_panel) + 1 - row_in_panel,
    y_pos = facet_id * 10 + row_rev * 1.8
  ) %>%
  ungroup()

df_a <- df_a %>%
  left_join(
    y_key_a,
    by = c("sleep_measure", "outcome_label")
  )

# Original trimester offsets
offset_a <- c(
  "avg" =  0.27,
  "T1"  =  0.09,
  "T2"  = -0.09,
  "T3"  = -0.27
)

df_a <- df_a %>%
  mutate(
    y_plot = y_pos +
      unname(offset_a[as.character(timepoint)])
  )

#
plot_behavior <- ggplot(
  df_a,
  aes(
    x = estimate,
    y = y_plot,
    color = timepoint
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  
  geom_segment(
    aes(
      x = conf_low,
      xend = conf_high,
      yend = y_plot
    ),
    linewidth = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  
  geom_point(
    shape = 21,
    fill = "white",
    size = 3.1,
    stroke = 1.0
  ) +
  
  # fdr sig
  geom_point(
    data = ~ dplyr::filter(.x, sig_a == "fdr_sig"),
    aes(fill = timepoint),
    shape = 21,
    size = 3.1,
    stroke = 1.0,
    show.legend = FALSE
  ) +
  
  # Nom sig
  geom_point(
    data = ~ dplyr::filter(.x, sig_a == "nominal_only"),
    shape = 16,
    size = 1.25,
    stroke = 0,
    show.legend = FALSE
  ) +
  
  # facet
  facet_wrap(
    ~ sleep_measure,
    ncol = 3,
    scales = "free_y"
  ) +
  
  scale_y_continuous(
    breaks = y_key_a$y_pos,
    labels = as.character(y_key_a$outcome_label),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  
  
  scale_x_continuous(
    limits = x_limits,
    breaks = x_breaks,
    labels = sprintf("%.1f", x_breaks)
  ) +
  
  scale_color_manual(
    values = trimester_colors,
    breaks = c("avg", "T1", "T2", "T3"),
    labels = c(
      "avg" = "Average",
      "T1" = "T1",
      "T2" = "T2",
      "T3" = "T3"
    ),
    name = "Trimester"
  ) +
  
  scale_fill_manual(
    values = trimester_colors,
    guide = "none"
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = "Infant behavior"
  ) +
  
  theme_classic(
    base_size = 15 * font_multiplier
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16 * font_multiplier,
      hjust = 0
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12 * font_multiplier,
      lineheight = 0.95,
      margin = margin(t = 5, b = 5)
    ),
    
    strip.background = element_rect(
      color = "black",
      fill = "white",
      linewidth = 0.8
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.text = element_text(
      color = "black",
      size = 13 * font_multiplier
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 12 * font_multiplier
    ),
    
    axis.title.x = element_text(
      size = 15 * font_multiplier
    ),
    
    axis.title.y = element_text(
      size = 15 * font_multiplier,
      margin = margin(r = 20)
    ),
    
    axis.ticks.y = element_blank(),
    
    panel.spacing.x = unit(0.8, "lines"),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 12 * font_multiplier,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12 * font_multiplier
    ),
    
    plot.margin = margin(5, 10, 2, 10)
  ) +
  
  guides(
    color = guide_legend(
      order = 1,
      nrow = 1,
      override.aes = list(
        shape = 21,
        size = 3.2,
        stroke = 1.0,
        fill = unname(
          trimester_colors[
            c("avg", "T1", "T2", "T3")
          ]
        )
      )
    )
  )

main_title <- paste(
  "Associations of maternal sleep problems during pregnancy\n",
  "with hypothesized infant behavioral phenotypes"
)

print(plot_behavior)
#--------------------------------------------------------------------
#bayley figure
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(NatParksPalettes)

df <- read.csv("~/Downloads/main_panel_results_bayleyv2.csv", stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    outcome = trimws(outcome),
    sleep_measure = trimws(sleep_measure),
    timepoint = trimws(timepoint),
    p_value = suppressWarnings(as.numeric(p_value))
  ) %>%
  mutate(
    outcome_label = case_when(
      outcome == "bayley_cog_comp" ~ "Cognition composite",
      TRUE ~ outcome
    ),
    sleep_measure = case_when(
      sleep_measure == "Sleep_qual" ~ "Maternal Sleep Quality",
      TRUE ~ sleep_measure
    ),
    sig = if_else(!is.na(p_value) & p_value < 0.05, "sig", "ns")
  )

# keep only outcome  that  exist in the csv
desired_outcome_order <- c(
  "Cognition composite",
  "Externalizing",
  "Internalizing",
  "Dysregulation",
  "Competence"
)

outcome_order <- intersect(desired_outcome_order, unique(df$outcome_label))
if (length(outcome_order) == 0) {
  outcome_order <- unique(df$outcome_label)
}

sleep_order <- unique(df$sleep_measure)

yellowstone_cols <- natparks.pals("Yellowstone", 4)

trimester_colors <- c(
  "avg" = "grey55",
  "T1"  = yellowstone_cols[1],
  "T2"  = yellowstone_cols[2],
  "T3"  = yellowstone_cols[3]
)

x_limits <- c(-1.0, 0.8)
x_breaks <- seq(-1.0, 0.5, by = 0.5)

# example panel a
df_a <- df %>%
  filter(timepoint %in% c("avg", "T1", "T2", "T3")) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("avg", "T1", "T2", "T3")),
    outcome_label = factor(outcome_label, levels = outcome_order),
    sleep_measure = factor(sleep_measure, levels = sleep_order),
    fill_a = if_else(sig == "sig", as.character(timepoint), "ns")
  ) %>%
  filter(!is.na(outcome_label), !is.na(sleep_measure))

facet_key_a <- df_a %>%
  distinct(sleep_measure) %>%
  arrange(sleep_measure) %>%
  mutate(facet_id = row_number())

y_key_a <- df_a %>%
  distinct(sleep_measure, outcome_label) %>%
  left_join(facet_key_a, by = "sleep_measure") %>%
  arrange(sleep_measure, outcome_label) %>%
  group_by(sleep_measure) %>%
  mutate(
    row_in_panel = dense_rank(outcome_label),
    row_rev = max(row_in_panel) + 1 - row_in_panel,
    y_pos = row_rev * 1.8
  ) %>%
  ungroup()

df_a <- df_a %>%
  left_join(y_key_a, by = c("sleep_measure", "outcome_label"))

offset_a <- c(
  "avg" =  0.27,
  "T1"  =  0.09,
  "T2"  = -0.09,
  "T3"  = -0.27
)

df_a <- df_a %>%
  mutate(
    y_plot = y_pos + unname(offset_a[as.character(timepoint)])
  )

p_development <- ggplot(
  df_a,
  aes(x = estimate, y = y_plot, color = timepoint)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  geom_segment(
    aes(x = conf_low, xend = conf_high, yend = y_plot),
    linewidth = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = fill_a),
    shape = 21,
    size = 3.1,
    stroke = 1.0
  ) +
  facet_wrap(~ sleep_measure, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    breaks = y_key_a$y_pos,
    labels = as.character(y_key_a$outcome_label),
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  scale_x_continuous(
    limits = x_limits,
    breaks = x_breaks
  ) +
  scale_color_manual(
    values = trimester_colors,
    breaks = c("avg", "T1", "T2", "T3"),
    labels = c(
      "avg" = "Average",
      "T1"  = "T1",
      "T2"  = "T2",
      "T3"  = "T3"
    ),
    name = "Trimester"
  ) +
  scale_fill_manual(
    values = c(trimester_colors, "ns" = "white"),
    guide = "none"
  ) +
  labs(
    title = "Trimester-specific associations",
    x = NULL,
    y = "Infant Developmental Domain"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(color = "black", size = 11),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14, margin = margin(r = 20)),
    axis.ticks.y = element_blank(),
    panel.spacing.x = unit(0.6, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(5, 10, 2, 10)
  ) +
  guides(
    color = guide_legend(
      order = 1,
      nrow = 1,
      override.aes = list(
        shape = 21,
        size = 3.2,
        stroke = 1.0,
        fill = unname(trimester_colors[c("avg", "T1", "T2", "T3")])
      )
    )
  )

print(p_development)




