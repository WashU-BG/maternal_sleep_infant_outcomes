
#figure 3
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(grid)
library(ragg)

# set the figure dimensions in advance
figure_width <- 28
figure_height <- 20
pathway_font_size <- 34

#load csv
med <- read.csv(
  "~/Downloads/mediation_resultsv2.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# organize and clean cols
med <- med %>%
  mutate(
    predictor = trimws(predictor),
    mediator = trimws(mediator),
    outcome = trimws(outcome),
    parameter = trimws(parameter),
    
    est = as.numeric(est),
    se = as.numeric(se),
    z = as.numeric(z),
    pvalue = as.numeric(pvalue),
    ci_lower = as.numeric(`ci lower`),
    ci_upper = as.numeric(`ci upper`),
    std.all = as.numeric(std.all),
    n_obs = as.numeric(n_obs)
  )

# drop all but indrect
indirect_df <- med %>%
  filter(parameter == "indirect") %>%
  mutate(
    
    predictor_label = case_when(
      predictor == "whiirs_avg" ~ "Insomnia avg",
      predictor == "whiirs_t1" ~ "Insomnia T1",
      predictor == "rls_avg" ~ "RLS avg",
      predictor == "rls_t1" ~ "RLS T1",
      predictor == "rls_t3" ~ "RLS T3",
      predictor == "sleepqual_t1" ~ "Sleep quality T1",
      predictor == "sleepqual_t3" ~ "Sleep quality T3",
      TRUE ~ predictor
    ),
    
    sleep_group = case_when(
      str_detect(predictor, "whiirs") ~
        "Maternal insomnia\nsymptoms",
      
      str_detect(predictor, "rls") ~
        "Maternal restless legs\nsymptoms",
      
      str_detect(predictor, "sleepqual") ~
        "Maternal sleep\nquality",
      
      TRUE ~ "Other"
    ),
    
    timepoint = case_when(
      str_detect(predictor, "avg") ~ "avg",
      str_detect(predictor, "t1") ~ "T1",
      str_detect(predictor, "t2") ~ "T2",
      str_detect(predictor, "t3") ~ "T3",
      TRUE ~ "avg"
    ),
    
    outcome_label = case_when(
      outcome == "itsea_competence_tscore" ~
        "Competence",
      
      outcome == "itsea_dysreg_tscore" ~
        "Dysregulation",
      
      outcome == "itsea_internal_tscore" ~
        "Internalizing",
      
      outcome == "itsea_external_tscore" ~
        "Externalizing",
      
      outcome == "bayley_cog_comp" ~
        "Cognition",
      
      TRUE ~ outcome
    ),
    
    mediator_label = case_when(
      
      mediator == "ctx_lh_temporalpole_gi" ~
        "L temporal pole GI",
      
      mediator == "Usp_vs_SAL_adult" ~
        "FC between unimodal sensory & salience network (A)",
      
      mediator == "ctx_rh_entorhinal_surface_area" ~
        "R entorhinal surface area",
      
      mediator == "DMN_vs_AUD_adult" ~
        "FC between default mode & auditory network (A)",
      
      mediator == "Usp_vs_FPN_adult" ~
        "FC between unimodal sensory and frontoparietal control network (A)",
      
      mediator == "ctx_rh_parsopercularis_gi" ~
        "R parsopercularis GI",
      
      mediator == "ctx_rh_entorhinal_volume" ~
        "R entorhinal volume",
      
      mediator == "FPN_vs_AUD_adult" ~
        "FC between frontoparietal control and auditory networks (A)",
      
      mediator == "Amyg_vs_Caud_bcp" ~
        "FC between amygdala and caudate nucleus (C)",
      
      mediator == "Palli_vs_USp_bcp" ~
        "FC between globus pallidus and unimodal sensory network (C)",
      
      mediator == "lt_fx_fa_fsl" ~
        "L fornix FA",
      
      mediator == "ifof_fa_fsl_TOTAL" ~
        "total-inferior fronto-occipital fasciculus FA",
      
      mediator == "ctx_lh_frontalpole_volume" ~
        "L frontal pole volume",
      
      mediator == "ctx_rh_frontalpole_volume" ~
        "R frontal pole volume",
      
      mediator == "DMN_vs_SAL_adult" ~
        "FC between default mode and salience networks (A)",
      
      mediator == "DAN_vs_SMH_adult" ~
        "FC between dorsal attention and somatomotor-hand networks (A)",
      
      mediator == "DAN_vs_Temporal_bcp" ~
        "FC between dorsal attention and temporal association networks (C)",
      
      mediator == "Amyg_vs_PCC_bcp" ~
        "FC between amygdala and posterior cingulate cortex (C)",
      
      mediator == "CO_vs_AUD_adult" ~
        "FC between cingulo-opercular and auditory networks (A)",
      
      mediator == "ctx_rh_parsorbitalis_gi" ~
        "R pars orbitalis GI",
      
      TRUE ~ mediator
    )
  )

# include CI
indirect_df <- indirect_df %>%
  mutate(
    ci_sig = if_else(
      ci_lower > 0 | ci_upper < 0,
      "CI excludes 0",
      "CI includes 0"
    )
  )


indirect_df <- indirect_df %>%
  mutate(
    pathway_label = paste0(
      mediator_label,
      " \u2192 ",
      outcome_label
    ),
    
    outcome_label = factor(
      outcome_label,
      levels = c(
        "Cognition",
        "Competence",
        "Dysregulation",
        "Externalizing",
        "Internalizing"
      )
    ),
    
    sleep_group = factor(
      sleep_group,
      levels = c(
        "Maternal insomnia\nsymptoms",
        "Maternal restless legs\nsymptoms",
        "Maternal sleep\nquality"
      )
    ),
    
    timepoint = factor(
      timepoint,
      levels = c("avg", "T1", "T2", "T3")
    )
  )

# order
indirect_df <- indirect_df %>%
  arrange(outcome_label, est) %>%
  group_by(outcome_label) %>%
  mutate(
    pathway_label = fct_inorder(pathway_label)
  ) %>%
  ungroup()

# ensure no wrapping and not figure cut offs
pathway_layout <- indirect_df %>%
  mutate(
    pathway_label = as.character(pathway_label)
  ) %>%
  distinct(
    outcome_label,
    pathway_label,
    mediator_label,
    .keep_all = TRUE
  ) %>%
  group_by(outcome_label) %>%
  mutate(
    original_order = row_number()
  ) %>%
  ungroup() %>%
  arrange(
    outcome_label,
    desc(original_order)
  ) %>%
  mutate(
    display_label = as.character(pathway_label),
    label_lines = 1L,
    row_height = 1.5
  ) %>%
  group_by(outcome_label) %>%
  mutate(
    y_local = sum(row_height) -
      cumsum(row_height) +
      row_height / 2,
    
    y_plot = 100 * (
      5 - as.integer(outcome_label)
    ) + y_local
  ) %>%
  ungroup()

# combine y axis positions 
indirect_df <- indirect_df %>%
  mutate(
    pathway_label = as.character(pathway_label)
  ) %>%
  left_join(
    pathway_layout %>%
      select(
        outcome_label,
        pathway_label,
        y_plot
      ),
    by = c(
      "outcome_label",
      "pathway_label"
    )
  )

if (anyNA(indirect_df$y_plot)) {
  stop(
    "Some pathways could not be matched to display positions."
  )
}


y_bounds <- pathway_layout %>%
  group_by(outcome_label) %>%
  summarise(
    y_bottom = min(
      y_plot - row_height / 2
    ) - 0.35,
    
    y_top = max(
      y_plot + row_height / 2
    ) + 0.35,
    
    .groups = "drop"
  )

bounds_long <- bind_rows(
  y_bounds %>%
    transmute(
      outcome_label,
      y_bound = y_bottom
    ),
  
  y_bounds %>%
    transmute(
      outcome_label,
      y_bound = y_top
    )
)

# consistent trimester cols
trimester_cols <- c(
  "avg" = "grey55",
  "T1" = "#1f78b4",
  "T2" = "#d9bf6b",
  "T3" = "#d67c1c"
)

#make figure
p_indirect_med <- ggplot(
  indirect_df,
  aes(
    x = est,
    y = y_plot,
    color = timepoint
  )
) +
  
  geom_blank(
    data = bounds_long,
    aes(
      x = 0,
      y = y_bound
    ),
    inherit.aes = FALSE
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  
  geom_segment(
    aes(
      x = ci_lower,
      xend = ci_upper,
      yend = y_plot
    ),
    linewidth = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  
  # open circles
  geom_point(
    shape = 21,
    fill = "white",
    size = 3.4,
    stroke = 1.0,
    show.legend = TRUE
  ) +
  
  # filled circles when CI excludes zero
  geom_point(
    data = ~ dplyr::filter(
      .x,
      ci_sig == "CI excludes 0"
    ),
    aes(fill = timepoint),
    shape = 21,
    size = 3.4,
    stroke = 1.0,
    show.legend = FALSE
  ) +
  
  facet_grid(
    outcome_label ~ sleep_group,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
scale_y_continuous(
  breaks = pathway_layout$y_plot,
  labels = pathway_layout$display_label,
  expand = expansion(
    mult = c(0, 0)
  )
) +
  
 #x axis
scale_x_continuous(
  breaks = function(lims) {
    
    b <- pretty(lims, n = 3)
    
    b <- b[
      b >= lims[1] &
        b <= lims[2]
    ]
    
    if (length(b) > 3) {
      b <- b[
        c(
          1,
          ceiling(length(b) / 2),
          length(b)
        )
      ]
    }
    
    unique(b)
  },
  
  labels = function(x) {
    sprintf("%.2f", x)
  }
) +
  
#color scheme
scale_color_manual(
  values = trimester_cols,
  breaks = c("avg", "T1", "T2", "T3"),
  labels = c(
    "Average", "T1", "T2", "T3"
  ),
  name = "Trimester",
  drop = FALSE
) +
  
  scale_fill_manual(
    values = trimester_cols,
    guide = "none",
    drop = FALSE
  ) +
  

  #titles of axis and main
labs(
  title = paste0(
    "Indirect associations of maternal sleep\n",
    "with infant outcomes via infant brain phenotypes"
  ),
  x = "Indirect effect (95% bootstrap CI)",
  y = NULL
) +
  
  theme_classic(
    base_size = 18
  ) +
  

theme(
  plot.title.position = "plot",
  
  plot.title = element_text(
    face = "bold",
    size = 32,
    hjust = 0.5,
    lineheight = 1.08,
    margin = margin(b = 14)
  ),
  
  strip.background = element_blank(),
  
  strip.text.x = element_text(
    face = "bold",
    size = 24,
    lineheight = 0.95,
    margin = margin(
      t = 6,
      b = 8
    )
  ),
  
  strip.text.y.left = element_text(
    face = "bold",
    size = 22,
    angle = 0
  ),
  
  axis.text.y = element_text(
    size = pathway_font_size,
    color = "black",
    lineheight = 0.94,
    hjust = 1,
    margin = margin(r = 8)
  ),
  
  axis.text.x = element_text(
    size = 18,
    color = "black"
  ),
  
  axis.title.x = element_text(
    size = 24,
    color = "black",
    margin = margin(t = 10)
  ),
  
  axis.ticks.y = element_line(
    color = "black"
  ),
  
  legend.position = "bottom",
  legend.justification = "center",
  legend.box = "horizontal",
  
  legend.title = element_text(
    size = 22,
    face = "bold"
  ),
  
  legend.text = element_text(
    size = 22
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
  
  panel.spacing.x = unit(
    1.2,
    "lines"
  ),
  
  panel.spacing.y = unit(
    0.8,
    "lines"
  ),
  
  plot.margin = margin(
    t = 18,
    r = 30,
    b = 15,
    l = 25
  )
) +
  
  # add legend symbols 
guides(
  color = guide_legend(
    nrow = 1,
    byrow = TRUE,
    
    override.aes = list(
      shape = rep(21, 4),
      size = rep(4.0, 4),
      stroke = rep(1.0, 4),
      
      fill = unname(
        trimester_cols[
          c("avg", "T1", "T2", "T3")
        ]
      )
    )
  )
)


print(p_indirect_med)
