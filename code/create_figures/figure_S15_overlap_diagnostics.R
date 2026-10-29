# --------------------------------------------------------------------------------
# Create overlap diagnostics figures for sensitivity analysis
# --------------------------------------------------------------------------------

here::i_am("code/create_figures/figure_S15_overlap_diagnostics.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

no_etiology <- readRDS(here::here(
  "results/no_etiology/ipd_no_etiology_guideline_abx_with_prop_vec.Rds"
))

case_control <- readRDS(here::here(
  "results/case_control/ipd_results_msd_case_control_tac_or_culture_guideline_abx_with_prop_vec.Rds"
))

# ------------------------------------------------------------------------
# Plot propensity score distributions
# ------------------------------------------------------------------------

plot_propensity_histograms <- function(prop_vector_list,
                                       bins = 30,
                                       case_control = FALSE) {
  
  abx_labels <- c(
    "abx_0" = "No guideline recommended antibiotics",
    "abx_1" = "Guideline recommended antibiotics"
  )
  
  # Consistent colors across analyses
  blue_cols <- c(
    "abx_0" = "#A6CEE3",
    "abx_1" = "#08519C"
  )
  
  red_cols <- c(
    "abx_0" = "#FCAE91",
    "abx_1" = "#CB181D"
  )
  
  control_col <- "#E57373"
  
  # Generate legend labels automatically
  treatment_cols <- c(
    setNames(
      blue_cols,
      paste0("Attributable: ", abx_labels[names(blue_cols)])
    ),
    setNames(
      red_cols,
      paste0("No etiology: ", abx_labels[names(red_cols)])
    ),
    "Control" = control_col
  )
  
  # Extract propensity matrices
  get_prop <- function(prefix) {
    prop_vector_list %>%
      select(starts_with(paste0(prefix, "."))) %>%
      rename_with(~ sub(paste0("^", prefix, "\\."), "", .x))
  }
  
  prop_1a <- get_prop("prop_vectors_1a")
  prop_3a <- get_prop("prop_vectors_3a")
  
  # ----------------------------------------------------------------------
  # Prepare histogram data
  # ----------------------------------------------------------------------
  
  if (!case_control) {
    
    prop_1b <- get_prop("prop_vectors_1b")
    prop_3b <- get_prop("prop_vectors_3b")
    
    # Antibiotic propensity
    df1 <- bind_rows(
      prop_1a %>%
        pivot_longer(
          everything(),
          names_to = "abx",
          values_to = "prob"
        ) %>%
        mutate(
          group = paste0("Attributable: ", abx_labels[abx])
        ),
      
      prop_1b %>%
        pivot_longer(
          everything(),
          names_to = "abx",
          values_to = "prob"
        ) %>%
        mutate(
          group = paste0("No etiology: ", abx_labels[abx])
        )
    )
    
    # Outcome observation probability
    df3 <- bind_rows(
      (1 - prop_3a) %>%
        pivot_longer(
          everything(),
          names_to = "abx",
          values_to = "prob"
        ) %>%
        mutate(
          group = paste0("Attributable: ", abx_labels[abx])
        ),
      
      (1 - prop_3b) %>%
        pivot_longer(
          everything(),
          names_to = "abx",
          values_to = "prob"
        ) %>%
        mutate(
          group = paste0("No etiology: ", abx_labels[abx])
        )
    )
    
  } else {
    
    # --------------------------------------------------------------------
    # Identify observations with valid case-specific antibiotic propensities
    #
    # The case-control vectors contain:
    #   - 2,783 rows where antibiotic probabilities sum to 1
    #   - 5,898 rows where both probabilities equal 0.01
    #
    # Exclude apparent placeholder rows from case-specific distributions.
    # --------------------------------------------------------------------
    
    case_rows <- abs(
      prop_vector_list$prop_vectors_1a.abx_0 +
        prop_vector_list$prop_vectors_1a.abx_1 - 1
    ) < 1e-6
    
    # Check number of retained observations
    message(
      "Case-control analysis: retaining ",
      sum(case_rows, na.rm = TRUE),
      " of ",
      length(case_rows),
      " observations for case-specific propensity distributions."
    )
    
    # Antibiotic propensity among cases only
    df1 <- prop_1a %>%
      filter(case_rows) %>%
      pivot_longer(
        everything(),
        names_to = "abx",
        values_to = "prob"
      ) %>%
      mutate(
        group = paste0("Attributable: ", abx_labels[abx])
      )
    
    # Outcome observation among cases and controls
    df3 <- bind_rows(
      (1 - prop_3a) %>%
        filter(case_rows) %>%
        pivot_longer(
          everything(),
          names_to = "abx",
          values_to = "prob"
        ) %>%
        mutate(
          group = paste0("Attributable: ", abx_labels[abx])
        ),
      
      data.frame(
        prob = 1 - prop_vector_list$control,
        group = "Control"
      )
    )
  }
  
  # ----------------------------------------------------------------------
  # Shared histogram styling
  # ----------------------------------------------------------------------
  
  hist_geom <- function() {
    geom_histogram(
      bins = bins,
      alpha = 0.45,
      color = "white",
      linewidth = 0.15,
      position = "identity",
      boundary = 0
    )
  }
  
  common_theme <- theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
    )
  
  # ----------------------------------------------------------------------
  # Plot 1: Antibiotic propensity
  # ----------------------------------------------------------------------
  
  p1 <- ggplot(df1, aes(x = prob, fill = group)) +
    hist_geom() +
    scale_fill_manual(
      values = treatment_cols,
      breaks = names(treatment_cols)[
        names(treatment_cols) %in% unique(df1$group)
      ],
      name = "Antibiotic treatment / group"
    ) +
    guides(
      fill = guide_legend(
        ncol = if (case_control) 1 else 2,
        byrow = TRUE
      )
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = if (case_control) {
        "Antibiotic Propensity: Case-Control Analysis"
      } else {
        "Antibiotic Propensity: Case-Only Analysis"
      },
      x = "Predicted Probability",
      y = "Count"
    ) +
    common_theme
  
  # ----------------------------------------------------------------------
  # Plot 2: Outcome observation propensity
  # ----------------------------------------------------------------------
  
  p3 <- ggplot(df3, aes(x = prob, fill = group)) +
    hist_geom() +
    scale_fill_manual(
      values = treatment_cols,
      breaks = names(treatment_cols)[
        names(treatment_cols) %in% unique(df3$group)
      ],
      name = "Antibiotic treatment / group"
    ) +
    guides(
      fill = guide_legend(
        ncol = if (case_control) 1 else 2,
        byrow = TRUE
      )
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = if (case_control) {
        "Outcome Observation Propensity: Case-Control Analysis"
      } else {
        "Outcome Observation Propensity: Case-Only Analysis"
      },
      x = "Predicted Probability",
      y = "Count"
    ) +
    common_theme
  
  return(list(p1 = p1, p3 = p3))
}

# ------------------------------------------------------------------------
# Case-only (no-etiology) analysis
# ------------------------------------------------------------------------

prop_vector_list <- no_etiology$aipw_est$aipw_models$prop_vector_list

plots <- plot_propensity_histograms(
  prop_vector_list,
  case_control = FALSE
)

# ------------------------------------------------------------------------
# Case-control analysis
# ------------------------------------------------------------------------

prop_vector_list_cc <- case_control$aipw_est$aipw_models$prop_vector_list

plots_cc <- plot_propensity_histograms(
  prop_vector_list_cc,
  case_control = TRUE
)

# ------------------------------------------------------------------------
# Combine into one figure, with legends under each subplot
#
# A: Antibiotic Propensity - Case-Only Analysis
# B: Outcome Observation Propensity - Case-Only Analysis
# C: Antibiotic Propensity - Case-Control Analysis
# D: Outcome Observation Propensity - Case-Control Analysis
# ------------------------------------------------------------------------

combined_plot <- (
  (plots$p1 | plots$p3) /
    (plots_cc$p1 | plots_cc$p3)
) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = grid::unit(0.4, "cm"),
    legend.background = element_blank(),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.tag = element_text(
      face = "bold",
      size = 13
    )
  )

combined_plot

# ------------------------------------------------------------------------
# Save figure
# ------------------------------------------------------------------------

ggsave(
  filename = here::here("figures/fig_S18_overlap_diagnostics.pdf"),
  plot = combined_plot,
  width = 13,
  height = 9,
  units = "in"
)
