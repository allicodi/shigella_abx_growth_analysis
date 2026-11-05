# --------------------------------------------------------------------------------
# Create overlap diagnostics figures for sensitivity analysis
# --------------------------------------------------------------------------------

here::i_am("code/create_figures/figure_S17_overlap_diagnostics.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

no_etiology <- readRDS(here::here(
  "results/no_etiology/ipd_no_etiology_guideline_abx_with_prop_vec.Rds"
))

# no etiology dataset -- use to get marginalization over shig_attr
data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology_bin_abx.Rds"))

# check
# nrow(data) == length(no_etiology$aipw_est$aipw_models$prop_vector_list$prop_vectors_1a.abx_1)
shig_attr_idx <- which(data$shig_attr == 1)

case_control <- readRDS(here::here(
  "results/case_control/ipd_results_msd_case_control_tac_or_culture_guideline_abx_with_prop_vec.Rds"
))

# case-control dataset -- use to get marginalization over cases -- FROM ONE HOT ENCODED VERSION (reordered)
ipd_data <- readRDS(here::here("data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig_bin_abx.Rds"))

# Impute missing values, cases and controls separate
case_ipd_data <- ipd_data[ipd_data$case == 1,]
control_ipd_data <- ipd_data[ipd_data$case == 0,]

imp_data_case <- impute_covariates(data = case_ipd_data,
                                   imp_covariates = c("sex",
                                                      "age",
                                                      "ses_quintile",
                                                      "enr_haz",
                                                      "edu_bin",
                                                      "site",
                                                      "I_followup_days",
                                                      "I_followup_days_x_followup_days",
                                                      "dysentery",
                                                      "any_vom",
                                                      "any_fev",
                                                      "dehyd_level",
                                                      "lsstools",
                                                      "duration_pre_enroll",
                                                      "shigella_new" ,
                                                      "rotavirus_new",
                                                      "adenovirus_new" ,
                                                      "etec_new" ,
                                                      "cryptosporidium_new"  ,
                                                      "astrovirus_new" ,
                                                      "norovirus_new" ,
                                                      "tepec_new",
                                                      "campylobacter_new" ,
                                                      "sapovirus_new" ,
                                                      "giardia_new"  ,
                                                      "e_bieneusi_new"  ,
                                                      "eaec_new"))

imp_data_control <- impute_covariates(data = control_ipd_data,
                                      imp_covariates = c("sex",
                                                         "age",
                                                         "ses_quintile",
                                                         "enr_haz",
                                                         "edu_bin",
                                                         "site",
                                                         "I_followup_days",
                                                         "I_followup_days_x_followup_days"))


imp_data_case <- imp_data_case[,colnames(imp_data_case) %in% colnames(ipd_data), drop = FALSE]
imp_data_control <- imp_data_control[,colnames(imp_data_control) %in% colnames(ipd_data), drop = FALSE]
imp_data <- rbind(imp_data_case, imp_data_control)

one_hot_data <- one_hot_encode(imp_data,
                               laz_var_name = "final_haz",
                               covariate_list = c("sex",
                                                  "age",
                                                  "ses_quintile",
                                                  "enr_haz",
                                                  "I_followup_days",
                                                  "I_followup_days_x_followup_days",
                                                  "site",
                                                  "edu_bin",
                                                  "num_hh_lt5",
                                                  "imp_water",
                                                  "imp_sanit"),
                               abx_var_name = "guideline_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("dysentery",
                                                 "any_vom",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

# check
# nrow(case_control_data) == length(case_control$aipw_est$aipw_models$prop_vector_list$prop_vectors_1a.abx_0)
case_idx <- which(one_hot_data$data$case == 1)

# ------------------------------------------------------------------------
# Plot propensity score distributions
# ------------------------------------------------------------------------

plot_propensity_histograms <- function(prop_vector_list,
                                       marg_idx,
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
  
  # Separate control color from no-etiology colors
  control_col <- "#E57373"
  
  # Generate treatment legend labels automatically
  treatment_cols <- c(
    setNames(
      unname(blue_cols),
      paste0("Attributable: ", abx_labels[names(blue_cols)])
    ),
    setNames(
      unname(red_cols),
      paste0("No etiology: ", abx_labels[names(red_cols)])
    ),
    "Control" = control_col
  )
  
  # Group colors for attribution and case-control propensities
  group_cols <- c(
    "Attributable" = blue_cols[["abx_0"]],
    "No etiology" = red_cols[["abx_0"]],
    "Control" = control_col
  )
  
  n <- nrow(prop_vector_list)
  
  stopifnot(
    length(marg_idx) > 0,
    !anyNA(marg_idx),
    all(marg_idx %in% seq_len(n))
  )
  
  # Extract propensity matrices for marginalization population
  get_prop <- function(prefix) {
    prop_vector_list %>%
      slice(marg_idx) %>%
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
    prop_2a <- get_prop("prop_vectors_2a")
    prop_3b <- get_prop("prop_vectors_3b")
    
    # Histogram 1: Antibiotic propensity
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
    
    # Histogram 2: Attribution propensity
    # Histogram 2: Attribution propensity - No etiology only
    df2 <- prop_2a %>%
      transmute(
        group = "No etiology",
        prob = no_attr
      )
    
    # Histogram 3: Outcome observation probability
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
    
    # Histogram 1: Antibiotic propensity among cases
    df1 <- prop_1a %>%
      pivot_longer(
        everything(),
        names_to = "abx",
        values_to = "prob"
      ) %>%
      mutate(
        group = paste0("Attributable: ", abx_labels[abx])
      )
    
    # Histogram 2: Case-control propensity
    df2 <- prop_vector_list %>%
      slice(marg_idx) %>%
      transmute(
        group = "Control",
        prob = 1 - case
      )
    
    # Histogram 3: Outcome observation probability
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
      
      prop_vector_list %>%
        slice(marg_idx) %>%
        transmute(
          prob = 1 - control,
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
        size = 11,
        face = "bold",
        hjust = 0.5
      )
    )
  
  # Retain only legend categories present in each plot
  get_breaks <- function(df, colors) {
    names(colors)[names(colors) %in% unique(df$group)]
  }
  
  # ----------------------------------------------------------------------
  # Plot 1: Antibiotic propensity (A and D)
  # ----------------------------------------------------------------------
  
  p1 <- ggplot(df1, aes(x = prob, fill = group)) +
    hist_geom() +
    scale_fill_manual(
      values = treatment_cols,
      breaks = get_breaks(df1, treatment_cols),
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
  # Plot 2: Attribution / case-control propensity (B and E)
  # ----------------------------------------------------------------------
  
  p2 <- ggplot(df2, aes(x = prob, fill = group)) +
    hist_geom() +
    scale_fill_manual(
      values = group_cols,
      breaks = get_breaks(df2, group_cols),
      name = "Attribution / comparison group"
    ) +
    guides(
      fill = guide_legend(nrow = 1)
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = if (case_control) {
        "Case-Control Propensity: Case-Control Analysis"
      } else {
        "Attribution Propensity: Case-Only Analysis"
      },
      x = "Predicted Probability",
      y = "Count"
    ) +
    common_theme
  
  # ----------------------------------------------------------------------
  # Plot 3: Outcome observation propensity (C and F)
  # ----------------------------------------------------------------------
  
  p3 <- ggplot(df3, aes(x = prob, fill = group)) +
    hist_geom() +
    scale_fill_manual(
      values = treatment_cols,
      breaks = get_breaks(df3, treatment_cols),
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
  
  return(list(p1 = p1, p2 = p2, p3 = p3))
}


# Case-only (no-etiology comparator) analysis
prop_vector_list <- no_etiology$aipw_est$aipw_models$prop_vector_list

stopifnot(nrow(data) == nrow(prop_vector_list))

plots <- plot_propensity_histograms(
  prop_vector_list,
  marg_idx = shig_attr_idx,
  case_control = FALSE
)

# Case-control analysis
prop_vector_list_cc <- case_control$aipw_est$aipw_models$prop_vector_list

stopifnot(nrow(ipd_data) == nrow(prop_vector_list_cc))

plots_cc <- plot_propensity_histograms(
  prop_vector_list_cc,
  marg_idx = case_idx,
  case_control = TRUE
)

# ------------------------------------------------------------------------
# Combine into six-panel figure (3 rows x 2 columns)
#
# A: Attribution Propensity - Case-Only Analysis
# B: Case-Control Propensity - Case-Control Analysis
#
# C: Antibiotic Propensity - Case-Only Analysis
# D: Antibiotic Propensity - Case-Control Analysis
#
# E: Outcome Observation Propensity - Case-Only Analysis
# F: Outcome Observation Propensity - Case-Control Analysis
# ------------------------------------------------------------------------

combined_plot <- (
  (plots$p2 | plots_cc$p2) /
    (plots$p1 | plots_cc$p1) /
    (plots$p3 | plots_cc$p3)
) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(0.45, "cm"),
    legend.spacing.x = grid::unit(0.2, "cm"),
    legend.background = element_blank(),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

combined_plot

ggsave(
  filename = here::here("figures/fig_S17_overlap_diagnostics.pdf"),
  plot = combined_plot,
  width = 15,
  height = 15,
  units = "in"
)

