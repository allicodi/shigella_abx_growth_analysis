# ---------------------------------------------------------------------------------
# Create supplementary figures for effect heterogeneity by study and baseline HAZ
# ---------------------------------------------------------------------------------

here::i_am("code/create_figures/figure_S14_effect_heterogeneity_site.R")

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------
# Load effect heterogeneity results
# ------------------------------------------------------------------------

study_case_only <- readRDS(here::here(
  "results/no_etiology/ipd_no_etiology_guideline_abx_site_no_models.Rds"
))

study_case_control <- readRDS(here::here(
  "results/case_control/ipd_results_msd_case_control_tac_or_culture_guideline_abx_study_no_models.Rds"
))

# ------------------------------------------------------------------------
# P-values for effect modification
# ------------------------------------------------------------------------

test_effect_mod <- function(object, ...) {
  
  msm_list <- object$aipw_est$results_object$aipw_msm
  
  stopifnot(
    !is.null(msm_list),
    all(c("msm_0", "msm_1") %in% names(msm_list))
  )
  
  # Combine influence functions from both antibiotic groups
  inf_fn_matrix <- do.call(
    cbind,
    lapply(msm_list[c("msm_0", "msm_1")], function(msm) {
      as.matrix(msm$inf_fns)
    })
  )
  
  # Covariance matrix of joint MSM coefficient estimates
  cov_matrix <- cov(inf_fn_matrix) / nrow(inf_fn_matrix)
  
  # Identify intercept positions in combined matrix
  n_coef_0 <- length(msm_list$msm_0$coefficients)
  n_coef_1 <- length(msm_list$msm_1$coefficients)
  
  intercept_idx <- c(1, n_coef_0 + 1)
  
  # Remove intercept rows and columns
  cov_beta <- cov_matrix[
    -intercept_idx,
    -intercept_idx,
    drop = FALSE
  ]
  
  # Combine study coefficients from both MSMs
  beta <- c(
    msm_list$msm_0$coefficients[-1],
    msm_list$msm_1$coefficients[-1]
  )
  
  # Joint Wald test
  wald_stat <- as.numeric(
    t(beta) %*% solve(cov_beta, beta)
  )
  
  # Degrees of freedom = total number of study coefficients
  df <- length(beta)
  
  pval <- pchisq(
    wald_stat,
    df = df,
    lower.tail = FALSE
  )
  
  return(pval)
}

# Calculate joint study heterogeneity p-values
p_study_case_only <- test_effect_mod(study_case_only)
p_study_case_control <- test_effect_mod(study_case_control)

print(p_study_case_only)
print(p_study_case_control)

# ------------------------------------------------------------------------
# Function to extract AIPW estimates
# ------------------------------------------------------------------------

extract_loo_results <- function(res, study_label) {
  
  res_df <- res$aipw_est$results_object$results_df
  se_vec <- res$aipw_est$results_object$se
  
  abx_levels <- res_df$abx_levels
  
  plot_data <- do.call(
    rbind,
    lapply(abx_levels, function(abx_level) {
      
      pt_est <- res_df$effect_inf_abx_level[
        res_df$abx_levels == abx_level
      ]
      
      se <- se_vec[paste0("effect_", abx_level)]
      
      data.frame(
        subgroup = as.character(abx_level),
        pt_est = as.numeric(pt_est),
        lower_ci = as.numeric(pt_est - 1.96 * se),
        upper_ci = as.numeric(pt_est + 1.96 * se),
        study_label = study_label
      )
    })
  )
  
  return(plot_data)
}

# ------------------------------------------------------------------------
# Function to create LOO forest plot
# ------------------------------------------------------------------------

plot_loo_results <- function(results,
                             study_labels,
                             heterogeneity_p = NULL,
                             title = NULL) {
  
  stopifnot(length(results) == length(study_labels))
  
  plot_data <- do.call(
    rbind,
    lapply(seq_along(results), function(i) {
      extract_loo_results(results[[i]], study_labels[i])
    })
  )
  
  # Relabel treatment groups
  plot_data$subgroup <- gsub(
    "_", " ", plot_data$subgroup
  )
  
  plot_data$subgroup <- gsub(
    "No azithromycin", "Placebo", plot_data$subgroup
  )
  
  plot_data$subgroup[plot_data$subgroup == "1"] <-
    "Guideline recommended antibiotics"
  
  plot_data$subgroup[plot_data$subgroup == "0"] <-
    "No guideline recommended antibiotics"
  
  subgroup_colors <- c(
    "Guideline recommended antibiotics" = "#00468BFF",
    "No guideline recommended antibiotics" = "#ED0000FF",
    "Azithromycin" = "#00468BFF",
    "Placebo" = "#ED0000FF"
  )
  
  plot_data$subgroup <- factor(
    plot_data$subgroup,
    levels = c(
      "No guideline recommended antibiotics",
      "Placebo",
      "Guideline recommended antibiotics",
      "Azithromycin"
    )
  )
  
  # Create estimate and confidence interval labels
  plot_data$estimate_label <- paste0(
    sprintf("%.2f", plot_data$pt_est),
    " (",
    sprintf("%.2f", plot_data$lower_ci),
    ", ",
    sprintf("%.2f", plot_data$upper_ci),
    ")"
  )
  
  # Format single joint heterogeneity p-value
  p_label <- NULL
  
  if (!is.null(heterogeneity_p)) {
    
    stopifnot(
      length(heterogeneity_p) == 1,
      is.numeric(heterogeneity_p)
    )
    
    if (!is.na(heterogeneity_p)) {
      p_label <- if (heterogeneity_p < 0.001) {
        "p < 0.001"
      } else {
        paste0("p = ", sprintf("%.3f", heterogeneity_p))
      }
    }
  }
  
  # X-axis scaling
  global_x_min <- min(plot_data$lower_ci, na.rm = TRUE)
  global_x_max <- max(plot_data$upper_ci, na.rm = TRUE)
  
  x_range <- global_x_max - global_x_min
  
  if (x_range == 0) {
    x_range <- 1
  }
  
  annotation_padding <- x_range * 0.5
  global_x_max <- global_x_max + annotation_padding
  x_loc <- global_x_max - annotation_padding / 2
  
  # Create individual panels
  plot_results <- lapply(study_labels, function(x) {
    
    sub_data <- plot_data[
      plot_data$study_label == x, ,
      drop = FALSE
    ]
    
    # Add joint heterogeneity p-value to Full IPD title only
    panel_title <- x
    
    if (x == "Full IPD" && !is.null(p_label)) {
      panel_title <- paste0(
        "Full IPD (Study heterogeneity ", p_label, ")"
      )
    }
    
    ggplot(
      data = sub_data,
      aes(x = pt_est, y = subgroup, color = subgroup)
    ) +
      geom_point(size = 5) +
      geom_errorbarh(
        aes(xmin = lower_ci, xmax = upper_ci),
        height = 0.3
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "#ED0000FF"
      ) +
      geom_text(
        aes(
          x = x_loc,
          label = estimate_label
        ),
        hjust = 0,
        vjust = 0.5,
        size = 4,
        color = "black"
      ) +
      labs(
        title = panel_title,
        x = "HAZ difference (95% CI)",
        y = "Treatment received"
      ) +
      scale_color_manual(
        values = subgroup_colors,
        drop = FALSE
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(
          size = 14,
          hjust = 0.5
        ),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      ) +
      coord_cartesian(
        xlim = c(global_x_min, global_x_max)
      )
  })
  
  # Combine plots vertically
  combined_plot <- patchwork::wrap_plots(
    plot_results,
    ncol = 1,
    axis_titles = "collect"
  ) &
    theme(legend.position = "none")
  
  if (!is.null(title)) {
    combined_plot <- combined_plot +
      patchwork::plot_annotation(title = title)
  }
  
  return(combined_plot)
}

# ------------------------------------------------------------------------
# 1. No etiology IPD - five studies
# ------------------------------------------------------------------------

studies <- c("ABCD", "MALED", "GEMS", "VIDA", "EFGH")

ipd_results <- readRDS(here::here(
  "results/no_etiology/ipd_no_etiology_guideline_abx_no_models.Rds"
))

loo_results <- lapply(studies, function(study) {
  
  readRDS(here::here(
    "results/no_etiology",
    paste0(
      "ipd_loo_",
      tolower(study),
      "_no_etiology_guideline_abx_no_models.Rds"
    )
  ))
})

results <- c(list(ipd_results), loo_results)

study_labels <- c(
  "Full IPD",
  "Excluding ABCD",
  "Excluding MAL-ED",
  "Excluding GEMS",
  "Excluding VIDA",
  "Excluding EFGH"
)

plot_no_etiology <- plot_loo_results(
  results = results,
  study_labels = study_labels,
  heterogeneity_p = p_study_case_only
)

ggsave(
  filename = here::here(
    "figures/fig_S15_sens_loo_no_etiology.png"
  ),
  plot = plot_no_etiology,
  width = 12,
  height = 12,
  dpi = 300
)

ggsave(
  filename = here::here(
    "figures/fig_S15_sens_loo_no_etiology.pdf"
  ),
  plot = plot_no_etiology,
  width = 12,
  height = 12
)

# ------------------------------------------------------------------------
# 2. Case-control IPD - three studies
# ------------------------------------------------------------------------

studies <- c("MALED", "GEMS", "VIDA")

ipd_results_cc <- readRDS(here::here(
  "results/case_control/ipd_results_msd_case_control_tac_or_culture_guideline_abx_no_models.Rds"
))

loo_results_cc <- lapply(studies, function(study) {
  
  readRDS(here::here(
    "results/case_control",
    paste0(
      "ipd_loo_",
      tolower(study),
      "_msd_case_control_tac_or_culture_guideline_abx_no_models.Rds"
    )
  ))
})

results_cc <- c(list(ipd_results_cc), loo_results_cc)

study_labels_cc <- c(
  "Full IPD",
  "Excluding MAL-ED",
  "Excluding GEMS",
  "Excluding VIDA"
)

plot_case_control <- plot_loo_results(
  results = results_cc,
  study_labels = study_labels_cc,
  heterogeneity_p = p_study_case_control
)

ggsave(
  filename = here::here(
    "figures/fig_S16_sens_loo_case_control.png"
  ),
  plot = plot_case_control,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  filename = here::here(
    "figures/fig_S16_sens_loo_case_control.pdf"
  ),
  plot = plot_case_control,
  width = 12,
  height = 9
)
