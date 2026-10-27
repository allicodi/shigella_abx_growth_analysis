
# ----------------------------------------------------------------
# Create effect modification by baseline HAZ figures for IPD
# ----------------------------------------------------------------

here::i_am("code/create_figures/figure_S15_effect_heterogeneity_haz.R")

library(ggplot2)
library(patchwork)
library(grid)
library(dplyr)

# ----------------------------------------------------------------
# Load results
# ----------------------------------------------------------------

no_etiology_results <- readRDS(here::here(
  "results/no_etiology/ipd_no_etiology_guideline_abx_enr_haz_no_models.Rds"
))

case_control_results <- readRDS(here::here(
  "results/case_control/ipd_results_msd_case_control_tac_or_culture_guideline_abx_enr_haz_no_models.Rds"
))

no_etiology_msm <- no_etiology_results$aipw_est$results_object$aipw_msm
case_control_msm <- case_control_results$aipw_est$results_object$aipw_msm

# Two-level antibiotic groups
who_no_etiology <- no_etiology_msm$msm_1
no_no_etiology <- no_etiology_msm$msm_0

who_case_control <- case_control_msm$msm_1
no_case_control <- case_control_msm$msm_0

# ----------------------------------------------------------------
# P-values for effect modification
# ----------------------------------------------------------------

test_effect_mod <- function(object) {
  
  stopifnot(!is.null(object$aipw_est$results_object$aipw_msm))
  
  pval <- sapply(
    object$aipw_est$results_object$aipw_msm,
    function(msm) {
      
      beta <- msm$coefficients[-1]
      cov_beta <- as.matrix(msm$cov)[-1, -1, drop = FALSE]
      
      wald_stat <- t(beta) %*% solve(cov_beta) %*% beta
      
      pchisq(
        as.numeric(wald_stat),
        df = length(beta),
        lower.tail = FALSE
      )
    }
  )
  
  return(pval)
}

p_no_etiology <- test_effect_mod(no_etiology_results)
p_case_control <- test_effect_mod(case_control_results)

print(p_no_etiology)
print(p_case_control)

# ----------------------------------------------------------------
# Function to calculate MSM predictions and confidence intervals
# ----------------------------------------------------------------

msm_haz_plot_data <- function(haz, coefs, cov_matrix, study_name) {
  
  X <- cbind(1, haz)
  
  est <- as.numeric(X %*% coefs)
  se <- sqrt(rowSums((X %*% cov_matrix) * X))
  
  lower_ci <- est - 1.96 * se
  upper_ci <- est + 1.96 * se
  
  data.frame(
    haz = haz,
    est = est,
    se = se,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    study_name = study_name
  )
}

# ----------------------------------------------------------------
# Create prediction data
# ----------------------------------------------------------------

haz_range <- seq(-4, 2, by = 0.05)

plot_who_ne <- msm_haz_plot_data(
  haz = haz_range,
  coefs = who_no_etiology$coefficients,
  cov_matrix = who_no_etiology$cov,
  study_name = "Guideline recommended antibiotics"
)

plot_no_ne <- msm_haz_plot_data(
  haz = haz_range,
  coefs = no_no_etiology$coefficients,
  cov_matrix = no_no_etiology$cov,
  study_name = "No guideline recommended antibiotics"
)

plot_who_cc <- msm_haz_plot_data(
  haz = haz_range,
  coefs = who_case_control$coefficients,
  cov_matrix = who_case_control$cov,
  study_name = "Guideline recommended antibiotics"
)

plot_no_cc <- msm_haz_plot_data(
  haz = haz_range,
  coefs = no_case_control$coefficients,
  cov_matrix = no_case_control$cov,
  study_name = "No guideline recommended antibiotics"
)

plot_ne <- bind_rows(plot_who_ne, plot_no_ne)
plot_cc <- bind_rows(plot_who_cc, plot_no_cc)

# ----------------------------------------------------------------
# Colors and factor levels
# ----------------------------------------------------------------

study_colors <- c(
  "Guideline recommended antibiotics" = "#00468BFF",
  "No guideline recommended antibiotics" = "#ED0000FF"
)

study_levels <- c(
  "Guideline recommended antibiotics",
  "No guideline recommended antibiotics"
)

plot_ne$study_name <- factor(
  plot_ne$study_name,
  levels = study_levels
)

plot_cc$study_name <- factor(
  plot_cc$study_name,
  levels = study_levels
)

# ----------------------------------------------------------------
# Format p-values for plot annotations
# ----------------------------------------------------------------

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

# Include effect modification p-values in each panel
p_labels_ne <- c(
  paste0(
    "Guideline recommended antibiotics: ",
    format_p(p_no_etiology["msm_1"])
  ),
  paste0(
    "No guideline recommended antibiotics: ",
    format_p(p_no_etiology["msm_0"])
  )
)

p_labels_cc <- c(
  paste0(
    "Guideline recommended antibiotics: ",
    format_p(p_case_control["msm_1"])
  ),
  paste0(
    "No guideline recommended antibiotics: ",
    format_p(p_case_control["msm_0"])
  )
)

# ----------------------------------------------------------------
# Plot theme
# ----------------------------------------------------------------

plot_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(size = 16, hjust = 0),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

# ----------------------------------------------------------------
# A. Case-only IPD
# ----------------------------------------------------------------

p1 <- ggplot(
  plot_ne,
  aes(x = haz, y = est, color = study_name, fill = study_name)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = study_colors) +
  scale_fill_manual(values = study_colors) +
  labs(
    title = expression(bold("A.") ~ "Diarrhea with no etiology"),
    x = "Baseline HAZ",
    y = "HAZ difference (95% CI)",
    color = "Antibiotics received",
    fill = "Antibiotics received"
  ) +
  coord_cartesian(xlim = c(-4, 2)) +
  plot_theme +
  theme(legend.position = "none")

# ----------------------------------------------------------------
# B. Case-control IPD
# ----------------------------------------------------------------

p2 <- ggplot(
  plot_cc,
  aes(x = haz, y = est, color = study_name, fill = study_name)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = study_colors) +
  scale_fill_manual(values = study_colors) +
  labs(
    title = expression(bold("B.") ~ "Non-diarrheal control"),
    x = "Baseline HAZ",
    y = "HAZ difference (95% CI)",
    color = "Antibiotics received",
    fill = "Antibiotics received"
  ) +
  coord_cartesian(xlim = c(-4, 2)) +
  plot_theme

# ----------------------------------------------------------------
# Combine plots with shared legend
# ----------------------------------------------------------------

p1_noleg <- p1 +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

p2_noleg <- p2 +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

extract_legend <- function(p) {
  
  grobs <- ggplotGrob(p)$grobs
  
  legend_index <- which(
    sapply(grobs, function(x) x$name) == "guide-box"
  )
  
  if (length(legend_index) == 0) return(NULL)
  
  grobs[[legend_index]]
}

legend_grob <- extract_legend(
  p2 + theme(legend.position = "bottom")
)

combined <- p1_noleg + p2_noleg +
  plot_layout(nrow = 1)

x_lab <- textGrob(
  "Baseline HAZ",
  gp = gpar(fontsize = 16),
  just = "center"
)

final_plot <- combined /
  patchwork::wrap_elements(x_lab) /
  patchwork::wrap_elements(legend_grob) +
  plot_layout(heights = c(10, 1, 1))

# ----------------------------------------------------------------
# Save figures
# ----------------------------------------------------------------

ggsave(
  here::here("figures/fig_S17_effect_mod_enr_haz.png"),
  final_plot,
  width = 12,
  height = 6,
  dpi = 300
)

ggsave(
  here::here("figures/fig_S17_effect_mod_enr_haz.pdf"),
  final_plot,
  width = 12,
  height = 6
)

# ----------------------------------------------------------------
# Print effect modification p-values
# ----------------------------------------------------------------

cat("\nCase-only IPD:\n")
print(p_labels_ne)

cat("\nCase-control IPD:\n")
print(p_labels_cc)
