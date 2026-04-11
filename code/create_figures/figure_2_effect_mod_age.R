# ----------------------------------------------------------------
# Create effect modification by age figures for IPD 
# ----------------------------------------------------------------

here::i_am("code/create_figures/figure_2_effect_mod_age.R")

library(ggplot2)
library(patchwork)
library(grid)

no_etiology_results <- readRDS(here::here("results/no_etiology/ipd_no_etiology_no_models.Rds"))
case_control_results <- readRDS(here::here("results/case_control/ipd_results_msd_case_control_no_models_tac_or_culture.Rds"))

no_etiology_msm <- no_etiology_results$aipw_est$results_object$aipw_msm
who_no_etiology <- no_etiology_msm$msm_Guideline_recommended_antibiotics
maybe_no_etiology <- no_etiology_msm$msm_Possibly_effective_antibiotics
no_no_etiology <- no_etiology_msm$msm_No_or_ineffective_antibiotics

case_control_msm <- case_control_results$aipw_est$results_object$aipw_msm
who_case_control <- case_control_msm$msm_WHO_recommended_antibiotics
maybe_case_control <- case_control_msm$msm_Maybe_effective_antibiotics
no_case_control <- case_control_msm$msm_No_or_ineffective_antibiotics

age_range <- seq(0,59)

msm_age_plot_data <- function(age, coefs, cov_matrix, study_name){
  X <- cbind(1, age)
  est <- X %*% coefs
  se <- sqrt(rowSums((X %*% cov_matrix) * X))
  lower_ci <- est - 1.96*se
  upper_ci <- est + 1.96*se
  data.frame(age, est, se, lower_ci, upper_ci, study_name)
}

plot_who_ne <- msm_age_plot_data(age = age_range,
                              coefs = who_no_etiology$coefficients,
                              cov_matrix = who_no_etiology$cov,
                              study_name = "Guideline recommended antibiotics")

plot_maybe_ne <- msm_age_plot_data(age = age_range,
                                coefs = maybe_no_etiology$coefficients,
                                cov_matrix = maybe_no_etiology$cov,
                                study_name = "Possibly effective antibiotics")

plot_no_ne <- msm_age_plot_data(age = age_range,
                             coefs = no_no_etiology$coefficients,
                             cov_matrix = no_no_etiology$cov,
                             study_name = "No or ineffective antibiotics")

plot_who_cc <- msm_age_plot_data(age = age_range,
                                 coefs = who_case_control$coefficients,
                                 cov_matrix = who_case_control$cov,
                                 study_name = "Guideline recommended antibiotics")

plot_maybe_cc <- msm_age_plot_data(age = age_range,
                                   coefs = maybe_case_control$coefficients,
                                   cov_matrix = maybe_case_control$cov,
                                   study_name = "Possibly effective antibiotics")


plot_no_cc <- msm_age_plot_data(age = age_range,
                                coefs = no_case_control$coefficients,
                                cov_matrix = no_case_control$cov,
                                study_name = "No or ineffective antibiotics")

plot_ne <- bind_rows(plot_who_ne, plot_maybe_ne, plot_no_ne)
plot_cc <- bind_rows(plot_who_cc, plot_maybe_cc, plot_no_cc)

study_colors <- c("Guideline recommended antibiotics" = '#42B540FF',
                  "Possibly effective antibiotics" = '#00468BFF',
                  "No or ineffective antibiotics" = "#ED0000FF")

# Make study_name an ordered factor to control legend order
plot_ne$study_name <- factor(plot_ne$study_name,
                             levels = c("Guideline recommended antibiotics",
                                        "Possibly effective antibiotics",
                                        "No or ineffective antibiotics"))

plot_cc$study_name <- factor(plot_cc$study_name,
                             levels = c("Guideline recommended antibiotics",
                                        "Possibly effective antibiotics",
                                        "No or ineffective antibiotics"))

theme_no_grid <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

p1 <- ggplot(plot_ne, aes(x = age, y = est, color = study_name, fill = study_name)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.15, color = NA) +
  scale_color_manual(values = study_colors) +
  scale_fill_manual(values = study_colors) +
  labs(
    title = expression(bold("A.")~"Diarrhea with no etiology"),
    x = "Age (months)",
    y = "HAZ difference (95% CI)",
    color = "Antibiotics received",
    fill = "Antibiotics received"
  ) +
  coord_cartesian(ylim = c(-0.20, 0.10)) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    #axis.title.y = element_blank(),  
    #axis.title.x = element_blank(),  
    axis.line = element_line(color = "black", size = 0.4),
    axis.ticks = element_line(color = "black", size = 0.3),
    plot.title = element_text(size = 16, hjust = 0),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


p2 <- ggplot(plot_cc, aes(x = age, y = est, color = study_name, fill = study_name)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.15, color = NA) +
  scale_color_manual(values = study_colors) +
  scale_fill_manual(values = study_colors) +
  labs(
    title = expression(bold("B.")~"Non-diarrheal control"),
    x = "Age (months)",
    y = "HAZ difference (95% CI)",
    color = "Antibiotics received",
    fill = "Antibiotics received"
  ) +
  coord_cartesian(ylim = c(-0.20, 0.10)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.4),
    axis.ticks = element_line(color = "black", size = 0.3),
    plot.title = element_text(size = 16, hjust = 0),
    #axis.title.y = element_blank(), 
    #axis.title.x = element_blank(),  
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

# Remove legends from both plots
p1_noleg <- p1 + theme(legend.position = "none", axis.title.x = element_blank())
p2_noleg <- p2 + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

# Function to extract legend grob from ggplot
extract_legend <- function(p) {
  grobs <- ggplotGrob(p)$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  if(length(legend_index) == 0) return(NULL)
  grobs[[legend_index]]
}

# Extract legend from p2 (which still has legend)
legend_grob <- extract_legend(p2 + theme(legend.position = "bottom"))

# Combine plots side-by-side without legends
combined <- p1_noleg + p2_noleg + plot_layout(nrow = 1)

# Create shared x-axis label grob
x_lab <- textGrob(
  "Age (months)",
  gp = gpar(fontsize = 16),
  just = "center"
)

# Stack plots / x-axis label / legend
final_plot <- combined / 
  patchwork::wrap_elements(x_lab) / 
  patchwork::wrap_elements(legend_grob) + 
  plot_layout(heights = c(10, 1, 1))

ggsave(here::here("figures/fig_2_effect_mod_ipd.png"), final_plot, width = 12, height = 6, dpi = 300)

# ------------------------------------------------------------------------
# P-values for effect modification
# ------------------------------------------------------------------------

test_effect_mod <- function(
    object, ...
){
  stopifnot(!is.null(object$aipw_est$results_object$aipw_msm))
  pval <- rep(NA, length(object$aipw_est$results_object$aipw_msm))
  for(i in 1:length(object$aipw_est$results_object$aipw_msm)){
    msm <- object$aipw_est$results_object$aipw_msm[[i]]
    beta <- msm$coefficients[-1]
    p <- length(beta)
    cov_beta <- msm$cov[-1, -1]
    wald_stat <- t(beta) %*% solve(cov_beta) %*% beta
    pval[i] <- pchisq(wald_stat, df = p, lower.tail = FALSE)
  }
  names(pval) <- names(object$aipw_est$results_object$aipw_msm)
  return(pval)
}

test_effect_mod(no_etiology_results)
# no etiology -- 0.199 guideline rec, 0.011 possibly, 0.063 no

test_effect_mod(case_control_results)
# case control -- 0.198 guideline rec, 0.001 maybe, 0.007 no 
