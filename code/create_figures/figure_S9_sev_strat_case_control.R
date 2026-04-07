# ---------------------------------------------------
# Severity stratified observed antibiotics use
# ---------------------------------------------------

here::i_am("code/create_figures/figure_S9_sev_strat_case_control.R")

no_abx <- readRDS("results/case_control/ipd_no_abx.Rds")
no_abx_lsd <- readRDS("results/case_control/ipd_no_abx_lsd.Rds")
no_abx_msd <- readRDS("results/case_control/ipd_no_abx_msd.Rds")

plot_data <- data.frame(subgroup = c("Less severe diarrhea", "All diarrhea", "Moderate to severe diarrhea"),
                        pt_est = c(no_abx_lsd$results_object$results_df$effect_inf_abx_level,
                                   no_abx$results_object$results_df$effect_inf_abx_level,
                                   no_abx_msd$results_object$results_df$effect_inf_abx_level),
                        lower_ci = c(no_abx_lsd$results_object$results_df$effect_inf_abx_level - 1.96*no_abx_lsd$results_object$se['effect_observed'],
                                     no_abx$results_object$results_df$effect_inf_abx_level - 1.96*no_abx$results_object$se['effect_observed'],
                                     no_abx_msd$results_object$results_df$effect_inf_abx_level - 1.96*no_abx_msd$results_object$se['effect_observed']),
                        upper_ci = c(no_abx_lsd$results_object$results_df$effect_inf_abx_level + 1.96*no_abx_lsd$results_object$se['effect_observed'],
                                     no_abx$results_object$results_df$effect_inf_abx_level + 1.96*no_abx$results_object$se['effect_observed'],
                                     no_abx_msd$results_object$results_df$effect_inf_abx_level + 1.96*no_abx_msd$results_object$se['effect_observed']))
# Set the desired order of subgroups
plot_data$subgroup <- factor(plot_data$subgroup,
                             levels = c("All diarrhea", "Moderate to severe diarrhea", "Less severe diarrhea"))

# Get xaxis scaling 
global_x_min <- min(plot_data$lower_ci)
global_x_max <- max(plot_data$upper_ci)

# Add padding for annotations
annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

# Define custom colors for each subgroup
subgroup_colors <- c(
  "Less severe diarrhea" = "#00468BFF" , # blue
  "Moderate to severe diarrhea" = "#ED0000FF",  # red
  "All diarrhea" = "#42B540FF"         # green
)

ggplot2::ggplot(data = plot_data, aes(x = pt_est, y = subgroup, color = subgroup)) +
  geom_point(size = 5) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(x = -0.0225,
                label = paste0(round(pt_est, 3),
                               " (", round(lower_ci, 3), ", ", round(upper_ci, 3), ")")),
            hjust = 0, vjust = 0.5, size = 4) +
  labs(x = "HAZ difference (95% CI)",
       y = "Shigella severity",
       color = "Shigella severity") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = subgroup_colors) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(here::here("figures/fig_S9_obs_abx_use_severity.png"), width = 14, height = 7)
