# ---------------------------------------------------
# Age stratified observed antibiotics use
# ---------------------------------------------------

here::i_am("code/create_figures/figure_S2_observed_abx_use.R")

ipd_0_12 <- readRDS("results/case_control/ipd_no_abx_0_11_msd.Rds")
ipd_12_24 <- readRDS("results/case_control/ipd_no_abx_12_23_msd.Rds")
ipd_24_59 <- readRDS("results/case_control/ipd_no_abx_24_59_msd.Rds")

plot_data <- data.frame(subgroup = c("0-11 months", "12-23 months", "24-59 months"),
                        pt_est = c(ipd_0_12$results_object$results_df$effect_inf_abx_level,
                                   ipd_12_24$results_object$results_df$effect_inf_abx_level,
                                   ipd_24_59$results_object$results_df$effect_inf_abx_level),
                        lower_ci = c(ipd_0_12$results_object$results_df$effect_inf_abx_level - 1.96*ipd_0_12$results_object$se['effect_observed'],
                                     ipd_12_24$results_object$results_df$effect_inf_abx_level - 1.96*ipd_12_24$results_object$se['effect_observed'],
                                     ipd_24_59$results_object$results_df$effect_inf_abx_level - 1.96*ipd_24_59$results_object$se['effect_observed']),
                        upper_ci = c(ipd_0_12$results_object$results_df$effect_inf_abx_level + 1.96*ipd_0_12$results_object$se['effect_observed'],
                                     ipd_12_24$results_object$results_df$effect_inf_abx_level + 1.96*ipd_12_24$results_object$se['effect_observed'],
                                     ipd_24_59$results_object$results_df$effect_inf_abx_level + 1.96*ipd_24_59$results_object$se['effect_observed']))
# Set the desired order of subgroups
plot_data$subgroup <- factor(plot_data$subgroup,
                             levels = c("24-59 months", "12-23 months", "0-11 months"))

# Get xaxis scaling 
global_x_min <- min(plot_data$lower_ci)
global_x_max <- max(plot_data$upper_ci)

# Add padding for annotations
annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

# Define custom colors for each subgroup
subgroup_colors <- c(
  "0-11 months" = "#00468BFF" , # blue
  "12-23 months" = "#ED0000FF",  # red
  "24-59 months" = "#42B540FF"         # green
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
       y = "Age strata",
       color = "Age strata") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = subgroup_colors) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(here::here("figures/fig_S2_obs_abx_use_age_stratified.png"), width = 14, height = 7)
