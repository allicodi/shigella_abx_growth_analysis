# ------------------------------------------------------------------------
# Create figure of the year (MAL-ED longitudinal figure)
# ------------------------------------------------------------------------

here::i_am("code/create_figures/figure_3_MALED_longitudinal.R")

# Vector of months
months <- 1:12

# Load all the RDS files into a list
mo_list <- lapply(months, function(mo) {
  readRDS(here::here(sprintf("results/case_control/longitudinal/case_control_no_abx_%d.Rds", mo)))
})

mo_list_severe <- lapply(months, function(mo) {
  readRDS(here::here(sprintf("results/case_control/longitudinal/case_control_no_abx_msd_%d.Rds", mo)))
})

mo_list_mild <- lapply(months, function(mo) {
  readRDS(here::here(sprintf("results/case_control/longitudinal/case_control_no_abx_lsd_%d.Rds", mo)))
})

# Create a data frame with point estimates and confidence intervals
plot_data <- data.frame(
  subgroup = paste0(months, " month"),
  pt_est = sapply(mo_list, function(mo) mo$results_object$results_df$effect_inf_abx_level),
  lower_ci = sapply(mo_list, function(mo) {
    est <- mo$results_object$results_df$effect_inf_abx_level
    se <- mo$results_object$se['effect_observed']
    est - 1.96 * se
  }),
  upper_ci = sapply(mo_list, function(mo) {
    est <- mo$results_object$results_df$effect_inf_abx_level
    se <- mo$results_object$se['effect_observed']
    est + 1.96 * se
  })
)

plot_data_severe <- data.frame(
  subgroup = paste0(months, " month"),
  pt_est = sapply(mo_list_severe, function(mo) mo$results_object$results_df$effect_inf_abx_level),
  lower_ci = sapply(mo_list_severe, function(mo) {
    est <- mo$results_object$results_df$effect_inf_abx_level
    se <- mo$results_object$se['effect_observed']
    est - 1.96 * se
  }),
  upper_ci = sapply(mo_list_severe, function(mo) {
    est <- mo$results_object$results_df$effect_inf_abx_level
    se <- mo$results_object$se['effect_observed']
    est + 1.96 * se
  })
)

plot_data_mild <- data.frame(
  subgroup = paste0(months, " month"),
  pt_est = sapply(mo_list_mild, function(mo) mo$results_object$results_df$effect_inf_abx_level),
  lower_ci = sapply(mo_list_mild, function(mo) {
    est <- mo$results_object$results_df$effect_inf_abx_level
    se <- mo$results_object$se['effect_observed']
    est - 1.96 * se
  }),
  upper_ci = sapply(mo_list_mild, function(mo) {
    est <- mo$results_object$results_df$effect_inf_abx_level
    se <- mo$results_object$se['effect_observed']
    est + 1.96 * se
  })
)

# Add numeric month column for plotting
plot_data$month_num <- months
plot_data_severe$month_num <- months
plot_data_mild$month_num <- months

plot_data$group <- "Any Shigella"
plot_data_severe$group <- "Moderate-to-Severe Shigella"
plot_data_mild$group <- "Less-severe Shigella"

# Combine the datasets
combined_plot_data <- rbind(plot_data_mild, plot_data, plot_data_severe)
combined_plot_data$group <- factor(combined_plot_data$group, levels = c("Less-severe Shigella", "Any Shigella", "Moderate-to-Severe Shigella"))

final_plot <- ggplot(combined_plot_data, aes(x = month_num, y = pt_est, color = group, fill = group)) +
  # CI ribbon for raw estimates
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
  # Raw line and points
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  # Smoothed piecewise linear fit
  # stat_smooth(
  #   method = "glm",
  #   formula = y ~ -1 + x + I(pmax(0, x - 3)) + I(pmax(0, x - 9)),
  #   size = 1.2,
  #   linetype = "dotted",
  #   se = FALSE
  # ) +
  # Axis and labels
  scale_x_continuous(breaks = months) +
  labs(
    x = "Months since diarrhea episode",
    y = "HAZ difference (95% CI)",
    # title = "Longitudinal Shigella growth effects in MAL-ED with observed antibiotic use",
    color = "Shigella severity",
    fill = "Shigella severity"
  ) +
  scale_color_manual(values = c(
    "Any Shigella" = "#00468BFF",
    "Moderate-to-Severe Shigella" = "#ED0000FF",
    "Less-severe Shigella" = "#42B540FF"
  )) +
  scale_fill_manual(values = c(
    "Any Shigella" = "#00468BFF",
    "Moderate-to-Severe Shigella" = "#ED0000FF",
    "Less-severe Shigella" = "#42B540FF"
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(filename = here::here("figures/fig_3_longitudinal.png"), final_plot, width = 10, height = 6)

