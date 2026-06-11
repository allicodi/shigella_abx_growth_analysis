# ---------------------------------------------------
# Dysentery stratified 
# ---------------------------------------------------

# 0 = no shigella, 1 = flex, 2 = sonnei

here::i_am("code/create_figures/figure_S11_serotype.R")

no_etiology <- readRDS("results/no_etiology/ipd_no_etiology_serotype.Rds")

plot_data <- expand.grid(
  abx = c(
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics",
    "Possibly effective antibiotics"
  ),
  subgroup = c("S. flexneri", "S. sonnei")
)

plot_data$pt_est <- c(
  no_etiology$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
  no_etiology$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level
)

plot_data$se <- no_etiology$aipw_est$results_object$se[10:15]

plot_data <- plot_data %>%
  mutate(
    lower_ci = pt_est - 1.96 * se,
    upper_ci = pt_est + 1.96 * se
  )

# Define custom colors for each antibiotic category
abx_colors <- c(
  "Guideline recommended antibiotics" = "#42B540FF",
  "Possibly effective antibiotics" = "#00468BFF",
  "No or ineffective antibiotics" = "#ED0000FF"
)

# Order antibiotics in legend
plot_data$abx <- factor(
  plot_data$abx,
  levels = c(
    "No or ineffective antibiotics",
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics"
  )
)

# Get x-axis scaling
global_x_min <- min(plot_data$lower_ci)
global_x_max <- max(plot_data$upper_ci)

# Add padding for annotations
annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

p <- ggplot2::ggplot(
  data = plot_data,
  aes(
    x = pt_est,
    y = subgroup,
    color = abx
  )
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 5
  ) +
  geom_errorbarh(
    aes(xmin = lower_ci, xmax = upper_ci),
    position = position_dodge(width = 0.5),
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
      label = paste0(
        round(pt_est, 2),
        " (",
        round(lower_ci, 2),
        ", ",
        round(upper_ci, 2),
        ")"
      )
    ),
    position = position_dodge(width = 0.5),
    hjust = 0,
    vjust = 0.5,
    size = 4
  ) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Serotype (No etiology)",
    color = "Treatment received"
  ) +
  scale_color_manual(values = abx_colors) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(global_x_min, global_x_max))

ggsave(
  here::here("figures/fig_S11_serotype_no_etiology.png"),
  p,
  width = 12,
  height = 6
)

# 2. Case control ------------------------------------------------------------

case_control <- readRDS("results/case_control/ipd_cc_serotype.Rds")

plot_data <- expand.grid(
  abx = c(
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics"
  ),
  subgroup = c("S. flexneri", "S. sonnei")
)

plot_data$pt_est <- c(
  case_control$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  case_control$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level
)

plot_data$se <- case_control$aipw_est$results_object$se[8:13]

plot_data <- plot_data %>%
  mutate(
    lower_ci = pt_est - 1.96 * se,
    upper_ci = pt_est + 1.96 * se
  )

# Define custom colors for each antibiotic category
abx_colors <- c(
  "Guideline recommended antibiotics" = "#42B540FF",
  "Possibly effective antibiotics" = "#00468BFF",
  "No or ineffective antibiotics" = "#ED0000FF"
)

# Order antibiotics in legend
plot_data$abx <- factor(
  plot_data$abx,
  levels = c(
    "No or ineffective antibiotics",
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics"
  )
)

# Get x-axis scaling
global_x_min <- min(plot_data$lower_ci)
global_x_max <- max(plot_data$upper_ci)

# Add padding for annotations
annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

p <- ggplot2::ggplot(
  data = plot_data,
  aes(
    x = pt_est,
    y = subgroup,
    color = abx
  )
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 5
  ) +
  geom_errorbarh(
    aes(xmin = lower_ci, xmax = upper_ci),
    position = position_dodge(width = 0.5),
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
      label = paste0(
        round(pt_est, 2),
        " (",
        round(lower_ci, 2),
        ", ",
        round(upper_ci, 2),
        ")"
      )
    ),
    position = position_dodge(width = 0.5),
    hjust = 0,
    vjust = 0.5,
    size = 4
  ) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Serotype (Case control)",
    color = "Treatment received"
  ) +
  scale_color_manual(values = abx_colors) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(global_x_min, global_x_max))

ggsave(
  here::here("figures/fig_S11_serotype_case_control.png"),
  p,
  width = 12,
  height = 6
)
