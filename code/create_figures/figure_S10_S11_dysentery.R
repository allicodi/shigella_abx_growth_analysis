# ---------------------------------------------------
# Dysentery stratified 
# ---------------------------------------------------

here::i_am("code/create_figures/figure_S10_S11_dysentery.R")

library(ggplot2)
library(tidyverse)

# 1. Case control  ------------------------------------------------------------

### Age specific results
abx_colors <- c("Guideline recommended antibiotics" = "#0072b2", 
                     "Possibly effective antibiotics" = "#cc79a7", 
                     "No or ineffective antibiotics" = "#d55e00")


res_0_11 <- readRDS(here::here("results/case_control/ipd_cc_dysentery_0_11.Rds"))
res_12_23 <- readRDS(here::here("results/case_control/ipd_cc_dysentery_12_23.Rds"))
res_24_59 <- readRDS(here::here("results/case_control/ipd_cc_dysentery_24_59.Rds")) # note the abx are all in same orders in results df here

plot_data <- expand.grid(
  abx = c(
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  age = c("0-11 months", "12-23 months", "24-59 months")
)

plot_data$pt_est <- c(
  res_0_11$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  res_0_11$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level,
  
  res_12_23$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  res_12_23$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level,
  
  res_24_59$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  res_24_59$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level
)

plot_data$se <- c(res_0_11$aipw_est$results_object$se[8:13],
                  res_12_23$aipw_est$results_object$se[8:13],
                  res_24_59$aipw_est$results_object$se[8:13])

plot_data <- plot_data %>%
  mutate(
    lower_ci = pt_est - 1.96 * se,
    upper_ci = pt_est + 1.96 * se,
    abx = factor(
      abx,
      levels = c(
        "No or ineffective antibiotics",
        "Possibly effective antibiotics",
        "Guideline recommended antibiotics"
      )
    ),
    age = factor(age, levels = c("0-11 months", "12-23 months", "24-59 months"))
  )

global_x_min <- min(plot_data$lower_ci, na.rm = TRUE)
global_x_max <- max(plot_data$upper_ci, na.rm = TRUE)

annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

p_age <- ggplot2::ggplot(
  data = plot_data,
  aes(
    x = pt_est,
    y = subgroup,
    color = abx
  )
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 4
  ) +
  geom_errorbarh(
    aes(xmin = lower_ci, xmax = upper_ci),
    position = position_dodge(width = 0.5),
    height = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#d55e00"
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
    size = 3.3
  ) +
  facet_wrap(~ age, ncol = 1) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Dysentery subgroup",
    color = "Treatment received"
  ) +
  scale_color_manual(values = abx_colors) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(global_x_min, global_x_max))

ggsave(
  here::here("figures/fig_S10_dysentery_age_specific_cc.png"),
  p_age,
  width = 12,
  height = 10
)

# 2. No etiology ------------------------------------------------------------

### Age specific results
abx_colors <- c("Guideline recommended antibiotics" = "#0072b2", 
                "Possibly effective antibiotics" = "#cc79a7", 
                "No or ineffective antibiotics" = "#d55e00")


res_0_11 <- readRDS(here::here("results/no_etiology/ipd_no_etiology_dysentery_0_11.Rds"))
res_12_23 <- readRDS(here::here("results/no_etiology/ipd_no_etiology_dysentery_12_23.Rds"))
res_24_59 <- readRDS(here::here("results/no_etiology/ipd_no_etiology_dysentery_24_59.Rds")) # note the abx are all in same orders in results df here

plot_data <- expand.grid(
  abx = c(
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics",
    "Possibly effective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  age = c("0-11 months", "12-23 months", "24-59 months")
)

plot_data$pt_est <- c(
  res_0_11$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
  res_0_11$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level,
  
  res_12_23$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
  res_12_23$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level,
  
  res_24_59$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
  res_24_59$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level
)

plot_data$se <- c(res_0_11$aipw_est$results_object$se[10:15],
                  res_12_23$aipw_est$results_object$se[10:15],
                  res_24_59$aipw_est$results_object$se[10:15])

plot_data <- plot_data %>%
  mutate(
    lower_ci = pt_est - 1.96 * se,
    upper_ci = pt_est + 1.96 * se,
    abx = factor(
      abx,
      levels = c(
        "No or ineffective antibiotics",
        "Possibly effective antibiotics",
        "Guideline recommended antibiotics"
      )
    ),
    age = factor(age, levels = c("0-11 months", "12-23 months", "24-59 months"))
  )

global_x_min <- min(plot_data$lower_ci, na.rm = TRUE)
global_x_max <- max(plot_data$upper_ci, na.rm = TRUE)

annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

p_age <- ggplot2::ggplot(
  data = plot_data,
  aes(
    x = pt_est,
    y = subgroup,
    color = abx
  )
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 4
  ) +
  geom_errorbarh(
    aes(xmin = lower_ci, xmax = upper_ci),
    position = position_dodge(width = 0.5),
    height = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#d55e00"
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
    size = 3.3
  ) +
  facet_wrap(~ age, ncol = 1) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Dysentery subgroup",
    color = "Treatment received"
  ) +
  scale_color_manual(values = abx_colors) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(global_x_min, global_x_max))

ggsave(
  here::here("figures/fig_S11_dysentery_age_specific_no_etiology.png"),
  p_age,
  width = 12,
  height = 10
)


