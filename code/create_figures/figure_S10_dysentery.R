# ---------------------------------------------------
# Dysentery stratified 
# ---------------------------------------------------

here::i_am("code/create_figures/figure_S10_dysentery.R")

library(ggplot2)
library(tidyverse)

# 1. No etiology ------------------------------------------------------------
no_etiology <- readRDS("results/no_etiology/ipd_no_etiology_dysentery.Rds")

plot_data <- expand.grid(
  abx = c(
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics",
    "Possibly effective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery")
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
    y = "Dysentery subgroup (No etiology)",
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
  here::here("figures/fig_S10_dysentery_no_etiology.png"),
  p,
  width = 12,
  height = 6
)

#######################################################
# Site-specific
#######################################################

maled_res <- readRDS(here::here("results/no_etiology/maled_no_etiology_dysentery.Rds"))
gems_res <- readRDS(here::here("results/no_etiology/gems_no_etiology_dysentery.Rds"))
vida_res <- readRDS(here::here("results/no_etiology/vida_no_etiology_dysentery.Rds"))
efgh_res <- readRDS(here::here("results/no_etiology/efgh_no_etiology_dysentery.Rds"))

plot_data_maled <- expand.grid(
  abx = c(
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics",
    "Possibly effective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("MALED")
)

plot_data_gems <- expand.grid(
  abx = c(
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("GEMS")
)

plot_data_vida <- expand.grid(
  abx = c(
    "No or ineffective antibiotics",
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("VIDA")
)

plot_data_efgh <- expand.grid(
  abx = c(
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics",
    "Possibly effective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("EFGH")
)

plot_data <- rbind(plot_data_maled, plot_data_gems, plot_data_vida, plot_data_efgh)

plot_data$pt_est <- c(
 maled_res$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
 maled_res$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level,
 
 gems_res$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
 gems_res$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level,
 
 vida_res$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
 vida_res$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level,
 
 efgh_res$aipw_est$results_object$results_df$effect_inf_subinfect1_abx_level,
 efgh_res$aipw_est$results_object$results_df$effect_inf_subinfect2_abx_level
)

plot_data$se <- c(maled_res$aipw_est$results_object$se[10:15],
                  gems_res$aipw_est$results_object$se[10:15],
                  vida_res$aipw_est$results_object$se[10:15],
                  efgh_res$aipw_est$results_object$se[10:15])

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
    study = factor(study, levels = c("MALED", "GEMS", "VIDA", "EFGH"))
  )

global_x_min <- min(plot_data$lower_ci, na.rm = TRUE)
global_x_max <- max(plot_data$upper_ci, na.rm = TRUE)

annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

p_study <- ggplot2::ggplot(
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
    size = 3.3
  ) +
  facet_wrap(~ study, ncol = 2) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Dysentery subgroup (No etiology)",
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
  here::here("figures/fig_S10_dysentery_study_specific_no_etiology.png"),
  p_study,
  width = 14,
  height = 9
)

### Age specific results

abx_colors <- c(
  "Guideline recommended antibiotics" = "#42B540FF",
  "Possibly effective antibiotics" = "#00468BFF",
  "No or ineffective antibiotics" = "#ED0000FF"
)

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
  here::here("figures/fig_S10_dysentery_age_specific_no_etiology.png"),
  p_age,
  width = 12,
  height = 8
)


# 1. Case control  ------------------------------------------------------------

case_control <- readRDS("results/case_control/ipd_cc_dysentery.Rds")

plot_data <- expand.grid(
  abx = c(
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery")
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
    y = "Dysentery subgroup (Case control)",
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
  here::here("figures/fig_S10_dysentery_case_control.png"),
  p,
  width = 12,
  height = 6
)

## Study specific results

abx_colors <- c(
  "Guideline recommended antibiotics" = "#42B540FF",
  "Possibly effective antibiotics" = "#00468BFF",
  "No or ineffective antibiotics" = "#ED0000FF"
)

maled_res <- readRDS(here::here("results/case_control/maled_cc_dysentery.Rds"))
gems_res <- readRDS(here::here("results/case_control/gems_cc_dysentery.Rds"))
vida_res <- readRDS(here::here("results/case_control/vida_cc_dysentery.Rds")) # note the abx are all in different orders in results df here

plot_data_maled <- expand.grid(
  abx = c(
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics",
    "Possibly effective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("MALED")
)

plot_data_gems <- expand.grid(
  abx = c(
    "Possibly effective antibiotics",
    "Guideline recommended antibiotics",
    "No or ineffective antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("GEMS")
)

plot_data_vida <- expand.grid(
  abx = c(
    "Possibly effective antibiotics",
    "No or ineffective antibiotics",
    "Guideline recommended antibiotics"
  ),
  subgroup = c("No dysentery", "Dysentery"),
  study = c("VIDA")
)

plot_data <- rbind(plot_data_maled, plot_data_gems, plot_data_vida)

plot_data$pt_est <- c(
  maled_res$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  maled_res$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level,
  
  gems_res$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  gems_res$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level,
  
  vida_res$aipw_est$results_object$results_df$effect_case_subinfect1_abx_level,
  vida_res$aipw_est$results_object$results_df$effect_case_subinfect2_abx_level
)

plot_data$se <- c(maled_res$aipw_est$results_object$se[8:13],
                  gems_res$aipw_est$results_object$se[8:13],
                  vida_res$aipw_est$results_object$se[8:13])

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
    study = factor(study, levels = c("MALED", "GEMS", "VIDA"))
  )

global_x_min <- min(plot_data$lower_ci, na.rm = TRUE)
global_x_max <- max(plot_data$upper_ci, na.rm = TRUE)

annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

p_study <- ggplot2::ggplot(
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
    size = 3.3
  ) +
  facet_wrap(~ study, ncol = 2) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Dysentery subgroup (Case control)",
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
  here::here("figures/fig_S10_dysentery_study_specific_cc.png"),
  p_study,
  width = 14,
  height = 9
)

### Age specific results

abx_colors <- c(
  "Guideline recommended antibiotics" = "#42B540FF",
  "Possibly effective antibiotics" = "#00468BFF",
  "No or ineffective antibiotics" = "#ED0000FF"
)

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
  height = 8
)
