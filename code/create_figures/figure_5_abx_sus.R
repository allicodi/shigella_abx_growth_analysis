# ---------------------------------------------------------------------------
# Create figure for antibiotic susceptibility in EFGH
# ---------------------------------------------------------------------------

here::i_am("code/create_figures/figure_5_abx_sus.R")

# iteration 1 same - we care abt no etiology, WHO approved
# iteration 2 - keep shigella attr same, no_etiology = I(resistant to at least one WHO approved thing == 1) -- abx var = (0,1,2) where (I/R, S, NA to the drug you got)
#               then take the 'no etiology' aka resistant shigella output, no etiology + abx = resistant
#               compare to no etiology, WHO approved (iteration 1)
# iteration 3 - keep shigella attr same, no_etiology = I(susceptible to at least one WHO approved thing == 1) -- abx var = (0,1,2) where (I/R, S, NA to the drug you got)

efgh_results <- readRDS(here::here("results/no_etiology/susceptibility/no_etiology.Rds"))
efgh_results_IR <- readRDS(here::here("results/no_etiology/susceptibility/resistant_results.Rds"))
efgh_results_S <- readRDS(here::here("results/no_etiology/susceptibility/susceptible_results.Rds"))

main_result_aipw <- efgh_results$aipw_est$results_object$results_df
no_etiology_who_drug <- main_result_aipw$abx_level_inf_0[main_result_aipw$abx_levels == "Guideline_recommended_antibiotics"]

main_result_eifs <- efgh_results$aipw_est$results_object$eif_matrix

# Contrast 1: Resistant vs no etiology
# Resistant Shigella (no etiology) & drug resistant to from iteration 2 vs No etiology + WHO approved abx from iteration 1

resistant_result <- efgh_results_IR$aipw_est$results_object$results_df
resistant_shigella_resistant_drug <- resistant_result$abx_level_inf_0[resistant_result$abx_levels == 2]

IR_results_eifs <- efgh_results_IR$aipw_est$results_object$eif_matrix

# Contrast 2: Susceptible vs no etiology
# Susceptible Shigella (no etiology) & drug susceptible to from iteration 2 vs No etiology + WHO approved abx from iteration 1

susceptible_result <- efgh_results_S$aipw_est$results_object$results_df
susceptible_shigella_susceptible_drug <- susceptible_result$abx_level_inf_0[susceptible_result$abx_levels == 1]

S_results_eifs <- efgh_results_S$aipw_est$results_object$eif_matrix

# Put EIFS together and get SE

eifs_effect_IR_S <- data.frame(matrix(ncol = 2, nrow = nrow(main_result_eifs)))
colnames(eifs_effect_IR_S) <- c("effect_IR", "effect_S")

# NOTE IR/S in order 0,2,1 
all_eifs <- cbind(IR_results_eifs, S_results_eifs, main_result_eifs)
colnames(all_eifs) <- c("inf_eif_0_IR","inf_eif_2_IR","inf_eif_1_IR", "no_attr_eif_0_IR","no_attr_2_IR","no_attr_1_IR", "effect_0_IR","effect_2_IR","effect_1_IR",
                        "inf_eif_0_S","inf_eif_2_S","inf_eif_1_S", "no_attr_eif_0_S","no_attr_2_S","no_attr_1_S", "effect_0_S","effect_2_S","effect_1_S",
                        "inf_eif_1_main","inf_eif_0_main","no_attr_eif_1_main","no_attr_eif_0_main","effect_1_main","effect_0_main")

gradient_IR <- c(0,0,0,0,1,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,
                 0,0,0,-1,0,0,0,0,0)

gradient_IR <- matrix(gradient_IR, ncol = 1)

eif_effect_IR <- as.numeric(as.matrix(all_eifs) %*% gradient_IR)
eifs_effect_IR_S[,1] <- eif_effect_IR

gradient_S <- c(0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,1,0,0,0,
                0,0,0,-1,0,0,0,0,0)

gradient_S <- matrix(gradient_S, ncol = 1)

eif_effect_S <- as.numeric(as.matrix(all_eifs) %*% gradient_S)
eifs_effect_IR_S[,2] <- eif_effect_S

all_eifs <- cbind(all_eifs, eifs_effect_IR_S)

cov_matrix <- stats::cov(all_eifs)
eif_hat <- sqrt(diag(cov_matrix) / nrow(all_eifs))

# Child with resistant Shigella to the drug they got vs child with no etiology with WHO recommended drug
ir_effect <- resistant_shigella_resistant_drug - no_etiology_who_drug
ir_se <- eif_hat['effect_IR']
ir_ci <- c(ir_effect - 1.96*ir_se, ir_effect + 1.96*ir_se)

# Child with susceptible Shigella to the drug they got vs child with no etiology with WHO recommended drug
s_effect <- susceptible_shigella_susceptible_drug - no_etiology_who_drug
s_se <- eif_hat['effect_S']
s_ci <- c(s_effect - 1.96*s_se, s_effect + 1.96*s_se)


# Create a dataframe with the results
subgroup_df <- data.frame(
  subgroup = c("Intermediate or Resistant", 
               "Susceptible"),
  pt_est = c(ir_effect, s_effect),
  lower_ci = c(ir_ci[1], s_ci[1]),
  upper_ci = c(ir_ci[2], s_ci[2]),
  x_loc = max(c(ir_ci[2], s_ci[2])) + 0.02  # Adjust rightward as needed
)

subgroup_df$subgroup <- factor(subgroup_df$subgroup, levels = c(
  "Susceptible", 
  "Intermediate or Resistant"
))

# Set color scheme for subgroups
subgroup_colors <- c(
  "Intermediate or Resistant" = "#ED0000FF",
  "Susceptible" = "#00468BFF"
)

# Define x-axis limits for consistent padding
global_x_min <- min(subgroup_df$lower_ci) - 0.05
global_x_max <- max(subgroup_df$upper_ci) + 0.1

# Create the plot
sir_plot <- ggplot(subgroup_df, aes(x = pt_est, y = subgroup, color = subgroup)) +
  geom_point(size = 5) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#ED0000FF") +
  geom_text(aes(x = x_loc,
                label = paste0(formatC(pt_est, format = "f", digits = 2),
                               " (", formatC(lower_ci, format = "f", digits = 2),
                               ", ", formatC(upper_ci, format = "f", digits = 2), ")")),
            hjust = 0, vjust = 0.5, size = 4) +
  labs(
    x = "HAZ difference (95% CI)",
    y = "Shigella susceptibility",
    color = "Shigella susceptibility"
  ) +
  scale_color_manual(values = subgroup_colors) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(global_x_min, global_x_max))

ggsave(here::here("figures/fig_5_abx_sus.png"), sir_plot, width = 12, height = 6, dpi = 300)

