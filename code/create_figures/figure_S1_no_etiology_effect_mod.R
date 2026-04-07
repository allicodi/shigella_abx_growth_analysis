# --------------------------------------------------------------------
# Figure S1: Study-specific heterogeneity by age for no etiology
# --------------------------------------------------------------------

here::i_am("code/create_figures/figure_S1_no_etiology_effect_mod.R")

abcd_results <- readRDS(here::here("results/no_etiology/abcd_results_no_models.Rds"))
efgh_results <- readRDS(here::here("results/no_etiology/efgh_results_no_models.Rds"))
maled_results <- readRDS(here::here("results/no_etiology/maled_results_MSD_no_models.Rds"))
gems_results <- readRDS(here::here("results/no_etiology/gems_results_no_models.Rds"))
vida_results <- readRDS(here::here("results/no_etiology/vida_results_no_models.Rds"))

abcd_msm <- abcd_results$aipw_est$results_object$aipw_msm
efgh_msm <- efgh_results$aipw_est$results_object$aipw_msm
gems_msm <- gems_results$aipw_est$results_object$aipw_msm
vida_msm <- vida_results$aipw_est$results_object$aipw_msm
maled_msm <- maled_results$aipw_est$results_object$aipw_msm

# ipd_msm <- ipd_results$aipw_est$results_object$aipw_msm

who_abcd <- abcd_msm$msm_Azithromycin
who_efgh <- efgh_msm$msm_Guideline_recommended_antibiotics
who_gems <- gems_msm$msm_Guideline_recommended_antibiotics
who_vida <- vida_msm$msm_Guideline_recommended_antibiotics
who_maled <- maled_msm$msm_Guideline_recommended_antibiotics

maybe_efgh <- efgh_msm$msm_Possibly_effective_antibiotics
maybe_gems <- gems_msm$msm_Possibly_effective_antibiotics
maybe_vida <- vida_msm$msm_Possibly_effective_antibiotics
maybe_maled <- maled_msm$msm_Possibly_effective_antibiotics

no_ineff_abcd <- abcd_msm$msm_No_azithromycin
no_ineff_efgh <- efgh_msm$msm_No_or_ineffective_antibiotics
no_ineff_gems <- gems_msm$msm_No_or_ineffective_antibiotics
no_ineff_vida <- vida_msm$msm_No_or_ineffective_antibiotics
no_ineff_maled <- maled_msm$msm_No_or_ineffective_antibiotics

age_range_abcd <- seq(2, 24)
age_range_efgh <- seq(5, 35)
age_range_gems <- seq(0, 59)
age_range_vida <- seq(0, 59)
age_range_maled <- seq(0, 24)

msm_age_plot_data <- function(age, coefs, cov_matrix, study_name){
  X <- cbind(1, age)
  est <- X %*% coefs
  se <- sqrt(rowSums((X %*% cov_matrix) * X))
  lower_ci <- est - 1.96*se
  upper_ci <- est + 1.96*se
  data.frame(age, est, se, lower_ci, upper_ci, study_name)
}

# WHO Approved Abx
plot_abcd <- msm_age_plot_data(age = age_range_abcd, 
                               coefs = who_abcd$coefficients, 
                               cov_matrix = who_abcd$cov,
                               study_name = "ABCD")
plot_efgh <- msm_age_plot_data(age = age_range_efgh,
                               coefs = who_efgh$coefficients,
                               cov_matrix = who_efgh$cov,
                               study_name = "EFGH")
plot_gems <- msm_age_plot_data(age = age_range_gems,
                               coefs = who_gems$coefficients,
                               cov_matrix = who_gems$cov,
                               study_name = "GEMS")
plot_vida <- msm_age_plot_data(age = age_range_vida,
                               coefs = who_vida$coefficients,
                               cov_matrix = who_vida$cov,
                               study_name = "VIDA")
plot_maled <- msm_age_plot_data(age = age_range_maled,
                                coefs = who_maled$coefficients,
                                cov_matrix = who_maled$cov,
                                study_name = "MALED")

plot_all_who <- rbind(plot_gems,
                      plot_maled,
                      plot_vida,
                      plot_abcd,
                      plot_efgh)

ggplot(plot_all_who, aes(x = age, y = est, color = study_name, fill = study_name)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.1) +
  labs(title = "MSM Results - WHO Approved Abx", y = "Estimate", x = "Age", color = "Study name", fill = "Study name") +
  theme_minimal()

# Maybe effective Abx
plot_efgh_maybe <- msm_age_plot_data(age = age_range_efgh,
                                     coefs = maybe_efgh$coefficients,
                                     cov_matrix = maybe_efgh$cov,
                                     study_name = "EFGH")
plot_gems_maybe <- msm_age_plot_data(age = age_range_gems,
                                     coefs = maybe_gems$coefficients,
                                     cov_matrix = maybe_gems$cov,
                                     study_name = "GEMS")
plot_vida_maybe <- msm_age_plot_data(age = age_range_vida,
                                     coefs = maybe_vida$coefficients,
                                     cov_matrix = maybe_vida$cov,
                                     study_name = "VIDA")
plot_maled_maybe <- msm_age_plot_data(age = age_range_maled,
                                      coefs = maybe_maled$coefficients,
                                      cov_matrix = maybe_maled$cov,
                                      study_name = "MALED")

plot_all_maybe <- rbind(plot_gems_maybe,
                        plot_maled_maybe,
                        plot_vida_maybe,
                        plot_efgh_maybe)

ggplot(plot_all_maybe, aes(x = age, y = est, color = study_name, fill = study_name)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.1) +
  labs(title = "MSM Results - Maybe Effective Abx", y = "Estimate", x = "Age", color = "Study name", fill = "Study name") +
  theme_minimal()

# No / Ineffective Abx
plot_abcd_no <- msm_age_plot_data(age = age_range_abcd, 
                                  coefs = no_ineff_abcd$coefficients, 
                                  cov_matrix = no_ineff_abcd$cov,
                                  study_name = "ABCD")
plot_efgh_no <- msm_age_plot_data(age = age_range_efgh,
                                  coefs = no_ineff_efgh$coefficients,
                                  cov_matrix = no_ineff_efgh$cov,
                                  study_name = "EFGH")
plot_gems_no <- msm_age_plot_data(age = age_range_gems,
                                  coefs = no_ineff_gems$coefficients,
                                  cov_matrix = no_ineff_gems$cov,
                                  study_name = "GEMS")
plot_vida_no <- msm_age_plot_data(age = age_range_vida,
                                  coefs = no_ineff_vida$coefficients,
                                  cov_matrix = no_ineff_vida$cov,
                                  study_name = "VIDA")
plot_maled_no <- msm_age_plot_data(age = age_range_maled,
                                   coefs = no_ineff_maled$coefficients,
                                   cov_matrix = no_ineff_maled$cov,
                                   study_name = "MALED")

plot_all_no <- rbind(plot_efgh_no,
                     plot_gems_no,
                     plot_maled_no,
                     plot_vida_no,
                     plot_abcd_no,
                     plot_efgh_no)

ggplot(plot_all_no, aes(x = age, y = est, color = study_name, fill = study_name)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.1) +
  labs(title = "MSM Results - No / Ineffective Abx", y = "Estimate", x = "Age", color = "Study name", fill = "Study name") +
  theme_minimal()

# Define a fixed color palette for all studies
# study_colors <- c("GEMS" = "#42B54099", 
#                   "MALED" = "#ED000099",
#                   "VIDA" = "#925E9F99", 
#                   "ABCD" = "#00468B99",
#                   "EFGH" = "#AD002A99")

study_colors <- c("GEMS" = "#42B540FF", 
                  "MALED" = "#ED0000FF",
                  "VIDA" = "#00468BFF", 
                  "ABCD" = "#0099B4FF",
                  "EFGH" = "#925E9FFF")


# Combine all datasets and add a category column
plot_all_who$category <- "Guideline recommended antibiotics"
plot_all_maybe$category <- "Possibly effective antibiotics"
plot_all_no$category <- "No or ineffective antibiotics"

plot_all <- bind_rows(plot_all_who, plot_all_maybe, plot_all_no)

# Set factor levels to control facet order
plot_all$category <- factor(plot_all$category, levels = c("No or ineffective antibiotics", 
                                                          "Possibly effective antibiotics", 
                                                          "Guideline recommended antibiotics"))

plot_all$study_name <- factor(plot_all$study_name, levels = names(study_colors))

msm_all_plot <- ggplot(plot_all, aes(x = age, y = est, color = study_name, fill = study_name)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.15, color = NA) +
  labs(y = "HAZ difference (95% CI)", x = "Age (months)", 
       color = "Study name", fill = "Study name") +
  facet_wrap(~category, scales = "free_x") +
  scale_color_manual(values = study_colors) +
  scale_fill_manual(values = study_colors) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                               # remove grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),  # add border to each panel
    axis.line = element_blank(),                                # don't rely on axis.line
    axis.ticks = element_line(color = "black", size = 0.2)      # add back thin tick marks
  )

ggsave(here::here("figures/fig_S1_no_etiology_effect_mod_by_study.png"), msm_all_plot,height = 6, width = 8, dpi = 300)
