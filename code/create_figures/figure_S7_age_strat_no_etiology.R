# ----------------------------------------------------------------------------
# Create figures for age-stratified no etiology analyses
# ----------------------------------------------------------------------------

library(ggplot2)

here::i_am("code/create_figures/figure_S7_age_strat_no_etiology.R")

no_etiology_0_11 <- readRDS(here::here("results/no_etiology/ipd_no_etiology_0_11_no_models.Rds"))
no_etiology_12_23 <- readRDS(here::here("results/no_etiology/ipd_no_etiology_12_23_no_models.Rds"))
no_etiology_24_59 <- readRDS(here::here("results/no_etiology/ipd_no_etiology_24_59_no_models.Rds"))

# ------------------------------------------------------------------------------
# No etiology plot
# ------------------------------------------------------------------------------

no_etiology_results <- list(no_etiology_0_11,
                            no_etiology_12_23,
                            no_etiology_24_59)

abx_labels = NULL
inf_labels = NULL
outcome_labels = NULL
age_labels <- c("0–11 months", "12–23 months", "24–59 months")

# Make dataset with results from all studies
plot_data <- lapply(1:length(no_etiology_results), function(i) {
  res <- no_etiology_results[[i]]
  
  # Get labels for plot if not passed in
  if(is.null(abx_labels)){
    abx_label <- res$parameters$abx_var_name
  } else {
    abx_label <- abx_labels[[i]]
  }
  
  if(is.null(inf_labels)) {
    inf_label <- res$parameters$infection_var_name
    # if no inf_var_name, use case_var_name (case control)
    if(is.null(inf_label)){
      inf_label <- res$parameters$case_var_name
    }
  } else{
    inf_label <- inf_labels[[i]]
  }
  
  if(is.null(outcome_labels)) { 
    outcome_label <- res$parameters$laz_var_name
  } else {
    outcome_label <- outcome_labels[[i]]
  }
  
  if(is.null(age_labels)){
    study_label <- i
  } else{
    study_label <- age_labels[[i]]
  }
  
  plot_label <- paste0("Antibiotics: ", abx_label, ", Infection: ", inf_label) 
  plot_data <- data.frame()
  
  # make dataframe for plotting
  if(class(res) == "agaipw_res"){
    # Plot AIPW 
    res_df <- res$aipw_est$results_object$results_df
    se_vec <- res$aipw_est$results_object$se
    
    abx_levels <- res_df$abx_levels
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      tmp <- data.frame(subgroup = abx_level,
                        pt_est = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level],
                        lower_ci = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level] - 1.96*se_vec[paste0("effect_", abx_level)],
                        upper_ci = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level] + 1.96*se_vec[paste0("effect_", abx_level)],
                        plot_label = plot_label,
                        study_label = study_label)
      
      plot_data <- rbind(plot_data, tmp)
    }
    
  } else{
    # Plot GComp
    res_df <- res$results$pt_est
    se_vec <- res$results$bootstrap_results
    
    abx_levels <- res_df$abx_levels
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      tmp <- data.frame(subgroup = abx_level,
                        pt_est = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level],
                        lower_ci = se_vec[[paste0('Abx = ', abx_level)]]$lower_ci_effect_inf_abx_level,
                        upper_ci = se_vec[[paste0('Abx = ', abx_level)]]$upper_ci_effect_inf_abx_level,
                        plot_label = plot_label,
                        study_label = study_label)
      plot_data <- rbind(plot_data, tmp)
    }
    
  }
  
  return(plot_data)
})

plot_data <- do.call(rbind, plot_data)

# Change _ to space for subgroup labels
plot_data$subgroup <- gsub("_", " ", plot_data$subgroup)

subgroup_colors <- c("Guideline recommended antibiotics" = "#42B540FF", 
                     "Possibly effective antibiotics" = "#00468BFF", 
                     "No or ineffective antibiotics" = "#ED0000FF")

plot_data$subgroup <- factor(plot_data$subgroup, levels = c("No or ineffective antibiotics",
                                                            "Possibly effective antibiotics",  
                                                            "Guideline recommended antibiotics"))
# Get xaxis scaling 
global_x_min <- min(plot_data$lower_ci)
global_x_max <- max(plot_data$upper_ci)

# Add padding for annotations
annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

plot_results <- lapply(unique(plot_data$study_label), function(x, plot_data){
  sub_data <- plot_data[plot_data$study_label == x,]
  
  ggplot2::ggplot(data = sub_data, aes(x = pt_est, y = subgroup, color = subgroup)) +
    geom_point(size = 5) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#ED0000FF") +
    geom_text(aes(x = x_loc,
                  label = paste0(round(pt_est, 2),
                                 " (", round(lower_ci, 2), ", ", round(upper_ci, 2), ")")),
              hjust = 0, vjust = 0.5, size = 4) +
    labs(title = sub_data$study_label[1],
         x = "HAZ difference (95% CI)",
         y = "Treatment received") +
    scale_color_manual(values = subgroup_colors) +
    theme_minimal(base_size = 14) +
    labs(color = "Treatment received") +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    coord_cartesian(xlim = c(global_x_min, global_x_max))
  
}, plot_data = plot_data)

# Combine plots with patchwork
combined_plot <- patchwork::wrap_plots(plot_results, ncol = 1, axis_titles = 'collect') & 
  theme(legend.position = "none") 

ggsave(filename = here::here("results/final/figures/ipd_no_etiology_results_age_stratified.png"), plot = combined_plot, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# Case-control plot
# ------------------------------------------------------------------------------

cc_results <- list(cc_0_11,
                   cc_12_23,
                   cc_24_59)

abx_labels = NULL
inf_labels = NULL
outcome_labels = NULL
age_labels <- c("0–11 months", "12–23 months", "24–59 months")

# Make dataset with results from all studies
plot_data <- lapply(1:length(cc_results), function(i) {
  res <- cc_results[[i]]
  
  # Get labels for plot if not passed in
  if(is.null(abx_labels)){
    abx_label <- res$parameters$abx_var_name
  } else {
    abx_label <- abx_labels[[i]]
  }
  
  if(is.null(inf_labels)) {
    inf_label <- res$parameters$infection_var_name
    # if no inf_var_name, use case_var_name (case control)
    if(is.null(inf_label)){
      inf_label <- res$parameters$case_var_name
    }
  } else{
    inf_label <- inf_labels[[i]]
  }
  
  if(is.null(outcome_labels)) { 
    outcome_label <- res$parameters$laz_var_name
  } else {
    outcome_label <- outcome_labels[[i]]
  }
  
  if(is.null(age_labels)){
    study_label <- i
  } else{
    study_label <- age_labels[[i]]
  }
  
  plot_label <- paste0("Antibiotics: ", abx_label, ", Infection: ", inf_label) 
  plot_data <- data.frame()
  
  # make dataframe for plotting
  if(class(res) == "agaipw_res"){
    # Plot AIPW 
    res_df <- res$aipw_est$results_object$results_df
    se_vec <- res$aipw_est$results_object$se
    
    abx_levels <- res_df$abx_levels
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      tmp <- data.frame(subgroup = abx_level,
                        pt_est = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level],
                        lower_ci = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level] - 1.96*se_vec[paste0("effect_", abx_level)],
                        upper_ci = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level] + 1.96*se_vec[paste0("effect_", abx_level)],
                        plot_label = plot_label,
                        study_label = study_label)
      
      plot_data <- rbind(plot_data, tmp)
    }
    
  } else{
    # Plot GComp
    res_df <- res$results$pt_est
    se_vec <- res$results$bootstrap_results
    
    abx_levels <- res_df$abx_levels
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      tmp <- data.frame(subgroup = abx_level,
                        pt_est = res_df$effect_inf_abx_level[res_df$abx_levels == abx_level],
                        lower_ci = se_vec[[paste0('Abx = ', abx_level)]]$lower_ci_effect_inf_abx_level,
                        upper_ci = se_vec[[paste0('Abx = ', abx_level)]]$upper_ci_effect_inf_abx_level,
                        plot_label = plot_label,
                        study_label = study_label)
      plot_data <- rbind(plot_data, tmp)
    }
    
  }
  
  return(plot_data)
})

plot_data <- do.call(rbind, plot_data)

# Change _ to space for subgroup labels
plot_data$subgroup <- gsub("_", " ", plot_data$subgroup)
plot_data$subgroup <- gsub("Maybe effective antibiotics", "Possibly effective antibiotics", plot_data$subgroup)
plot_data$subgroup <- gsub("WHO recommended antibiotics", "Guideline recommended antibiotics", plot_data$subgroup)

subgroup_colors <- c("Guideline recommended antibiotics" = "#42B540FF", 
                     "Possibly effective antibiotics" = "#00468BFF", 
                     "No or ineffective antibiotics" = "#ED0000FF")

plot_data$subgroup <- factor(plot_data$subgroup, levels = c("No or ineffective antibiotics",
                                                            "Possibly effective antibiotics",  
                                                            "Guideline recommended antibiotics"))
# Get xaxis scaling 
global_x_min <- min(plot_data$lower_ci)
global_x_max <- max(plot_data$upper_ci)

# Add padding for annotations
annotation_padding <- (global_x_max - global_x_min) * 0.5
global_x_max <- global_x_max + annotation_padding
x_loc <- global_x_max - annotation_padding / 2

plot_results <- lapply(unique(plot_data$study_label), function(x, plot_data){
  sub_data <- plot_data[plot_data$study_label == x,]
  
  ggplot2::ggplot(data = sub_data, aes(x = pt_est, y = subgroup, color = subgroup)) +
    geom_point(size = 5) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#ED0000FF") +
    geom_text(aes(x = x_loc,
                  label = paste0(round(pt_est, 2),
                                 " (", round(lower_ci, 2), ", ", round(upper_ci, 2), ")")),
              hjust = 0, vjust = 0.5, size = 4) +
    labs(title = sub_data$study_label[1],
         x = "HAZ difference (95% CI)",
         y = "Treatment received") +
    scale_color_manual(values = subgroup_colors) +
    theme_minimal(base_size = 14) +
    labs(color = "Treatment received") +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    coord_cartesian(xlim = c(global_x_min, global_x_max))
  
}, plot_data = plot_data)

# Combine plots with patchwork
combined_plot <- patchwork::wrap_plots(plot_results, ncol = 1, axis_titles = 'collect') & 
  theme(legend.position = "none") 

ggsave(filename = here::here("figures/fig_S7_no_etiology_age_stratified.png"), plot = combined_plot, width = 14, height = 8, dpi = 300)
