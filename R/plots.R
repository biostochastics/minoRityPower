#' Plot Power Analysis Results
#' 
#' @description
#' Creates visualizations of power analysis results from both trial-level
#' and participant-level analyses.
#' 
#' @param results_dir Directory containing power analysis results
#' @param analysis_type Character, either "trial" or "participant"
#' @param alpha_levels Significance levels to analyze
#' @return List of ggplot objects
#' @export
plot_power_results <- function(results_dir,
                             analysis_type = c("trial", "participant"),
                             alpha_levels = c(0.05, 0.01, 0.005)) {
  
  analysis_type <- match.arg(analysis_type)
  
  # Find all result files based on analysis type
  pattern <- if(analysis_type == "trial") {
    "p_values_.*\\.txt$"
  } else {
    "bootstrap_p_values_.*\\.txt$"
  }
  
  filenames <- list.files(results_dir, pattern = pattern, full.names = TRUE)
  
  if(length(filenames) == 0) {
    stop("No power analysis results found in directory: ", results_dir)
  }
  
  library(dplyr)
  library(ggplot2)
  library(ggsci)  # For scale_color_jco()
  
  calc_avg_pval <- function(filename) {
    datme <- read.table(filename) %>% filter(is.na(V1) == F)
    
    # Extract parameters from filename
    parameters <- basename(filename) %>%
      gsub(ifelse(analysis_type == "trial", "p_values_", "bootstrap_p_values_"), "", .) %>%
      gsub(".txt", "", .) %>%
      strsplit("_", fixed = TRUE) %>%
      unlist()
    
    effect_size1 <- as.numeric(parameters[1])
    effect_size2 <- as.numeric(parameters[2])
    N <- as.numeric(parameters[3])
    
    p_values_05 <- datme %>% pull(V1) < 0.05
    p_values_01 <- datme %>% pull(V1) < 0.01
    p_values_005 <- datme %>% pull(V1) < 0.005
    
    return(data.frame(
      effect_size1 = effect_size1,
      effect_size2 = effect_size2,
      N = N,
      avg_p_val_05 = mean(p_values_05, na.rm = TRUE),
      avg_p_val_01 = mean(p_values_01, na.rm = TRUE),
      avg_p_val_005 = mean(p_values_005, na.rm = TRUE)
    ))
  }
  
  # Apply the function to each file and bind the results into a data frame
  results <- do.call(rbind, lapply(filenames, calc_avg_pval))
  
  df <- results
  df$effect_size1_percent <- paste0((df$effect_size1 - 1) * 100, "%")
  df$effect_size2_percent <- paste0((df$effect_size2 - 1) * 100, "%")
  df <- df %>% filter(effect_size2 < 4)
  df$effect_size1_percent <- factor(df$effect_size1_percent,
                                  levels = c("1%","5%","10%","25%","50%","100%"))
  df$effect_size2_percent <- factor(df$effect_size2_percent,
                                  levels = c("1%","5%","10%","25%","50%","100%"))
  
  # Melt the dataframe for plotting
  df_melt <- reshape2::melt(df, 
                           id.vars = c("effect_size1", "effect_size2", "N", 
                                     "effect_size1_percent", "effect_size2_percent"),
                           variable.name = "significance_level", 
                           value.name = "avg_p_val")
  
  # Convert significance_level to factor with nice labels
  df_melt$significance_level <- factor(df_melt$significance_level,
                                     levels = c("avg_p_val_05", "avg_p_val_01", "avg_p_val_005"),
                                     labels = c("p<0.05", "p<0.01", "p<0.005"))
  
  # Create main plot
  title_prefix <- if(analysis_type == "trial") {
    "Trial-level Analysis"
  } else {
    "Participant-level Analysis"
  }
  
  x_label <- if(analysis_type == "trial") {
    "N trials per arm"
  } else {
    "N participants per trial"
  }
  
  main_plot <- ggplot(df_melt, aes(x = N, y = avg_p_val, color = significance_level)) +
    geom_line() + 
    geom_point(alpha = 0.9, size = 0.85) +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.5, color = 'black', alpha = 0.7) +
    coord_cartesian(ylim = c(0.4, 1), expand = F) +
    scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1)) +
    labs(title = sprintf('%s: Power to Detect CARE Intervention Effects', title_prefix),
         x = x_label,
         y = "Power") +
    theme_bw() +   
    theme(
      title = element_text(size = 8),
      text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      strip.text = element_text(size = 8)
    ) +
    facet_grid(rows = vars(effect_size2_percent), 
               cols = vars(effect_size1_percent), 
               scales = "free",
               labeller = label_both) +
    scale_color_jco() +
    guides(color = guide_legend(
      title = "Type I error rate (P)",
      override.aes = list(size = 5, shape = 15),
      keywidth = 1, 
      keyheight = 1
    ))
  
  # Return plot
  return(main_plot)
}