
boxPlot2 <- function(data) {
  
  ##Box Plot 
  long_data <- data |>
    
    pivot_longer(everything(), names_to = 'sample', values_to = 'intensity') |>
    mutate(sample = str_replace(sample, 'Reporter.intensity.', ''),
           intensity = intensity + 1)
  
  box_plot2 <- ggplot(long_data, aes(x = sample, y = intensity, group = sample)) +
    geom_boxplot() +
    scale_y_log10() +  
    theme_minimal() +
    labs(title = "Box Plot 2", x = "Sample", y = "Intensity (Log Scale)") +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold"))
  
  return(box_plot2)
}

