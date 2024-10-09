
barPlot <- function(data) {
  
  ##Box Plot 
  long_data <- data |>
    
    pivot_longer(everything(), names_to = 'sample', values_to = 'intensity') |>
    mutate(sample = str_replace(sample, 'Reporter.intensity.corrected.', ''),
           missing = is.na (intensity))
  
  bar_plot <- ggplot(long_data, aes(x = sample, fill = missing)) +
    geom_bar() +
    theme_minimal() + 
    labs(title = "Bar Plot", x = "Sample", y = "Count") +
    scale_fill_brewer(palette = 'Paired') +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold"))
  
  return(bar_plot)
  
}
