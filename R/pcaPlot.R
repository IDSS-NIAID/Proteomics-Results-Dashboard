#' pcaPlot
#' Generate a PCA plot for visual inspection of clusters
#' 
#' @param data A data frame containing the data to be plotted. All columns will be passed to `prcomp` and should be numeric.
#' @param group A character vector containing the group information
#' 
#' @return A ggplot object
#' @export
#' @importFrom stats prcomp var
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal scale_color_brewer
pcaPlot <- function(data, group) {
  
  # for those pesky no visible binding warnings
  if(FALSE)
    PC1 <- PC2 <- NULL
  
  # transpose if necessary
  if(nrow(data) != length(group)) {
    data <- t(data)
    if(nrow(data) != length(group)) {
      stop("Number of samples and length of group are not equal.")
    }
  }
  
  # remove 0 variance columns
  colvars <- apply(data, 2, var)
  
  data <- data[, colvars > 0]
  
  # Perform PCA
  pca <- prcomp(data, scale = TRUE)
  
  # Extract the PCA scores
  pca_scores <- as.data.frame(pca$x) |>
    mutate(group = group)
  
  # Plot the PCA
  pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2")
  
  return(pca_plot)
}