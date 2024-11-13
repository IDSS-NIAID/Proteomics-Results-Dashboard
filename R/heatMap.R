

        heatMap <- function(data_matrix) {
        
        # Color function
        col_fun <- colorRamp2(c(min(data_matrix), 
                                median(data_matrix), max(data_matrix)), c("blue", "white", "red"))
        
        # Draw heatmap
        Heatmap(
          data_matrix,
          # show_row_names = TRUE,
          # show_column_names = TRUE,
          
          name = "Intensity",
          column_title = "Heatmap",
          row_title = "Protein intensity",
          col = col_fun,
          
          #cluster_rows = TRUE,
          #cluster_columns = TRUE

        )
        }

