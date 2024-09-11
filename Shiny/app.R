#Libraries
library(shiny)
library(DT)
library(bslib)
library(RSQLite)
library(dplyr)
library(dbplyr)
library(ProtResDash)

# move default data to sqlite file
extdata <- system.file("extdata", package = 'ProtResDash')

import_raw ( file.path(extdata, "peptides.txt"), 
             file.path(extdata, "proteinGroups.txt"), 
             "data/protein_peptidedb.sqlite")

con <- dbConnect(SQLite(), dbname = "data/protein_peptidedb.sqlite")
peptides <- tbl(con, "peptides")
proteins <- tbl(con, "proteins")


#bslib theme
base_theme <- bs_theme(bootswatch = "quartz")

protein_theme <- bs_add_rules(base_theme,
                              ".protein-theme { 
    background-color: #f8f9fa; 
    padding: 10px; 
    border-radius: 5px; 
    box-shadow: 0 4px 8px rgba(0,0,0,0.1); 
  }"
)

peptide_theme <- bs_add_rules(base_theme,
                              ".peptide-theme { 
    background-color: #e9ecef; 
    padding: 10px; 
    border-radius: 5px; 
    box-shadow: 0 4px 8px rgba(0,0,0,0.1); 
  }"
)

#UI
ui <- div(
  style = "overflow-x: scroll; width:100%; max-width:100%;",
  navbarPage(
  theme = base_theme,
  title = "Protein and Peptide Viewer",
  
  tabPanel(
    "Home",
    fluidRow(
      column(
        width = 6,
        h4("Select a Protein:"),
        tags$div(
          class = "protein-theme",
          DTOutput("proteinTable"),
          tags$style(HTML("
            .dataTables_filter input {
              width: 50% !important;
              overflow: scroll;
            }
           .shiny-input-container{
           color: #474747;
           }
          "))
        )
      ),
      column(
        width = 6,
        h3("Peptide Description:"),
        tags$div(
          class = "peptide-theme",
          DTOutput("peptideTable"),
         
        )
      )
    )
  ),
  
  # New QC Tab
  tabPanel(
    "QC",
    fluidRow(
      column(
        width = 6,
        h3("PCA Plot"),
        plotOutput("pcaPlot")
      )
     )
    )
  )
)

server <- function(input, output, session) {
  
  output$proteinTable <- renderDT({
    req(proteins) 
    datatable(collect(proteins), selection = 'single', options = list(pageLength = 4))
  })
  
  output$peptideTable <- renderDT({
    req(input$proteinTable_rows_selected)
    
    selectedProteinID <- proteins %>%
      slice(row_number() == input$proteinTable_rows_selected) %>% 
      pull(Protein.group.IDs) #Extract the Protein.group.IDs from the selected row. 
    
    selectedPeptides <- peptides %>% 
      filter(Protein.group.IDs == !!selectedProteinID[1]) %>%
      collect()
    
    datatable(selectedPeptides, selection = 'single', options = list(pageLength = 4))
  })
  
  output$pcaPlot <- renderPlot({
    pca_data <- proteins %>%
      select(starts_with('Reporter.intensity.corrected')) %>%
      collect() %>%
      na.omit()
    
    # Check if there is enough data for PCA (more than one row and column)
    if (nrow(pca_data) > 1 && ncol(pca_data) > 1) {
      
      pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
      pca_df <- as.data.frame(pca_res$x[, 1:2]) # Extract the first two principal components and convert to a data frame
      
      # Create a group vector with the correct length
      pca_df$Group <- rep(c('Group 1', 'Group 2'), length.out = nrow(pca_df))
      
      # Scatter plot
      plot(pca_df$PC1, pca_df$PC2, col = as.factor(pca_df$Group),
           xlab = "Group 1",
           ylab = "Group 2",
           main = "PCA Plot",
           pch = 19)
      legend("topright", legend = levels(as.factor(pca_df$Group)), 
             col = 1:2, pch = 19)
      
    } 
    
  }, width = "auto", height = "auto", res = 72)
}

shinyApp(ui, server)





