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
          DTOutput("peptideTable")
         
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
      ),
      column(
        width = 6,
        h3("Box plot"),
        plotOutput("Box Plot")
      ),
      column(
        width = 6,
        h3("Bar Plot"),
        plotOutput("Bar Plot")
      )
     )
    )
  )
)

server <- function(input, output, session) {
  
  output$proteinTable <- renderDT({
     req(proteins)
    
    datatable(collect(select(proteins, c('Sequence','Mass','Proteins', 
                                         'Reporter.intensity.corrected.1',
                                         'Reporter.intensity.corrected.2',
                                         'Reporter.intensity.corrected.3',
                                         'Reporter.intensity.corrected.4',
                                         'Reporter.intensity.corrected.5',
                                         'Reporter.intensity.corrected.6',
                                         'Reporter.intensity.corrected.7',
                                         'Reporter.intensity.corrected.8',
                                         'Reporter.intensity.corrected.9',
                                         'Reporter.intensity.corrected.10'))),  selection = 'single', options = list(pageLength = 4))
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
  
  #PCA Plot
  if (nrow(data) > 1 && ncol(data) > 1) 
    
    {
      output$pcaPlot <- renderPlot ({
        data <- proteins %>%
          select(starts_with('Reporter.intensity.corrected')) %>%
          select(1:10) |>
          collect() 
        
        pca_res <- prcomp(data, center = TRUE, scale. = TRUE)
        pca_df <- as.data.frame(pca_res$x[, 1:2]) # Extract the first two principal components and convert to a data frame
        
        pcaPlot(data,
                rep(c('Control', 'IgM'), each = 5))
      
        
        } ) 
      
      #Box Plot
      {
        output$boxPlot <- renderPlot ({
          data <- proteins %>%
            select(starts_with('Reporter.intensity.corrected')) %>%
            select(1:10) |>
            collect() 
          
          box_res <- prcomp(data, center = TRUE, scale. = TRUE)
          box_df <- as.data.frame(box_res$x[, 1:2]) # Extract the first two principal components and convert to a data frame
          
          boxPlot(data,
                  rep(c('Control', 'IgM'), each = 5))
          
        } ) 
  
  }
  }
    
   else {
    plot(NA, NA, xlim = 0:1, ylim = 0:1, type = "n", xlab = "",
         main = "Not enough data for PCA")
      }
} 

shinyApp(ui, server) 
  






  