#Libraries
library(shiny)
library(DT)
library(bslib)
library(RSQLite)
library(dplyr)
library(dbplyr)
library(ProtResDash)
library(tidyr)
library(stringr)

library(ComplexHeatmap)
library(colorRamp2)

library(tidyr)
library(stringr)
library(ggplot2)

#Box plot
file.path(here::here(), "R", "boxPlot.R") |>
  source()
#Box plot
file.path(here::here(), "R", "boxPlot2.R") |>
  source()

#Bar plot
 file.path(here::here(), "R", "barPlot.R") |>
   source()
 
#Heatmap
file.path(here::here(), "R", "heatMap.R") |>
 source()

# move default data to sqlite file
extdata <- system.file("extdata", package = 'ProtResDash')

if(!file.exists("data/protein_peptidedb.sqlite")) {
  import_raw ( file.path(extdata, "peptides.txt"), 
               file.path(extdata, "proteinGroups.txt"), 
               "data/protein_peptidedb.sqlite")
}

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
      )),
    fluidRow(
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
  
  #QC Tab with plots: PCA, Box, and Bar
  tabPanel(
    "QC",
    fluidRow(
      column(
        width = 6,
        h3("Box plot: Normalized"),
        plotOutput("boxPlot")
      ),
      column(
        width = 6,
        h3("Box plot 2: No-normalized"),
        plotOutput("boxPlot2")
      ),
      column(
        width = 6,
        h3("PCA Plot"),
        plotOutput("pcaPlot")
      ),
      column(
        width = 6,
        h3("Bar Plot"),
        plotOutput("barPlot")
      ),
       column(
         width = 6,
         h3("Heatmap"),
         plotOutput("heatMap")
       )
     )
    )
  )
)

#Server
server <- function(input, output, session) {
  
  #Protein table
  output$proteinTable <- renderDT({
     req(proteins)
    
    #remember to generalize this part
    protein_data <- collect(select(proteins, c('Sequence','Mass','Proteins', 
                                         'Reporter.intensity.corrected.1',
                                         'Reporter.intensity.corrected.2',
                                         'Reporter.intensity.corrected.3',
                                         'Reporter.intensity.corrected.4',
                                         'Reporter.intensity.corrected.5',
                                         'Reporter.intensity.corrected.6',
                                         'Reporter.intensity.corrected.7',
                                         'Reporter.intensity.corrected.8',
                                         'Reporter.intensity.corrected.9',
                                         'Reporter.intensity.corrected.10')))  
    datatable(protein_data, selection = 'single', options = list(pageLength = 4))
  })
  
  
  
  #Peptide table
  output$peptideTable <- renderDT({
    req(input$proteinTable_rows_selected)
    
    
      selected_row <- input$proteinTable_rows_selected
      
      selectedProteinID <- proteins %>%
      filter(row_number() == selected_row) %>% 
      pull(Protein.group.IDs) #Extract the Protein.group.IDs from the selected row.
      
    
    selectedPeptides <- peptides %>%  #remember to generalize this part
      filter(Protein.group.IDs == !!selectedProteinID[1]) %>%
      collect() #Collect the data as a df
      datatable(select(selectedPeptides, c('Protein.group.IDs',
                                           'Reporter.intensity.corrected.1',
                                           'Reporter.intensity.corrected.2',
                                           'Reporter.intensity.corrected.3',
                                           'Reporter.intensity.corrected.4',
                                           'Reporter.intensity.corrected.5',
                                           'Reporter.intensity.corrected.6',
                                           'Reporter.intensity.corrected.7',
                                           'Reporter.intensity.corrected.8',
                                           'Reporter.intensity.corrected.9',
                                           'Reporter.intensity.corrected.10',
                                           'Mod..peptide.IDs')), 
                
  selection = 'single', options = list(pageLength = 4))
  })
  
  #PCA Plot#
  output$pcaPlot <- renderPlot ({
    data <- proteins %>%
      select(starts_with('Reporter.intensity.corrected')) %>%
      select(1:10) |>
      collect()  #generalize 
        
    if (nrow(data) > 1 && ncol(data) > 1) 
    {
      pca_res <- prcomp(data, center = TRUE, scale. = TRUE)
      pca_df <- as.data.frame(pca_res$x[, 1:2]) # Extract the first two principal components and convert to a data frame
        
      pcaPlot(data,
              rep(c('Control', 'IgM'), each = 5))
    }else{
      plot(NA, NA, xlim = 0:1, ylim = 0:1, type = "n", xlab = "",
           main = "Not enough data for PCA")
    }
  }) 
      
    #Box Plot 2: No-normalized
  output$boxPlot2 <- renderPlot ({
    data <- proteins %>%
      select(starts_with('Reporter.intensity')) %>%
      select(11:20) |>
      collect() 
          
    boxPlot2(data)
  } ) 
  
  #Box Plot: Normalized
  output$boxPlot <- renderPlot ({
    data <- proteins %>%
      select(starts_with('Reporter.intensity.corrected')) %>%
      select(1:10) |>
      collect() 
    
    boxPlot(data)
  } ) 
  
  #Bar Plot#
  output$barPlot <- renderPlot ({
    data <- proteins %>%
      select(starts_with('Reporter.intensity.corrected')) %>%
      select(1:10) |>
      collect() 
    
    barPlot(data)
  } ) 
  
  #ComplexHeatmap
  output$heatMap <- renderPlot({
   data_matrix <- proteins %>%
    select(starts_with('Reporter.intensity.corrected')) %>%
    select(1:10) |>
    collect() |> 
    as.matrix()
   
   sub_mat <- data_matrix[1:10, 1:10] #take first 10 rows & col from data_matrix
   #we convert the data -> matrix, then we create a sub_matrix -> display heatmap of sub_matrix
    
    heatMap(sub_mat)  
  })

} 

shinyApp(ui, server) 
  



  