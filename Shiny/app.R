

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
ui <- navbarPage(
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
           plotOutput("PCA Plot")
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
       slice(row_number()== input$proteinTable_rows_selected) %>% 
       pull(Protein.group.IDs) # Extract the Protein.group.IDs from the selected row
    
    selectedPeptides <- peptides %>% 
      filter(Protein.group.IDs == !!selectedProteinID[1]) |>
      collect() # Convert the filtered result to a data frame
     
    
    datatable(selectedPeptides, selection = 'single', options = list(pageLength = 4))
  })
  
  
}

shinyApp(ui, server)
