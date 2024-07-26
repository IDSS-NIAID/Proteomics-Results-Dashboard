

#Libraries
library(shiny)
library(DT)
library(bslib)
library(RSQLite)

library(ProtResDash)

# move default data to sqlite file
extdata <- system.file("extdata", package = 'ProtResDash')

import_raw ( file.path(extdata, "peptides.txt"), 
             file.path(extdata, "proteinGroups.txt"), 
             "data/protein_peptidedb.sqlite")


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
  )
)


server <- function(input, output, session) {
  
  output$proteinTable <- renderDT({
    datatable(proteins, selection = 'single', options = list(pageLength = 4))
  })
  
  output$peptideTable <- renderDT({
    req(input$proteinTable_rows_selected)
    selectedProteinID <- proteins$ProteinID[input$proteinTable_rows_selected]
    selectedPeptides <- peptides[peptides$ProteinID == selectedProteinID, ]
    datatable(selectedPeptides, selection = 'single', options = list(pageLength = 4))
  })
  
}

shinyApp(ui, server)

