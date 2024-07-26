

#Libraries
library(shiny)
library(DT)
library(bslib)
library(RSQLite)

# Sample data
#proteins <- data.frame(
#  ProteinID = 1:4,
 # ProteinName = paste("Protein", 1:4)
#)

#peptides <- data.frame(
#  ProteinID = rep(1:4, each = 4),
 # PeptideName = paste("Peptide", 1:4)
#)

#data.frames 
import_raw<- function(peptide_file, protein_file, db_file) {
  peptide <- read.table (peptide_file, header = TRUE, sep = "\t")
  protein <- read.table(protein_file, header = TRUE, sep = "\t")
  
  #new SQLites database
  con <- dbConnect(SQLite(), dbname = db_file)
  
  #write the data to the database 
  dbWriteTable(con, 'peptide', peptide, overwrite = TRUE)
  dbWriteTable(con, 'protein', peptide, overwrite = TRUE)
  
  #Close the connection
  dbDisconnect(con)
}


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
    datatable(protein, selection = 'single', options = list(pageLength = 4))
  })
  
  output$peptideTable <- renderDT({
    req(input$proteinTable_rows_selected)
    selectedProteinID <- proteins$ProteinID[input$proteinTable_rows_selected]
    selectedPeptides <- peptides[peptides$ProteinID == selectedProteinID, ]
    datatable(selectedPeptides, selection = 'single', options = list(pageLength = 4))
  })
  
}

shinyApp(ui, server)

