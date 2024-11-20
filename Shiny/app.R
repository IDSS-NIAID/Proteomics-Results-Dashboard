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
library(tidyr)
library(stringr)
library(ggplot2)
library(yaml)


#Box plot
file.path(here::here(), "R", "boxPlot.R") |> source()
#Box plot
file.path(here::here(), "R", "boxPlot2.R") |> source()
#Bar plot
file.path(here::here(), "R", "barPlot.R") |> source()

config <- yaml::read_yaml("config.yaml", eval.expr = TRUE)


# move default data to sqlite file
extdata <- system.file("extdata", package = 'ProtResDash')


#config
#if(!file.exists("data/protein_peptidedb.sqlite")) {
if(!file.exists(config$database$file)) {
  
  import_raw ( file.path(extdata, "peptides.txt"), 
               file.path(extdata, "proteinGroups.txt"), 
               config$database$file)
               #"data/protein_peptidedb.sqlite")
}
#config
#con <- dbConnect(SQLite(), dbname = "data/protein_peptidedb.sqlite")
con <- dbConnect(SQLite(), dbname = config$database$file)
peptides <- tbl(con, "peptides")
proteins <- tbl(con, "proteins")

#UI bslib theme
#config
#base_theme <- bs_theme(bootswatch = "quartz")
base_theme <- bs_theme(bootswatch = config$themes$base_theme)

protein_theme <- bs_add_rules(base_theme,
                              ".protein-theme { 
      background-color: ", config$themes$protein_theme$background_color, "; 
      padding: ", config$themes$protein_theme$padding, "; 
      border-radius: ", config$themes$protein_theme$border_radius, "; 
      box-shadow: ", config$themes$protein_theme$box_shadow, "; 
    }"
)

peptide_theme <- bs_add_rules(base_theme,
                              ".peptide-theme { 
      background-color: ", config$themes$peptide_theme$background_color, "; 
      padding: ", config$themes$peptide_theme$padding, "; 
      border-radius: ", config$themes$peptide_theme$border_radius, "; 
      box-shadow: ", config$themes$peptide_theme$box_shadow, "; 
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
  
  #QC Tab with plots
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
     )
    )
  )
)

#Server
server <- function(input, output, session) {
  
  #Protein table
  output$proteinTable <- renderDT({
     req(proteins)
    
  #   protein_data <- collect(select(proteins, c("Sequence", "Mass", "Proteins",
  #               reporter_intensity_columns("Reporter.intensity.corrected", 1:10))))
  #   datatable(protein_data, selection = 'single', options = list(pageLength = 4))
  # })
    
    #config
    protein_data <- collect(select(proteins, config$columns$protein_table))
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
    
  #config
  #   datatable(select(selectedPeptides, c("Protein.group.IDs", "Mod..peptide.IDs"),
  #           reporter_intensity_columns("Reporter.intensity.corrected", 1:10)),
  #           selection = 'single', options = list(pageLength = 4))
  # })
    datatable(select(selectedPeptides, config$columns$peptide_table), 
      selection = 'single', options = list(pageLength = 4))
  })
    
    
  #PCA Plot
  output$pcaPlot <- renderPlot ({
    data <- proteins %>%
      select(starts_with('Reporter.intensity.corrected')) %>%
      select(1:10) |>
      collect()

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
} 

shinyApp(ui, server) 
  



  