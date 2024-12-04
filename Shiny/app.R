
library(shiny)
library(DT)
library(bslib)
library(RSQLite)
library(dplyr)
library(dbplyr)
library(ProtResDash)
library(tidyr)
library(stringr)
library(ggplot2)
library(yaml)

# YAML configuration
config <- yaml::read_yaml("config.yaml")

# Box plot
file.path(here::here(), "R", "boxPlot.R") |> source()
# Box plot
file.path(here::here(), "R", "boxPlot2.R") |> source()
# Bar plot
file.path(here::here(), "R", "barPlot.R") |> source()

# Move default data to SQLite file if it doesn't exist
extdata <- system.file("extdata", package = config$database$extdata$package)

if (!file.exists(config$database$path)) {
  import_raw(file.path(extdata, config$database$extdata$peptides_file),
             file.path(extdata, config$database$extdata$protein_groups_file),
             config$database$path)
}

# Database connection
con <- dbConnect(SQLite(), dbname = config$database$path)
peptides <- tbl(con, "peptides")
proteins <- tbl(con, "proteins")

# bslib theme
base_theme <- bs_theme(bootswatch = config$app$theme)

protein_theme <- bs_add_rules(base_theme, paste0(
  ".protein-theme {",
  paste(
    paste0(names(config$app$base_theme_rules$protein), ": ", 
                 config$app$base_theme_rules$protein, ";"),
    collapse = " "
  ),
  "}"
))

peptide_theme <- bs_add_rules(base_theme, paste0(
  ".peptide-theme {",
  paste(
    paste0(names(config$app$base_theme_rules$peptide), ": ", 
           config$app$base_theme_rules$peptide, ";"),
    collapse = " "
  ),
  "}"
))

# UI
ui <- div(
  style = "overflow-x: scroll; width:100%; max-width:100%;",
  navbarPage(
    theme = base_theme,
    title = config$app$title,
    
    # Home Tab
    tabPanel(
      config$ui$tabs[[1]]$name,
      fluidRow(
        column(
          width = config$ui$tabs[[1]]$layout$fluidRow[[1]]$column,
          h4("Select a Protein:"),
          tags$div(
            class = "protein-theme",
            DTOutput("proteinTable")
          )
        )
      ),
      fluidRow(
        column(
          width = config$ui$tabs[[1]]$layout$fluidRow[[2]]$column,
          h3("Peptide Description:"),
          tags$div(
            class = "peptide-theme",
            DTOutput("peptideTable")
          )
        )
      )
    ),
    
    # QC Tab
    tabPanel(
      config$ui$tabs[[2]]$name,
      fluidRow(
        lapply(config$plots, function(plot) {
          column(
            width = 6,
            h3(plot$type),
            plotOutput(plot$output)
          )
        })
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Protein Table
  output$proteinTable <- renderDT({
    req(proteins)
    
    protein_data <- collect(select(proteins, all_of(config$tables$protein_table$columns)))
    datatable(protein_data, selection = 'single', options = list(pageLength = config$tables$protein_table$page_length))
  })
  
  # Peptide Table
  output$peptideTable <- renderDT({
    req(input$proteinTable_rows_selected)
    
    selected_row <- input$proteinTable_rows_selected
    selectedProteinID <- proteins %>%
    filter(row_number() == selected_row) %>%
    pull(Protein.group.IDs) # Extract the Protein.group.IDs from the selected row
    
    selectedPeptides <- peptides %>%
    filter(Protein.group.IDs == !!selectedProteinID[1]) %>%
    collect() # Collect the data as a df
    
    datatable(select(selectedPeptides, all_of(config$tables$peptide_table$columns)),
    selection = 'single', options = list(pageLength = config$tables$peptide_table$page_length))
  })
  
  # PCA Plot
  output$pcaPlot <- renderPlot({
    data <- proteins %>%
      select(starts_with(config$plots$pca_plot$starts_with)) %>%
      select(config$plots$pca_plot$range[1]:config$plots$pca_plot$range[2]) |>
      collect()
    
    if (nrow(data) > 1 && ncol(data) > 1) {
      print("Not enough data for PCA")
      pca_res <- prcomp(data, center = TRUE, scale. = TRUE)
      pca_df <- as.data.frame(pca_res$x[, 1:2]) # Here it extract the first 2 principal components and convert it to a df
      
      pcaPlot(data, 
              rep(c("Control", "IgM"), each = 5))
    } else {
      print("Enough data for PCA")
      plot(NA, NA, xlim = 0:1, ylim = 0:1, type = "n", xlab = "",
           main = "Not enough data for PCA")
    }
  })
  
  # Box Plot: Normalized
  output$boxPlot <- renderPlot({
    data <- proteins %>%
      select(starts_with(config$plots$box_plot$starts_with)) %>%
      select(config$plots$box_plot$range[1]:config$plots$box_plot$range[2]) |>
      collect()
    
    boxPlot(data)
  })
  
  # Box Plot 2: No-normalized 
  output$boxPlot2 <- renderPlot({
    data <- proteins %>%
      select(starts_with(config$plots$box_plot2$starts_with)) %>%
      select(config$plots$box_plot2$range[1]:config$plots$box_plot2$range[2]) |>
      collect()
    
    boxPlot2(data)
  })
  
  # Bar Plot
  output$barPlot <- renderPlot({
    data <- proteins %>%
      select(starts_with(config$plots$bar_plot$starts_with)) %>%
      select(config$plots$bar_plot$range[1]:config$plots$bar_plot$range[2]) |>
      collect()
    
    barPlot(data)
  })
}

shinyApp(ui, server)





