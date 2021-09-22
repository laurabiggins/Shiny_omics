library(shiny)
library(tidyverse)
library(plotly)
library(DT)


dataset <- readRDS("data/meta.rds") #used for the DT - searching
data_long <- readRDS("data/acid_long.rds")
acid_pval_fc <- readRDS("data/acid_pval_fc.rds")
genes_long <- readRDS("data/genes_long.rds")
#acid_pval_fc <- readRDS("data/pval_check_temp.rds")

volcano_dataset <- acid_pval_fc

table_data <- dataset %>%
    rowid_to_column()


plot_colours <- c("red3", "blue3", "green3", "orange3")

conditions <- c("Naive", "Naive+PRC2i", "Primed", "Primed+PRCi")

acid_plot_height <- "250px"

datatypes <- c(
  "Gene expression" = "gene_expr", 
  "Acid extractome protein abundance" = "acid_protein",
  "Chromatin-associated protein abundance" = "chr_protein",
  "Histone modification abundance" = "histones"
)

ui <- fluidPage(

  shinyalert::useShinyalert(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "omics.css")
  ),
  
  # Application title
  wellPanel(id = "title_section", 
    titlePanel("Chromatin and histone proteomics of human pluripotent states"),
    h3("Integrated Multi-Omics Analyses Reveal Polycomb Repressive Complex 2
  Restricts Naive Human Pluripotent Stem Cell to Trophoblast Fate Induction"),
    h4("Zijlmans, Talon, Verhelst, Bendall et al."),
    br()
  ),
  wellPanel(id = "table_volcano_panel",
    fluidRow(
      column(
        width = 7,
        wellPanel(
          id = "table_panel", 
          actionButton("clear_table", label = "Clear table selections"),
          DT::dataTableOutput("pp_table"))
      ),
      column(
        width = 5,
        wellPanel(
          id = "volcano_panel",
          br(),
          plotlyOutput("volcano", height = "400px")
        )
      )
    )
  ),
  br(),
  wellPanel(id = "plot_type_selections",
    fluidRow(
      column(
        width = 2,
        checkboxInput(inputId = "select_all_plots", label = strong("Show all data types"), value = TRUE)
      ),
      column(
        width = 10,
        checkboxGroupInput(
          inputId = "plot_panels_to_display",
          label = NULL,
          choices = datatypes,
          inline = TRUE
        )
      )
    )
  ),
  br(),
  fluidRow(
    conditionalPanel(
      condition = "input.plot_panels_to_display.includes('gene_expr')",
      column(width = 6, uiOutput(outputId = "gene_expr", class = "plot_box"))
    ),
    conditionalPanel(
        condition = "input.plot_panels_to_display.includes('acid_protein')",
        column(width = 6, uiOutput("protein1", class = "plot_box"))
    ),
    conditionalPanel(
        condition = "input.plot_panels_to_display.includes('chr_protein')",
        column(width = 6, uiOutput("protein_second", class = "plot_box"))
    ),
    conditionalPanel(
        condition = "input.plot_panels_to_display.includes('histones')",
        column(width = 6, uiOutput("histones", class = "plot_box"))
    ) 
 ),
  checkboxGroupInput(
      inputId = "conditions_to_display",
      label = "",
      choices = conditions,
      selected = c("Naive", "Primed"),
      inline = TRUE
  ),
  br(),
  actionButton("browser", "browser")
)

key <- row.names(volcano_dataset)

server <- function(input, output, session) {
    
  table_proxy <- dataTableProxy("pp_table")
  
  selected_ids <- reactiveValues()
  
  output$pp_table <- DT::renderDataTable({
    dt_setup(
      dataset,
      #lineHeight = "50%",
      selection = "multiple",
      dt_options = list(
        dom = "ftip", 
        columnDefs = list(
          list(
            targets = 2,
            width = "600px",
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 70 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 70) + '...</span>' : data;",
              "}"
            )
          )
        )  
      )
    )
  })
  
  # dataset filters ----
  
  filtered_acid_dataset <- reactive({
    data_long %>%
      filter(condition %in% input$conditions_to_display)
  })
  
  filtered_gene_dataset <- reactive({
    genes_long %>%
      filter(condition %in% input$conditions_to_display)
  })
  
  
  # volcano plot ----
  
  output$volcano <- renderPlotly({

    p <- volcano_dataset %>%
          ggplot(aes(x = log2fc_naive_primed, y = -log10(Naive_vs_Primed), key = key)) + #colour = species, fill = species)) +
            geom_point(shape = 21) +
            theme(legend.position="none")
      
      if(! is.null(selected_ids$protein_acid)) {
          
          selected_subset <- volcano_dataset %>%
             filter(Accession %in% selected_ids$protein_acid)
          #colours_needed <- plot_colours[1:nrow(selected_subset)]

          p <- p + geom_point(
                    data = selected_subset,
                    aes(key = NULL),
                    shape = 21,
                    colour = "red",
                    # fill = colours_needed, # this works with ggplot but not plotly
                    stroke = 1,
                    size = 3
                )
        }
            
        ggplotly(p)
        
    })
    
  # observeEvents ----
    
  observeEvent(input$clear_table, {
      selectRows(table_proxy, selected = NULL)
    selected_ids$protein_acid <- NULL
  })
  
  observeEvent(input$select_all_plots, {
      
    if(input$select_all_plots) {
      updateCheckboxGroupInput(session, "plot_panels_to_display", selected = datatypes)
    } 
  })
   
  observeEvent(input$plot_panels_to_display, {
    
    if(length(input$plot_panels_to_display) != 4) {
      updateCheckboxInput(session, "select_all_plots", value = FALSE)
    }
  })
   
  ## table row selections ----
  observeEvent(input$pp_table_rows_selected, {
      
    row_numbers <- as.numeric(input$pp_table_rows_selected)
    
    #print("updating selected rows")
    
    if(length(row_numbers) > 4){
      shinyalert::shinyalert(title = "", text = "Maximum of 4 rows of data will be shown.")
      row_numbers <- row_numbers[1:4]
    } 
    selected_ids$protein_acid <- table_data %>%
      slice(row_numbers) %>%
      pull(Accession)
    
    selected_ids$gene <- table_data %>%
      slice(row_numbers) %>%
      pull(Gene)
  })
    
  ## update highlighted rows on table ----
  observeEvent(selected_ids$protein_acid, ignoreNULL = FALSE, {
  
    if(length(selected_ids$protein_acid) != length(input$pp_table_rows_selected)){
      print("discrepancy in selections!!")
  
      if(is.null(selected_ids$protein_acid)) {
        selectRows(table_proxy, selected = NULL)
      } else {
        selected_table_rows <- table_data %>%
          filter(Accession %in% selected_ids$protein_acid) %>%
          pull(rowid)
        selectRows(table_proxy, selected = selected_table_rows)
      }
    }
  })
    
    
  # plotly events ----
    
  observeEvent(event_data("plotly_doubleclick"), {
      selected_ids$protein_acid <- NULL
  })
  
  observeEvent(event_data("plotly_click"), {

    d <- event_data("plotly_click")
    req(d)
    row_no <- as.numeric(d$key)
    
    selected_accessions <- volcano_dataset %>%
      slice(row_no) %>%
      pull(Accession)
    
    if(length(selected_ids$protein_acid) < 4) { # add up to 4 datasets
      selected_accessions <- c(selected_ids$protein_acid, selected_accessions)
    }    
    selected_ids$protein_acid <- selected_accessions
  })
  
  observeEvent(event_data("plotly_selected"), {
   
    d <- event_data("plotly_selected")
    req(d)
    row_numbers <- as.numeric(d$key)
    
    if(length(row_numbers) > 4){
      shinyalert::shinyalert(
        title = "",
        text = "Maximum of 4 data points will be shown. \n
        Try zooming in if points are too close together. Double click to deselect all."
      )
      row_numbers <- row_numbers[1:4]
    }
    selected_ids$protein_acid <- slice(volcano_dataset, row_numbers) %>%
      pull(Accession)
  })

      
  observeEvent(event_data("plotly_deselect"), {
    print("none selected")
    selected_ids$protein_acid <- NULL
    #selectRows(table_proxy, selected = NULL)
  })

  # plot panels -----
  output$gene_expr <- renderUI({
    
    req(selected_ids$gene)
    gene_exprUI <- mod_plotsUI("gene_expr_panel")
    mod_plotsServer("gene_expr_panel", filtered_gene_dataset,  selected_ids, id_type = "gene", plot_colours)
    
    wellPanel(
      id = "gene_expr_panel", 
      class = "plot_panel",
      h2("Gene expression", class = "panel_title"),
      gene_exprUI
    )
  }) 
    
  output$protein1 <- renderUI({
      
    req(selected_ids$protein_acid)
    protein1UI <- mod_plotsUI("protein1_panel")
    mod_plotsServer("protein1_panel", filtered_acid_dataset,  selected_ids, id_type = "protein_acid", plot_colours)
  
    wellPanel(
      id = "prot_acid_panel", 
      class = "plot_panel",
      h2("Acid extractome protein abundance", class = "panel_title"),
      protein1UI
    )
  })
    
  output$protein_second <- renderUI({
    
    req(selected_ids$protein_acid)
    protein2UI <- mod_plotsUI("protein2_panel")
    mod_plotsServer("protein2_panel", filtered_acid_dataset, selected_ids, id_type = "protein_acid", plot_colours)
    
    wellPanel(
      id = "prot_chromatin_panel", 
      class = "plot_panel",
      h2("Chromatin protein abundance", class = "panel_title"),
      protein2UI
    )
  })
    
  output$histones <- renderUI({
    
    req(selected_ids$protein_acid)
    histoneUI <- mod_plotsUI("histone_panel")
    mod_plotsServer("histone_panel", filtered_acid_dataset,  selected_ids, id_type = "protein_acid", plot_colours)
    
    wellPanel(
      id = "histone_panel", 
      class = "plot_panel",
      h2("Histone abundance", class = "panel_title"),
      histoneUI
    )
    
  })
        
  observeEvent(input$browser, browser())
    
}

# Run the application 
shinyApp(ui = ui, server = server)
