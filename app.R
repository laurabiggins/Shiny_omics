library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(htmlwidgets)
# 
# # TODO: sort out main data table
# # tool tip text on volcano plot
# 
dataset <- readRDS("data/meta.rds") #used for the DT - searching
data_long <- readRDS("data/acid_long.rds")
acid_pval_fc <- readRDS("data/acid_pval_fc.rds")
genes_long <- readRDS("data/genes_long.rds")
histone_data <- readRDS("data/histone_data.rds")
histone_pval_fc <- readRDS("data/histone_pval_fc.rds")
chep_data <- readRDS("data/chep_data.rds")
chep_pval_fc <- readRDS("data/chep_pval_fc.rds")
#acid_pval_fc <- readRDS("data/pval_check_temp.rds")

table_data <- dataset %>%
 select(-Protein.names, -Description) %>%
 #select(-rowid, -Protein.names, -Description) %>%
 relocate(Gene_id) %>%
 arrange(desc(Accession)) %>%
 rowid_to_column()


#plot_colours <- c("red3", "blue3", "green3", "orange3")

conditions <- c("Naive", "Naive+PRC2i", "Primed", "Primed+PRC2i")
histone_media <- c("PXGL", "ENHSM", "t2iLGo")

acid_plot_height <- "250px"

datatypes <- c(
 "Gene expression" = "gene_expr", 
 "Acid extractome protein abundance" = "acid_protein",
 "Chromatin-associated protein abundance" = "chr_protein",
 "Histone modification abundance" = "histones"
)

# UI ----
ui <- fluidPage(

  shinyalert::useShinyalert(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "omics.css")
  ),
  
  ## Title panel ----
  wellPanel(id = "title_section", 
    titlePanel("Chromatin and histone proteomics of human pluripotent states"),
    h3("Integrated Multi-Omics Analyses Reveal Polycomb Repressive Complex 2
  Restricts Naive Human Pluripotent Stem Cell to Trophoblast Fate Induction"),
    h4("Zijlmans, Talon, Verhelst, Bendall et al."),
    br()
  ),
  ## table and volcano panel ----
  wellPanel(id = "table_volcano_panel",
    fluidRow(
      column(class = "table_padding",
        width = 7,
        wellPanel(
          id = "table_panel", 
          actionButton("clear_table", label = "Clear table selections"),
          tabsetPanel(
            tabPanel("all", DT::dataTableOutput("pp_table")),
            tabPanel("selected", br(), DT::dataTableOutput("selected_table"))
          )
        )
      ),
      column(class = "volcano_padding",
        width = 5,
        wellPanel(
          id = "volcano_panel",
          br(),
          plotlyOutput("volcano", height = "400px"),
          br(),
          radioButtons(
            inputId = "volcano_type", 
            label = NULL, 
            choices = datatypes[c(2,4,3)], 
            inline = TRUE
          ),
          radioButtons(
            inputId = "volcano_condition_type", 
            label = NULL, 
            choices = list("naive vs primed" = "naive_primed", 
                          "naive vs naive+PRC2i" = "naive_naive2i",
                           "primed vs primed+PRC2i" = "primed_primed2i"),
            inline = TRUE
          )
        )
      )
    )
  ),
  ## data type checkboxes ----
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
  ## plot panels ----
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
 ## condition type checkboxes ----
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


# server ----
server <- function(input, output, session) {
   
  library(DT)
   
  table_proxy <- dataTableProxy("pp_table")
  
  hideCols(table_proxy, hide = 0)
  
  ## filtered datatable ----
  filtered_meta <- reactiveVal()
  
  selected_ids <- reactiveValues()
  
  # key for data table ----
  key <- reactive({
    row.names(volcano_dataset())
  })
  
  ## main datatable ----
  output$pp_table <- DT::renderDataTable({
     
    table_data <- table_data %>% replace(is.na(.), "")
    
    datatable(
      table_data,
      rownames = FALSE,
      options = list(
        dom = "ftip",
        columnDefs = list(
          list(
            targets = 2,
            render = htmlwidgets::JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 10 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
              "}")
          ),
          list(
            targets = 3,
            render = htmlwidgets::JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 15 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
              "}")
          )
        )
      )
    )
                  
  })
  
  ## render filtered meta ----
  output$selected_table <- DT::renderDataTable({
    req(filtered_meta())
    dt_setup(filtered_meta())
  })
  
  
  ## dataset filters ----
  filtered_acid_dataset <- reactive({
    data_long %>%
      filter(condition %in% input$conditions_to_display)
  })
  
  filtered_chep_dataset <- reactive({
    chep_data %>%
      filter(condition %in% input$conditions_to_display)
  })
  
  filtered_gene_dataset <- reactive({
    genes_long %>%
      filter(condition %in% input$conditions_to_display)
  })
  
  filtered_histone_dataset <- reactive({
    histone_data %>%
      filter(condition %in% input$conditions_to_display) %>%
      filter(medium %in% input$histone_media)
  })
  
  
  volcano_dataset <- reactiveVal()
  
  volcano_dataset <- reactive({
    
    volcano_ds <- switch(input$volcano_type,
                         histones = histone_pval_fc,
                         acid_protein = acid_pval_fc,
                         chr_protein = chep_pval_fc)
    
   volcano_ds %>%
      filter(condition == input$volcano_condition_type) %>%
      drop_na() 
  })
  
  # observeEvent(input$volcano_type, {
  #   
  #   volcano_ds <- switch(input$volcano_type,
  #          histones = histone_pval_fc,
  #          acid_protein = acid_pval_fc,
  #          chr_protein = chep_pval_fc)
  #   
  #   volcano_ds <- volcano_ds %>%
  #     filter(condition == "naive_primed") %>%
  #     drop_na() 
  #   
  #   volcano_dataset(volcano_ds)
  # })
  
  ## volcano plot ----
  output$volcano <- renderPlotly({
    
    first_col <- colnames(volcano_dataset())[1]
    
    p <- volcano_dataset() %>%
      ggplot(aes(x = log2fc, y = -log10(pval), key = key(), label = .data[[first_col]])) + #colour = species, fill = species)) +
      geom_point(shape = 20) +
      geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
      geom_vline(xintercept = c(-1,1), linetype = "dashed", col = "darkgray") +
      #geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
      theme(legend.position="none")
    
      if(! is.null(filtered_meta())) {

        selected_subset <- left_join(filtered_meta(), volcano_dataset())
        
        if(! is.null(selected_subset)){
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
      }
            
      ggplotly(p, tooltip = "label")
        
    })
    
  ## observeEvents ----
  ### clear selections ----  
  observeEvent(input$clear_table, {
    selectRows(table_proxy, selected = NULL)
    set_ids_to_null()
  })
  
  ### select datatypes ----
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
   
  ### table row selections ----
  observeEvent(input$pp_table_rows_selected, ignoreNULL = FALSE, {
    
    
      
    row_numbers <- as.numeric(input$pp_table_rows_selected)
    
    print("updating selected rows")
    
    if(length(row_numbers) == 0){
      set_ids_to_null()
    } else {
        if(length(row_numbers) > 4){
          shinyalert::shinyalert(title = "", text = "Maximum of 4 rows of data will be shown.")
          row_numbers <- row_numbers[1:4]
        } 
      
        selected_rows <- slice(table_data, row_numbers)
        
        filtered_meta(selected_rows)
        
    }
  })
    
  ### update highlighted rows on table ----
  #observeEvent(filtered_meta(), ignoreNULL = FALSE, {
  observeEvent(filtered_meta(), {

   # browser()
    
    if(is.null(filtered_meta())){
      selectRows(table_proxy, selected = NULL)
    } else if (length(filtered_meta()$rowid) != length(input$pp_table_rows_selected)){
      selectRows(table_proxy, selected = filtered_meta()$rowid)
      
    } else if (any(filtered_meta()$rowid != input$pp_table_rows_selected)){
        print("discrepancy in selections!!")
        # we go with the filtered_meta() object as this gets updated from other sources i.e. the volcano plot
        selectRows(table_proxy, selected = filtered_meta()$rowid)
    }
  })
   
  set_ids_to_null <- function(){
    print("setting ids to null")
    filtered_meta(NULL)
    selectRows(table_proxy, selected = NULL)
  } 
    
  ## plotly events ----
    
  observeEvent(event_data("plotly_doubleclick"), {
      set_ids_to_null()
  })
  
  ### click - selecting one point at a time ----
  ### update the selected rows on table
  observeEvent(event_data("plotly_click"), {

    d <- event_data("plotly_click")
    req(d)
    row_no <- as.numeric(d$key)
    
    selected_accessions <- volcano_dataset() %>%
      slice(row_no) %>%
      pull(1)

    id_type <- colnames(volcano_dataset())[1]
    row_to_add <- filter(table_data, .data[[id_type]] == selected_accessions)
    
    #row_to_add <- left_join(selected_accessions, table_data)
    
    if(is.null(filtered_meta())) {
      filtered_meta(row_to_add)
    } else if (nrow(filtered_meta()) < 4) {
      filtered_meta(bind_rows(filtered_meta(), row_to_add))
    }  

  })

  # this wipes out other selections at the moment  
  observeEvent(event_data("plotly_selected"), {
   
    d <- event_data("plotly_selected")
    
    req(d)
    filtered_meta(NULL)
    row_numbers <- as.numeric(d$key)
    
    if(length(row_numbers) > 4){
      shinyalert::shinyalert(
        title = "",
        text = "Maximum of 4 data points will be shown. \n
        Try zooming in if points are too close together. Double click to deselect all."
      )
      row_numbers <- row_numbers[1:4]
    }
    
    selected_accessions <- volcano_dataset() %>%
      slice(row_numbers) %>%
      pull(1)
    
    id_type <- colnames(volcano_dataset())[1]
    selected_rows <- filter(table_data, .data[[id_type]] %in% selected_accessions)

    filtered_meta(selected_rows)
    
  })

      
  observeEvent(event_data("plotly_deselect"), {
    print("none selected")
    filtered_meta(NULL)
    #selectRows(table_proxy, selected = NULL)
  })

  ## plot panels -----
  ### gene expr ----
  output$gene_expr <- renderUI({
    
    if(!is.null(filtered_meta()[["Gene_expr_id"]])){
      if(!isTruthy(filtered_meta()[["Gene_expr_id"]])){
        wellPanel(
          id = "gene_expr_panel", 
          class = "plot_panel",
          h2("Gene expression", class = "panel_title"),
          p(class = "no_data", "No data for selected genes")
        )
      } else {
      
        req(filtered_meta()$Gene_expr_id)
        gene_exprUI <- mod_plotsUI("gene_expr_panel")
        mod_plotsServer("gene_expr_panel", filtered_gene_dataset,  filtered_meta, id_type = "Gene_expr_id")
        
        wellPanel(
          id = "gene_expr_panel", 
          class = "plot_panel",
          h2("Gene expression", class = "panel_title"),
          gene_exprUI
        )
      }
    }  
  }) 
   
  ### protein acid ----
   
  output$protein1 <- renderUI({
      
    req(filtered_meta()[["Accession"]])
    protein1UI <- mod_plotsUI("protein1_panel")
    mod_plotsServer("protein1_panel", filtered_acid_dataset,  filtered_meta, id_type = "Accession")
  
    wellPanel(
      id = "prot_acid_panel", 
      class = "plot_panel",
      h2("Acid extractome protein abundance", class = "panel_title"),
      protein1UI
    )
  })
   
  ### protein chr ----
   
  output$protein_second <- renderUI({
    
    req(filtered_meta()[["Majority.protein.IDs"]])
    protein2UI <- mod_plotsUI("protein2_panel")
    mod_plotsServer("protein2_panel", filtered_chep_dataset, filtered_meta, id_type = "Majority.protein.IDs", accession_col = "Majority.protein.IDs")
    
    wellPanel(
      id = "prot_chromatin_panel", 
      class = "plot_panel",
      h2("Chromatin protein abundance", class = "panel_title"),
      protein2UI
    )
  })
   
  ### histones ----
   
  output$histones <- renderUI({
    
    req(filtered_meta()[["histone_mark"]])
    histoneUI <- mod_plotsUI("histone_panel")
    mod_plotsServer("histone_panel", filtered_histone_dataset,  filtered_meta, id_type = "histone_mark", accession_col = "histone_mark", second_factor = "medium")
    
    wellPanel(
      id = "histone_panel", 
      class = "plot_panel",
      h2("Histone abundance", class = "panel_title"),
      histoneUI,
      checkboxGroupInput(inputId = "histone_media", label = "", inline = TRUE, 
                         choices = histone_media, selected = histone_media[1])
    )
    
  })
        
  observeEvent(input$browser, browser())
    
}

# Run the application 
shinyApp(ui = ui, server = server)
