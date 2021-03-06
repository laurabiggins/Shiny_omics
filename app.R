library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(htmlwidgets)
# 
# # TODO: sort out main data table
# # tool tip text on volcano plot
# 
browser_buttons <- TRUE
# 
gene_id_table <- readRDS("data/gene_id_table.rds") #%>%
#  rename(`Acid extract` = `acid extractome`, `Gene expr` = `gene expr`, ChEP = chep, Histone = histone) %>%
#  select(rowid, Gene_id, `Gene expr`, ChEP, `Acid extract`, Histone, everything())

#used for the DT - searching
table_data <- readRDS("data/table_data.rds")
data_long <- readRDS("data/acid_long.rds")
acid_pval_fc <- readRDS("data/acid_pval_fc.rds")
genes_long <- readRDS("data/genes_long.rds") #%>%
#  mutate(Gene_id = Gene_expr_id)
gene_pval_fc <- readRDS("data/gene_pval_fc.rds") 
histone_data <- readRDS("data/histone_data.rds")
all_histone_links <- readRDS("data/all_histone_links.rds")
histone_pval_fc <- readRDS("data/histone_pval_fc.rds")
chep_data <- readRDS("data/chep_data.rds")
chep_pval_fc <- readRDS("data/chep_pval_fc.rds")

conditions <- c("Naive", "Naive+PRC2i", "Primed", "Primed+PRC2i")
histone_media <- c("PXGL", "ENHSM", "t2iLGo")

datatypes <- c(
 "Gene expression" = "gene_expr", 
 "Acid extractome protein abundance" = "acid_protein",
 "Chromatin-associated protein abundance" = "chr_protein",
 "Histone modifications" = "histones"
)

volcano_highlights <- scale_color_manual(
  name = "selected", 
  values = c(not = "#ffffff00", selected = "red" )
)

# UI ----
ui <- fluidPage(

  shinyalert::useShinyalert(),
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "omics.css")
  ),
  
  ## Title panel ----
  wellPanel(id = "title_section", 
    titlePanel("Chromatin and histone proteomics of human pluripotent states"),
    h3("Integrated Multi-Omics Reveal Polycomb Repressive Complex 2 Restricts Human Trophoblast Induction"),
    h4("Zijlmans, Talon, Verhelst, Bendall et al.")
  ),
  ## table and volcano panel ----
  wellPanel(id = "table_volcano_panel",
            h5(" Explore the datasets associated with this paper by searching the table or browsing the points in the volcano plot. To search the table, enter a gene name, protein name or ID, or histone mark. Select up to 4 IDs to view the data in plots below."),
            h5("Selections can also be made by clicking on points in the volcano plot. The subset of selected IDs can be viewed by switching to the 'selected' tab below."),
    fluidRow(
      column(class = "table_padding",
        width = 6,
        wellPanel(
          id = "table_panel", 
          actionButton("clear_table", label = "Clear all table selections"),
          tabsetPanel(
            tabPanel("all", 
                     br(),
                     DT::dataTableOutput("pp_table")),
            tabPanel("selected", 
                     br(), 
                     DT::dataTableOutput("selected_table"),
                     br(),
                     actionButton("remove_selected", label = "remove selected"))
          )
        )
      ),
      column(class = "volcano_padding",
        width = 6,
        wellPanel(
          id = "volcano_panel",
          br(),
          radioButtons(
            inputId = "volcano_type", 
            label = NULL, 
            choices = datatypes[c(2,3,1,4)], 
            inline = TRUE
          ),
          radioButtons(
            inputId = "volcano_condition_type", 
            label = NULL, 
            choices = list("naive vs primed" = "naive_primed", 
                           "naive vs naive+PRC2i" = "naive_naive2i",
                           "primed vs primed+PRC2i" = "primed_primed2i"),
            inline = TRUE
          ),
          br(),
          plotlyOutput("volcano", height = "400px")
        )
      )
    )
  ),
  ## data type checkboxes ----
  wellPanel(id = "plot_type_selections", class = "bordered_panel",
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
    ),
    fluidRow(
      column(3, strong("Select conditions to display in plots")),
      column(9, 
             checkboxGroupInput(
               inputId = "conditions_to_display",
               label = NULL,
               choices = conditions,
               selected = c("Naive", "Primed"),
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
        column(width = 10, uiOutput("all_histones", class = "plot_box"))
    ),
    conditionalPanel(
      condition = "input.plot_panels_to_display.includes('histones')",
      column(width = 2, uiOutput("histone_media", class = "plot_box"))
    )
  ),
  br(),
  p("Any problems please email laura.biggins@babraham.ac.uk", 
    style = "font-size:14px", align = "right"),
  if(browser_buttons) actionButton("browser", "browser")
)


# server ----
server <- function(input, output, session) {
   
  library(DT)
   
  table_proxy <- dataTableProxy("pp_table")
  
  hideCols(table_proxy, hide = c(0, 6:10))
  
  ## filtered datatable ----
  filtered_meta <- reactiveVal()
  
  selected_histone_media <- reactiveVal("PXGL")
  
  shinyjs::hide("remove_selected")
  
  # key for data table ----
  key <- reactive({
    row.names(volcano_dataset())
  })
  
  ## main datatable ----
  output$pp_table <- DT::renderDataTable({

    datatable(
      gene_id_table, 
      rownames = FALSE,
      options = list(dom = "ftip", pageLength = 13, scrollX = TRUE, autoWidth = FALSE)
    ) %>% 
      formatStyle(
        c('Gene expr', 'ChEP', 'Acid extract', 'Histone'),
        backgroundColor = styleEqual(levels = "no data", values = '#CFCCC9', default = '#B5E3D9'),
        `font-size` = '60%'
      ) %>% formatStyle(0, target = 'row', lineHeight = '50%')
      
  })
  
  selected_datatable <- reactive({
    
    req(filtered_meta())
    selected_table <- filtered_meta() %>%
      select(Gene_id, Description, Gene_expr_id, Majority.protein.IDs, Accession, histone_mark) %>%
      rename(`Gene id` = Gene_id,
             `Acid id` = Accession, 
             `ChEP id` = Majority.protein.IDs, 
             `Gene expr id`=Gene_expr_id, 
             `Histone mark` = histone_mark) %>%
      replace(is.na(.), "no data")
    
    datatable(
      selected_table,
      rownames = FALSE,
      options = list(
        dom = "t", 
        scrollX = TRUE, 
        autoWidth = FALSE,
        columnDefs = list(
          list(
            targets = c(3,4),
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 15 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
              "}")
          ),
          list(
            targets = 1,
            width = "900px",
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 30 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
              "}")
          )
        )
      )  
    ) %>%
      formatStyle(
        c("Acid id", "ChEP id", "Gene expr id", "Histone mark"),
        backgroundColor = styleEqual(levels = "no data",  values = '#CFCCC9', default = '#B5E3D9')
      ) %>% 
      formatStyle(0, target = 'row', `font-size` = '70%') %>%
      formatStyle(1, `font-size` = '120%')
    
  })
  
  
  ## render filtered meta ----
  output$selected_table <- DT::renderDataTable({
    
    req(selected_datatable())
    selected_datatable()
  })
  
  
  ## dataset filters ----
  filtered_acid_dataset <- reactive({
    
    data_long %>%
      filter(condition %in% input$conditions_to_display) %>%
      filter(Accession %in% filtered_meta()$Accession)
  })
  
  filtered_chep_dataset <- reactive({
    
    chep_data %>%
      filter(condition %in% input$conditions_to_display) %>%
      filter(Majority.protein.IDs %in% filtered_meta()$Majority.protein.IDs)
    
  })
  
  filtered_gene_dataset <- reactive({

    genes_long %>%
      filter(condition %in% input$conditions_to_display) %>%
      filter(Gene_expr_id %in% filtered_meta()$Gene_expr_id)
      
  })
  
  filtered_histone_dataset <- reactive({

    #req(filtered_meta()$histone_mark)
    selected_histones <- filtered_meta()$histone_mark
    
    if(any(str_detect(selected_histones, ","), na.rm = TRUE)){
      multi_links <- which(str_detect(selected_histones, ","))
      multi_genes <- filtered_meta()$Gene_id[multi_links]
      
      histones_to_add <- all_histone_links %>%
        filter(Gene_id %in% multi_genes) %>%
        pull(histone_mark)
      
      selected_histones <- selected_histones[-multi_links]
      selected_histones <- c(selected_histones, histones_to_add)
    }
    
    selected_histones <- unique(selected_histones[!is.na(selected_histones)])

    histone_data %>%
      filter(condition %in% input$conditions_to_display) %>%
      #filter(medium %in% input$histone_media) %>%
      filter(medium %in% selected_histone_media()) %>%
      filter(histone_mark %in% selected_histones)
    
  })

  volcano_dataset <- reactive({
    
    volcano_ds <- switch(
      input$volcano_type,
        histones = mutate(
          histone_pval_fc, 
          selected = if_else(histone_mark %in% filtered_histone_dataset()$histone_mark, "selected", "not")
        ),
        acid_protein = mutate(
          acid_pval_fc,
          selected = if_else(Accession %in% filtered_acid_dataset()$Accession, "selected", "not")
        ),
        chr_protein  = mutate(
          chep_pval_fc,
          selected = if_else(Majority.protein.IDs %in% filtered_chep_dataset()$Majority.protein.IDs , "selected", "not")
        ),
        gene_expr = mutate(
          gene_pval_fc, 
          selected = if_else(Gene_expr_id %in% filtered_gene_dataset()$Gene_expr_id, "selected", "not")
        )
    )
    
    volcano_ds %>%
      filter(condition == input$volcano_condition_type)
  })
  
  volcano_title <- reactive({
    if_else(input$volcano_type == "histones", "PXGL", "")
  })
  
  
  ### volcano plot ----
  output$volcano <- renderPlotly({
    
    first_col <- colnames(volcano_dataset())[1]
    
    p <- volcano_dataset() %>%
      ggplot(aes(x = log2fc, y = -log10(pval), key = key(), color = selected, label = Gene_id, text = substr(.data[[first_col]], 1, 15))) +
      geom_point(shape = 21, stroke = 1, size = 2, fill = "black") +
      geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
      geom_vline(xintercept = c(-1,1), linetype = "dashed", col = "darkgray") +
      theme(legend.position = "none") +
      ggtitle(volcano_title()) +
      volcano_highlights 
    
    #p <- p + geom_text(nudge_x = 0.5, nudge_y = 0.5, size = 4, fontface = "bold")
    ggplotly(p, tooltip = c("label", "text"))
        
    })
    
  ## observeEvents ----
  ### clear selections ---- 
  #### clear all ----  
  observeEvent(input$clear_table, {
    selectRows(table_proxy, selected = NULL)
    set_ids_to_null()
  })
  
  #### remove selected ----
  observeEvent(input$remove_selected, {
    req(input$selected_table_rows_selected)
    #filtered_meta(slice(filtered_meta(), -input$selected_table_rows_selected))
    new_selections <- slice(filtered_meta(), -input$selected_table_rows_selected) %>%
      pull(rowid)
    selectRows(table_proxy, selected = new_selections)
  })
  
  observeEvent(input$histone_media, {
    #print("selected media = ")
    #print(input$histone_media)
    selected_histone_media(input$histone_media)
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
          shinyalert::shinyalert(title = "", text = "Maximum of 4 rows of data will be shown. \n Clear row selections using the button above the table or by double clicking on an empty area of the plot.")
          row_numbers <- row_numbers[1:4]
          selectRows(table_proxy, selected = row_numbers)
        } 
        
        selected_rows <- slice(table_data, row_numbers)
        filtered_meta(selected_rows)
        shinyjs::show("remove_selected")
        
    }
  })
    
  ### update highlighted rows on table ----
  #observeEvent(filtered_meta(), ignoreNULL = FALSE, {
  observeEvent(filtered_meta(), {

   # browser()
    
    if(is.null(filtered_meta())){
      selectRows(table_proxy, selected = NULL)
      shinyjs::hide("remove_selected")
    } else if (length(filtered_meta()$rowid) != length(input$pp_table_rows_selected)){
      selectRows(table_proxy, selected = filtered_meta()$rowid)
      shinyjs::show("remove_selected")
    } else if (any(filtered_meta()$rowid != input$pp_table_rows_selected)){
        print("discrepancy in selections!!")
        # we go with the filtered_meta() object as this gets updated from other sources i.e. the volcano plot
        selectRows(table_proxy, selected = filtered_meta()$rowid)
        shinyjs::show("remove_selected")
    }
  })
   
  set_ids_to_null <- function(){
    print("setting ids to null")
    filtered_meta(NULL)
    selectRows(table_proxy, selected = NULL)
    shinyjs::hide("remove_selected")
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
    
    if(input$volcano_type == "histones"){
      shinyalert::shinyalert(
        title = "",
        text = "Choose a different type of data to enable selections."
      )
    }
    else {
      row_no <- as.numeric(d$key)
      
      selected_accessions <- volcano_dataset() %>%
        slice(row_no) %>%
        pull(1)
      
      id_type <- colnames(volcano_dataset())[1]
      row_to_add <- filter(table_data, .data[[id_type]] == selected_accessions)
  
      if(is.null(filtered_meta())) {
        filtered_meta(row_to_add)
      } else if (nrow(filtered_meta()) < 4) {
        filtered_meta(bind_rows(filtered_meta(), row_to_add))
      } else if (nrow(filtered_meta()) == 4){
        shinyalert::shinyalert(
          title = "",
          text = "Maximum of 4 data points will be shown. \n
           Double click on an empty area of the plot to deselect all."
        )
      } 
    }
  })

  # this wipes out other selections at the moment  
  observeEvent(event_data("plotly_selected"), {
   
    d <- event_data("plotly_selected")
    
    req(d)
    
    if(input$volcano_type == "histones"){
      shinyalert::shinyalert(
        title = "",
        text = "Select a different type of data to enable selections."
      )
    } else {
    
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
    }
  })

      
  observeEvent(event_data("plotly_deselect"), {
    print("none selected")
    filtered_meta(NULL)
  })

  ## plot panels -----
  ### gene expr ----
  output$gene_expr <- renderUI({
    
    if(!is.null(filtered_meta()[["Gene_expr_id"]])){
      if(!isTruthy(filtered_meta()[["Gene_expr_id"]])){
        welltags <- p(class = "no_data", "No data for selected genes")
      } else {
        req(filtered_meta()$Gene_expr_id)
        gene_exprUI <- mod_plotsUI("gene_expr_panel")
        welltags <- gene_exprUI
        mod_plotsServer(
          "gene_expr_panel", 
          filtered_gene_dataset, 
          selected_ids = filtered_meta()$Gene_expr_id,
          panel_name = "gene_expression", 
          ylabel = "log2 normalised counts"
        )
      }  
      wellPanel(
        id = "gene_expr_panel", 
        class = "plot_panel",
        h2("Gene expression", class = "panel_title"),
        welltags
      )
    }  
  }) 
 
  # This works but it stops the flow of the fluidRow elements. I think to make it work better, the javascript for the conditional panel should include an argument triggered by a close button. 
#  observeEvent(input$close_gene_panel, {
#    shinyjs::hideElement("gene_expr")
#  })
  
   
  ### protein acid ----
   
  output$protein1 <- renderUI({
     
    if(!is.null(filtered_meta()[["Accession"]])){
      if(!isTruthy(filtered_meta()[["Accession"]])){
        welltags <- p(class = "no_data", "No data for selected genes")
      } else {
        req(filtered_meta()[["Accession"]])
        protein1UI <- mod_plotsUI("protein1_panel")
        welltags <- protein1UI
        mod_plotsServer(
          "protein1_panel", 
          filtered_acid_dataset,  
          selected_ids = filtered_meta()$Accession, 
          panel_name = "acid_extractome", 
          ylabel = "normalised abundance"
        )
      }
      wellPanel(
        id = "prot_acid_panel", 
        class = "plot_panel",
        h2("Acid extractome protein abundance", class = "panel_title"),
        welltags
      )
    }
  })
   
  ### protein chr ----
   
  output$protein_second <- renderUI({
    
    if(!is.null(filtered_meta()[["Majority.protein.IDs"]])){
      if(!isTruthy(filtered_meta()[["Majority.protein.IDs"]])){
        welltags <- p(class = "no_data", "No data for selected genes")
      } else {
        req(filtered_meta()[["Majority.protein.IDs"]])
        protein2UI <- mod_plotsUI("protein2_panel")
        welltags <- protein2UI
        mod_plotsServer(
          "protein2_panel", 
          filtered_chep_dataset, 
          selected_ids = filtered_meta()$Majority.protein.IDs, 
          panel_name = "chromatin_associated_protein", 
          ylabel = "LFQ intensity"
        )
      }   
      wellPanel(
        id = "prot_chromatin_panel", 
        class = "plot_panel",
        h2("Chromatin protein abundance", class = "panel_title"),
        welltags
      )
    }
  })
  
  output$histone_media <- renderUI({
    wellPanel(
      #id = "all_histone_panel", 
      class = "plot_panel",
      h2("Histone media", class = "panel_title"),
      checkboxGroupInput(inputId = "histone_media", label = NULL, 
                         choices = histone_media, selected = histone_media[1]),
      p("Data for different histone modifications may be available for one or more of the media.")
    )
  })
  
   
  ### all histones ----
  output$all_histones <- renderUI({
    
    if(!is.null(filtered_meta()[["histone_mark"]])){
      if(!isTruthy(filtered_histone_dataset()) | nrow(filtered_histone_dataset()) == 0){
          welltags <- p(class = "no_data", "No histone data for current selections")
      } else {
        req(filtered_histone_dataset())
        histoneUI <- mod_histone_plotsUI("all_histone_panel")
        mod_histone_plotsServer(
          "all_histone_panel", 
          data_long = filtered_histone_dataset, 
          panel_name = "histone_mark", 
          ylabel = "normalised abundance"
        )
        welltags <- histoneUI
      }
      wellPanel(
        #id = "all_histone_panel", 
        class = "plot_panel",
        h2("Histone modifications", class = "panel_title"),
        welltags#,
        #br(),
       # p("Histone data may only be available for a different medium."),
        # checkboxGroupInput(inputId = "histone_media", label = NULL, inline = TRUE, 
        #                    choices = histone_media, selected = selected_histone_media())
      )
    }
  })
  
  
        
  observeEvent(input$browser, browser())
    
}

# Run the application 
shinyApp(ui = ui, server = server)
