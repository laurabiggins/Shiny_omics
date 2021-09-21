library(shiny)
library(tidyverse)
library(plotly)
library(DT)


dataset <- readRDS("data/acid_meta.rds") #used for the DT - searching
data_long <- readRDS("data/acid_long.rds")
#acid_pval_fc <- readRDS("data/acid_pval_fc.rds")
acid_pval_fc <- readRDS("data/pval_check_temp.rds")

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
    titlePanel("Chromatin and histone proteomics of human pluripotent states"),
    h3("Integrated Multi-Omics Analyses Reveal Polycomb Repressive Complex 2
Restricts Naive Human Pluripotent Stem Cell to Trophoblast Fate Induction"),
    h4("Zijlmans, Talon, Verhelst, Bendall et al."),
    br(),
    br(),

    actionButton("clear_table", label = "Clear table selections"),
    fluidRow(
        column(
            width = 7,
            DT::dataTableOutput("pp_table")
        ),
        column(
            width = 5,
            br(),
            plotlyOutput("volcano")
        )
    ),
    br(),
    checkboxInput(inputId = "select_all_plots", label = "Show all data types", value = TRUE),
    checkboxGroupInput(
        inputId = "plot_panels_to_display",
        label = "",
        choices = datatypes,
        inline = TRUE
    ),
    
    fluidRow(
        conditionalPanel(
            condition = "input.plot_panels_to_display.includes('gene_expr')",
            column(
                width = 6,
                uiOutput(outputId = "gene_expr", class = "plot_box")
            )
        ),
        conditionalPanel(
            condition = "input.plot_panels_to_display.includes('acid_protein')",
            column(
                width = 6, 
                uiOutput("protein1", class = "plot_box")
            )
        )
    ),
    fluidRow(
        conditionalPanel(
            condition = "input.plot_panels_to_display.includes('chr_protein')",
            column(
               width = 6,
               uiOutput("protein_second", class = "plot_box")
           )
        ),
        conditionalPanel(
            condition = "input.plot_panels_to_display.includes('histones')",
            column(
               width = 6, 
               uiOutput("histones", class = "plot_box")
            )
        )    
   ),
    br(),
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

server <- function(input, output) {
    
    table_proxy <- dataTableProxy("pp_table")
    
    output$pp_table <- DT::renderDataTable(dataset)
    
    selected_ids <- reactiveVal()
    
    filtered_dataset <- reactive({
        data_long %>%
            filter(condition %in% input$conditions_to_display)
    })
    
    output$volcano <- renderPlotly({
  
        p <- volcano_dataset %>%
            ggplot(aes(x = log2fc_naive_primed, y = -log10(Naive_vs_Primed), key = key)) + #colour = species, fill = species)) +
                geom_point(shape = 21) +
                theme(legend.position="none")
        
        if(! is.null(selected_ids())) {
            
            selected_subset <- volcano_dataset %>%
               filter(Accession %in% selected_ids())
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
    
    observeEvent(input$clear_table, {
        selectRows(table_proxy, selected = NULL)
        selected_ids(NULL)
    })
    
    observeEvent(input$select_all_plots, {
        
        print("observed")
        
        if(input$select_all_plots) {
            updateCheckboxGroupInput(
                inputId = "plot_panels_to_display",
                selected = datatypes
            )
        } else {
            print("FALSE, do something")
            updateCheckboxGroupInput(
                inputId = "plot_panels_to_display",
                selected = ""
            )
        }
        
    })
    
    observeEvent(input$pp_table_rows_selected, {
        
        row_numbers <- as.numeric(input$pp_table_rows_selected)
        
        print("updating selected rows")
        
        if(length(row_numbers) > 4){
            shinyalert::shinyalert(
                title = "",
                text = "Maximum of 4 rows of data will be shown.")
            row_numbers <- row_numbers[1:4]
        } 
        selected_ids(
            table_data %>%
                slice(row_numbers) %>%
                pull(Accession)
        )
    })
    
    # update highlighted rows on table
    observeEvent(selected_ids(), ignoreNULL = FALSE, {
        
        if(length(selected_ids()) != length(input$pp_table_rows_selected)){
            print("discrepancy in selections!!")
        
            if(is.null(selected_ids())) {
                selectRows(table_proxy, selected = NULL)
            } else {
                selected_table_rows <- table_data %>%
                    filter(Accession %in% selected_ids()) %>%
                    pull(rowid)
                selectRows(table_proxy, selected = selected_table_rows)
            }
        }
    })
    
    observeEvent(event_data("plotly_doubleclick"), {
        print("double_clicked")
        selected_ids(NULL)
    })
    
    observeEvent(event_data("plotly_click"), {
        print("clicked")
        d <- event_data("plotly_click")
        req(d)
        row_no <- as.numeric(d$key)
        
        selected_accessions <- volcano_dataset %>%
            slice(row_no) %>%
            pull(Accession)
        
        if(length(selected_ids()) < 4) { # add up to 4 datasets
            selected_accessions <- c(selected_ids(), selected_accessions)
        }    
        selected_ids(selected_accessions)
    })
    
    observeEvent(event_data("plotly_selected"), {
        print("selected")
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
        selected_ids(
            volcano_dataset %>%
                slice(row_numbers) %>%
                pull(Accession)
        )
    })

        
    observeEvent(event_data("plotly_deselect"), {
        print("none selected")
        selected_ids(NULL)
        #selectRows(table_proxy, selected = NULL)
    })
    
    
    output$protein1 <- renderUI({
            
        if(is.null(selected_ids())) tags <- NULL
        
        req(selected_ids())
        
        if(length(selected_ids()) > 4){
            tags <- NULL
        }
        
        if(length(selected_ids()) == 1){
            tags <- tagList(
                plotOutput("prot1", height = acid_plot_height)
            )
        } else if (length(selected_ids()) == 2){
            tags <- tagList(
                splitLayout(
                    plotOutput("prot1", height = acid_plot_height),
                    plotOutput("prot2", height = acid_plot_height)
                )
            )
        } else if (length(selected_ids()) == 3){
            tags <- tagList(
                fluidRow(
                    column(width = 6, plotOutput("prot1", height = acid_plot_height)),
                    column(width = 6, plotOutput("prot2", height = acid_plot_height))
                ),
                br(),
                fluidRow(
                    column(width = 6, plotOutput("prot3", height = acid_plot_height))
                )
            )
        } else if (length(selected_ids()) == 4){
            tags <- tagList(
                fluidRow(
                    column(width = 6, plotOutput("prot1", height = acid_plot_height)),
                    column(width = 6, plotOutput("prot2", height = acid_plot_height))
                ),
                br(),
                fluidRow(
                    column(width = 6, plotOutput("prot3", height = acid_plot_height)),
                    column(width = 6, plotOutput("prot4", height = acid_plot_height))
                )
            )
        }
        wellPanel(id = "prot_acid_panel", 
                  class = "plot_panel",
                  h2("Acid extractome protein abundance", class = "panel_title"),
                  tags)
    })
    
    
    output$protein_second <- renderUI({
        
        req(selected_ids())
        
        protein2UI <- mod_plotsUI("protein2_panel")
        
        mod_plotsServer(
            "protein2_panel",
            filtered_dataset, 
            selected_ids,
            plot_colours
        )
        
        wellPanel(
            id = "prot_chromatin_panel", 
            class = "plot_panel",
            h2("Chromatin protein abundance", class = "panel_title"),
            protein2UI
        )
        
    })
    
       
    output$gene_expr <- renderUI({
        
        req(selected_ids())
        
        gene_exprUI <- mod_plotsUI("gene_expr_panel")
        
        mod_plotsServer(
            "gene_expr_panel",
            filtered_dataset, 
            selected_ids,
            plot_colours
        )
        
        wellPanel(
            id = "gene_expr_panel", 
            class = "plot_panel",
            h2("Gene expression", class = "panel_title"),
            gene_exprUI
        )
    })
    
    output$histones <- renderUI({
        
        req(selected_ids())
        
        histoneUI <- mod_plotsUI("histone_panel")
        
        mod_plotsServer(
            "histone_panel",
            filtered_dataset, 
            selected_ids,
            plot_colours
        )
        
        wellPanel(
            id = "histone_panel", 
            class = "plot_panel",
            h2("Histone abundance", class = "panel_title"),
            histoneUI
        )
        
    })
        

    ## protein acid plots ----
    acid_boxplot <- function(data, title, box_colour){
        
        ggplot(data, aes(x = condition, y = value)) +
            geom_boxplot(fill = box_colour) +
            xlab("") +
            ggtitle(title)
        
    }
    
    output$prot1 <- renderPlot({
        id <- selected_ids()[1]
        req(id)
        
        data_long %>%
            filter(Accession == id) %>%
            filter(condition %in% input$conditions_to_display) %>%
            ggplot(aes(x = condition, y = value)) +
                geom_boxplot(fill = plot_colours[1])
    })
    
    output$prot2 <- renderPlot({
        
        req(selected_ids()[2])
        #id <- dataset[["Accession"]][selected_ids()[2]]
        id <- selected_ids()[2]
        
        data_filt <- data_long %>%
            filter(Accession == id) %>%
            filter(condition %in% input$conditions_to_display) # could extract this but leave for now
        
        acid_boxplot(data_filt, id, plot_colours[2])
        
    })
    
    output$prot3 <- renderPlot({
        
        req(selected_ids()[3])
        id <- selected_ids()[3]
        
        data_filt <- data_long %>%
            filter(Accession == id) %>%
            filter(condition %in% input$conditions_to_display) # could extract this but leave for now
        
        acid_boxplot(data_filt, id, plot_colours[3])
    })
    
    output$prot4 <- renderPlot({
        
        req(selected_ids()[4])
        id <- selected_ids()[4]
        
        data_filt <- data_long %>%
            filter(Accession == id) %>%
            filter(condition %in% input$conditions_to_display) # could extract this but leave for now
        
        acid_boxplot(data_filt, id, plot_colours[4])
    })
    
    observeEvent(input$browser, browser())
    
}

# Run the application 
shinyApp(ui = ui, server = server)
