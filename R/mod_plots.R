
mod_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
      uiOutput(ns("plot_panel"), class = "plot_box"),
  )  
      
}

mod_plotsServer <- function(id, data_long, selected_ids, plot_colours, plot_height = "200px") {
  moduleServer(id, function(input, output, session) {
    
    ns_server <- NS(id)

    observeEvent(input$browser, browser())
    
   # data_initial_filt <- reactive(filter(data_long(), Accession %in% selected_ids()))

    
    output$plot_panel <- renderUI({
      
      ids <- selected_ids()
      
      if(is.null(ids)) tags <- NULL
      
      req(ids)
      
      if(length(ids) > 4){
        tags <- NULL
      }
      
      if(length(ids) == 1){
        tags <- tagList(
          plotOutput(ns_server("plot1"), height = plot_height)
        )
      } else if (length(ids) == 2){
        tags <- tagList(
          splitLayout(
            plotOutput(ns_server("plot1"), height = plot_height),
            plotOutput(ns_server("plot2"), height = plot_height)
          )
        )
      } else if (length(ids) == 3){
        tags <- tagList(
          fluidRow(
            column(width = 6, plotOutput(ns_server("plot1"), height = plot_height)),
            column(width = 6, plotOutput(ns_server("plot2"), height = plot_height))
          ),
          br(),
          fluidRow(
            column(width = 6, plotOutput(ns_server("plot3"), height = plot_height))
          )
        )
      } else if (length(selected_ids()) == 4){
        tags <- tagList(
          fluidRow(
            column(width = 6, plotOutput(ns_server("plot1"), height = plot_height)),
            column(width = 6, plotOutput(ns_server("plot2"), height = plot_height))
          ),
          br(),
          fluidRow(
            column(width = 6, plotOutput(ns_server("plot3"), height = plot_height)),
            column(width = 6, plotOutput(ns_server("plot4"), height = plot_height))
          )
        )
      }
      tags
    })
    
    ## protein acid plots ----
    acid_boxplot <- function(data, title, box_colour){

      ggplot(data, aes(x = condition, y = value)) +
        geom_boxplot(fill = box_colour) +
        xlab("") +
        ggtitle(title)
    }
    
    output$plot1 <- renderPlot({
      
      id <- selected_ids()[1]
      req(id)

      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[1])
        
    }) %>% bindCache(selected_ids()[1])#, data_long()) # could probably pass the name of the dataset - would be 
    # more efficient to compare than the whole dataset
    
    output$plot2 <- renderPlot({
      
      req(selected_ids()[2])
      id <- selected_ids()[2]
      
      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[2])
      
    }) %>% bindCache(selected_ids()[2], data_long())
    
    output$plot3 <- renderPlot({
      
      req(selected_ids()[3])
      id <- selected_ids()[3]
      
      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[3])
    }) %>% bindCache(selected_ids()[3], data_long())
    
    output$plot4 <- renderPlot({
      
      req(selected_ids()[4])
      id <- selected_ids()[4]
      
      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[4])
    }) %>% bindCache(selected_ids()[4], data_long())
  })   
}            