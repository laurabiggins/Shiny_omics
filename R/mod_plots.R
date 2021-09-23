
## Try just having one layout!!

mod_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
      uiOutput(ns("plot_panel"), class = "plot_box"),
  )  
      
}

mod_plotsServer <- function(id, data_long, selected_ids, id_type, plot_colours, plot_height = "200px") {
  moduleServer(id, function(input, output, session) {
    
    ns_server <- NS(id)

    observeEvent(input$browser, browser())
    
   # data_initial_filt <- reactive(filter(data_long(), Accession %in% selected_ids()))
    ids <- reactive(selected_ids[[id_type]])
    
    tags_plot <- function(id, plot_name, plot_height = "200px"){
      if(id == "") {
        tags <- tagList(p(class = "no_data", "No data available"))
      } else {
        tags <- tagList(
          plotOutput(ns_server(plot_name), height = plot_height)
        )
      }
     tags 
    }
    
    
    
    output$plot_panel <- renderUI({
      
      #ids <- selected_ids[[id_type]]
      
      if(is.null(ids())) tags <- NULL
      
      req(ids())
      
      if(length(ids()) > 4){
        tags <- NULL
      } 
      
      if(length(ids()) == 1){
        tags <- tags_plot(ids()[1], plot_name = "plot1")
      } else if (length(ids()) == 2){
        tags <- tagList(
          fluidRow(
            column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
            column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
          )  
        )
      } else if (length(ids()) == 3){
        tags <- tagList(
          fluidRow(
            column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
            column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
          ),
          br(),
          fluidRow(
            column(width = 6, tags_plot(ids()[3], plot_name = "plot3"))
          )
        )
      } else if (length(ids()) == 4){
        tags <- tagList(
          fluidRow(
            column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
            column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
          ),
          br(),
          fluidRow(
            column(width = 6, tags_plot(ids()[3], plot_name = "plot3")),
            column(width = 6, tags_plot(ids()[4], plot_name = "plot4"))
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
      
      id <- ids()[1]
      req(id)

      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[1])
        
    }) %>% bindCache(ids()[1], data_long()) # could probably pass the name of the dataset - would be 
    # more efficient to compare than the whole dataset
    
    output$plot2 <- renderPlot({
      
      req(ids()[2])
      id <- ids()[2]
      
      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[2])
      
    }) %>% bindCache(ids()[2], data_long())
    
    output$plot3 <- renderPlot({
      
      req(ids()[3])
      id <- ids()[3]
      
      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[3])
    }) %>% bindCache(ids()[3], data_long())
    
    output$plot4 <- renderPlot({
      
      req(ids()[4])
      id <- ids()[4]
      
      data_filt <- data_long() %>%
        filter(Accession == id)
      
      acid_boxplot(data_filt, id, plot_colours[4])
    }) %>% bindCache(ids()[4], data_long())
  })   
}            