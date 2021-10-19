
## Try just having one layout!!

mod_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
      uiOutput(ns("plot_panel"), class = "plot_box"),
      actionButton(ns("browser"), "browser")
  )  
      
}

mod_plotsServer <- function(id, data_long, selected_ids, id_type, plot_colours = c("#7EC247","#53A2DA"), second_factor = FALSE, accession_col = "Accession", plot_height = "200px") {
  moduleServer(id, function(input, output, session) {
    
    ns_server <- NS(id)

    observeEvent(input$browser, browser())
    
    ids <- reactive(selected_ids()[[id_type]])
    
    tags_plot <- function(id, plot_name, plot_height = "200px"){
      if(!isTruthy(id)) {
        tags <- tagList(p(class = "no_data", "No data available", height = plot_height))
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
          fluidRow(class = "plotRow",
            column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
            column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
          )  
        )
      } else if (length(ids()) == 3){
        tags <- tagList(
          fluidRow(class = "plotRow",
            column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
            column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
          ),
          br(),
          fluidRow(class = "plotRow",
            column(width = 6, tags_plot(ids()[3], plot_name = "plot3"))
          )
        )
      } else if (length(ids()) == 4){
        tags <- tagList(
          fluidRow(class = "plotRow",
            column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
            column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
          ),
          br(),
          fluidRow(class = "plotRow",
            column(width = 6, tags_plot(ids()[3], plot_name = "plot3")),
            column(width = 6, tags_plot(ids()[4], plot_name = "plot4"))
          )
        )
      }
      tags
    })
    
    ## protein acid plots ----
    acid_boxplot <- function(data, title, box_colour, second_factor = FALSE){

      if(isTruthy(second_factor)){
        ggplot(data, aes(x = condition, y = value, fill = condition, colour = .data[[second_factor]])) +
          geom_boxplot(lwd = 1.2, fatten = 0.5) +
          xlab("") +
          scale_fill_manual(values = c("#7EC247", "#53A2DA")) +
          scale_colour_manual(values = c("black", "red4", "blue4")) +
          ggtitle(title)
      } else {
        ggplot(data, aes(x = condition, y = value, fill = condition)) +
          geom_boxplot() +
          xlab("") +
          scale_fill_manual(values = c("#7EC247", "#53A2DA")) +
          ggtitle(title)
      }
    }
    
    output$plot1 <- renderPlot({
      
      id <- ids()[1]
      req(id)

      data_filt <- data_long() %>%
        filter(.data[[accession_col]] == id)
      
      acid_boxplot(data_filt, id, plot_colours, second_factor = second_factor)
        
    }) %>% bindCache(ids()[1], data_long()) # could probably pass the name of the dataset - would be 
    # more efficient to compare than the whole dataset
    
    output$plot2 <- renderPlot({
      
      req(ids()[2])
      id <- ids()[2]
      
      data_filt <- data_long() %>%
        filter(.data[[accession_col]] == id)
      
      acid_boxplot(data_filt, id, plot_colours, second_factor = second_factor)
      
    }) %>% bindCache(ids()[2], data_long())
    
    output$plot3 <- renderPlot({
      
      req(ids()[3])
      id <- ids()[3]
      
      data_filt <- data_long() %>%
        filter(.data[[accession_col]] == id)
      
      acid_boxplot(data_filt, id, plot_colours, second_factor = second_factor)
    }) %>% bindCache(ids()[3], data_long())
    
    output$plot4 <- renderPlot({
      
      req(ids()[4])
      id <- ids()[4]
      
      data_filt <- data_long() %>%
        filter(.data[[accession_col]] == id)
      
      acid_boxplot(data_filt, id, plot_colours, second_factor = second_factor)
    }) %>% bindCache(ids()[4], data_long())
  })   
}            