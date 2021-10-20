browser_buttons <- TRUE
## Try just having one layout!!

mod_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
      uiOutput(ns("plot_panel"), class = "plot_box"),
      if(browser_buttons) actionButton(ns("browser"), "browser")
  )  
      
}

mod_plotsServer <- function(id, data_long, selected_ids, id_type, panel_name, title_id = "Gene_id",  plot_colours = c("#7EC247","#53A2DA"), second_factor = FALSE, accession_col = "Accession", plot_height = "200px") {
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
      tagList(tags, downloadButton(outputId = ns_server("download_plots"), label = "Download pdf"))
    })
    
    output$download_plots <- downloadHandler(
      
      filename = function() {
        all_ids <- unique(pull(filtered_data(), .data[[title_id]]))
        id_text <- paste0(all_ids, collapse = "_")
        paste0(id_text, "_", panel_name, ".pdf")
      },
      content = function(file) {
        n_plots <- length(plot_list())
        ggsave(file, 
               gridExtra::marrangeGrob(
                 grobs = plot_list(), 
                 nrow = ceiling(n_plots / 2), 
                 ncol = if_else(n_plots == 1, 1, 2)
                )
               )
      }
    )
    
    # we get an error if we call one of the gg_plot objects if they don't exist,
    # so using a rather inelegant solution here
    plot_list <- reactive({
      req(ids())
      no_ids <- length(ids())
      
      x <- switch(no_ids,
             list(gg_plot1()),
             list(gg_plot1(), gg_plot2()),
             list(gg_plot1(), gg_plot2(), gg_plot3()),
             list(gg_plot1(), gg_plot2(), gg_plot3(), gg_plot4())
            )
      Filter(Negate(is.null), x)
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
    
    filtered_data <- reactive({
      req(ids())
      filter(data_long(), .data[[accession_col]] %in% ids())
    })
    
    gg_plot1 <- reactive({
      if(!isTruthy(ids()[1])) return (NULL)
      id <- ids()[1]
      req(id)
      data_filt <- filter(filtered_data(), .data[[accession_col]] == id)
      plot_title <- pull(data_filt, .data[[title_id]])[1]
      p <- acid_boxplot(data_filt, title = plot_title, plot_colours, second_factor = second_factor)
      p
    })
    
    gg_plot2 <- reactive({
      if(!isTruthy(ids()[2])) return (NULL)
      req(ids()[2])
      id <- ids()[2]
      data_filt <- filter(data_long(), .data[[accession_col]] == id)
      plot_title <- pull(data_filt, .data[[title_id]])[1]
      acid_boxplot(data_filt, title = plot_title, plot_colours, second_factor = second_factor)
    }) 
    
    gg_plot3 <- reactive({
      if(!isTruthy(ids()[3])) return (NULL)
      else {
        req(ids()[3])
        id <- ids()[3]
        data_filt <- filter(data_long(), .data[[accession_col]] == id)
        plot_title <- pull(data_filt, .data[[title_id]])[1]
        acid_boxplot(data_filt, title = plot_title, plot_colours, second_factor = second_factor)
      }
    }) 
    
    gg_plot4 <- reactive({
      if(!isTruthy(ids()[4])) return (NULL)
      req(ids()[4])
      id <- ids()[4]
      data_filt <- filter(data_long(), .data[[accession_col]] == id)
      plot_title <- pull(data_filt, .data[[title_id]])[1]
      acid_boxplot(data_filt, title = plot_title, plot_colours, second_factor = second_factor)
    }) 
     
    output$plot1 <- renderPlot(gg_plot1()) %>% bindCache(ids()[1], data_long()) # could probably pass the name of the dataset - would be 
    # more efficient to compare than the whole dataset
    
    output$plot2 <- renderPlot(gg_plot2()) %>% bindCache(ids()[2], data_long())
    
    output$plot3 <- renderPlot(gg_plot3()) %>% bindCache(ids()[3], data_long())
    
    output$plot4 <- renderPlot(gg_plot4()) %>% bindCache(ids()[4], data_long())
  })   
}            