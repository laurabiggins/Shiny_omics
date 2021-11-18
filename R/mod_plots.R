browser_buttons <- FALSE

# plot colours
col_scale <- ggplot2::scale_fill_manual(
  name = "condition", 
  values = c(
    Naive = "#7EC247", 
    Primed = "#53A2DA", 
    `Naive+PRC2i` = "#C8E5B0", 
    `Primed+PRC2i` = "#B5D7EF"
  )
)


mod_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
      uiOutput(ns("plot_panel"), class = "plot_box"),
      if(browser_buttons) actionButton(ns("browser"), "browser")
  )  
      
}

#' Title
#'
#' @param id 
#' @param filtered_data dataset to plot
#' @param selected_ids ids that have been selected
#' @param ylabel label to use for y axis
#' @param panel_name name of panel used in output filename for download plots
#'
#' @return
#' @export
#'
#' @examples
mod_plotsServer <- function(
  id, 
  filtered_data,
  selected_ids, 
  panel_name, 
  ylabel
) {
  moduleServer(id, function(input, output, session) {
    
    ns_server <- NS(id)

    observeEvent(input$browser, browser())
    
    ## renderUI plot panel ----
    output$plot_panel <- renderUI({

      tags <- tagList(
        fluidRow(purrr::map(1:n_plots(), individual_plot))
      )
      tagList(
        wellPanel(class = "plot_box", tags, br(),
          fluidRow(
            column(6, downloadButton(outputId = ns_server("download_plots"), label = "Download pdf")),
            column(6, checkboxInput(inputId = ns_server("common_scale"), label = "common scale")        )
          )
        )
      )
    })
    
    ## renderPlots ----
    output$plot1 <- renderPlot(gg_plot1()) %>% 
      bindCache(filtered_data(), selected_ids[1], input$common_scale)
    output$plot2 <- renderPlot(gg_plot2()) %>% 
      bindCache(filtered_data(), selected_ids[2], input$common_scale)
    output$plot3 <- renderPlot(gg_plot3()) %>% 
      bindCache(filtered_data(), selected_ids[3], input$common_scale)
    output$plot4 <- renderPlot(gg_plot4()) %>% 
      bindCache(filtered_data(), selected_ids[4], input$common_scale)
    
    ## reactives ----
    ### n_plots ----
    n_plots <- reactive(length(selected_ids))
    
    ### filtered data range ----
    filtered_data_min_max <- reactive({
      req(filtered_data())
      range(filtered_data()$value)
    })
    
    ### gg_plot reactives ----
    gg_plot1 <- reactive(ggplot_object(1, filtered_data()))
    gg_plot2 <- reactive(ggplot_object(2, filtered_data()))
    gg_plot3 <- reactive(ggplot_object(3, filtered_data()))
    gg_plot4 <- reactive(ggplot_object(4, filtered_data()))
    
    ## functions ----
    ### UI plot functions ----
    tags_plot <- function(id, plot_name, plot_height = "200px"){
      if(!isTruthy(id)) {
        tags <- tagList(div(class = "no_data", "No data available", height = plot_height))
      } else {
        tags <- tagList(
          plotOutput(ns_server(plot_name), height = plot_height)
        )
      }
      tags 
    }
    
    individual_plot <- function(i){
      column(width = 6,  class = "plotArea", tags_plot(selected_ids[i], plot_name = paste0("plot", i)))
    }
    
    ### server plot functions ----
    custom_boxplot <- function(data, title, ylabel = "custom label"){
      p <- ggplot(data, aes(x = condition, y = value, fill = condition)) +
        geom_boxplot() +
        theme(legend.position = "none") + 
        xlab("") +
        col_scale +
        ggtitle(title) +
        ylab(ylabel)
      
      if(input$common_scale){
        p + ylim(filtered_data_min_max())
      } else p
    }
    
    ggplot_object <- function(index, dataset){
      data_filt <- dataset %>% filter(.[[1]] == selected_ids[index])
      plot_title <- if_else(is.na(data_filt$Gene_id[1]), data_filt[[1]][1], data_filt$Gene_id[1])  
      custom_boxplot(data_filt, title = plot_title, ylabel = ylabel)
    }
    
    # download plots ----
    output$download_plots <- downloadHandler(
      
      filename = function() {
        #id_text <- paste0(plot_names, collapse = "_")
        paste0(panel_name, ".pdf")
      },
      content = function(file) {
        pdf(file, onefile = FALSE)
        print(arranged_plots())
        dev.off()
      }
    )
    
    # this is for the plot download
    # we get an error if we call one of the gg_plot objects if they don't exist,
    # so using a rather inelegant solution here
    arranged_plots <- reactive({
      x <- switch(n_plots(),
                  list(gg_plot1()),
                  list(gg_plot1(), gg_plot2()),
                  list(gg_plot1(), gg_plot2(), gg_plot3()),
                  list(gg_plot1(), gg_plot2(), gg_plot3(), gg_plot4())
      )
      plot_list <- Filter(Negate(is.null), x)
      
      gridExtra::marrangeGrob(
        grobs = plot_list, 
        nrow = 2,
        #nrow = ceiling(n_plots / 2), 
        ncol = 2, #if_else(n_plots == 1, 1, 2),
        layout_matrix = t(matrix(1:4, nrow= 2, ncol = 2))
      )
    })
    
  })   
}        
