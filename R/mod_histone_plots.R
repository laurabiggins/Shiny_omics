## Try just having one layout!!

mod_histone_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("plot_panel"), class = "plot_box"),
    #br(),
    #actionButton(ns("browser"), "browser")
  )  
  
}

#' Title
#'
#' @param id 
#' @param data_long   long data
#' @param selected_ids ids that have been selected
#' @param id_type 
#' @param panel_name name of panel used in output filename for download plots
#' @param title_id   id type to use for individual plots - most are the gene ids, except for histones
#' @param plot_legend  TRUE or FALSE whether to display legend for plots
#' @param plot_colours  plot colours
#' @param second_factor  default FALSE - for histones, a vector of media is supplied
#' @param accession_col  column name that contains the accessions 
#' @param plot_height  default 200px plot height 
#'
#' @return
#' @export
#'
#' @examples
mod_histone_plotsServer <- function(id, data_long, panel_name, ylabel) {
  moduleServer(id, function(input, output, session) {
    
    ns_server <- NS(id)
    
    custom_plot_colours <- c(
      Naive          = "#7EC247", 
      Primed         = "#53A2DA", 
      `Naive+PRC2i`  = "#C8E5B0", 
      `Primed+PRC2i` = "#B5D7EF"
    )
    col_scale <- ggplot2::scale_fill_manual(name = "condition", values = custom_plot_colours)
    
    observeEvent(input$browser, browser())
    
    plot_height <- reactive({
      n <- dplyr::n_distinct(data_long()$histone_mark)
      if(n == 3) 200
      else {
        width <- ceiling(sqrt(n))
        ceiling(n/width)*200
      }
    })

    plot_width <- reactive({
      n <- dplyr::n_distinct(data_long()$histone_mark)
      if(n == 1) 350
      else "auto"
    })
    
    y_scale <- reactive(dplyr::if_else(input$common_scale, "fixed", "free"))
    
    ## renderUI plot panel ----
    output$plot_panel <- renderUI({
      
      tagList(
        wellPanel(
          class = "plot_box", 
          plotOutput(ns_server("facet_plot"), height = plot_height()),
          br(),
          fluidRow(
            column(6,
                   downloadButton(outputId = ns_server("download_plots"), label = "Download pdf")
            ),
            column(6,
                   checkboxInput(inputId = ns_server("common_scale"), label = "common scale")
            )
          )
        ) 
      )
    }) 
 
    plot_object <- reactive({
      if(nrow(data_long()) == 0) return (NULL)
      data_long() %>%
        ggplot2::ggplot(ggplot2::aes(x = condition, y = value, fill = condition, colour = medium)) +
        ggplot2::geom_boxplot(lwd = 1.2, fatten = 0.5) +
        ggplot2::scale_colour_manual(values = c("black", "red4", "blue4")) +
        col_scale +
        ggplot2::ylab("normalised_abundance") +
        ggplot2::guides(fill = "none") +
        ggplot2::facet_wrap(~histone_mark, scales = y_scale())
    })
    
    output$facet_plot <- renderPlot({
      req(plot_object())
      plot_object()
    }, height = plot_height(), width = plot_width())
    
    output$download_plots <- downloadHandler(

      filename = function() {
        all_ids <- unique(data_long()$histone_mark)
        id_text <- paste0(all_ids, collapse = "_")
        paste0("histones_", id_text, ".pdf")
      },
      content = function(file) {
        pdf(file, onefile = FALSE)
        print(plot_object())
        dev.off()
      }
    )
  })
}
            