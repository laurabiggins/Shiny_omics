browser_buttons <- TRUE
## Try just having one layout!!

mod_histone_plotsUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("plot_panel"), class = "plot_box"),
    br(),
    if(browser_buttons) actionButton(ns("browser"), "browser")
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
    col_scale <- scale_fill_manual(name = "condition", values = custom_plot_colours)
    
    observeEvent(input$browser, browser())
    
    plot_height <- reactive({
      n <- n_distinct(data_long()$histone_mark)
      if(n == 3) 200
      else {
        width <- ceiling(sqrt(n))
        ceiling(n/width)*200
      }
    })

    plot_width <- reactive({
      n <- n_distinct(data_long()$histone_mark)
      if(n == 1) 350
      else "auto"
    })
    
    
    y_scale <- reactive(dplyr::if_else(input$common_scale, "fixed", "free"))

    
    # ## renderUI plot panel ----
    output$plot_panel <- renderUI({
      
      tagList(
        wellPanel(
          class = "plot_box", 
          plotOutput(ns_server("facet_plot"), height = plot_height()),
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
      data_long() %>%
        ggplot(aes(x = condition, y = value, fill = condition, colour = medium)) +
        geom_boxplot(lwd = 1.2, fatten = 0.5) +
        scale_colour_manual(values = c("black", "red4", "blue4")) +
        col_scale +
        ylab("normalised_abundance") +
        guides(fill = "none") +
        facet_wrap(~histone_mark, scales = y_scale())
    })
    
    output$facet_plot <- renderPlot(plot_object(), height = plot_height(), width = plot_width())
    
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
    #   
    #   #ids <- selected_ids[[id_type]]
    #   
    #   if(is.null(ids())) tags <- NULL
    #   
    #   req(ids())
    #   
    #   if(length(ids()) > 4){
    #     tags <- NULL
    #   } 
    #   
    #   if(length(ids()) == 1){
    #     tags <- tags_plot(ids()[1], plot_name = "plot1")
    #   } else if (length(ids()) == 2){
    #     tags <- tagList(
    #       fluidRow(
    #         column(width = 6, tags_plot(ids()[1], plot_name = "plot1")),
    #         column(width = 6, tags_plot(ids()[2], plot_name = "plot2"))
    #       )  
    #     )
    #   } else if (length(ids()) == 3){
    #     tags <- tagList(
    #       fluidRow(
    #         column(width = 6, class = "plotArea", tags_plot(ids()[1], plot_name = "plot1")),
    #         column(width = 6, class = "plotArea", tags_plot(ids()[2], plot_name = "plot2")),
    #         column(width = 6, class = "plotArea", tags_plot(ids()[3], plot_name = "plot3"))
    #       )
    #     )
    #   } else if (length(ids()) == 4){
    #     tags <- tagList(
    #       fluidRow(
    #         column(width = 6, class = "plotArea", tags_plot(ids()[1], plot_name = "plot1")),
    #         column(width = 6, class = "plotArea", tags_plot(ids()[2], plot_name = "plot2")),
    #         column(width = 6, class = "plotArea", tags_plot(ids()[3], plot_name = "plot3")),
    #         column(width = 6, class = "plotArea", tags_plot(ids()[4], plot_name = "plot4"))
    #       )
    #     )
    #   }
    #   tagList(
    #     wellPanel(
    #       class = "plot_box", 
    #       tags,
    #       br(),
    #       fluidRow(
    #         column(6, 
    #                downloadButton(outputId = ns_server("download_plots"), label = "Download pdf")
    #         ),
    #         column(6, 
    #                checkboxInput(inputId = ns_server("common_scale"), label = "common scale")
    #         )
    #       )
    #     )
    #   )
    #   
    #   #tagList(tags, downloadButton(outputId = ns_server("download_plots"), label = "Download pdf"))
    # })
    # 
    # 
    # 
    # 
    # # we get an error if we call one of the gg_plot objects if they don't exist,
    # # so using a rather inelegant solution here
    # arranged_plots <- reactive({
    #   req(ids())
    #   no_ids <- length(ids())
    #   
    #   x <- switch(no_ids,
    #               list(gg_plot1()),
    #               list(gg_plot1(), gg_plot2()),
    #               list(gg_plot1(), gg_plot2(), gg_plot3()),
    #               list(gg_plot1(), gg_plot2(), gg_plot3(), gg_plot4())
    #   )
    #   plot_list <- Filter(Negate(is.null), x)
    #   
    #   n_plots <- length(plot_list)
    #   
    #   gridExtra::marrangeGrob(
    #     grobs = plot_list, 
    #     nrow = 2,
    #     #nrow = ceiling(n_plots / 2), 
    #     ncol = 2, #if_else(n_plots == 1, 1, 2),
    #     layout_matrix = t(matrix(1:4, nrow= 2, ncol = 2))
    #   )
    # })
    # 
    # # ## box plot function ----
    # custom_boxplot <- function(data, title, box_colour, ylabel = "custom label", second_factor = FALSE){
    # 
    #   if(isTruthy(second_factor)){
    #     p <- ggplot(data, aes(x = condition, y = value, fill = condition, colour = .data[[second_factor]])) +
    #       geom_boxplot(lwd = 1.2, fatten = 0.5) +
    #       scale_colour_manual(values = c("black", "red4", "blue4"))
    #   } else {
    #     p <- ggplot(data, aes(x = condition, y = value, fill = condition)) +
    #       geom_boxplot() +
    #       theme(legend.position = "none")
    #   }
    #   p <- p +
    #     xlab("") +
    #     col_scale +
    #     ggtitle(title) +
    #     ylab(ylabel)
    # 
    #   if(input$common_scale){
    #     p + ylim(filtered_data_min_max())
    #   } else p
    # 
    # }
    # 

      
    
    # 
    # 
    # ## filtered data ----
    # 
    # filtered_data <- reactive({
    #   req(ids())
    #   filter(data_long(), .data[[accession_col]] %in% ids())
    # })
    # 
    # ## filtered data range ----
    # filtered_data_min_max <- reactive({
    #   req(filtered_data)
    #   range(filtered_data()$value)
    # })
    # 
    # 
    # ## plot reactives ----
    # gg_plot1 <- reactive({
    #   if(!isTruthy(ids()[1])) return (NULL)
    #   id <- ids()[1]
    #   req(id)
    #   data_filt <- filter(filtered_data(), .data[[accession_col]] == id)
    #   plot_title <- pull(data_filt, .data[[title_id]])[1]
    #   custom_boxplot(data_filt, title = plot_title, ylabel = ylabel, plot_colours, second_factor = second_factor)
    # })
    # 
    # gg_plot2 <- reactive({
    #   if(!isTruthy(ids()[2])) return (NULL)
    #   req(ids()[2])
    #   id <- ids()[2]
    #   data_filt <- filter(data_long(), .data[[accession_col]] == id)
    #   plot_title <- pull(data_filt, .data[[title_id]])[1]
    #   custom_boxplot(data_filt, title = plot_title, ylabel = ylabel, plot_colours, second_factor = second_factor)
    # }) 
    # 
    # gg_plot3 <- reactive({
    #   if(!isTruthy(ids()[3])) return (NULL)
    #   else {
    #     req(ids()[3])
    #     id <- ids()[3]
    #     data_filt <- filter(data_long(), .data[[accession_col]] == id)
    #     plot_title <- pull(data_filt, .data[[title_id]])[1]
    #     custom_boxplot(data_filt, title = plot_title, ylabel = ylabel, plot_colours, second_factor = second_factor)
    #   }
    # }) 
    # 
    # gg_plot4 <- reactive({
    #   if(!isTruthy(ids()[4])) return (NULL)
    #   req(ids()[4])
    #   id <- ids()[4]
    #   data_filt <- filter(data_long(), .data[[accession_col]] == id)
    #   plot_title <- pull(data_filt, .data[[title_id]])[1]
    #   custom_boxplot(data_filt, title = plot_title, ylabel = ylabel, plot_colours, second_factor = second_factor)
    # }) 
    # 
    # output$plot1 <- renderPlot(gg_plot1()) %>% bindCache(ids()[1], input$common_scale)#, data_long()) # could probably pass the name of the dataset - would be 
    # # more efficient to compare than the whole dataset
    # 
    # output$plot2 <- renderPlot(gg_plot2()) %>% bindCache(ids()[2], input$common_scale)
    # 
    # output$plot3 <- renderPlot(gg_plot3()) %>% bindCache(ids()[3], input$common_scale)
    # 
    # output$plot4 <- renderPlot(gg_plot4()) %>% bindCache(ids()[4], input$common_scale)
            