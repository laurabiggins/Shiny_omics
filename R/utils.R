#' Set DT::datatable options
#'
#' @param data  A tibble or dataframe. 
#' @param n_rows The number of rows to display (integer).
#' @param lineHeight The height of each row in table (percentage). Default "80%". 
#' @param dom_opt dom options for table components: 
#' * l - length changing input control
#' * f - filtering input
#' * t - table
#' * i - table info 
#' * p - pagination     
#'  Any combination of these options should be supplied as a single character 
#'  string and will be displayed in the order supplied. 
#'  Default: "tip" (table, info, pagination).
#' @param cols_to_round columns to round (vector of column numbers). 
#'  See `sig_digits`. Default: NULL 
#' @param dt_options list of options to supply to `datatable()`.
#'  The default is NULL, meaning options are taken from the relevant supplied parameters, 
#'  i.e. `dom`, `pageLength = n_rows`, `lengthMenu = table_lengths`. 
#'  This argument is for if further customisation is required.
#' @param sig_digits  The number of significant digits to round to if `cols_to_round` 
#'  is not NULL. Default: 3
#' @param regex If `TRUE`, this allows regex searchng of the table e.g. gene1 | gene2.
#' @param selection Whether to allow multiple rows to be selected at once, or just a 
#'  single row. Either `"multiple"` or `"single"`. Default: `"single"`
#' @param table_lengths Integer vector of available table lengths. 
#'  Default: c(10,20,50,100). This will only show if dom_opt includes "l".
#' @param filter_pos Location of column filters. One of c("none", "bottom", "top"). 
#'  Default: "none"
#' @param show_rownames Whether to show rownames or not. Default: FALSE
#' @return DT::datatable() object
#' @export
#' @md
#' @examples
#' dt_setup(iris)
dt_setup <-  function(data, 
                      n_rows = 10, 
                      lineHeight = "50%", 
                      font_size = "90%",
                      dom_opt = "tip", 
                      cols_to_round = NULL, 
                      dt_options = NULL, 
                      sig_digits = 3, 
                      regex = FALSE, 
                      selection = "single", 
                      table_lengths = c(10,20,50,100), 
                      filter_pos = "none",
                      show_rownames = FALSE,
                      style = "bootstrap4") {
  
  assertthat::assert_that(
    tibble::is_tibble(data) | base::is.data.frame(data) | base::is.matrix(data), 
    msg = "data supplied to dt_setup must be a tibble or data frame"
  )
  
  if (is.null(dt_options)) {
    dt_options = list(dom = dom_opt, 
                      pageLength = n_rows,
                      lengthMenu = table_lengths)
  }
  
  if (regex) {
    dt_options[["search"]] <- list(regex = TRUE, caseInsensitive = TRUE)
  }
  
  dt_table <- DT::datatable(
    data,
    style = style,
    rownames = show_rownames,
    escape   = FALSE,
    filter   = filter_pos,
    options  = dt_options,
    selection = selection
    
  ) %>%
    DT::formatStyle(0, target = 'row', lineHeight = lineHeight, fontSize = font_size)
  
  if (!is.null(cols_to_round)) {
    dt_table <- DT::formatRound(dt_table, cols_to_round, sig_digits)
  }  
  dt_table   
}


dt_options <- list(
  dom = 'fltip',
  lengthMenu = c(5, 10, 20, 50),
  pageLength = 10,
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
      targets = 5,
      width = "800px",
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.length > 30 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
        "}")
    ),
    list(
      targets = c(3,4),
      width = "300px"
    ),
    list(
      targets = 5,
      width = "50px"
    )
  )
)


