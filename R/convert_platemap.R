#' convert_platemap
#'
#' function to convert excel plate map to long data frame with map locations.
#' uses \code{readxl} to read in excel file so additional arguments
#' are passed to \code{readxl::read_excel}
#'
#' @param map_file excel file containing plate map
#' @param sheet name or number of excel sheet. defaults to first sheet
#' @param map_range string. cell range containing plate map format
#'                  using excel notation (e.g. 'A1:H12')
#' @param drop_empty boolean. default TRUE. drop unlabeled wells
#' @param ... Arguments passed to \code{readxl::read_excel}
#' @returns datafame of platemap in long form
#'

convert_platemap <- function(map_file, sheet=NULL, map_range, drop_empty = TRUE,
                             ...) {

  # required libraries
  require(readxl)
  require(dplyr)
  require(tibble)

  # check inputs----------------------------------------------------------------
  if(!file.exists(map_file)) {
    stop("map_file not found")
  }

  # read in file----------------------------------------------------------------
  raw <- read_excel(map_file, sheet = sheet, range=map_range, col_names = FALSE, ...)

  # report plate size
  msg <- sprintf("Reading in plate map for %s well plate", ncol(raw)*nrow(raw))
  message(msg)

  # add well names
  colnames(raw) <- 1:ncol(raw)
  raw$row <- LETTERS[1:nrow(raw)]

  # reshape to long
  map_df <- raw %>%
    gather('col','sample_name',-row) %>%
    mutate(well_id = paste0(row, col))

  if(drop_empty) {
    map_df <- map_df %>%
      filter(!is.na(sample_name))
  }

  return(map_df)

}