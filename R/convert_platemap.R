#' convert_platemap
#'
#' @description
#' function to convert excel plate map to long data frame with map locations.
#' uses \code{readxl} to read in excel file so additional arguments
#' are passed to \code{readxl::read_excel}
#'
#' @details
#' When `from_file=FALSE`, plate rows (e.g. A, B, C etc.) are set as rownames.
#'
#' @param from_file logical. Set to `TRUE` if supplying platemap from excel file.
#'                  Set to `FALSE` if supplying a `data.frame`. Default TRUE.
#' @param plate_map excel file containing plate map
#' @param sheet name or number of excel sheet. defaults to first sheet
#' @param map_range string. cell range containing plate map format
#'                  using excel notation (e.g. 'A1:H12')
#' @param drop_empty boolean. default TRUE. drop unlabeled wells
#' @param ... Arguments passed to \code{readxl::read_excel}
#' @returns datafame of platemap in long form
#' @export
#' @import readxl
#' @import dplyr
#' @import tibble

convert_platemap <- function(from_file = TRUE, plate_map, sheet=NULL, map_range, drop_empty = TRUE,
                             ...) {

  if(from_file) {
    # check inputs----------------------------------------------------------------
    if(!file.exists(plate_map)) {
      stop("plate_map not found")
    }

    # read in file----------------------------------------------------------------
    raw <- read_excel(plate_map, sheet = sheet, range=map_range,
                      col_names = FALSE, ...)

    # report plate size
    msg <- sprintf("Reading in plate map for %s well plate", ncol(raw)*nrow(raw))
    message(msg)

    # add well names
    colnames(raw) <- 1:ncol(raw)
    raw$row <- LETTERS[1:nrow(raw)]

  } else {
    raw <- plate_map %>%
      rownames_to_column('row')
  }

  # reshape to long
  map_df <- raw %>%
    gather('col','well_value',-row) %>%
    mutate(well_id = paste0(row, col))

  if(drop_empty) {
    map_df <- map_df %>%
      filter(!is.na(well_value))
  }

  return(map_df)

}