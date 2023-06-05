#' nsample_by_var
#'
#' @description
#' Tallies the number of samples per identifier for metadata variables
#' by counting the number of non-NA values for each metadata variable
#' for a given identifier.
#' This can be helpful in time course data where you want to check the
#' number of samples per individual. This can also be used to check how complete
#' your metadata is.
#'
#' @param ddata dataframe. samples in rows. with metadata variables in columns
#' @param id character. sample identifier
#' @param var character vector. metadata variables to be tallied
#'
#' @returns dataframe. tally of number of id per variable
#'
#' @import tibble
#' @import dplyr
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' # time course data: checking number of samples per patient for each metadata
#' # variable
#' df <- data.frame(sample_id = paste0("sample", 1:100),
#'                  patient_id = rep(LETTERS[1:25], 4),
#'                  var1 = sample(c(rnorm(30, 10, 0.5), rnorm(40, 25, 2),
#'                                  rep(NA, 30)), 100),
#'                  var2 = sample(c(rnorm(65, 0.5, 0.01),
#'                                  rep(0, 20),
#' rep(NA, 15)), 100),
#'                  var3 = sample(c(letters[1:5], NA), 100, replace=TRUE))
#'
#' nsample_by_var(df, 'patient_id', c('var1','var2','var3'))

nsample_by_var <- function(ddata, id, var) {

  # only working with specified metadata variables
  df <- ddata[c(id, var)]

  out <- distinct(ddata, !!sym(id))
  # work on one variable at a time
  for( v in var) {
    entry <- df[,c(id, v)] %>%
      group_by(!!sym(id)) %>%
      reframe(n_sample = sum(!is.na(!!sym(v))))
    colnames(entry)[2] <- v

    out <- left_join(out, entry, id)
  }
  return(out)
}