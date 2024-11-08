#' get_shortnames
#' return the shortnames of taxonomy labels
#'
#' @param longnames vector; vector of long taxonomic labels
#' @export
#' @return vector
#' @examples
#' get_shortnames(rownames(abundance_matrix))

get_shortnames <- function(longnames){
  return(trimws(longnames, whitespace = ".*__"))
}
