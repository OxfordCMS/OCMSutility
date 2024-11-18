#' get_shortnames
#' return the shortnames of taxonomy labels
#'
#' @param longnames vector; vector of long taxonomic labels
#' @export
#' @return vector
#' @examples
#' longnames <- c("k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus",
#'                "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides")
#' get_shortnames(longnames)

get_shortnames <- function(longnames){
  return(trimws(longnames, whitespace = ".*__"))
}
