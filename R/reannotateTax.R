#' reannotateTax
#'
#' @description
#' Reannotate taxonomy
#' `r lifecycle::badge("deprecated")`
#'
#' `reannotateTax` was renamed to `reannotate_tax` for more consistent
#' function naming nomenclature in `OCMSutility`. Usage is still exactly the same
#' @keywords internal
#' Reannotates taxonomy table so that "unclassfied" assignments include higher
#' level classifications. This helps preserve the biological meaning of an unclassfied genus (as it could be classfied at the Family level).
#'
#' @param taxonomy taxonomy table. Must have \code{c('Kingdom','Phylum','Class','Order','Family','Genus','Species')} as columns. unclassified entries are denoted as NA
#'
#' @import dplyr
#' @export
#' @examples
#' data(asv_example)
#' # adding Kingdom column; removing sequence column because don't need asv IDs in this example
#' old_tax <- tax_example
#' colnames(old_tax)[1] <- 'Kingdom'
#' old_tax$Kingdom <- 'Bacteria'
#' knitr::kable(head(old_tax))
#'
#' new_tax <- reannotateTax(old_tax)
#' # ->
#' new_tax <- reannotate_tax(old_tax)
#' knitr::kable(head(new_tax))
#' @return updated taxonomy table with "unclassified" prepended with higher level classifications

reannotateTax <- function(taxonomy) {

  reannotate_tax(taxonomy)
}
