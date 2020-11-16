#' ocms_reannotateTax
#'
#' Reannotate taxonomy
#'
#' Reannotates taxonomy table so that "unclassfied" assignments include higher
#' level classifications. This helps preserve the biological meaning of an unclassfied genus (as it could be classfied at the Family level).
#'
#' @param taxonomy taxonomy table. Must have \code{c('Kingdom','Phylum','Class','Order','Family','Genus','Species')} as columns. unclassified entries are denoted as NA
#'
#' @import dplyr
#' @export
#' @example
#' @return updated taxonomy table with "unclassified" prepended with higher level classifications
#'

ocms_reannotateTax <- function(taxonomy) {

  tax_level <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

  # check for columns ----------------------------------------------------------
  if(any(!tax_level %in% colnames(taxonomy))) {
    stop("taxonomy table must have c('Kingdom','Phylum','Class','Order','Family','Genus','Species') as columns")
  }

  # make sure all values are character
  taxonomy <- mutate_all(taxonomy, as.character)

  # change NA taxonomy in "[level-up-taxonomy]_unclassified"--------------------
  # work with one column at a time -- not checking Kingdom level
  for(i in 2:length(tax_level)) {

    # find row with na in current tax_level
    na_ind <- which(is.na(taxonomy[tax_level[i]]))

    if(length(na_ind) != 0) {

      # look at column before
      curr <- taxonomy[, c(tax_level[i-1], tax_level[i])]

      # make updated taxonomy labels
      curr <- curr %>%
        mutate(updated = ifelse(is.na(.data[[tax_level[i]]]), # if current taxon is NA
                                # prefix with prev level
                                paste(.data[[tax_level[i-1]]], 'unclassified',
                                      sep = '_'),
                                .data[[tax_level[i]]])) # else keep as current taxaon

      # update entire row taxonomy table
      taxonomy[na_ind, tax_level[i:length(tax_level)]] <- curr$updated[na_ind]
    }

  }

  return(taxonomy)
}
