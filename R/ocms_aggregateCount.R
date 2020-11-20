#' ocms_aggregateCount
#'
#' Aggregate read counts based on taxonomic level
#'
#' @param count_df dataframe. count table with samples in columns and
#'     ASV in rows. First column is 'featureID'
#' @param tax_df dataframe. featureID must match \code{count_df}.
#'     has columns \code{'featureID','Kingdom','Phylum',
#'     'Class','Order','Family','Genus','Species'}
#' @param aggregate_by Aggregate counts by taxonomic level.
#'     Set to \code{NULL} to keep reads at ASV level. default \code{NULL} .
#'     Must be one of \code{c('Kingdom','Phylum','Class','Order','Family',
#'     'Genus','Species')}.
#'
#' @export
#' @import stringr
#' @return list: aggregated count table (\code{count_df}) and
#'     updated taxonomy table (\code{tax_df})
#'
#' @details
#' If \code{tax_df} has sequences column and features are aggregated,
#' the sequence column is set to NA. If tax_df has Taxon column with
#' notation of p__phylum;c__class;o__order... then those are truncated
#' at the aggregation level.
#'
#' All taxa downstream of aggregation level are set to NA
#'
#' The aggregated taxonomy table has an additional column, \code{n_collapse},
#' which is the number of ASVs that were aggregated.
#'

ocms_aggregateCount <-  function(count_df, tax_df, aggregate_by = NULL) {

  tax_level <- c('Kingdom','Phylum','Class','Order', 'Family','Genus',
                 'Species','featureID')

  # check inputs----------------------------------------------------------------
  # count_df must be dataframe
  if(class(count_df) != 'dataframe') {
    stop("count_df must be dataframe")
  }

  # check order of count_df is same as order of tax_df
  if(!identical(count_df$featureID, tax_df$featureID)) {
    stop("Order of features in count_df and tax_df must be identical")
  }

  # aggregate_by must be one of "Kingdom","Phylum","Class","Order",
  # "Family","Genus","Species", NULL
  if(is.null(aggregate_by)) {
    aggregate_by <- 'featureID'
  } else {
    if(!aggregate_by %in% tax_level) {
      stop("Reads can only be aggregated by 'Kingdom','Phylum','Class','Family','Genus','Species'.")
    }
  }

  # aggregate counts -----------------------------------------------------------

  sampleID <- colnames(count_df)
  sampleID <- sampleID[sampleID != 'featureID']

  # set featureID in count_df to aggregation level
  count_df$featureID <- tax_df[,aggregate_by]

  # build formula
  yvar <- paste(sampleID, collapse=',')
  f <- sprintf("cbind(%s) ~ featureID", yvar)

  # perform aggregation
  aggregated <- aggregate(formula = formula(f), data = count_df, FUN = sum)

  # update taxonomy table-------------------------------------------------------

  if(aggregate_by != 'featureID') {
    tax_column <- colnames(tax_df)

    # set sequence column to NA
    if('sequences' %in% tax_column) {
      tax_df$sequence <- NA
    }

    # set aggregated level and all lower levels to NA
    ind <- which(tax_level == aggregate_by) + 1
    to_na <- tax_level[ind:length(tax_level)]

    tax_df[, to_na] <- NA

    # update Taxon column
    if('Taxon' %in% tax_column) {
      ind <- which(tax_level == aggregate_by)
      Taxon <- stringr::str_split(tax_df$Taxon, ";", simplify = TRUE)
      if(ind == 1) {
        Taxon <- Taxon[,ind]
      } else {
        Taxon <- apply(Taxon[,1:ind], 1, paste, collapse = ";")
      }
      tax_df$Taxon <- Taxon
    }

    # record how many ASVs aggregated
    n_collapse <- tax_df %>%
      group_by(.data[[aggregate_by]]) %>%
      mutate(n_collapse = n(),
             featureID = .data[[aggregate_by]])

    # remove redundant entries
    tax_df <- unique(tax_df)
  }

  return(list('count_df'=aggregated, 'tax_df' = tax_df))

}