#' aggregateCount
#'
#' Aggregate read counts based on taxonomic level
#'
#' @param count_df dataframe. count table with samples in columns and
#'     ASV in rows. feature IDs in row names.
#' @param tax_df dataframe. featureID column should
#'     match rownames of \code{count_df}.
#'     Has columns \code{'Kingdom','Phylum',
#'     'Class','Order','Family','Genus','Species'}.
#'     Can also have \code{Taxon} column (see details).
#' @param aggregate_by Aggregate counts by taxonomic level.
#'     Set to \code{NULL} to keep reads at ASV level. default \code{NULL} .
#'     Must be one of \code{c('Kingdom','Phylum','Class','Order','Family',
#'     'Genus','Species')}.
#'
#' @export
#' @import stringr
#' @import tidyr
#' @import tibble
#' @import dplyr
#'
#' @returns
#' returns a list of aggregated count table
#'     (\code{count_df}) and updated taxonomy table (\code{tax_df})
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
#'
#' @examples
#' data(dss_example)
#' # featureID should be row names
#' feature_count <- dss_example$merged_abundance_id %>%
#'   tibble::column_to_rownames('featureID')
#'
#' # cleanup sample names
#' colnames(feature_count) <- paste0('id', colnames(feature_count))
#' feature_tax <- dss_example$merged_taxonomy
#'
#' # set row order of count and tax tables to be the same
#' feature_count <- feature_count[feature_tax$featureID,]
#' aggregated_list <- aggregateCount(feature_count, feature_tax,
#'                                   aggregate_by = "Family")
#'
#' summary(aggregated_list)

aggregateCount <-  function(count_df, tax_df=NULL, aggregate_by = NULL) {

  tax_level <- c('featureID','Kingdom','Phylum','Class',
                 'Order','Family','Genus','Species')

  # check inputs----------------------------------------------------------------
  # count_df must be dataframe
  if(class(count_df) != 'data.frame') {
    stop("count_df must be data.frame")
  }

  if(class(tax_df) != 'data.frame') {
    stop("tax_df must be a data.frame")
  }

  # aggregate_by must be one of "Kingdom","Phylum","Class","Order",
  # "Family","Genus","Species"
  if(!aggregate_by %in% tax_level) {
    stop("Reads can only be aggregated by 'Kingdom','Phylum','Order','Class','Family','Genus','Species'.")
  }

  # check order of count_df is same as order of tax_df
  if(!identical(rownames(count_df), tax_df$featureID)) {
    stop("Order of features in count_df and tax_df must be identical")
  }

  # samples can't start with number
  if(any(grepl("^[0-9]", colnames(count_df)))) {
    stop('Samples cannot start with a number')
  }
  # aggregate counts -----------------------------------------------------------

  sampleID <- colnames(count_df)

  # set featureID column in count_df to aggregation level
  count_df$featureID <- tax_df[,aggregate_by]

  # build formula
  yvar <- paste(sampleID, collapse=',')
  f <- sprintf("cbind(%s) ~ featureID", yvar)
  # perform aggregation
  aggregated <- stats::aggregate(formula = formula(f),
                                 data = count_df, FUN = sum)

  aggregated <- aggregated %>%
    column_to_rownames('featureID')

  # update taxonomy table-------------------------------------------------------
  tax_column <- colnames(tax_df)

  # set sequence column to NA
  if('sequence' %in% tax_column) {
    tax_df$sequence <- NA
  }

  # set levels lower than aggregated level to NA
  ind <- which(tax_level == aggregate_by) + 1
  to_na <- tax_level[ind:length(tax_level)]

  tax_df[, to_na] <- NA

  # update Taxon column
  if('Taxon' %in% tax_column) {
    ind <- which(tax_level == aggregate_by)
    Taxon <- stringr::str_split(tax_df$Taxon, ";", simplify = TRUE)
    Taxon <- apply(Taxon[,1:ind], 1, paste, collapse = ";")
    tax_df$Taxon <- Taxon
  }

  # record how many ASVs aggregated
  tax_agg <- tax_df %>%
    group_by(.data[[aggregate_by]]) %>%
    mutate(n_collapse = n(),
              featureID = .data[[aggregate_by]]) %>%
    ungroup() %>%
    distinct()

  return(list('count_df' = aggregated, 'tax_df' = as.data.frame(tax_agg)))

}