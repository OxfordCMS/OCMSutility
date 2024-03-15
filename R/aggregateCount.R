#' aggregateCount
#'
#' @description
#' Aggregate read counts based on taxonomic level
#' `r lifecycle::badge("deprecated")`r
#'
#' `aggregateCount` was renamed to `aggregate_count` for more consistent
#' function naming nomenclature in `OCMSutility`. Usage is still exactly the same
#' @keywords internal
#' @export
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
#'
#' aggregated_list <- aggregateCount(feature_count, feature_tax,
#'                                   aggregate_by = "Family")
#' # ->
#' aggregated_list <- aggregate_count(feature_count, feature_tax,
#'                                   aggregate_by = "Family")
#'
#' summary(aggregated_list)

aggregateCount <-  function(count_df, tax_df=NULL, aggregate_by = NULL) {
  lifecycle::deprecate_warn("0.3.0", "aggregateCount()", "aggregate_count()")
  aggregate_count(count_df, tax_df, aggregate_by)
}