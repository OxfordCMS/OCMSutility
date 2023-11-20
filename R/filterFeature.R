#' filterFeature.R
#'
#' @description
#' filter out reads based on cutoff threshold and asv prevalence across samples
#' `r lifecycle::badge("deprecated")`
#'
#' `filterFeature` was renamed to `filter_feature` for more consistent
#' function naming nomenclature in `OCMSutility`. Usage is still exactly the same
#' @keywords internal
#' @export
#' @examples
#' data(dss_example)
#'
#' # put featureID as rownames
#' tax_df <- dss_example$merged_taxonomy
#' count_df <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#' # set features in count tax to be in same order
#' count_df <- count_df[tax_df$featureID,]
#'
#' filtered_ls <- filterFeature(count_df, tax_df, 'percent_sample', 0.001, 2)
#' # ->
#' filtered_ls <- filter_eature(count_df, tax_df, 'percent_sample', 0.001, 2)
#' summary(filtered_ls)
#' filtered_count <- filtered_ls$filtered
#' dim(filtered_count)
#' head(filtered_count)

filterFeature <- function(count_df, tax_df,
                           filter_method = 'abs_count',
                           asv_cutoff = 1, prev_cutoff = 2) {

  filter_feature(count_df, tax_df, filter_method, asv_cutoff, prev_cutoff)
}
