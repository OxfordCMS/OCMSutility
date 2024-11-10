#' prevalence_abundance
#' return a list that contains a dataframe of the prevalence and abundance (mean across samples)
#' of each taxon (from relative abundance matrix) and a plot.
#'
#' @param abundance_matrix dataframe; samples are columns and features are rows (relative abundance)
#' @import reshape2
#' @export
#' @return list
#' @examples
#' feature_count <- dss_example$merged_abundance_id %>%
#' tibble::column_to_rownames('featureID')
#'
#' colnames(feature_count) <- paste0('id', colnames(feature_count))
#' feature_tax <- dss_example$merged_taxonomy
#' feature_count <- feature_count[feature_tax$featureID,]
#' aggregated_list <- aggregate_count(feature_count, feature_tax,
#'                                    aggregate_by = "Genus")
#'
#' abundance_matrix <- relab(aggregated_list[['count_df']])
#' pa <- prevalence_abundance(abundance_matrix)
#' pa$data
#' pa$plot

prevalence_abundance <- function(abundance_matrix){
  
  prevalence <- rowSums(abundance_matrix > 0)
  ave_abund <- rowMeans(abundance_matrix)
  prev_abund <- data.frame(Taxon = get_shortnames(rownames(abundance_matrix)),
                           Prevalence = prevalence,
                           Abundance = ave_abund)
  prev_plt <- ggplot(prev_abund, aes(x=Abundance, y=Prevalence)) +
              geom_point() +
              scale_x_log10() +
              theme_classic() +
              xlab("Mean(relative abundance)")
  
  return(list(data = prev_abund, plot = prev_plt))
}


