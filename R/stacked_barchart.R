#' stacked_barchart
#' stacked barchart of the top_n features in a relative abundance table.
#'
#' @param abundance_matrix dataframe; samples are columns and features are rows (relative abundance)
#' @param top_n integer; number of taxa to plot in the stacked barchart. Default is 10.
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
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
#' st <- stacked_barchart(abundance_matrix)
#' st$data
#' st$plot

stacked_barchart <- function(abundance_matrix, top_n = 10){
  
  # The top_n is defined by the average relative abundance
  # across samples in the matrix
  ave_abund <- rowMeans(log10(abundance_matrix + 1e-6))
  relabund <- abundance_matrix[order(ave_abund, decreasing=TRUE),]
  
  # take the top n
  toppers <- rownames(relabund[1:top_n,])
  agg_type <- ifelse(rownames(relabund) %in% toppers, rownames(relabund), "Other")
  relabund <- aggregate(relabund, by = list(agg_type), FUN = "sum")
  relabund$Taxon <- relabund$Group.1 

  # make long
  relabund_m <- reshape2::melt(relabund)
  colnames(relabund_m) <- c("Taxon_", "Taxon", "Sample", "Relab") 
  
  # shortnames for taxa: visual purposes
  relabund_m$Taxon <- get_shortnames(relabund_m$Taxon)
  st_bar <- ggplot(relabund_m, aes(x = Sample, y = Relab, fill = Taxon)) +
    geom_bar(stat="identity") +
    theme_classic() +
    scale_fill_manual(values = getPalette(n = top_n+1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  ret <- list(data = relabund_m, plot = st_bar)
  return(ret)
}

