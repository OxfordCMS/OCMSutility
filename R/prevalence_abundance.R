#' prevalence_abundance
#' return a dataframe that gives the prevalence and abundance (mean across samples)
#' of each taxon (from relative abundance matrix)
#'
#' @param abundance_matrix dataframe; samples are columns and features are rows (relative abundance)
#' @import reshape2
#' @export
#' @return dataframe
#' @examples
#' prevalence_abundance(abundance_matrix)

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


