#' truePosRate
#'
#' @description
#' Calculate rate of true positives in positive control standards. Used in OCMS_zymobioimcs report.
#' `r lifecycle::badge("deprecated")`
#'
#' `truePosRate` was renamed to `true_pos_rate` for more consistent
#' function naming nomenclature in `OCMSutility`. Usage is still exactly the same
#' @keywords internal
#' @param relab relative abundance data frame with samples in columns, features in rows
#' @param annotations reference annotations for true positives
#' @param level taxonomy level at which analysis is done. can be one of
#'    \code{'species','genus','family'}. default \code{'species'}
#' @param cutoff threshold for read relative abundance in the sample. default is
#'    \code{0.01} (0.01 %)
#' @return dataframe ranking features as true positive or false positive in each sample.
#' @export
#'
#' @examples
#'  # this would be better exemplified with actual std data rather than the example smaples
#'  data("dss_example")
#' data(zymobiomics)
#'
#' # set count feature ids as rownames
#' count_df <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#'
#' # clean up some sample names
#' colnames(count_df) <- paste0('id', colnames(count_df))
#' tax_df <- dss_example$merged_taxonomy
#'
#' # aggregate counts
#' agg_gen <- aggregateCount(count_df[tax_df$featureID,], tax_df, "Genus")
#' genus_relab <- relab(agg_gen$count_df)
#'
#' true_pos_result <- truePosRate(relab=genus_relab,
#'                                annotations=zymobiomics$anno_ncbi_16s,
#'                                level='genus', cutoff=0.01)
#' # ->
#' true_pos_result <- true_pos_rate(relab=genus_relab,
#'                                annotations=zymobiomics$anno_ncbi_16s,
#'                                level='genus', cutoff=0.01)
#' # plot true pos rate
#' p <- ggplot(true_pos_result,
#'             aes(x=rank, y=true.pos.rate, colour=label, group=sample)) +
#'  geom_point() +
#'   theme_bw() +
#'   ylab("TP / (TP + FP)") +
#'  scale_colour_manual(values=c("grey", "purple")) +
#'  facet_wrap(~sample, scale="free")
#'
#' p

truePosRate <- function(relab, annotations, level="species", cutoff=0.01){

  true_pos_rate(relab, annotations, level, cutoff)
}