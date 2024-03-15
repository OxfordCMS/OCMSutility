#' plotPCoA
#' `r lifecycle::badge("deprecated")`
#' This is a simple PCoA function that colours all points by one
#' metadata variable. It can be helpful to visualise metadata variables
#' independently when assessing potential confounding metadtaa factors
#'
#' @param relab dataframe. relative abundance data with features in rows
#'              and samples in columns. feature IDs in rowname.
#' @param met dataframe. metadata with samples in rows. sample IDs in rowname
#' @param colour string. defulat NULL. metadata variable to colour points by
#' @param shape string. default NULL. metadata variable to set shape of points by
#' @param CI numeric. Default 0.95. Confidence interval used to draw ellipse
#'           around colour variable. set to NULL to omit drawing ellipse
#' @import vegan
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' data(dss_example)
#' met_df <- dss_example$metadata
#'
#' count_df <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#' count_df <- count_df[,met_df$sampleID]
#' relab <- relab(count_df)
#'
#' iter_var <- c('Genotype','Phenotype')
#' for(i in iter_var) {
#'   plotPCoA(relab, met_df, colour = i)
#' }


plotPCoA <- function(relab, met, colour=NULL, shape=NULL, CI=0.95) {

  lifecycle::deprecate_warn("0.3.0", "plotPCoA()", "plot_pcoa()")
  plot_pcoa(relab, met, colour, shape, CI)
}