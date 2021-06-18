#' plotPCoA
#'
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

  require(vegan)
  require(dplyr)
  require(ggplot2)

  # check inputs----------------------------------------------------------------
  # check data is in relative abundance format
  if(all(relab %% 1 == 0)) {
    stop("All values are integers. Values should be in relative abundance for Bray Curtis distance used in PCoA")
  }

  # check same sampleIDs found in both relab and met
  if(!identical(sort(rownames(met)), sort(colnames(relab)))) {
    stop("Sample IDs in met and relab do not match")
  }

  # check colour and shape variables are found in met
  if(!is.null(colour)){
    if(!colour %in% colnames(met)) {
      stop(sprintf("colour variable '%s' not a column in met"))
    }
  }
  if(!is.null(shape)) {
    if(!shape %in% colnames(met)) {
      stop(sprintf("shape variable '%s' not a column in met"))
    }
  }

  # check CI is between 0 and 1
  if(!is.null(CI)) {
    if(CI > 1 || CI < 0) {
      stop("CI must be between 0 and 1 or set to NULL")
    }
  }

  # set rownames as sampleID in met
  if(!'sampleID' %in% colnames(met)) {
    met$sampleID <- rownames(met)
  }

  # make pcoa-------------------------------------------------------------------
  # calculate Bray Curtis distance
  bc_dist <- vegan::vegdist(t(relab), method = "bray")

  # calculate pcoa
  pcoa <- cmdscale(bc_dist, k = 4)
  eig <- cmdscale(bc_dist, k=2, eig = TRUE)
  var_explained <- round(eig$eig*100/sum(eig$eig),1)

  pcoa <- as.data.frame(pcoa)
  colnames(pcoa) <- paste0('PC', 1:ncol(pcoa))
  pcoa$sampleID <- rownames(pcoa)

  pdata <- pcoa %>%
    left_join(met, 'sampleID')

  p <- ggplot(pdata, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(colour=colour, shape=shape), size = 2, alpha = 0.6)

  if(!is.null(CI)) {
    p <- p +
      stat_ellipse(aes_string(colour=colour), level=CI)
  }

  p <- p +
    xlab(sprintf('PC1 (%s%%)', var_explained[1])) +
    ylab(sprintf('PC2 (%s%%)', var_explained[2])) +
    theme_classic(10)
  p
  return(p)
}