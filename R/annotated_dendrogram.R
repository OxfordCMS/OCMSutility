#' annotated_dendrogram
#'
#' @description Creates an annotated dendrogram from a distance matrix.
#'              Dendrogram annotations are based on metadata supplied
#' @param dist class \code{dist}. distance matrix
#' @param met dataframe. samples in rows, metadata variables in columns.
#' @param id string. sample identifier. must be a column in met
#' @param method string. agglomeration method used by \code{hclust}.
#' @param coord numeric vector. optional. allows you to specify coordinates
#'         of annotations. Vector must be same length as number of variables supplied.
#'         Leave as \code{NULL} for default coordinates.
#' @param pal unnamed vector. Colours to be used for each variable. see details section.
#'            Set to \code{NULL} to use default (RColorBrewer) colours.
#' @details Colour values needs to match the number variable values \code{pal}
#' @returns ggplot
#' @export
#' @import ggplot2
#' @import ggnewscale
#' @import ggdendro
#' @import dplyr
#' @import RColorBrewer
#' @examples
#' set.seed(1)
#' # get relative abundance data
#' data(dss_example)
#' ddata <- dss_example$merged_abundance_id[,2:26]
#' rownames(ddata) <- dss_example$merged_abundance_id[,1]
#' ddata <- t(OCMSutility::relab(ddata))
#' # distance matrix
#' mydist <- vegan::vegdist(ddata, method='bray')
#'
#' # metdata variable
#' mdata <- dss_example$metadata
#' mdata <- mdata[,c('sampleID','Genotype','Phenotype')]
#' annotated_dendrogram(mydist, mdata, 'sampleID')
#' # custom colours
#' col_geno <- RColorBrewer::brewer.pal(9, "Paired")[1:2]
#' col_phen <- RColorBrewer::brewer.pal(9, "Paired")[3:4]
#' annotated_dendrogram(mydist, mdata, 'sampleID', pal=c(col_geno, col_phen))

annotated_dendrogram <- function(dist, met, id, method='complete',
                                 coord=NULL, pal=NULL){

  if(!id %in% colnames(met)) {
    stop("identifier supplied not found as a column in met")
  }

  if(class(dist) != 'dist') {
    stop("need distance matrix of class dist")
  }
  if(!is.null(coord) && !is.numeric(coord)) {
    stop("coord should be NULL or a vector of numbers equal to the number of metadata variables")
  }
  coord_df <- data.frame(var=met_var, coord=coord)

  pdata <- ggdendro::label(dend_data) %>%
    left_join(met, c(label = id))  %>%
    gather('var','value', -x, -y, -label) %>%
    left_join(coord_df, 'var') %>%
    mutate(fill = paste(var, value, sep=":"))

  if (!is.null(pal) && length(pal) != (length(unique(pdata$fill)))) {
    stop("need a colour palette for each variable")
  }
  # annotation locations
  if(is.null(coord)) {
    coord <- seq(from=-0.5, by=-0.8, length.out=ncol(met)-1)
  }

  # colour palettes
  if (is.null(pal)) {
    pal_pool <- c(brewer.pal(12, "Paired"), brewer.pal(12,
                                                       "Set3"))
    pal <- pal_pool[1:length(unique(pdata$fill))]
  } else {
    names(pal) <- NULL
  }

  text_df <- data.frame(y=coord,
                        x=round(max(dend_data$segments$x) * 1.02),
                        label = met_var)

  p <- ggplot(ggdendro::segment(dend_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = text_df,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5, angle = 90) +
    geom_tile(data = pdata,
              aes(x = x, y = coord, fill = fill),
              height = abs(coord[1] - coord[2])) +
    scale_fill_manual(values = pal, name = NULL) +
    expand_limits(x = c(0, round(max(dend_data$segments$x) * 1.3))) +
    coord_flip() +
    scale_y_reverse() +
    ylab("Distance") +
    theme_classic(14) +
    theme(axis.title.y = element_blank())

  return(p)

}#' annotated_dendrogram
#'
#' @description Creates an annotated dendrogram from a distance matrix.
#'              Dendrogram annotations are based on metadata supplied
#' @param dist class \code{dist}. distance matrix
#' @param met dataframe. samples in rows, metadata variables in columns.
#' @param id string. sample identifier. must be a column in met
#' @param method string. agglomeration method used by \code{hclust}.
#' @param coord numeric vector. optional. allows you to specify coordinates
#'         of annotations. Vector must be same length as number of variables supplied.
#'         Leave as \code{NULL} for default coordinates.
#' @param pal unnamed vector. Colours to be used for each variable. see details section.
#'            Set to \code{NULL} to use default (RColorBrewer) colours.
#' @details Colour values needs to match the number variable values \code{pal}
#' @returns ggplot
#' @export
#' @import ggplot2
#' @import ggnewscale
#' @import ggdendro
#' @import dplyr
#' @import RColorBrewer
#' @examples
#' set.seed(1)
#' # get relative abundance data
#' data(dss_example)
#' ddata <- dss_example$merged_abundance_id[,2:26]
#' rownames(ddata) <- dss_example$merged_abundance_id[,1]
#' ddata <- t(OCMSutility::relab(ddata))
#' # distance matrix
#' mydist <- vegan::vegdist(ddata, method='bray')
#'
#' # metdata variable
#' mdata <- dss_example$metadata
#' mdata <- mdata[,c('sampleID','Genotype','Phenotype')]
#' annotated_dendrogram(mydist, mdata, 'sampleID')
#' # custom colours
#' col_geno <- RColorBrewer::brewer.pal(9, "Paired")[1:2]
#' col_phen <- RColorBrewer::brewer.pal(9, "Paired")[3:4]
#' annotated_dendrogram(mydist, mdata, 'sampleID', pal=c(col_geno, col_phen))

annotated_dendrogram <- function(dist, met, id, method='complete',
                                 coord=NULL, pal=NULL){

  if(!id %in% colnames(met)) {
    stop("identifier supplied not found as a column in met")
  }

  if(class(dist) != 'dist') {
    stop("need distance matrix of class dist")
  }
  if(!is.null(coord) && !is.numeric(coord)) {
    stop("coord should be NULL or a vector of numbers equal to the number of metadata variables")
  }
  # annotation locations
  if(is.null(coord)) {
    coord <- seq(from=-0.5, by=-0.8, length.out=ncol(met)-1)
  }

  met_var <- colnames(met)
  met_var <- met_var[met_var != id]
  dend <- as.dendrogram(hclust(dist, method = method))
  dend_data <- ggdendro::dendro_data(dend, type = "rectangle")
  coord_df <- data.frame(var=met_var, coord=coord)

  pdata <- ggdendro::label(dend_data) %>%
    left_join(met, c(label = id))  %>%
    gather('var','value', -x, -y, -label) %>%
    left_join(coord_df, 'var') %>%
    mutate(fill = paste(var, value, sep=":"))

  if (!is.null(pal) && length(pal) != (length(unique(pdata$fill)))) {
    stop("need a colour palette for each variable")
  }


  # colour palettes
  if (is.null(pal)) {
    pal_pool <- c(brewer.pal(12, "Paired"), brewer.pal(12,
                                                       "Set3"))
    pal <- pal_pool[1:length(unique(pdata$fill))]
  } else {
    names(pal) <- NULL
  }

  text_df <- data.frame(y=coord,
                        x=round(max(dend_data$segments$x) * 1.02),
                        label = met_var)

  p <- ggplot(ggdendro::segment(dend_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = text_df,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5, angle = 90) +
    geom_tile(data = pdata,
              aes(x = x, y = coord, fill = fill),
              height = abs(coord[1] - coord[2])) +
    scale_fill_manual(values = pal, name = NULL) +
    expand_limits(x = c(0, round(max(dend_data$segments$x) * 1.3))) +
    coord_flip() +
    scale_y_reverse() +
    ylab("Distance") +
    theme_classic(14) +
    theme(axis.title.y = element_blank())

  return(p)

}