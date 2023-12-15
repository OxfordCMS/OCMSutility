#' pcoa_by_var
#'
#' @description
#' Overlay numeric metadata variables onto a PCoA score plot.
#' Useful during exploratory analysis
#'
#' @param ddata dataframe. samples in rows. sample ids in rownames.
#'              rows should match mdata
#' @param mdata dataframe. samples in rows. sample ids in rownames.
#'              rows should match ddata. all variables can be numeric,
#'              characters, or factors
#' @param method string. distance method. Default \code('bray').
#'               Passed to \code(vegan::vegdist)
#' @param PC numeric vector of length 2. specifies the PCs to plot.
#'           default is `c(1,2)`
#' @param CI numeric. Default 0.95. Confidence interval used to draw ellipse
#'           around categorical variable. set to NULL to omit drawing ellipse
#' @param score_label string. label data ponts on score/biplot.
#'                    should be column name in mdata.
#'                    When set to `FALSE`, data ponts are not labelled.
#'                    Default `FALSE`
#'
#' @returns named list of plots. First plot in list is main_pca.
#'          Subsequent plots in list correspond to metadata variables
#'          (columns in metadata)
#'
#' @import RColorBrewer
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @import vegan
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' data(dss_example)
#' ddata <- dss_example$merged_abundance_id[,2:26]
#' rownames(ddata) <- dss_example$merged_abundance_id[,1]
#' ddata <- as.data.frame(t(relab(ddata)))
#'
#' mdata <- dss_example$metadata
#' mdata <- mdata[match(rownames(ddata), mdata$sampleID),]
#'
#' # creating some dummy metadata variable
#' mdata$var1 <- rnorm(25, 0.5, 3)
#' mdata$var2 <- rep(LETTERS[21:25], 5)
#' mdata$var3 <- as.factor(rep(letters[1:5], each=5))
#' mdata <- mdata[,c('Phenotype','var1','var2','var3')]
#' p_list <- pcoa_by_var(ddata, mdata, method='bray')
#'
#' # pcoa
#' p_list$main_pcoa
#' # pcoa with metadata variables overlayed. no ellipses draw when variables are numeric
#' p_list$Phenotype
#' p_list$var1
#' p_list$var2
#' p_list$var3
#'
#' # can use cowplot::plot_grid to put all plots into one
#' cowplot::plot_grid(plotlist=list(p_list$Phenotype, p_list$var1, p_list$var2, p_list$var3))

pcoa_by_var <- function(ddata, mdata, method='bray', PC=c(1,2),
                        CI = 0.95, score_label=FALSE) {

  # check inputs----------------------------------------------------------------
  # check data is in relative abundance format
  if(method == 'bray' && all(ddata %% 1 == 0)) {
    stop("All values are integers. Values should be in relative abundance for Bray Curtis distance used in PCoA")
  }
  if(!is.data.frame(ddata)) {
    stop("ddata should be data.frame with samples in rows and sample ids in rownames")
  }
  if(!is.data.frame(mdata)) {
    stop("mdata should be data.frame with samples in rows and sample ids in rownames")
  }
  if(!identical(rownames(ddata), rownames(mdata))) {
    stop("Rows of ddata and mdata (samples) are not in the same order")
  }
  if(length(PC) != 2 | !is.numeric(PC)) {
    stop("PC must be numeric vector of length 2")
  }
  if(!is.logical(score_label)) {
    stop("score_label should be TRUE or FALSE")
  }
  # check CI is between 0 and 1
  if(!is.null(CI)) {
    if(CI > 1 || CI < 0) {
      stop("CI must be between 0 and 1 or set to NULL")
    }
  }

  # make pcoa-------------------------------------------------------------------
  # calculate Bray Curtis distance
  ddist <- vegan::vegdist(ddata, method = method)

  # calculate pcoa
  pcoa <- cmdscale(ddist, k = 4)
  eig <- cmdscale(ddist, k=2, eig = TRUE)
  var_explained <- round(eig$eig*100/sum(eig$eig),1)

  pcoa <- as.data.frame(pcoa)
  colnames(pcoa) <- paste0('PC', 1:ncol(pcoa))
  pcoa$new_id <- rownames(pcoa)

  pdata <- pcoa %>%
    left_join(mdata %>% rownames_to_column('new_id'), 'new_id')

  x = sprintf("PC%d", PC[1])
  y = sprintf("PC%d", PC[2])
  p_main <- ggplot(pdata, aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(size = 2, alpha = 0.6) +
    xlab(sprintf('PC1 (%s%%)', var_explained[1])) +
    ylab(sprintf('PC2 (%s%%)', var_explained[2])) +
    theme_bw(14)

  pal <- c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Pastel2'))

  p_list <- list()
  p_list$main_pcoa <- p_main
  iter <- colnames(mdata)
  for(i in iter) {
    p <- ggplot(pdata, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(aes(fill=!!sym(i)), shape=21, size = 4, alpha = 0.6)


    if(!is.null(CI)) {
      if(is.numeric(mdata[,i])) {
        p <- p +
          scale_fill_gradientn(colours=brewer.pal(11,'Spectral'))
      } else  {
        # only draw ellipse when CI is not null and variable is not numeric
        p <- p +
          stat_ellipse(aes(colour=!!sym(i)), level=CI)

        # customize colour pallete if possible
        if(length(unique(mdata[,i])) <= length(pal)) {
          p <- p +
            scale_fill_manual(values=pal) +
            scale_colour_manual(values=pal)
        }
      }
    }

    p <- p +
      theme_bw(14) +
      theme(legend.title=element_blank(),
            axis.title = element_blank()) +
      labs(subtitle=i)

    p_list[[i]] <- p
  }

  return(p_list)
}