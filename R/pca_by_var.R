#' pca_by_var
#'
#' @description
#' Overlay numeric metadata variables onto a PCA score plot.
#' Useful during exploratory analysis
#'
#' @param ddata dataframe. samples in rows. sample ids in rownames.
#'              rows should match mdata
#' @param mdata dataframe. samples in rows. sample ids in rownames.
#'              rows should match ddata. all variables should be numeric.
#' @param PC numeric vector of length 2. specifies the PCs to plot.
#'           default is `c(1,2)`
#' @param biplot logical. show PCA as score plot or biplot
#' @param score_colour string. colour data points on score/biplot.
#'                    should be column name in mdata.
#'                    When set to `FALSE`, data points are black.
#'                    Default `FALSE`
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
#' @import ggfortify
#'
#' @export
#'
#' @examples
#'
#' data(iris)
#' dd <- iris
#'
#' dd$species_num <- NA
#' dd$species_num[dd$Species == 'setosa'] <- 1
#' dd$species_num[dd$Species == 'versicolor'] <- 2
#' dd$species_num[dd$Species == 'virginica'] <- 3
#'
#' set.seed(1)
#' dd$var1 <- rep(rnorm(15, 25, 3), each=10)
#'
#' set.seed(1)
#' dd$var2 <- rep(rnorm(10, 3, 0.5), 15)
#' rownames(dd) <- paste0('sample',1:nrow(dd))
#' ddata <- dd[,1:4]
#' mdata <- dd[,6:8]
#'
#' p_list <- pca_by_var(ddata, mdata)
#'
#' # biplot
#' p_list$main_pca
#'
#' # panel by metadata variables
#' cowplot::plot_grid(plotlist=list(p_list$species_num, p_list$var1, p_list$var2))

pca_by_var <- function(ddata, mdata, PC=c(1,2), biplot=TRUE,
                       score_colour = FALSE, score_label=FALSE) {

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
  if(!is.logical(biplot)) {
    stop("biplot should be TRUE or FALSE")
  }
  if(!is.logical(score_label)) {
    stop("score_label should be TRUE or FALSE")
  }
  if(score_colour != FALSE && !score_colour %in% colnames(mdata)) {
    stop("score_colour FALSE or a column in mdata")
  }

  # pca
  pcx <- prcomp(ddata, center=TRUE, scale=TRUE)

  # score plot/biplot
  if(biplot) {
    if(score_colour == FALSE) {
      p_biplot <- autoplot(pcx, x=PC[1], y=PC[2], data=mdata,
                           label=score_label,
                           loadings=TRUE,
                           loadings.label=TRUE, loadings.colour=NA,
                           loadings.label.colour='purple')
    } else {
      p_biplot <- autoplot(pcx, x=PC[1], y=PC[2], data=mdata,
                           label=score_label,
                           colour=score_colour,
                           loadings=TRUE,
                           loadings.label=TRUE, loadings.colour=NA,
                           loadings.label.colour='black')
    }

  } else {
    if(score_colour == FALSE) {
      p_biplot <- autoplot(pcx, data=mdata,
                           label=score_label,
                           loadings=FALSE)
    } else {
      p_biplot <- autoplot(pcx, data=mdata,
                           label=score_label, colour = score_colour,
                           loadings=FALSE)
    }
  }

  p_biplot <- p_biplot + theme_bw(14)

  if(score_label) {
    p_biplot$layers[[1]] <- NULL
  }

  # putting score plot/biplot into output list
  out_list <- list('main_pca'=p_biplot)

  # getting score data
  score_data <- as.data.frame(pcx$x) %>%
    rownames_to_column('new_id')

  # adding metadata
  pdata <- score_data %>%
    left_join(mdata %>% rownames_to_column('new_id'), 'new_id')

  # iterating through metadata variables
  iter <- colnames(mdata)
  for(i in iter) {
    x = sprintf("PC%d", PC[1])
    y = sprintf("PC%d", PC[2])
    p_curr <- ggplot(pdata, aes(x=!!sym(x), y=!!sym(y), fill=!!sym(i))) +
      geom_point(shape=21, size=4) +
      theme_bw(14) +
      theme(legend.title=element_blank(),
            axis.title = element_blank()) +
      scale_fill_gradientn(colours=brewer.pal(11,'Spectral')) +
      labs(subtitle=i)

    out_list[[i]] <- p_curr
  }

  return(out_list)
}