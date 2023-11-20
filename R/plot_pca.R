#' getVE
#'
#' Calculate variance explained in pca. used in plotPrincipalComponents
#'
#' @param pc prcomp object
#' @param component principal component
#'
#' @return matrix of variance explained for a principal component
#'
getVE <- function(pc, component="PC1"){

  pve <- summary(pc)$importance[,component][[2]]
  return (pve)
}


#' plot_pca
#'
#' plot PCA
#'
#' @param pc prcomp object
#' @param metadata metadata dataframe
#' @param colourby colour by. default "none"
#' @param shapeby shape by. default "none"
#' @param group group aes in ggplot. default "none"
#' @param continuous is colourby a continuous variable? default FALSE
#' @param pcs principal components to plot. default \code{c("PC1", "PC2")}
#'
#' @returns list of ggplot pca plot and corresponding dataframe used to make the plot
#' @export
#' @examples
#' pca_result <- prcomp(USArrests, scale = TRUE)
#' state_data <- data.frame(abb = state.abb, area = state.area,
#'                          center = state.center, region = state.region,
#'                          division = state.division)
#' rownames(state_data) <- state.name
#' plot_pca(pca_result, state_data, colourby='division', shapeby='region')

plot_pca <- function(
  pc, metadata, colourby="none", shapeby="none", group="none",
  continuous=FALSE,  pcs=c("PC1", "PC2")) {

  # check inputs----------------------------------------------------------------
  # pc must be prcomp object
  if(class(pc) != 'prcomp') {
    stop("pc must be of class prcomp")
  }
  # metadata must have same IDs as pc
  if(!identical(sort(rownames(pc$x)), sort(rownames(metadata)))) {
    stop("rownames between pc$x and metadata do not match")
  }
  # colourby must be "none" or in metadata
  if(!colourby %in% c('none', colnames(metadata))) {
    stop("colourby must be 'none' or a column name of metadata")
  }
  # group must be "none" or in metadata
  if(!group %in% c('none', colnames(metadata))) {
    stop("group must be 'none' or a column name of metadata")
  }
  # shapeby must be "none" or in metadata
  if(!shapeby %in% c('none', colnames(metadata))) {
    stop("shapeby must be 'none' or a column name of metadata")
  }
  # continuous must be logical
  if(class(continuous) != 'logical') {
    stop("continuous must be either TRUE or FALSE")
  }
  # pcs must be found in promp model (pc$x)
  if(any(!pcs %in% colnames(pc$x))) {
    stop("pcs specified not found in pc model")
  }

  # calculate covariates--------------------------------------------------------
  # covariate must be in same order as pc rownames

  # get variance explained for each component
  ve1 <- getVE(pc, component=pcs[1])
  ve2 <- getVE(pc, component=pcs[2])

  ve1 <- round(ve1, 2)*100
  ve2 <- round(ve2, 2)*100

  # get data frame of components
  pca <- data.frame(pc$x)

  # add conditions
  if (colourby == "none"){
    pca$condition <- "none"}else{
      pca$condition <- metadata[,colourby]}

  # add shape
  if (shapeby == "none"){
    pca$shape <- "none"}else{
      pca$shape <- metadata[,shapeby]}

  if (group == "none"){
    pca$group <- "none"}else{
      pca$group <- metadata[,group]}

  if (continuous==FALSE){
    pca$condition <- factor(pca$condition, levels=unique(pca$condition))
  }

  # plot
  pc1 <- pcs[1]
  pc2 <- pcs[2]

  # labels
  xlabel <- paste(pc1, ve1, sep=" (")
  xlabel <- paste(xlabel, "%", sep="")
  xlabel <- paste(xlabel, ")", sep="")
  ylabel <- paste(pc2, ve2, sep=" (")
  ylabel <- paste(ylabel, "%", sep="")
  ylabel <- paste(ylabel, ")", sep="")

  n <- length(unique(pca$condition))
  colours <- rainbow(n, s=0.7, v=0.6)

  plot1 <- ggplot(pca, aes_string(x=pc1, y=pc2, group="group",
                                  colour="condition", shape="shape"))
  plot2 <- plot1 + geom_point(size=5, alpha=0.8)
  plot3 <- plot2 + theme_bw()
  plot4 <- plot3 + xlab(xlabel) + ylab(ylabel)
  if (continuous==TRUE){
    plot4 <- plot4 + scale_colour_gradient()}
  else{
    plot4 <- plot4 + scale_colour_manual(values=colours)
  }
  return(list('p' = plot4, 'pdata' = pca))
}
