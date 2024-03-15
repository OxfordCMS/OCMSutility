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


#' plotPCA
#'
#' @description
#' plot PCA
#' `r lifecycle::badge("deprecated")`
#'
#' `plotPCA` was renamed to `plot_pca` for more consistent
#' function naming nomenclature in `OCMSutility`. Usage is still exactly the same
#' @keywords internal
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
#' plotPCA(pca_result, state_data, colourby='division', shapeby='region')
#' # ->
#' plot_pca(pca_result, state_data, colourby='division', shapeby='region')
#'
plotPCA <- function(
  pc, metadata, colourby="none", shapeby="none", group="none",
  continuous=FALSE,  pcs=c("PC1", "PC2")) {
  lifecycle::deprecate_warn("0.3.0", "plotPCA()", "plot_pca()")
  plot_pca(pc, metadata, colourby, shapeby, group, continuous, pcs)
}
