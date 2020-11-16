#' ocms_clr
#'
#' use aldex2 to perform clr transform.  \code{aldex.clr} performs Monte Carlo sampling from a Dirichlet distribution for each sample. Then clr transformation each value generated. This function uses the median of these Monte-Carlo Dirichlet instances as the clr-transformed value, and returns as a dataframe be default.
#'
#'  @param count_dataframe dataframe with samples in columns, ASVs in rows
#'  @param condition passed into \code{conds} argument of \code{aldex.clr}.
#'      When set to \code{NULL} (default), each sample is assigned a different condition.
#'      Alternatively, can supply a vector containing a descriptor for the
#'      samples, allowing them to be grouped and compared. This setting is
#'      arbitrary, unless \code{return_as_dataframe} is set to \code{FALSE} and
#'      performing diversity analysis with \code{ALDEx2}.
#' @param return_as_dataframe default \code{TRUE}. Set to \code{FALSE} to
#'      return as \code{aldex.clr} object
#'
#' @return dataframe of clr-tranfsormed counts or aldex.clr object with
#'      clr-transformed Monte-Carlo Dirichlet instances of counts
#' @import ALDEx2
#' @export
#' @examples
#' library(ALDEx2)
#' data(selex)
#' subset for efficiency
#' selex <- selex[1201:1600,]
#' df <- ocms_clr(count_dataframe = selex, condition = NULL, return_as_dataframe = TRUE)

ocms_clr <- function(count_dataframe, condition = NULL, return_as_dataframe = TRUE) {

  # check input classes
  if(!is.logical(return_as_dataframe)) {
    stop("return_as_dataframe must be logical value")
  }
  if(!class(count_dataframe) %in% c('data.frame', 'matrix')) {
    stop("count_dataframe must be of class dataframe or matrix")
  }

  # set condition to each sample
  if(is.null(condition)) {
    condition <- as.character(1:ncol(count_dataframe))
  }

  asv_clr <- ALDEx2::aldex.clr(count_dataframe, conds = condition, useMC = TRUE)

  if(return_as_dataframe) {
    clr_instance <- lapply(ALDEx2::getMonteCarloInstances(asv_clr),
                           function(m){t(apply(m,1,median))})
    ## samples in columns
    clr_df <- data.frame(matrix(unlist(clr_instance),
                                ncol = length(clr_instance),
                                byrow = FALSE,
                                dimnames = list(colnames(clr_instance[[1]]),
                                                names(clr_instance))),
                         stringsAsFactors=FALSE)

    out <- clr_df
  }
  else {
    out <- asv_clr
  }

  return(out)
}
