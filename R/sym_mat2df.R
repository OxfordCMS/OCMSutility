#' sym_mat2df
#'
#' Convert symmetrical matrix to long dataframe.
#' Helpful for correlation or distance matrices
#'
#' @param mat symmetrical matrix
#' @returns long dataframe with columns x, y, and value
#' @export
#' @examples
#' # load example data
#' data(dss_example)
#'
#' # subset features, features in columns
#' feat_mat <- dss_example$merged_abundance_id[1:6,2:26]
#' rownames(feat_mat) <- dss_example$merged_abundance_id[1:6,1]
#' feat_mat <- t(feat_mat)
#'
#' # correlation matrix
#' corr_result <- cor(feat_mat)
#'
#' sym_mat2df(corr_result)
#'
sym_mat2df <- function(mat) {
  # convert upper triangle matrix to long df
  xy <- t(combn(colnames(mat), 2))
  out <- data.frame(xy, value=mat[xy])

  return(out)
}