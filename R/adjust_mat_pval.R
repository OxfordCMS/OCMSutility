#' adjust_mat_pval
#'
#' Adjust matrix of p-values for multiple correction. Helpful for correlation matrices.
#' @param mat numeric matrix
#' @param method multiple testing correction method. default 'BH'
#' @param out_type 'matrix' or 'dataframe'. matrix returns symmetrical matrix.
#'                  dataframe returns long dataframe
#' @returns
#' Adjusted p-values as either matrix (default) or long data.frame
#' @export
#'
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
#' corr_result <- psych::corr.test(feat_mat)
#'
#' adjust_mat_pval(corr_result$p.value)
#' adjust_mat_pval(corr_result$p.value, out_type='dataframe')
adjust_mat_pval <- function(mat, method='BH', out_type='matrix') {

  tri <- mat[upper.tri(mat)]
  padj <- p.adjust(tri, method=method)

  # put adjusted p.values in symmetrical matrix
  adj_mat <- matrix(1, nrow(mat), ncol(mat))
  rownames(adj_mat) <- rownames(mat)
  colnames(adj_mat) <- colnames(mat)
  adj_mat[upper.tri(adj_mat)] <- padj

  adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  if(out_type == 'matrix') {
    out <- adj_mat
  }
  if(out_type == 'dataframe') {
    # convert upper triangle matrix to long df
    xy <- t(combn(colnames(adj_mat), 2))
    out <- data.frame(xy, padjust=adj_mat[xy])
  }
  return(out)
}