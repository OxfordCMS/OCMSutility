#' compare_cor_ci
#'
#' Performs pairwise correlations of features with adjusted p-values.
#' Correlations and confidence intervals calculated for each sample group.
#' @param mat matrix with features in columns.
#'            correlations calculated between feature pairs for each group
#' @param group group assignment of samples corresponding to rows in mat.
#' @param method correlation method pearson, spearman, kendall
#' @param adjust multiple testing correction method.
#'               BH, BY, fdr, hom, hochberg, hommel, bonferroni, none
#' @param alternative two.sided, greater, less
#' @param conf_level confidence level for the returned confidence interval.
#' @returns Returns a dataframe.
#' x = x variable
#' y = y variable
#' group = sample group
#' n = number of samples in group
#' r = correlation coefficient
#' p = p-value
#' p.adj = adjusted p-value
#' lower_ci = lower confidence interval
#' upper_ci = upper confidence interval
#'
#' @export
#' @importFrom psych corr.test
#' @importFrom DescTools CorCI
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
#' # metadata in same order
#' met_df <- dss_example$metadata
#' met_df <- met_df[match(rownames(feat_mat), met_df$sampleID),]
#'
#' compare_cor_ci(feat_mat, met_df$Phenotype)

compare_cor_ci <- function(mat, group, method='pearson', adjust='BH',
                           alternative='two.sided', conf_level=0.95) {
  # needs sym_mat2df function from OMCS_Sandbox/R_utility/sym_mat2df.R
  # confidence interval calculation (Fieller) for spearman and kendall correlation
  # obtained from https://rpubs.com/seriousstats/616206

  df <- cbind(as.data.frame(mat), group=as.character(group))
  df_ls <- split(df, as.character(group))

  out <- c()
  for(i in 1:length(df_ls)) {
    # correlation
    curr_cor <- corr.test(df_ls[[i]][,1:ncol(mat)],
                                 method=method, adjust=adjust)

    # correlation coefficient
    curr_r <- sym_mat2df(curr_cor$r)
    colnames(curr_r) <- c('x','y','r')

    # pvalue
    curr_p <- curr_cor$p
    curr_p[upper.tri(curr_p)] <- t(curr_p[lower.tri(curr_p)])
    curr_p <- sym_mat2df(curr_p)
    colnames(curr_p) <- c('x','y','p')

    # adjusted pvalue
    curr_padj <- curr_cor$p
    curr_padj[lower.tri(curr_padj)] <- NA
    curr_padj <- sym_mat2df(curr_padj)
    colnames(curr_padj) <- c('x','y','p.adj')

    # confidence intervals
    if(method == 'pearson') {
      curr_ci <- t(sapply(curr_r$r, CorCI, n=nrow(df_ls[[i]]),
                          alternative=alternative, conf.level=conf_level))

      curr_ci <- as.data.frame(curr_ci)
      colnames(curr_ci) <- c('r','lower_ci','upper_ci')
      curr_ci$x <- curr_r$x
      curr_ci$y <- curr_r$y

    } else if(method == 'spearman') {
      zrs.se <- (1.06/(nrow(df_ls[[i]]) - 3))^0.5
      moe <- qnorm(1 - (1 - conf_level)/2) * zrs.se
      zu <- atanh(curr_r$r) + moe
      zl <- atanh(curr_r$r) - moe
      vec_ci <- tanh(c(zl, zu))
      curr_ci <- data.frame(x=curr_r$x, y=curr_r$y,
                            lower_ci = vec_ci[1:length(vec_ci), 2],
                            upper_ci = vec_ci[2:length(vec_ci), 2])
    } else if(method == 'kendall') {
      tau.se <- (0.437/(nrow(df_ls[[i]]) - 4))^0.5
      moe <- qnorm(1 - (1 - conf_level)/2) * tau.se
      zu <- atanh(curr_r$r) + moe
      zl <- atanh(curr_r$r) - moe
      vec_ci <- tanh(c(zl, zu))
      curr_ci <- data.frame(x=curr_r$x, y=curr_r$y,
                            lower_ci = vec_ci[seq_along(vec_ci) %% 2 == 1],
                            upper_ci = vec_ci[seq_along(vec_ci) %% 2 == 0])
    }

    entry <- data.frame(x = curr_r$x,
                        y = curr_r$y,
                        group = unique(df_ls[[i]]$group),
                        n = nrow(df_ls[[i]]),
                        r = curr_r$r,
                        p = curr_p$p,
                        p.adj = curr_padj$p.adj,
                        lower_ci = curr_ci$lower_ci,
                        upper_ci = curr_ci$upper_ci)

    out <- rbind(out, entry)
  }
  return(out)
}
