#' bcdissimilarity
#'
#' Calculate Bray-Curtis Dissimilarity for samples within and between groups.
#'
#' @param relab_mat data.frame or matrix. relative abundance table with samples in
#'         columns and features in rows; feature IDs as rownames
#' @param met dataframe. metadata table
#' @param sample_id string. Sample identifier (column name in `met`)
#' @param var string. Metadata variable (column name in `met`) that specifies groups
#'
#' @examples
#'
#' data(dss_example)
#'
#' count_table <- dss_example$merged_abundance_id %>% column_to_rownames('featureID')
#' relab_table <- relab(count_table)
#'
#' met_table <- dss_example$metadata
#'
#' # bray curtis for one metadata variable
#' bc_result <- bcdissimilarity(relab_table, met_table, 'sampleID','Phenotype')
#'
#' pdata <- bc_result$bc_df
#' p_phen <- ggplot(pdata, aes(x=value, y=dist)) +
#'    geom_violin() +
#'    theme_bw(14) +
#'    ylab('Bray-Curtis Dissimilarity') +
#'    xlab('Phenotype')
#'
#' for multiple metadata variables
#' bc_data <- c()
#' for(var in c('Phenotype','Genotype')) {
#'     bc_result <- bcdissimilarity(relab_table, met_table, 'sampleID', var)
#'     bc_data <- rbind(bc_data, bc_result$bc_df)
#' }
#'
#' p_bc <- ggplot(bc_data, aes(x=comparison, y=dist)) +
#'    geom_violin() +
#'    theme_bw(14) +
#'    facet_wrap(~met_var) +
#'    ylab('Bray-Curtis Dissimilarity') +
#'    theme(axis.title.x=element_blank())
#'
#' @returns list.
#'          `bc_df` is long dataframe of samples being compared and corresponding distances.
#'          `bc_dist` matrix of bray curtis distances

bcdissimilarity <- function(relab_mat, met, sample_id, var) {
  curr_dist <- vegan::vegdist(t(relab_mat), method = "bray")
  curr_dist <- as.matrix(curr_dist)

  # convert upper triangle matrix to long df
  xy <- t(combn(colnames(curr_dist), 2))
  bc_df <- data.frame(xy, dist=curr_dist[xy])

  bc_df <- bc_df %>%
    left_join(met %>% select(!!sym(sample_id), !!sym(var)), c('X1'=sample_id))

  col_name <- colnames(bc_df)
  col_name[which(col_name == sample_id)] <- paste(sample_id, 'X1', sep='.')
  col_name[which(col_name == var)] <- paste(var, 'X1', sep='.')

  colnames(bc_df) <- col_name

  bc_df <- bc_df %>%
    left_join(met %>% select(!!sym(sample_id), !!sym(var)), c('X2'=sample_id))

  col_name <- colnames(bc_df)
  col_name[which(col_name == sample_id)] <- paste(sample_id, 'X2', sep='.')
  col_name[which(col_name == var)] <- paste(var, 'X2', sep='.')

  colnames(bc_df) <- col_name

  bc_df <- bc_df %>%
    mutate(comparison = ifelse(!!sym(sprintf("%s.X1", var)) ==
                                 !!sym(sprintf("%s.X2", var)),
                               'within', 'between'),
           met_var = var,
           value = ifelse(comparison == 'within', !!sym(sprintf('%s.X1', var)),
                          'between'))

  bc_df <- bc_df[,!grepl(var, colnames(bc_df))]

  return(list(bc_df=bc_df, bc_dist = dist(curr_dist)))
}