#' top_taxa
#'
#' Arrange top taxa in descending order. Assessed based on cumulative sum,
#' so only the top (80% for example) taxa are returned.
#'
#' @param relab dataframe or matrix. relative abundance matrix with taxa in
#' rows (as rownames) and samples in columns
#' @param cutoff numeric. between 0 and 1. Cut-off abundance for what is
#' considered as "top taxa." Default 0.8. Taxa with abundances less than the
#' cut-off are lumped together as "Other"
#' @examples
#' # put asv relative abundance features in same order as taxonomy table
#' # get relative abundance data
#' asv_mat <- dss_example$merged_abundance_id[,2:ncol(dss_example$merged_abundance_id)]
#' rownames(asv_mat) <- dss_example$merged_abundance_id$featureID
#' relab_data <- relab(asv_mat)
#' # match row order of relative abundance and taxonomy data
#' relab_data <- relab_data[match(dss_example$merged_taxonomy$featureID,
#'                                rownames(relab_data)),]
#' # prepend X to column names
#' colnames(relab_data) <- sprintf("X%s", colnames(relab_data))
#' # get genus-level relative abundance
#' genus_df <- aggregate_count(relab_data, dss_example$merged_taxonomy, 'Genus')$count_df
#'
#' # get top taxa
#' top80 <- top_taxa(genus_df, cutoff=0.8)
#' print(top80 %>% filter(sample_id == unique(top80$sample_id)[1]))
#'
#' @return dataframe of taxa, sample_id, relative abundance, and
#' cumulative sum of relative abundance.
#' @export

top_taxa <- function(relab_df, cutoff) {

  # parameter check
  if(!class(relab_df) %in% c('data.frame','matrix')) {
    stop("relab must be of class data.frame or matrix")
  }

  if(class(cutoff) != 'numeric') {
    stop("cutoff must be a number between 0 and 1")
  }

  if(cutoff < 0 | cutoff > 1) {
    stop("cutoff must be a number between 0 and 1")
  }

  top <- relab_df %>%
    rownames_to_column('taxa') %>%
    gather('sample_id','relab',-taxa) %>%
    group_by(sample_id) %>%
    arrange(desc(relab), .by_group = TRUE) %>%
    mutate(csum = round(cumsum(relab), 2),
           keep1 = ifelse(min(csum) > cutoff & relab > cutoff, 1, 0),
           keep2 = ifelse(keep1 == 0 & csum <= cutoff, 1, 0))

  out <- top %>%
    filter(keep1 == 1 | keep2 == 1) %>%
    arrange(desc(relab)) %>%
    select(taxa, sample_id, relab, csum)

  return(out)
}