#' filter_metagenomics.R
#'
#' Filter metagenomics data based on cutoff threshold and feature prevalence across samples.
#'
#' @param counts_table dataframe. A count table with samples in columns and features (e.g., genes or taxa) in rows. Feature ID in rownames.
#' @param filter_method default \code{"abs_count"} must be one of \code{c("abs_count", "percent_sample", "percent_dataset")}.
#'     This determines how the cutoff for feature abundance is applied:
#'     \code{"abs_count"} filters based on an absolute count threshold.
#'     \code{"percent_sample"} filters based on the relative abundance in each sample.
#'     \code{"percent_dataset"} filters based on the relative abundance across the entire dataset.
#' @param feature_cutoff cutoff used to filter features. Features are kept when their abundance is greater than this cutoff.
#'     If \code{filter_method = "abs_count"}, this is the minimum read count threshold.
#'     If \code{filter_method = "percent_sample" or "percent_dataset"}, this is the minimum relative abundance threshold (percentage).
#' @param prev_cutoff prevalence cutoff. Features must be present in at least this many samples to be kept.
#'     \code{prev_cutoff} is the minimum number of samples a feature must appear in (based on \code{feature_cutoff}).
#'     Default value is 2, meaning a feature must appear in at least two samples to be kept.
#' @return
#'     A list containing:
#'     \code{filtered}: A filtered count table containing features that passed the filtering criteria.
#'     \code{p_agg}: A plot of features removed/kept based on their relative abundance vs. feature prevalence (aggregated across samples).
#'     \code{p_exp}: An expanded view plot showing the relative abundance for each feature across all samples.
#'     \code{feature_keep}: A vector of features remaining after filtering.
#'     \code{feature_remove}: A vector of features removed during filtering.
#'     \code{msg}: A summary message detailing the number of features kept and removed.
#'
#' @details
#' Filtering is performed based on feature abundance and sample prevalence:
#'     - \code{abs_count}: Features are kept if their read count is greater than \code{feature_cutoff}.
#'     - \code{percent_sample}: Features are kept if their relative abundance in each sample is greater than \code{feature_cutoff}.
#'     - \code{percent_dataset}: Features are kept if their relative abundance across the entire dataset is greater than \code{feature_cutoff}.
#'     Prevalence filtering is applied by \code{prev_cutoff}, ensuring that features appear in at least this many samples after the abundance cutoff.
#'
#'     For \code{feature_cutoff}:
#'     - When \code{filter_method = 'abs_count'}, \code{feature_cutoff} should be an integer representing the minimum read count threshold.
#'     - When \code{filter_method = 'percent_sample' or 'percent_dataset'}, \code{feature_cutoff} should be a percentage value between 0 and 100, representing the minimum relative abundance (percentage).
#'
#'     \code{prev_cutoff} is the minimum number of samples a feature must appear in to be retained. The default value of 2 means that a feature must be present in at least two samples to pass the filter.
#'
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @export
#' @examples
#' data(dss_example)
#'
#' # Put featureID as rownames in the counts_table
#' counts_table <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#'
#' # Filter the metagenomics data by relative abundance in each sample
#' filtered_results <- filter_metagenomics(counts_table, filter_method = 'percent_sample', feature_cutoff = 0.01, prev_cutoff = 2)
#'
#' # View the filtered results
#' filtered_results$filtered
#' summary(filtered_results)
#' filtered_count <- filtered_results$filtered
#' dim(filtered_count)
#' head(filtered_count)

filter_metagenomics <- function(counts_table, 
                                filter_method = 'abs_count', 
                                feature_cutoff = 1, 
                                prev_cutoff = 2) {

  if(class(counts_table) != 'data.frame') stop("counts_table must be a dataframe")
  if(!filter_method %in% c('abs_count', 'percent_sample', 'percent_dataset')) stop("filter_method must be 'abs_count', 'percent_sample', or 'percent_dataset'")
  if(prev_cutoff > ncol(counts_table)) stop("prev_cutoff cannot exceed the number of samples")
  
  if(filter_method == 'abs_count') {
    if(feature_cutoff %% 1 != 0) stop("feature_cutoff must be an integer for 'abs_count' method")
    if(feature_cutoff > max(counts_table)) stop("feature_cutoff cannot exceed the max read count")
  } else if(feature_cutoff < 0 | feature_cutoff > 100) {
    stop("feature_cutoff must be between 0-100 for 'percent_sample' or 'percent_dataset'")
  }

  # Calculate relative abundance (by sample and by dataset)
  relab_by_sample <- apply(counts_table, 2, function(x) x / sum(x))
  relab_by_data <- counts_table / sum(counts_table)
  prev_data <- rowSums(counts_table != 0)

  # Perform filtering based on the chosen method
  filtered_count <- switch(
    filter_method,
    'abs_count' = counts_table[rowSums(counts_table >= feature_cutoff) >= prev_cutoff, ],
    'percent_sample' = relab_by_sample[rowSums(relab_by_sample >= feature_cutoff) >= prev_cutoff, ],
    'percent_dataset' = relab_by_data[rowSums(relab_by_data >= feature_cutoff) >= prev_cutoff, ]
  )

  if(is.null(filtered_count)) stop("No features remain after filtering")

  # Check for samples with no features remaining after filtering
  sampleID <- colnames(filtered_count)
  empty_sample_ind <- which(colSums(filtered_count) == 0)
  if(length(empty_sample_ind) > 0) {
    filtered_count <- filtered_count[, sampleID[-empty_sample_ind]]
    message(sprintf("%s samples contained 0 reads after feature filtering", length(empty_sample_ind)))
  }

  # Identify and record features that remain and those that were removed
  feature_keep <- rownames(filtered_count)
  feature_remove <- setdiff(rownames(counts_table), feature_keep)

  # Prepare data for plotting using pivot_longer instead of gather
  pdata <- counts_table %>%
    tibble::rownames_to_column('featureID') %>%
    tidyr::pivot_longer(cols = -featureID, names_to = 'sampleID', values_to = 'value') %>%
    left_join(enframe(prev_data, 'featureID', 'prevalence'), 'featureID') %>%
    mutate(colour = ifelse(featureID %in% feature_keep, 'to keep', 'to remove'))

  pdata_agg <- pdata %>%
    group_by(featureID) %>%
    mutate(y = mean(value)) %>%
    distinct(prevalence, y, colour)

  p_agg <- ggplot(pdata_agg, aes(x = prevalence, y = y, colour = colour)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab('Feature prevalence (# of samples)') +
    ylab("Mean Relative Abundance") +
    scale_x_continuous(labels = scales::label_number(), limits = c(0, ncol(counts_table))) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  p_exp <- ggplot(pdata, aes(x = prevalence, y = value, colour = colour)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab('Feature prevalence (# of samples)') +
    ylab("Relative Abundance") +
    scale_x_continuous(labels = scales::label_number(), limits = c(0, length(sampleID))) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  # Summary message
  msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s in at least %s/%s (%s%%) samples', 
                 length(feature_keep), nrow(counts_table), round(length(feature_keep) / nrow(counts_table) * 100, 2), 
                 feature_cutoff, prev_cutoff, ncol(counts_table), round(prev_cutoff / ncol(counts_table) * 100, 1))

  message(msg)

  return(list('filtered' = filtered_count,
              'p_agg' = p_agg, 'p_exp' = p_exp, 'feature_remove' = feature_remove,
              'feature_keep' = feature_keep, 'msg' = msg))
}
