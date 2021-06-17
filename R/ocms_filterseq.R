#' filterFeature.R
#'
#' filter out reads based on cutoff threshold and asv prevalence across samples
#'
#' @param count_df dataframe. count table with samples in columns and
#'     ASV in rows. feature ID in rownames.
#' @param tax_df dataframe. featureID must match rownames \code{count_df}.
#'     has columns \code{'featureID','Kingdom','Phylum',
#'     'Class','Order','Family','Genus','Species'}
#' @param filter_method default \code{"abs_count"} must be one of
#'     \code{c("abs_count", "percent_sample", "percent_dataset")}.
#'     Therefore, ASVs must reach a certain percentage of the entire dataset
#' @param asv_cutoff cutoff used to filter sequences. features are kept when
#'     they are greater than this cutoff
#' @param prev_cutoff  prevalence cutoff. ASVs must reach the \code{asv_cutoff}
#'     in at least this many samples to be kept.
#' @return
#'     list of:
#'     filtered_table - filtered
#'     Also returns list of:
#'     \code{p_agg} - plot of sequences removed/kept based on relative abundance
#'     vs asv prevalence in aggregated (mean ASV relative abundance)
#'     \code{p_exp} - expanded view (ASV relative abundance for every sample shown).
#'     \code{feat_keep} - vector of ASVs remaining after filtering
#'     \code{feat_remove} - vector of ASVs removed during filtering
#'
#' @details
#'
#' Filtering is performed based on read count and sample prevalence.
#'     ASVs are kept if they pass the ASV count cut-off OR if they pass
#'     the sample prevalence cut-off.
#' \code{asv_cutoff = 'abs_count'} uses a read count as a threshold cutoff.
#'     recommended default of \code{1}
#' When \code{asv_cutoff} is set to \code{'percent_sample'} uses percent of
#'     sample total read count as the threshold cutoff.
#'     Therefore, ASVs must reach a certain percentage of a given sample.
#'     Recommended default of \code{0.01} for 0.01% of each sample
#' When \code{asv_cutoff} is set to \code{'percent_dataset'} uses percent
#'     of dataset total read count as the threshold cutoff.
#'     Recommended default of \code{0.01} for 0.01% of entire dataset
#' \code{prev_cutoff} has minimum value of 1 (sequence must reach cutof in
#'     at least 1 sample, which would not filter out any sequences).
#'     Default value is set to 2, which is the most relaxed cutoff
#'     A recommended default is the number of samples to make up 5% of total number of samples.
#'
#' @import tibble
#' @import tidyr
#' @import dplyr
#'
#' @examples
#' data(dss_example)
#' # put featureID as rownames
#' tax_df <- dss_example$merged_taxonomy
#' count_df <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#' # set features in count tax to be in same order
#' count_df <- count_df[tax_df$featureID,]
#' filtered_ls <- filterFeature(count_df, tax_df, 'abs_count', 1, 2)
#' summary(filtered_ls)
#' filtered_count <- filtered_ls$filtered

filterFeature <- function(count_df, tax_df,
                           filter_method = 'abs_count',
                           asv_cutoff = 1, prev_cutoff = 2) {

  tax_level <- c('Kingdom','Phylum','Class','Order', 'Family','Genus',
                 'Species','featureID')

  # set x axis scale to 1 decimal
  scaleFUN <- function(x) sprintf("%.0f", x)

  # input checks ---------------------------------------------------------------
  # count_df must be dataframe
  if(class(count_df) != 'data.frame') {
    stop("count_df must be dataframe")
  }

  # filter_method must be one of "abs_count" "percent_sample" "percent_dataset
  if(!filter_method %in% c('abs_count','percent_sample','percent_dataset')) {
    stop("filter_method must be 'abs_count','percent_sample', or 'percent_dataset'")
  }

  # prev_cutoff must be integer between 0 and max number of samples
  if(prev_cutoff > ncol(count_df)) {
    stop("prev_cutoff cannot exceed the number of samples")
  }

  if(filter_method == 'abs_count') {
    # asv_cutoff must be integer when filter_method == 'abs_count'
    if(asv_cutoff %% 1 != 0) {
      stop("asv_cutoff must be an integer when filter_method = 'abs_count'")
    }
    # asv_cutoff must not be greater than the max number of reads when filter_method = 'abs_count'
    if(asv_cutoff > max(count_df[,2:ncol(count_df)])) {
      stop("asv_cutoff cannot be greater than the maximum number of reads when filter_method = 'abs_count'")
    }
  } else {
    # asv_cutoff must be between 0-100 when filter_method relies on percent
    if(asv_cutoff < 0 | asv_cutoff > 100) {
      stop("asv_cutoff must be between 0-100 when filter_method is 'percent_sample' or 'percent_dataset'")
    }
  }

  # check order of count_df is same as order of tax_df
  if(!identical(rownames(count_df), tax_df$featureID)) {
    stop("Order of features in count_df and tax_df must be identical")
  }

  # calculate different relative abundances-------------------------------------
  relab_by_sample <- apply(count_df, 2, function(x) x/sum(x))
  relab_by_data <- count_df / sum(count_df)

  # calculate prevalence based on cutoff set and perform filtering-------------
  binary <- (count_df != 0) * 1 # true means present
  prev_data <- rowSums(binary)

  # perform filtering-----------------------------------------------------------
  filtered_count <- switch(
    filter_method,
    # count cut-off
    abs_count = count_df[rowSums(count_df >= asv_cutoff) >= prev_data,],
    # cut-off based on percent of sample total
    percent_sample = relab_by_sample[rowSums(relab_by_sample >= asv_cutoff) >= prev_data,],
    # cut-off based on percent of dataset total
    percent_total = relab_by_data[rowSums(relab_by_data >= asv_cutoff) >= prev_data,]
  )

  if(is.null(filtered_count)) {
    stop("No features remain after filtering. Please select different filter cutoff.")
  }
  # secondary check for samples with no reads-----------------------------------
  sampleID <- colnames(filtered_count)
  empty_sample_ind <- which(colSums(filtered_count)==0)

  if(length(empty_sample_ind) > 0) {
    filtered_count <- filtered_count[sampleID[!empty_sample_ind],]

    msg <- sprintf("%s samples contained 0 reads after ASV filtering. The following samples have been removed: %s", length(empty_sample_ind), paste(sampleID[empty_sample_ind], collapse="\n"))

    message(msg)
  }

  # secondary check for empty asvs----------------------------------------------
  # identify empty asvs
  featID <- rownames(filtered_count)
  empty_feat_ind <- which(rowSums(filtered_count)==0)

  if(length(empty_feat_ind) > 0) {

    filtered_count <- filtered_count[,featID[!empty_which_ind]]
    msg <- sprintf("%s asvs contained 0 reads in all samples after ASV filtering. The following asvs have been removed: %s", length(empty_feat_ind), paste(featID[empty_feat_ind], collapse="\n"))

    message(msg)

  }

  # recording asvs kept and removed---------------------------------------------
  feat_keep <- rownames(filtered_count)
  feat_remove <- setdiff(rownames(count_df), rownames(filtered_count))

  # filter taxonomy table-------------------------------------------------------
  tax_filtered <- tax_df %>%
    filter(featureID %in% feat_keep)

  # preparing plot data---------------------------------------------------------

  if(filter_method == 'abs_count') {
    pdata <-  count_df %>%
      tibble::rownames_to_column('featureID') %>%
      gather('sampleID','value', -featureID)
    ylab <- 'Read count'
  }
  if(filter_method == 'percent_total') {
    pdata <- relab_by_data %>%
      as.data.frame() %>%
      tibble::rownames_to_column('featureID') %>%
      gather('sampleID','value', -featureID)
    ylab <- 'Relative abundance\n(% of total reads)'
  }
  if(filter_method == 'percent_sample') {
    pdata <- relab_by_sample %>%
      as.data.frame() %>%
      tibble::rownames_to_column('featureID') %>%
      gather('sampleID','value', -featureID)
    ylab <- 'Relative abundance\n(% of sample)'
  }

  pdata <- pdata %>%
    left_join(enframe(prev_data, 'featureID','prevalence'), 'featureID') %>%
    mutate(colour = ifelse(featureID %in% feat_keep,
                           'to keep', 'to remove'))

  # plot aggregated view--------------------------------------------------------
  pdata_agg <- pdata %>%
    group_by(featureID) %>%
    summarise(y = mean(value))

  p_agg <- ggplot(pdata_agg, aes(x = prevalence, y = y, colour = colour)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab('ASV prevalence (# of samples)') +
    ylab(sprintf("Mean %s", ylab)) +
    scale_x_continuous(labels = scaleFUN, limits = c(0, ncol(count_df))) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  # plot expanded view----------------------------------------------------------
  p_exp <- ggplot(pdata,
                  aes(x = prevalence, y = value, colour = colour)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab('Feature prevalence (# of samples)') +
    ylab(ylab) +
    scale_x_continuous(labels = scaleFUN, limits = c(0, length(sampleID))) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  # filter message--------------------------------------------------------------
  msg1 <- length(feat_keep)
  msg2 <- nrow(count_df)
  msg3 <- round(msg1 / msg2 * 100, 2)
  msg4 <- asv_cutoff
  msg5 <- prev_cutoff
  msg6 <- ncol(count_df)
  msg7 <- round(msg5 / msg6 * 100, 1)

  if(filter_method == 'abs_count') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s with total read count in >= %s/%s (%s%%) samples', msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  } else if(filter_method == 'percent_dataset') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s%% with dataset total read count in >= %s/%s (%s%%) samples', msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  } else if(filter_method == 'percent_sample') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s%% with sample total read count in >= %s/%s (%s%%) samples', msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  }

  message(msg)


  return(list('taxonomy' = tax_filtered, 'filtered' = filtered_count,
              'p_agg'=p_agg, 'p_exp'=p_exp, 'feat_remove'=feat_remove,
              'feat_keep' = feat_keep, 'msg' = msg))
}
