#' ocms_filterseq.R
#'
#' filter out reads based on cutoff threshold and asv prevalence across samples
#'
#' @param count_df dataframe. count table with samples in columns and
#'     ASV in rows. First column is 'featureID'
#' @param tax_df dataframe. featureID must match \code{count_df}.
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
#' @import stringr
#' @import tidyr
#' @import dplyr

ocms_filterseq <- function(count_df, tax_df,
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
  if(!identical(count_df$featureID, tax_df$featureID)) {
    stop("Order of features in count_df and tax_df must be identical")
  }

  # calculate different relative abundances-------------------------------------
  working_count <- count_df %>%
    gather('sampleID','read_count', -featureID) %>%
    left_join(tax_df, 'featureID')

  nsample <- ncol(count_df) -1
  dataset_total <- sum(working_count$read_count)

  working_count <- working_count %>%
    # calculate relative abundance as dataset total
    mutate(ds_rel_abund = read_count / dataset_total * 100) %>%
    # calculate relative abundance as sample total
    group_by(sampleID) %>%
    mutate(sample_total = sum(read_count)) %>%
    group_by(sampleID, featureID) %>%
    mutate(samp_rel_abund = read_count / sample_total * 100) %>%
    ungroup()

  # calculate prevalence based on cutoff set and perform filtering-------------
  # count cut-off
  if(filter_method == 'abs_count') {

    filtered_count <- working_count %>%
      group_by(featureID) %>%
      mutate(prev_observed = sum(read_count >= asv_cutoff)) %>%
      filter(read_count >= asv_cutoff || prev_observed >= prev_cutoff)
  }
  # cut-off based on percent of sample total
  if(filter_method == 'percent_sample') {
    filtered_count <- working_count %>%
      group_by(featureID) %>%
      mutate(prev_observed = sum(samp_rel_abund >= asv_cutoff)) %>%
      filter(samp_rel_abund >= asv_cutoff || prev_observed >= prev_cutoff)

  }
  # cut-off based on percent of dataset total
  if(filter_method == 'percent_dataset') {
    filtered_count <- working_count %>%
      group_by(featureID) %>%
      mutate(prev_observed = sum(ds_rel_abund >= asv_cutoff)) %>%
      filter(ds_rel_abund >= asv_cutoff || prev_observed >= prev_cutoff)
  }

  # secondary check for samples with no reads-----------------------------------
  empty_sample <- filtered_count %>%
    group_by(sampleID) %>%
    summarise(sample_total = sum(read_count)) %>%
    filter(sample_total == 0)

  if(length(unique(empty_sample$sampleID)) > 0) {
    filtered_count <- filtered_count %>%
      filter(!sampleID %in% empty_sample)

    msg <- sprintf("%s samples contained 0 reads after ASV filtering. The following samples have been removed: %s", length(empty_sample), paste(empty_sample, collapse="\n"))

    message(msg)
  }

  # secondary check for empty asvs----------------------------------------------
  # identify empty asvs
  empty_feat <- filtered_count %>%
    group_by(featureID) %>%
    summarise(feat_total = sum(read_count)) %>%
    filter(feat_total == 0)

  if(length(unique(empty_feat$featureID)) > 0) {

    filtered_count <- filtered_count %>%
      filter(!featureID %in% empty_feat)
    msg <- sprintf("%s asvs contained 0 reads in all samples after ASV filtering. The following asvs have been removed: %s", length(empty_feat), paste(empty_feat, collapse="\n"))

    message(msg)

  }

  # put filtered dataframe back as wide count table-----------------------------
  filtered <- filtered_count %>%
    select(featureID, sampleID, read_count) %>%
    spread(sampleID, read_count)

  # recording asvs kept and removed---------------------------------------------
  feat_keep <- unique(filtered_count$featureID)
  feat_remove <- working_count %>% filter(!featureID %in% feat_keep)
  feat_remove <- unique(feat_remove$featureID)

  # filter taxonomy table-------------------------------------------------------
  tax_filtered <- tax_df %>%
    filter(featureID %in% feat_keep)

  # preparing plot data---------------------------------------------------------
  pdata_agg <- working_count %>%
    group_by(featureID) %>%
    # number of samples in which asv meets read cutoff
    mutate(prev_observed = sum(read_count >= asv_cutoff)) %>%
    group_by(featureID) %>%
    summarise(prev_observed = prev_observed,
              agg_count = sum(read_count),
              agg_samp_rel = sum(samp_rel_abund),
              agg_ds_rel = sum(ds_rel_abund)) %>%
    mutate(colour = ifelse(featureID %in% feat_keep, 'to keep', 'to remove')) %>%
    distinct()

  if(filter_method == 'abs_count') {
    y_agg <- 'agg_count'
    ylab <- 'Aggregated Read count'
  }
  if(filter_method == 'percent_dataset') {
    y_agg <- 'agg_ds_rel'
    ylab <- 'Aggregated Relative abundance\n(% of total reads)'
  }
  if(filter_method == 'percent_sample') {
    y_agg <- 'agg_samp_rel'
    ylab <- 'Aggregated Relative abundance\n(% of sample)'
  }

  # plot aggregated view--------------------------------------------------------
  p_agg <- ggplot(pdata_agg, aes_string(x = 'prev_observed', y = y_agg, colour = 'colour')) +
    geom_point(aes(text=sprintf("featureID: %s", featureID)), alpha = 0.5, size = 2) +
    xlab('ASV prevalence (# of samples)') +
    ylab(ylab) +
    scale_x_continuous(labels = scaleFUN, limits = c(0, nsample)) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  # plot expanded view----------------------------------------------------------
  if(filter_method == 'abs_count') {
    y_exp <- 'read_count'
    ylab <- 'Read count'
  }
  if(filter_method == 'percent_dataset') {
    y_exp <- 'ds_rel_abund'
    ylab <- 'Relative abundance\n(% of total reads)'
  }
  if(filter_method == 'percent_sample') {
    y_exp <- 'samp_rel_abund'
    ylab <- 'Relative abundance\n(% of sample)'
  }

  pdata_exp <- working_count %>%
    group_by(featureID) %>%
    # number of samples in which asv meets read cutoff
    mutate(prev_observed = sum(read_count >= asv_cutoff),
           colour = ifelse(featureID %in% feat_keep, 'to keep', 'to remove')) %>%
    distinct()

  p_exp <- ggplot(pdata_exp, aes_string(x = 'prev_observed', y = y_exp, colour = 'colour')) +
    geom_point(aes(text=sprintf("featureID: %s<br>sampleID: %s",
                                featureID, sampleID)), alpha = 0.5, size = 2) +
    xlab('Feature prevalence (# of samples)') +
    ylab(ylab) +
    scale_x_continuous(labels = scaleFUN, limits = c(0, nsample)) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  # filter message--------------------------------------------------------------
  msg1 <- length(feat_keep)
  msg2 <- length(unique(working_count$featureID))
  msg3 <- round(msg1 / msg2 * 100, 2)
  msg4 <- asv_cutoff
  msg5 <- prev_cutoff
  msg6 <- nsample
  msg7 <- round(msg5 / msg6 * 100, 1)

  if(filter_method == 'abs_count') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s with total read count in >= %s/%s (%s%%) samples', msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  } else if(filter_method == 'percent_dataset') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s%% with dataset total read count in >= %s/%s (%s%%) samples', msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  } else if(filter_method == 'percent_sample') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s%% with sample total read count in >= %s/%s (%s%%) samples', msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  }

  message(msg)


  return(list('taxonomy' = tax_filtered, 'filtered' = filtered,
              'p_agg'=p_agg, 'p_exp'=p_exp, 'feat_remove'=feat_remove,
              'feat_keep' = feat_keep, 'msg' = msg))
}
