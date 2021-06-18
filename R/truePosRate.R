#' truePosRate
#'
#' Calculate rate of true positives in positive control standards. Used in OCMS_zymobioimcs report.
#'
#' @param relab relative abundance data frame with samples in columns, features in rows
#' @param annotations reference annotations for true positives
#' @param level taxonomy level at which analysis is done. can be one of
#'    \code{'species','genus','family'}. default \code{'species'}
#' @param cutoff threshold for read relative abundance in the sample. default is
#'    \code{0.01} (0.01 %)
#' @return dataframe ranking features as true positive or false positive in each sample.
#' @export
#'
#' @examples
#'  # this would be better exemplified with actual std data rather than the example smaples
#'  data("dss_example")
#' data(zymobiomics)
#'
#' # set count feature ids as rownames
#' count_df <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#'
#' # clean up some sample names
#' colnames(count_df) <- paste0('id', colnames(count_df))
#' tax_df <- dss_example$merged_taxonomy
#'
#' # aggregate counts
#' agg_gen <- aggregateCount(count_df[tax_df$featureID,], tax_df, "Genus")
#' genus_relab <- relab(agg_gen$count_df)
#'
#' true_pos_result <- truePosRate(relab=genus_relab,
#'                                annotations=zymobiomics$anno_ncbi_16s,
#'                                level='genus', cutoff=0.01)
#'
#' # plot true pos rate
#' p <- ggplot(true_pos_result,
#'             aes(x=rank, y=true.pos.rate, colour=label, group=sample)) +
#'  geom_point() +
#'   theme_bw() +
#'   ylab("TP / (TP + FP)") +
#'  scale_colour_manual(values=c("grey", "purple")) +
#'  facet_wrap(~sample, scale="free")
#'
#' p

truePosRate <- function(relab, annotations, level="species", cutoff=0.01){


  # check inputs----------------------------------------------------------------
  # count_df must be dataframe
  if(class(relab) != 'data.frame') {
    stop("count_df must be data.frame")
  }

  if(class(annotations) != 'data.frame') {
    stop("tax_df must be a data.frame")
  }

  # level must be either species or genus
  if(!level %in% c('species','genus','family')) {
    stop("True positive rates can be calculated by 'species','genus', or 'family'")
  }

  # cutoff must be between 0 and 100
  if(cutoff > 100 | cutoff < 0) {
    stop("cutoff must be between 0 and 100")
  }

  # check for empty samples
  if(any(colSums(relab) == 0)) {
    stop("Samples with no reads detected. Please remove empty samples")
  }

  # check for empty features
  if(any(rowSums(relab) == 0)) {
    stop("Features with no reads detected. Please remove empty features")
  }

  # calculate true positive rate------------------------------------------------
  result <- list()
  for (i in 1:ncol(relab)){
    ab <- data.frame(taxon=rownames(relab), abundance = relab[,i], stringsAsFactors = FALSE)
    ab <- ab[ab$abundance >= cutoff,]
    ab <- ab[order(ab$abundance, decreasing=TRUE),]
    trues <- 0
    falses <- 0
    tps <- c()
    labels <- c()
    for (j in 1:nrow(ab)){
      if (ab[j,]$taxon %in% annotations[,level]){
        trues = trues + 1
        labels <- append(labels, "TP")
      } else{
        falses = falses + 1
        labels <- append(labels, "FP")
      }
      tp <- trues/(trues + falses)
      tps <- append(tps, tp)
    }
    sample <- rep(colnames(relab)[i], length(tps))
    res <- data.frame(rank = seq(1:length(tps)),
                      true.pos.rate = tps,
                      sample = sample,
                      label = labels
    )
    result[[i]] <- res
  }
  result <- bind_rows(result)
  return(result)
}