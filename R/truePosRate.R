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
#'
#' @export
#'
#' @examples
#'  # this would be better exemplified with actual std data rather than the example smaples
#' # this would be better exemplified with actual std data rather than the example samples
#' data(asv_example)
#' data(zymobiomics)
#'
#' # set feature IDs in rownames
#' count_df <- asv_example %>%
#'   tibble::column_to_rownames('sequence')
#' count_df <- count_df[tax_df$sequence, ]
#'
#' # convert count to relative abundance
#' relab <- ocms_relab(count_df)
#'
#' # aggregate on genus level
#' agg_ls <- ocms_aggregateCount(relab, tax_df, 'Genus')
#' genus_relab <- agg_ls$count_df
#'
#' # examine subset of samples
#' genus_relab <- genus_relab[,1:10]
#'
#' true_pos_result <- truePosRate(relab=genus_relab,
#'                                     annotations=zymobiomics$anno_ncbi_16s,
#'                                     level='genus', cutoff=0.01)
#'
#' # plot true pos rate
#' p <- ggplot(true_pos_result,
#'             aes(x=rank, y=true.pos.rate, colour=label, group=sample)) +
#'   geom_point() +
#'   theme_bw() +
#'   ylab("TP / (TP + FP)") +
#'   scale_colour_manual(values=c("grey", "purple")) +
#'   facet_wrap(~sample, scale="free")

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
}