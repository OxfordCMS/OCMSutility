#' rarefaction
#'
#' produce rarefaction plot.
#' adapted from MetaSequencingSnake enumeration_report.Rmd
#'
#' @param df dataframe with samples in columns, features in rows
#' @return list of ggplot rarefaction curve and associated dataframe
#' @export
#' @import foreach
#' @import vegan
#' @examples
#' count_data <-  cbind(matrix(rnorm(100*5,mean=50000,sd=10000), 100, 5),
#'                      matrix(rnorm(100*5,mean=20000,sd=5000), 100, 5),
#'                      matrix(rnorm(100*5,mean=5000,sd=1000), 100, 5))
#' rarefaction(as.data.frame(count_data))

rarefaction <- function(df) {

  if(!class(df) %in% c('data.frame', 'matrix')) {
    stop('input data must be a dataframe or matrix with samples in columns, ASVs in rows, with respective column and row names')
  }
  `%dopar%` <- foreach::`%dopar%`
  # getting sample read depth
  sampmax = colSums(df)
  raredepths = round(c(seq(from=1, to=max(sampmax),
                         by=(max(sampmax)-1)/20)))

  # initiate rarefaction values matrix
  vals = matrix(nrow=length(raredepths), ncol=ncol(df))

  # calculate rarefaction
  forres = foreach::foreach(i = 1:length(raredepths)) %dopar% {
    depth = raredepths[i]
    res = suppressWarnings({
      vegan::rarefy(round(t(df)), depth)
    })
    res[sampmax < depth] = NA
    return(res)
  }

  # populate matrix with rarefaction values
  for(i in 1:length(forres)){
    vals[i,] = unlist(forres[i])
  }

  colnames(vals)=colnames(df)
  rownames(vals)=raredepths

  raremelt = reshape2::melt(vals)
  colnames(raremelt) = c("Depth","Sample","Richness")
  raremelt = raremelt[!is.na(raremelt$Richness),]

  labels = data.frame(Sample = colnames(df),
                      Depth = sampmax,
                      Rich = apply(vals, 2, function(x){max(x,na.rm = T)}))

  # # only label lowest 10 samples
  # labels$Sample = ifelse(rank(labels$Depth)<11,
  #                        as.character(labels$Sample), NA)

  rplot = ggplot(raremelt, aes(x=Depth, y=Richness, color=Sample)) +
    geom_line()+
    guides(color = 'none')+
    ggrepel::geom_label_repel(data = labels,
                              aes(x=Depth, y=Rich, label=Sample),
                              fill = alpha(c("white"),0.2))

  return(list('rare_df' = raremelt, 'rare_p' = rplot))
}