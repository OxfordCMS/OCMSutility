#' ocms_relab
#' convert counts table to relative abundance (%)
#'
#' @param counts dataframe; samples are columns and features are rows
#' @export
#' @return dataframe
#' @examples
#' counts <- data.frame(matrix(rnorm(1000), ncol=20, nrow=50))
#' rel_abundance <- ocms_relab(counts)

relab <- function(counts){

    relab <- (sweep(counts, 2, colSums(counts), "/"))*100
    return(relab)
    }
