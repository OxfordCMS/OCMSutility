#' ocms_palette
#' generates colour palette from RColorBrewer
#'
#' @param n  numeric; default NULL. number of colours; cannot exceed 335. when NULL, gives all colours (in palette, if specified).
#' @param palette character; allows you to specify a RColorBrewer palette by name
#'        Diverging
#'        BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
#'
#'        Qualitative
#'        Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
#'
#'        Sequential
#'        Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu,
#'        PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd
#' @param preview logic; default FALSE; when set, gives preview of colour
#' @export
#' @return vector of HEX colours, if full set to TRUE, gives all 335 colours, else gives n numbers
#' @examples
#' cms_palette(n=5, palette = 'Set3')
#' cms_palette(palette = c('Set1','Dark2'), preview = TRUE)


ocms_palette <- function(n = NULL, palette = NULL, preview = FALSE) {

  # pulling colours from RColorBrewer-------------------------------------------

  # putting all RColorBrewer colours into dataframe
  full_palette <- RColorBrewer::brewer.pal.info
  full_palette$palID <- rownames(full_palette)

  ordered_pal <- rbind(dplyr::filter(full_palette, category == 'qual'),
                       dplyr::filter(full_palette, category == 'div'),
                       dplyr::filter(full_palette, category == 'seq'))

  # subsetting based on number of colours
  if(!is.null(palette)) {
    ordered_pal <- dplyr::filter(ordered_pal, palID %in% palette)
    if(is.null(n)) {
      n <- sum(ordered_pal$maxcolors  )
    }
  }
  else {
    if(is.null(n)) {
      n <- 335
    }
  }

  col_vector = unlist(mapply(RColorBrewer::brewer.pal, ordered_pal$maxcolors,
                             ordered_pal$palID))

  if(n > length(col_vector)) {
    stop("Number of colours requested exceeds the number of colours in the specified palette(s). Specify more palettes using the palette argument.")
  }
  col_samp <- col_vector[1:n]

  out <- col_samp

  if(preview==TRUE) pie(rep(1,n), col=out)

  return(out)
}
