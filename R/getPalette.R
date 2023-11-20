#' getPalette
#'
#' @description
#' generates colour palette from RColorBrewer
#' `r lifecycle::badge("deprecated")`
#'
#' `getPalette` was renamed to `get_palette` for more consistent
#' function naming nomenclature in `OCMSutility`. Usage is still exactly the same
#' @keywords internal
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
#' getPalette(n=5, palette = 'Set3')
#' # ->
#' get_palette(palette = c('Set1','Dark2'), preview = TRUE)


getPalette <- function(n = NULL, palette = NULL, preview = FALSE) {
  get_palette(n, palette, preview)
}
