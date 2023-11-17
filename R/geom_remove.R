#' remove_geom
#'
#' Remove geom layer from ggplot. Based on https://stackoverflow.com/questions/13407236/remove-a-layer-from-a-ggplot2-chart
#' @param ggplot2_object ggplot2 object
#' @param geom_type
#'
#' @export
#'
#' @examples
#' d <- data.frame(
#' x = runif(10),
#' y = runif(10),
#' label = sprintf("label%s", 1:10)
#'
#' # ggplot with geom_text_repel from ggrepel
#' p <- ggplot(d, aes(x, y, label = label)) + geom_point() + geom_text_repel()
#'
#' # Remove the labels added by ggrepel.
#' p <- remove_geom(p, "GeomTextRepel")

remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  return(ggplot2_object)
}