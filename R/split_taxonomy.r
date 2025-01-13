#' split_taxonomy
#'
#' Splits taxonomy information from row names of a counts table and aggregates data by taxonomy levels.
#'
#' @param counts_table A data frame where row names represent taxonomy information (in the format "k__;p__;c__;o__;f__;g__;s__") and columns represent sample counts. The row names should contain taxonomy information for each species, with levels separated by semicolons.
#'
#' @details This function takes a counts table with taxonomy information encoded in the row names. It splits the taxonomy into separate levels (kingdom, phylum, class, order, family, genus, species) and aggregates the count data at each of these levels.
#'
#' @import stats
#'
#' @returns A list of data frames, each containing aggregated counts data at a specific taxonomy level (kingdom, phylum, class, order, family, genus, species).
#'
#' @export
#'
#' @examples
#' # Example usage:
#' result <- split_taxonomy(counts_table)
#' result$kingdom
#' result$phylum
#' result$class
#' result$order
#' result$family
#' result$genus
#' result$species
#'
split_taxonomy <- function(counts_table) {
  tax_split <- strsplit(rownames(counts_table), ";")
  
  tax_levels <- data.frame(
    kingdom = gsub("k__", "", sapply(tax_split, `[`, 1)),
    phylum = gsub("p__", "", sapply(tax_split, `[`, 2)),
    class = gsub("c__", "", sapply(tax_split, `[`, 3)),
    order = gsub("o__", "", sapply(tax_split, `[`, 4)),
    family = gsub("f__", "", sapply(tax_split, `[`, 5)),
    genus = gsub("g__", "", sapply(tax_split, `[`, 6)),
    species = gsub("s__", "", sapply(tax_split, `[`, 7))
  )
  
  species_names <- tax_levels$species
  counts_data <- as.data.frame(lapply(counts_table, as.numeric))
  
  process_level <- function(data, group) {
    result <- aggregate(data, by = list(group), sum, na.rm = TRUE)
    rownames(result) <- result[,1]
    result[,-1]
  }
  
  list(
    kingdom = process_level(counts_data, tax_levels$kingdom),
    phylum = process_level(counts_data, tax_levels$phylum),
    class = process_level(counts_data, tax_levels$class),
    order = process_level(counts_data, tax_levels$order),
    family = process_level(counts_data, tax_levels$family),
    genus = process_level(counts_data, tax_levels$genus),
    species = process_level(counts_data, species_names)
  )
}
