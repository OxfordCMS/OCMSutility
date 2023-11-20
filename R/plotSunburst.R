#' plotSunburst
#'
#' Creates interactive sunburst plot based on taxonomy. The sunburst plot can show areas based on relative abundance or based on the number of taxa at a given taxonomic level.
#'
#' @param relab dataframe. default NULL. relative abundance data with samples
#'              in columns and features in rows. Feature IDs in the rownames.
#'
#' @param tax dataframe. taxonomy table with featureID, Phylum, Class, Order,
#'            Family, Genus as columns.
#'            See Details about NAs in the taxonomy table
#'
#' @param palettes named vector. default NULL. Can specify a palette for each
#'            Phylum, where values are the colour palette to use and name
#'            is the corresponding phylum (e.g.
#'            \code{c("Bacteroidetes" = "Oranges", "Firmicutes" = "Greens")}).
#'            Palettes should be from rColorBrewer. If the number of palettes
#'            specified doesn't include all phyla in the tax table, only the
#'            specified ones will be coloured and the rest will be in grey.
#'            If palettes is set to NULL, the default colours selected by
#'            \code{sunbrustR} will be used
#'
#' @param highlight named list. default NULL. This can be used to
#'            highlight all taxa relevant to a specific taxon and any taxa that
#'            are not specified will be coloured as grey.
#'            e.g. \code{list("Family"=c("Enterococcaceae","Ruminacoccaceae")}.
#'            This is applied after palettes is used to colour by phylum
#'            if palettes argument is specified.
#' @param ... additional arguments for \code{sunburstR::sunburst()}
#'
#' @details
#' When \code{relab} is set to NULL, the sunburst plot will show the number of
#' taxa observed at each taxonomic level. When relative abundance is supplied,
#' sunburst leaves reflect mean relative abundance of each taxon
#' across all samples
#'
#' Note NAs in the taxonomy table cause colouring to be assigned in unexpected
#' order so it is best to use \code{ocms_reannotateTax} to apply a taxonomy
#' roll-down and remove all NAs. \code{sunburstR} uses hyphens (\code{-}) to
#' distinguish taxonomic levels so any hyphens in the taxonomy name will be
#' interpreted as two separate levels. Therefore, all hyphens are silently and
#' automatically removed from taxonomy names.
#'
#' @return interactive sunburst plot that can be included in rmarkdown or shiny.
#' @import tibble
#' @import dplyr
#' @import RColorBrewer
#' @import sunburstR
#' @examples
#' data("dss_example")
#' #' set count feature ids as rownames
#' count_df <- dss_example$merged_abundance_id %>%
#'   column_to_rownames('featureID')
#'
#' #' clean up some sample names
#' colnames(count_df) <- paste0('id', colnames(count_df))
#' tax_df <- dss_example$merged_taxonomy
#'
#' #' aggregate counts
#' agg_gen <- aggregateCount(count_df[tax_df$featureID,], tax_df, "Genus")
#' count_genus <- agg_gen$count_df
#'
#' #' reannotate taxonomy
#' tax_genus <- reannotateTax(agg_gen$tax_df)
#'
#' relab <- relab(count_genus)
#' #' color specific phyla
#' plotSunburst(relab = NULL, tax = tax_genus,
#'              palettes = c("Proteobacteria" = "Oranges",
#'                           "Bacteroidetes" = "Greens"))
#' #' color specific phyla taking into account of relative abundance
#' plotSunburst(relab = relab, tax = tax_genus,
#'              palettes = c("Proteobacteria" = "Oranges", "Bacteroidetes" = "Greens"))
#'
#' #' highlight specific genera
#' plotSunburst(relab = relab, tax = tax_genus,
#'              palettes = c("Bacteroidetes" = "Greens",'Firmicutes'='Blues'),
#'              highlight = list("Genus" = c("Bacteroides",'Clostridium XlVa')))


plotSunburst <- function(tax, relab=NULL, palettes=NULL, highlight=NULL, ...) {

  # check inputs----------------------------------------------------------------
  # check data is in relative abundance format
  if(!is.null(relab)) {
    if(all(relab %% 1 == 0)) {
      stop("All values are integers. Values should be in relative abundance format")
    }
  }
  # tax table has required columns
  if(any(!c('featureID','Phylum','Class','Order','Family','Genus') %in% colnames(tax))) {
    stop("tax must have columns featureID','Phylum','Class','Order','Family','Genus'")
  }
  # named phyla are found in tax table
  if(!is.null(palettes)) {
    phyla <- names(palettes)
    if(any(!phyla %in% tax$Phylum)) {
      not_found <- which(!phyla %in% tax$Phylum)
      msg <- sprintf("%s not a Phylum in taxonomy table", phyla[not_found])
      stop(msg)
    }
  }
  # highlight is named list
  if(!is.null(highlight)) {
    if(class(highlight) != 'list') {
      stop("highlight should be NULL or a list")
    }
    # named list in correct format
    if(any(!names(highlight) %in% colnames(tax))) {
      stop("List items should be named as Phylum, Class, Order, Family, or Genus")
    }
    if(unique(sapply(highlight, class)) != 'character') {
      stop("List items should be a character vector")
    }
    # named taxa are found in tax table
    for(i in 1:length(highlight)) {
      not_found <- highlight[[i]][!highlight[[i]] %in% pull(tax, names(highlight[i]))]
      if(length(not_found) != 0) {
        msg <- sprintf("%s not a %s in taxonomy table", not_found, names(highlight[i]))
        stop(msg)
      }
    }
  }

  # format taxonomy table-------------------------------------------------------
  tax <- tax %>% select(featureID, Phylum, Class, Order, Family, Genus)

  if(is.null(relab)) {
    pdata <- tax %>%
      # get number of taxa within each level
      gather('tax_level', 'taxonomy', -featureID) %>%
      group_by(tax_level, taxonomy) %>%
      summarise(featureID = featureID, n_taxa = n()) %>%
      filter(tax_level == 'Genus') %>%
      ungroup() %>%
      select(featureID, n_taxa) %>%
      right_join(tax, 'featureID') %>%
      mutate(path = paste(Phylum, Class, Order, Family, Genus, sep='-')) %>%
      select(path, n_taxa)
  } else {
    pdata <- relab %>%
      rownames_to_column('featureID') %>%
      gather('sampleID','relab', -featureID) %>%
      # aggregate all samples together
      group_by(featureID) %>%
      summarise(avg_relab = mean(relab, na.rm = TRUE)) %>%
      # add taxonomy data
      right_join(tax, 'featureID') %>%
      # format so tax levels are separated by "-"
      mutate(path = paste(Phylum, Class, Order, Family, Genus, sep='-')) %>%
      select(path, avg_relab)

    print(head(pdata))
  }

  # specify colours-------------------------------------------------------------
  # get unique leaf for each phylum
  phyla <- unique(tax$Phylum)
  leaf <- list()
  for(i in phyla) {
    curr_data <- tax %>%
      select(-featureID) %>%
      filter(Phylum == i)
    iter <- c(curr_data)
    leaf_entry <- c()
    # get unique taxa while preserving order
    for(j in 1:length(iter)) {
      entry <- iter[[j]][which(!iter[[j]] %in% leaf_entry)]
      entry <- gsub("-","", entry)
      leaf_entry <- unique(c(leaf_entry, as.character(entry)))
    }
    leaf[[i]] <- leaf_entry
  }

  # get colours from one palette for each phylum
  # from https://stackoverflow.com/questions/49993198/how-to-specify-the-colors-and-toggle-labels-for-each-category-in-r-sunburst
  if(!is.null(palettes)) {
    col_df <- c()
    # go through one phylum at a time
    for(i in 1:length(leaf)) {
      curr_phylum <- names(leaf[i])
      curr_leaf <- leaf[[i]]

      if(curr_phylum %in% names(palettes)) {
        # colours from specified palette
        get_col <- colorRampPalette(brewer.pal(9, palettes[curr_phylum]))
        entry <- data.frame(range = get_col(length(leaf[[i]])),
                            domain = leaf[[i]])

        if(any(is.na(entry$domain))) {
          msg <- sprintf('NAs detected in Phylum %s. Features assigned as NA are coloured as grey', curr_phylum)
          message(msg)
          entry$range[is.na(entry$range)] <- '#D9D9D9'
        }

      } else {
        # set colour to grey
        entry <- data.frame(range = rep('#D9D9D9', length(leaf[[i]])),
                            domain = leaf[[i]])
      }
      col_df <- rbind(col_df, entry)
    }
  }


  if(!is.null(highlight)) {
    tax_highlight <- c()
    for(i in 1:length(highlight)) {
      entry <- tax %>%
        select(-featureID) %>%
        filter(.data[[names(highlight[i])]] %in% highlight[[i]]) %>%
        gather('tax_level','taxonomy')
      tax_highlight <- rbind(tax_highlight, entry)
    }

    tax_highlight <- unique(tax_highlight)

    # check palette was specified for highlight taxa
    highlight_phyla <- tax_highlight$taxonomy[tax_highlight$tax_level == 'Phylum']
    no_palette <- highlight_phyla[!highlight_phyla %in% names(palettes)]
    if(length(no_palette) != 0) {
      msg <- sprintf("No palette specified for the phylum %s\n", no_palette)
      stop(msg)
    }

    # get index of values to keep coloured
    col_ind <- which(col_df$domain %in% tax_highlight$taxonomy)

    col_df$range[-col_ind] <- '#D9D9D9'

  }


  if(is.null(palettes) & is.null(highlight)) {
    p_sun <- sunburst(pdata, ...)
  } else {

    p_sun <- sunburst(pdata,
                      colors = list(range=col_df$range,
                                    domain=col_df$domain),
                      ...)
  }

  return(p_sun)
}