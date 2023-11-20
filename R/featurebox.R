#' featurebox
#' boxplot expression/abundance values of a feature(s) of interest from a matrix
#'
#' @param abundance_matrix dataframe; samples are columns and features are rows
#' @param metadata dataframe; samples are rows and rownames should match those in abundance_matrix
#' @param features vector; character vector of features to plot i.e. a feature(s) in rownames of abundance_matrix
#' @param group_by character; grouping variable for boxplot. If set then variable must be present in metadata.
#' If NULL each sample is plotted separately.
#' @import ggplot2
#' @import reshape2
#' @import gtools
#' @import RColorBrewer
#' @export
#' @return grob
#' @examples
#' abundance_matrix <- data.frame(matrix(rnorm(1000), ncol=20, nrow=50))
#' colnames(abundance_matrix) <- LETTERS[1:20]
#' metadata <- data.frame(group = c(rep("group1", 5), rep("group2", 5),
#'                                  rep("group3", 5), rep("group4", 5)))
#' rownames(metadata) <- LETTERS[1:20]
#' featurebox(abundance_matrix, metadata=metadata, features="10", group_by="group")
#'
#' # multiple features
#' featurebox(abundance_matrix, metadata=metadata, features=c("10", "20", "30", "40"))

featurebox <- function(abundance_matrix, metadata=NULL, features=NULL, group_by=NULL){

    # must have a feature specified
    if (is.null(features)){
        stop("must specify features to plot")
    }else{
        # check features are in the data
	for (feature in features){
	    if (!(feature %in% rownames(abundance_matrix))){
	        stop(paste0("feature ", feature, " not found in abundance_table"))
            }
	}
    }

    # check if the group_by variable is NULL
    if (is.null(group_by)){
        warning("no group_by variable specified...will plot each sample separately")
        metadata$Sample <- rownames(metadata)
	group_by <- "Sample"
    }else{
        # check that the group_by variable is in the metadata
        if (!(group_by %in% colnames(metadata))){
            stop(paste0("variable ", group_by, " not in metadata"))
	}
    }

    # check metadata and abundance matrix have matching samples
    if (!(isTRUE(all.equal(colnames(abundance_matrix), rownames(metadata))))){
        stop("samples between abundance_matrix and metadata do not match - do they need sorting?")
    }

    # add the feature as variable
    mat <- abundance_matrix[features,]
    mat$feature <- rownames(mat)

    # reshape
    mat.m <- reshape2::melt(mat)

    # add metadata
    mat.m$covariate <- metadata[mat.m$variable, group_by]

    # add enough colours
    getcolours <- colorRampPalette(brewer.pal(8,"Set2"))
    colours <- getcolours(length(unique(metadata[,group_by])))

    # plotting
    featureplot <- ggplot2::ggplot(mat.m, aes(x=factor(covariate, levels=gtools::mixedsort(unique(mat.m$covariate))), y=value, colour=covariate)) +
                            geom_boxplot(outlier.alpha=0) +
                            geom_jitter(height=0, width=0.15) +
                            theme_bw() +
                            ylab("Abundance") +
                            scale_colour_manual(values=colours) +
			    xlab("") +
			    theme(axis.text.x=element_text(angle=90)) +
			    facet_wrap(~feature)
    return(featureplot)
}


