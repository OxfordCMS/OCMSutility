#' ocms_featurebox
#' boxplot expression/abundance values of a feature of interest from a matrix
#'
#' @param abundance_matrix dataframe; samples are columns and features are rows
#' @param metadata dataframe; samples are rows and rownames should match those in abundance_matrix
#' @param feature character; feature to plot i.e. a feature in rownames of abundance_matrix
#' @param group_by character; grouping variable for boxplot. Must be present in metadata
#' @import ggplot2
#' @import reshape
#' @export
#' @return grob
#' @examples
#' abundance_matrix <- data.frame(matrix(rnorm(1000), ncol=20, nrow=50))
#' colnames(abundance_matrix) <- LETTERS[1:20]
#' metadata <- data.frame(group = c(rep("group1", 5), rep("group2", 5), rep("group3", 5), rep("group4", 5)))
#' rownames(metadata) <- LETTERS[1:20]
#' ocms_featurebox(abundance_matrix, metadata=metadata, feature="10", group_by="group")
#'
#' # multiple features
#' grobs.list <- lapply(rownames(abundance_matrix)[1:4], function(x) ocms_featurebox(abundance_matrix, metadata=metadata, feature=x, group_by="group"))
#' gridExtra::grid.arrange(grobs=grobs.list, ncol=2)

ocms_featurebox <- function(abundance_matrix, metadata=NULL, feature=NULL, group_by=NULL){

    # must have a feature specified
    if (is.null(feature)){
        stop("must specify a feature to plot")
    }else{
        if (!(feature %in% rownames(abundance_matrix))){
	    stop("feature not present in abundance_matrix - check rownames")
	}
    }

    # check if the group_by variable is NULL
    if (is.null(group_by)){
        metadata$Sample <- rownames(metadata)
	group_by <- "Sample"
    }else{
        # check that the group_by variable is in the metadata
        if (!(group_by %in% colnames(metadata))){
            stop(paste0("variable ", group_by, " not in metadata"))
	}
    }

    # check metadata and abundance matrix have matching samples
    if (!(all.equal(colnames(abundance_matrix), rownames(metadata)))){
        stop("samples between abundance_matrix and metadata do not match - do they need sorting?")
    }

    # add the test_id as variable
    mat <- abundance_matrix[feature,]
    mat$test_id <- rownames(mat)

    # reshape
    mat.m <- reshape::melt(mat)

    # add metadata
    mat.m$covariate <- metadata[mat.m$variable, group_by]

    colours <- rainbow(length(unique(metadata[,group_by])), s=0.7, v=0.6)

    # plotting
    featureplot <- ggplot2::ggplot(mat.m, aes(x=factor(covariate, levels=mixedsort(unique(mat.m$covariate))), y=value, colour=covariate)) +
                            geom_boxplot(outlier.alpha=0) +
                            geom_jitter(height=0, width=0.15) +
                            theme_bw() +
                            ggtitle(feature) +
                            ylab("Abundance") +
                            scale_colour_manual(values=colours) + xlab("") + theme(axis.text.x=element_text(angle=90))
    return(featureplot)
}


