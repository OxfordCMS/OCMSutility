#' dissimilarity
#'
#' Calculate either within-individual or between-individual Bray-Curtis dissimilarity.
#'
#' @param abundance_matrix dataframe; samples are columns and features are rows
#' @param metadata dataframe; samples are rows and rownames should match those in abundance_matrix
#' @param individual_variable character; name of the variable in metadata that specifies the individual
#' @param method character; choose 'within' or 'between' depending on whether there are paired samples available
#' @import phyloseq
#' @import dplyr
#' @export
#' @return data frame
#' @examples
#' abundance_matrix <- data.frame(matrix(rnorm(1000, mean=100, sd=1), ncol=8, nrow=50))
#' colnames(abundance_matrix) <- LETTERS[1:8]
#' abundance_matrix <- ocms_relab(abundance_matrix)
#' metadata <- data.frame(individual = c(rep("1", 2), rep("2", 2), rep("3", 2), rep("4", 2)))
#' metadata$other_variable <- rnorm(8)
#' rownames(metadata) <- LETTERS[1:8]
#' dissimilarity(abundance_matrix, metadata=metadata, individual_variable="individual", method="within")


dissimilarity <- function(abundance_matrix=NULL, metadata=NULL, individual_variable=NULL, method="within"){

    if (!(method %in% c("within", "between"))){
        stop("method must be one of 'within' or 'between'")
    }

    if (is.null(abundance_matrix)){
        stop("Must provide relative abundance data frame")
    }
    if (is.null(metadata)){
        stop("Must provide metadata data frame")
    }


    # check metadata and abundance matrix have matching samples
    if (!(isTRUE(all.equal(colnames(abundance_matrix), rownames(metadata))))){
        stop("samples between abundance_matrix and metadata do not match \n do they need sorting or making rownames the sample names in metadata")
    }

    sample.ids <- rownames(metadata)
    features <- abundance_matrix[,sample.ids]

    if (method == "within"){

        # make sure individual variable is specified
        if (is.null(individual_variable)){
            stop("Must provide the variable that specifies the individual column in the metadata")
        }
        # make sure there are multiple samples per individual
        if (length(metadata[,individual_variable]) == length(unique(metadata[,individual_variable]))){
            stop("could not find multiple samples for each individual - cannot compute within dissimilarity")
        }

        # get the patient ids to iterate over
        individual.ids <- unique(metadata[,individual_variable])

        # iterate over idividuals and calculate
        # the bray.curtis distance
        result.p <- list()
        for (i in 1:length(individual.ids)){
            metadata.p <- metadata[metadata[,individual_variable] == individual.ids[i],]

            # don't compute on single samples
            if (nrow(metadata.p) == 1){
	        cat(paste0("won't calculate for ", individual.ids[i], " as there is only one sample\n"))
	        next}

            features.p <- data.frame(features[,rownames(metadata.p)])
            #rownames(features.p) <- rownames(features)
            #colnames(features.p) <- rownames(metadata.p)
            features.p <- phyloseq::otu_table(features.p, taxa_are_rows=TRUE)
            metadata.p <- phyloseq::sample_data(metadata.p)
            phyob.p <- phyloseq::merge_phyloseq(features.p, metadata.p)

            diss <- phyloseq::distance(phyob.p, method="bray")
            diss <- as.data.frame(as.matrix(diss))
            res <- data.frame(dissimilarity = diss[2:nrow(diss),1])
            res$method <- "Within-individual"
            result.p[[i]] <- res
        }
        res.all <- bind_rows(result.p)
        res <- res.all
    }

    else if (method=="between"){
        cat("using method=between, make sure there is only one sample per individual")
        features <- phyloseq::otu_table(abundance_matrix, taxa_are_rows=TRUE)
        metadata <- phyloseq::sample_data(metadata)
        phyob <- phyloseq::merge_phyloseq(features, metadata)
        diss <- phyloseq::distance(phyob, method="bray")
        res <- data.frame(dissimilarity = as.vector(diss))
        res$method <- "Between-individual"
    }

    return(res)
}
