## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6, fig.height = 4,
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(OCMSutility)
library(ggplot2)
library(dplyr)

## ----ocms_clr-----------------------------------------------------------------

# example counts
count_dataframe <- as.data.frame(matrix(sample(0:1000, size = 10000, replace=TRUE), ncol=50, nrow=200))
clr_transformed <- ocms_clr(count_dataframe = count_dataframe, return_as_dataframe = TRUE)

## ----ocms_clr_object----------------------------------------------------------
clr_transformed <- ocms_clr(count_dataframe = count_dataframe, return_as_dataframe = FALSE)

## ----ocms_relab---------------------------------------------------------------
rel_abundance <- ocms_relab(count_dataframe)

## ----ocms_palette-------------------------------------------------------------

ocms_palette(n=10, palette="Set3", preview=TRUE)


## ----ocms_featurebox----------------------------------------------------------

# get example data
#data(asv_example)

# hardcoded during testing
load("/gfs/devel/nilott/OCMS_Sandbox/R_utility/OCMSutility/data/asv_example.RData")

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

# for plotting purposes we would transform the data e.g. clr
asv_clr <- ocms_clr(asv_counts)

# generate some random metadata for the 295 samples - 5 groups for example
metadata <- data.frame(Group = c(rep("Group 1", 59),
                                 rep("Group 2", 59),
                                 rep("Group 3", 59),
                                 rep("Group 4", 59),
                                 rep("Group 5", 59)),
                                 row.names=colnames(asv_clr))

# produce boxplot of random 4 features as an example grouping by Group variable
features <- sample(rownames(asv_clr), size=4)
ocms_featurebox(abundance_matrix=asv_clr, metadata=metadata, features=features, group_by="Group")

## ----ocms_featurebox_colour---------------------------------------------------
ocms_featurebox(abundance_matrix=asv_clr, metadata=metadata, features=features, group_by="Group") +
  scale_colour_manual(values=ocms_palette(n=5, palette="Set1"))


## ----ocms_dissimilarity-------------------------------------------------------

# get example data
#data(asv_example)

# hardcoded during testing
load("/gfs/devel/nilott/OCMS_Sandbox/R_utility/OCMSutility/data/asv_example.RData")

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

asv_relab <- ocms_relab(asv_counts)

# generate some random metadata for the 295 samples - 5 time points with each individual
# having a data point at each time point
metadata <- data.frame(Timepoint = c(rep("Time 1", 59),
                                     rep("Time 2", 59),
                                     rep("Time 3", 59),
                                     rep("Time 4", 59),
                                     rep("Time 5", 59)),
                       Individual = as.character(c(rep(c(1:59), 5))),
                       row.names=colnames(asv_relab),
                       stringsAsFactors = FALSE)

# remove samples with NA
asv_relab <- asv_relab[,!(is.na(colSums(asv_relab)))]

# make sure they are the same
metadata <- metadata[colnames(asv_relab),]

# ask the question - Are individuals more similar to each other than samples are within timepoints?

# within-individual dissimilarity
within_diss <- ocms_dissimilarity(asv_relab, metadata=metadata, individual_variable = "Individual", method="within")

knitr::kable(head(within_diss))

# between-individual dissimilarity at timpoint 1
metadata_t1 <- metadata[metadata$Timepoint == "Time 1",]
asv_relab_t1 <- asv_relab[,rownames(metadata_t1)]
between_diss <- ocms_dissimilarity(asv_relab_t1, metadata=metadata_t1, method="between")

knitr::kable(head(between_diss))

# we can then combine and plot them
diss <- bind_rows(within_diss, between_diss)
ggplot(diss, aes(x=method, y=dissimilarity)) +
  geom_violin() +
  xlab("Dissimilarity type") +
  ylab("Bray-Curtis dissimilarity") +
  theme_bw()


