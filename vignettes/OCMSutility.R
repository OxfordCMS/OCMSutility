## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6, fig.height = 4,
  echo = TRUE,
  message = TRUE,
  collapse = TRUE,
  result = 'asis',
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(OCMSutility)
library(ggplot2)
library(dplyr)

## -----------------------------------------------------------------------------
data(asv_example)

# 49 ASVs and 296 samples
dim(asv_example)

# 49 ASVs, with taxonomic levels in columns 
dim(tax_example)
colnames(tax_example)

## ----ocms_clr-----------------------------------------------------------------
# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

clr_transformed <- ocms_clr(count_dataframe = asv_counts, return_as_dataframe = TRUE)

# returns data frame with transformed abundance estamtes with imputed zeroes
class(clr_transformed)
dim(clr_transformed)

## ----ocms_clr_object----------------------------------------------------------
clr_transformed <- ocms_clr(count_dataframe = asv_counts, return_as_dataframe = FALSE)

# returns ALDEx2 object
class(clr_transformed)

## ----ocms_relab---------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

rel_abundance <- ocms_relab(asv_counts)

## ----ocms_palette-------------------------------------------------------------
ocms_palette(n=10, palette="Set3", preview=TRUE)

## ----ocms_featurebox----------------------------------------------------------
# get example data
data(asv_example)

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
data(asv_example)

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


## ----ocms_plotPCA-------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

asv_transformed <- ocms_clr(count_dataframe = asv_counts, return_as_dataframe = TRUE)

# generate some random metadata for the 295 samples - 5 time points with each individual
# having a data point at each time point
metadata <- data.frame(Timepoint = c(rep("Time 1", 59),
                                     rep("Time 2", 59),
                                     rep("Time 3", 59),
                                     rep("Time 4", 59),
                                     rep("Time 5", 59)),
                       Individual = as.character(c(rep(c(1:59), 5))),
                       row.names=colnames(asv_transformed),
                       stringsAsFactors = FALSE)
metadata$ID <- rownames(metadata)

pca_result <- prcomp(t(asv_transformed), scale = TRUE)
plot_data <- ocms_plotPCA(pca_result, metadata, colourby='Timepoint')

plot_data$p

# modify default plot
add_meta <- merge(plot_data$pdata, metadata, by = 'row.names' )
col_val <- ocms_palette(5, "Set3")
p <- plot_data$p +
  scale_colour_manual(values = col_val) + # pick own colours
  scale_shape_manual(values=21, guide = FALSE) + # change shape and remove from legend
  geom_text(data = add_meta, aes(x = PC1, y = PC2, label = ID)) # add text label
p

## ----ocms_rarefy--------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

rarefaction <- ocms_rarefy(asv_counts)

# default plot
p <- rarefaction$rare_p
p

# modify default plot -- remove geom_label_repel layer
p$layers[[2]] <- NULL
p

## ----reannotateTax-example1---------------------------------------------------
ex1 <- data.frame(ASV = paste0("ASV", 1:5),
                  Order = "order1",
                  Family = c(paste0("family", c(1,1,2,3)), 'unclassified'),
                  Genus = c("unclassified", 'genus1','unclassified','genus2',
                            "unclassified"),
                  read_count = 10)

knitr::kable(ex1)

## ----reannotateTax-example2---------------------------------------------------
ex2 <- ex1[,c('ASV','Order')]
ex2$Family <- c(paste0("family", c(1,1,2,3)), 'order1_unclassified')
ex2$Genus <- c('family1_unclassified','genus1','family2_unclassified','genus2',
               'order1_unclassified')
ex2$read_count <- 10

knitr::kable(ex2)

## ----ocms_reannotateTax-------------------------------------------------------
# showing the dummy example
old_tax <- ex1[,2:4]
old_tax$Kingdom <- 'kingdom1'
old_tax$Phylum <- 'phylum1'
old_tax$Class <- 'class1'
old_tax$Species <- 'unclassified'

old_tax <- old_tax[, c('Kingdom','Phylum','Class','Order','Family','Genus','Species')]
old_tax[old_tax == 'unclassified'] <- NA
knitr::kable(old_tax)

new_tax <- ocms_reannotateTax(old_tax)
knitr::kable(new_tax)

# try with example data
data(asv_example)

# adding Kingdom column; removing sequence column because don't need asv IDs in this example
old_tax <- tax_example
colnames(old_tax)[1] <- 'Kingdom'
old_tax$Kingdom <- 'Bacteria'
knitr::kable(head(old_tax))

new_tax <- ocms_reannotateTax(old_tax)
knitr::kable(head(new_tax))

