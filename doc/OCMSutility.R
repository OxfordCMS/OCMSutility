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
library(knitr)
library(OCMSutility)
library(ggplot2)
library(dplyr)
library(tibble)

## -----------------------------------------------------------------------------
data(asv_example)

# 49 ASVs and 296 samples
dim(asv_example)

# 49 ASVs, with taxonomic levels in columns 
dim(tax_example)
colnames(tax_example)

## ----clr----------------------------------------------------------------------
# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

clr_transformed <- clr(count_dataframe = asv_counts, return_as_dataframe = TRUE)

# returns data frame with transformed abundance estamtes with imputed zeroes
class(clr_transformed)
dim(clr_transformed)

## ----clr_object---------------------------------------------------------------
clr_transformed <- clr(count_dataframe = asv_counts, return_as_dataframe = FALSE)

# returns ALDEx2 object
class(clr_transformed)

## ----relab--------------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

rel_abundance <- relab(asv_counts)

## ----getPalette---------------------------------------------------------------
getPalette(n=10, palette="Set3", preview=TRUE)

## ----featurebox---------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

# for plotting purposes we would transform the data e.g. clr
asv_clr <- clr(asv_counts)

# generate some random metadata for the 295 samples - 5 groups for example
metadata <- data.frame(Group = c(rep("Group 1", 59),
                                 rep("Group 2", 59),
                                 rep("Group 3", 59),
                                 rep("Group 4", 59),
                                 rep("Group 5", 59)),
                                 row.names=colnames(asv_clr))

# produce boxplot of random 4 features as an example grouping by Group variable
features <- sample(rownames(asv_clr), size=4)
featurebox(abundance_matrix=asv_clr, metadata=metadata, features=features, group_by="Group")

## ----featurebox_colour--------------------------------------------------------
featurebox(abundance_matrix=asv_clr, metadata=metadata, features=features, group_by="Group") +
  scale_colour_manual(values=getPalette(n=5, palette="Set1"))


## ----dissimilarity------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

asv_relab <- relab(asv_counts)

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
within_diss <- dissimilarity(asv_relab, metadata=metadata, individual_variable = "Individual", method="within")

knitr::kable(head(within_diss))

# between-individual dissimilarity at timpoint 1
metadata_t1 <- metadata[metadata$Timepoint == "Time 1",]
asv_relab_t1 <- asv_relab[,rownames(metadata_t1)]
between_diss <- dissimilarity(asv_relab_t1, metadata=metadata_t1, method="between")

knitr::kable(head(between_diss))

# we can then combine and plot them
diss <- bind_rows(within_diss, between_diss)
ggplot(diss, aes(x=method, y=dissimilarity)) +
  geom_violin() +
  xlab("Dissimilarity type") +
  ylab("Bray-Curtis dissimilarity") +
  theme_bw()


## ----plotPCA------------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

asv_transformed <- clr(count_dataframe = asv_counts, return_as_dataframe = TRUE)

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
plot_data <- plotPCA(pca_result, metadata, colourby='Timepoint')

plot_data$p

# modify default plot
add_meta <- merge(plot_data$pdata, metadata, by = 'row.names' )
col_val <- getPalette(5, "Set3")
p <- plot_data$p +
  scale_colour_manual(values = col_val) + # pick own colours
  scale_shape_manual(values=21, guide = FALSE) + # change shape and remove from legend
  geom_text(data = add_meta, aes(x = PC1, y = PC2, label = ID)) # add text label
p

## ----rarefaction--------------------------------------------------------------
# get example data
data(asv_example)

# rownames have to be features
asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)

rarefaction <- rarefaction(asv_counts)

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

## ----reannotateTax------------------------------------------------------------
# showing the dummy example
old_tax <- ex1[,2:4]
old_tax$Kingdom <- 'kingdom1'
old_tax$Phylum <- 'phylum1'
old_tax$Class <- 'class1'
old_tax$Species <- 'unclassified'

old_tax <- old_tax[, c('Kingdom','Phylum','Class','Order','Family','Genus','Species')]
old_tax[old_tax == 'unclassified'] <- NA
knitr::kable(old_tax)

new_tax <- reannotateTax(old_tax)
knitr::kable(new_tax)

# try with example data
data(asv_example)

# adding Kingdom column; removing sequence column because don't need asv IDs in this example
old_tax <- tax_example
colnames(old_tax)[1] <- 'Kingdom'
old_tax$Kingdom <- 'Bacteria'
knitr::kable(head(old_tax))

new_tax <- reannotateTax(old_tax)
knitr::kable(head(new_tax))

## ----aggregateCount-----------------------------------------------------------
data(dss_example)
# featureID should be row names
feature_count <- dss_example$merged_abundance_id %>%
   tibble::column_to_rownames('featureID')

# cleanup sample names
colnames(feature_count) <- paste0('id', colnames(feature_count))
# taxonomy table must have columns 'Kingdom','Phylum',
# 'Class','Order','Family','Genus','Species'
# and feature IDs in rownames
feature_tax <- dss_example$merged_taxonomy

# set row order of count and tax tables to be the same
feature_count <- feature_count[feature_tax$featureID,]
aggregated_list <- aggregateCount(feature_count, feature_tax,
                                      aggregate_by = "Family")

summary(aggregated_list)
knitr::kable(head(aggregated_list[['count_df']][,1:5]))
knitr::kable(head(aggregated_list[['tax_df']]))

## ----plotPCoA-----------------------------------------------------------------
data(dss_example)
met_df <- dss_example$metadata

count_df <- dss_example$merged_abundance_id %>%
  column_to_rownames('featureID')
count_df <- count_df[,met_df$sampleID]
relab <- relab(count_df)

iter_var <- c('Genotype','Phenotype')
for(i in iter_var) {
  plotPCoA(relab, met_df, colour = i)
}

## ----plotSunburst-------------------------------------------------------------
data("dss_example")
# set count feature ids as rownames
count_df <- dss_example$merged_abundance_id %>%
  column_to_rownames('featureID')

# clean up some sample names
colnames(count_df) <- paste0('id', colnames(count_df))
tax_df <- dss_example$merged_taxonomy

# aggregate counts
agg_gen <- aggregateCount(count_df[tax_df$featureID,], tax_df, "Genus")
count_genus <- agg_gen$count_df

# reannotate taxonomy
tax_genus <- reannotateTax(agg_gen$tax_df)

relab <- relab(count_genus)
# color specific phyla
plotSunburst(relab = NULL, tax = tax_genus,  
             palettes = c("Proteobacteria" = "Oranges",
                          "Bacteroidetes" = "Greens"))
# color specific phyla taking into account of relative abundance
plotSunburst(relab = relab, tax = tax_genus,  
             palettes = c("Proteobacteria" = "Oranges", "Bacteroidetes" = "Greens"))

# highlight specific genera
plotSunburst(relab = relab, tax = tax_genus, 
             palettes = c("Bacteroidetes" = "Greens",'Firmicutes'='Blues'), 
             highlight = list("Genus" = c("Bacteroides",'Clostridium XlVa')))

## ----convert-platemap, eval=FALSE---------------------------------------------
#  plate_map <- convert_platemap(map_file = "my96wellplate.xlsx",
#                                map_range = 'A1:H12')

## ----convert-platemap-demo, echo = FALSE--------------------------------------
# example output of convert_platemap
data.frame(well_id = paste0(rep(LETTERS[1:8], each = 12), rep(1:12, 8)),
           col = rep(1:12, 8),
           row = rep(LETTERS[1:8], each = 12),
           sample_name = paste0('sample_name', 1:96)) 

## ----truPosRate---------------------------------------------------------------
# this would be better exemplified with actual std data rather than the example samples
data("dss_example")
data(zymobiomics)

# set count feature ids as rownames
count_df <- dss_example$merged_abundance_id %>%
  column_to_rownames('featureID')

# clean up some sample names
colnames(count_df) <- paste0('id', colnames(count_df))
tax_df <- dss_example$merged_taxonomy

# aggregate counts
agg_gen <- aggregateCount(count_df[tax_df$featureID,], tax_df, "Genus")
genus_relab <- relab(agg_gen$count_df)

true_pos_result <- truePosRate(relab=genus_relab,
                                    annotations=zymobiomics$anno_ncbi_16s,
                                    level='genus', cutoff=0.01)

# plot true pos rate
p <- ggplot(true_pos_result,
            aes(x=rank, y=true.pos.rate, colour=label, group=sample)) +
  geom_point() +
  theme_bw() +
  ylab("TP / (TP + FP)") +
  scale_colour_manual(values=c("grey", "purple")) +
  facet_wrap(~sample, scale="free")

p

## ----filterFeature------------------------------------------------------------
data(dss_example)

# put featureID as rownames
tax_df <- dss_example$merged_taxonomy
count_df <- dss_example$merged_abundance_id %>%
  column_to_rownames('featureID')
# set features in count tax to be in same order
count_df <- count_df[tax_df$featureID,]

filtered_ls <- filterFeature(count_df, tax_df, 'percent_sample', 0.001, 2)
summary(filtered_ls)
filtered_count <- filtered_ls$filtered
dim(filtered_count)
kable(head(filtered_count[,1:4]))

## ----metfile_init, eval=FALSE-------------------------------------------------
#  db_file <- "/path/to/db/file"
#  met <- metfile_init(db_file, dummy = "Group")
#  

