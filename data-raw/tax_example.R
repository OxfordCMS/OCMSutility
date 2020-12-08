# make taxonomy table from asv_example

data(asv_example)

library(stringr)
featureID <- str_extract(asv_example$sequence, "ASV[0-9]+")
taxon <- str_extract(asv_example$sequence, "(?<=:).*")
tax_split <- str_split_fixed(taxon, ";", 6)

colnames(tax_split) <- c('Phylum','Class','Order','Family','Genus','Species')

# remove taxon prefix
tax_split <- apply(tax_split, 2, str_extract, "(?<=[a-z]{1}__).*")

# convert "NA" from string to NA
tax_split[tax_split=='NA'] <- NA

# add asv IDs as column "sequence" to match asv_example
tax_example <- cbind(asv_example$sequence, tax_split)
colnames(tax_example)[1] <- 'sequence'

tax_example <- as.data.frame(tax_example)
# write RData file
save(asv_example, tax_example, file="./data/asv_example.RData")
