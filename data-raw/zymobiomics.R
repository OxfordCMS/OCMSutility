# zymoBIOMICS standard sample data
# pulled from OCMS_zymoBIOMICS

#####################################################################
# NCBI data for 16S analysis
#####################################################################
anno_ncbi_16s <- data.frame(species = c("Pseudomonas_aeruginosa",
                                        "Escherichia_coli",
                                        "Salmonella_enterica",
                                        "Lactobacillus_fermentum",
                                        "Enterococcus_faecalis",
                                        "Staphylococcus_aureus",
                                        "Listeria_monocytogenes",
                                        "Bacillus_subtilis"),
                            genus = c("Pseudomonas",
                                      "Escherichia/Shigella",
                                      "Salmonella",
                                      "Lactobacillus",
                                      "Enterococcus",
                                      "Staphylococcus",
                                      "Listeria",
                                      "Bacillus"),
                            family = c("Pseudomonadaceae",
                                       "Enterobacteriaceae",
                                       "Enterobacteriaceae",
                                       "Lactobacillaceae",
                                       "Enterococcaceae",
                                       "Staphylococcaceae",
                                       "Listeriaceae",
                                       "Bacillaceae"),
                            Expected = c(4.2,
                                         10.1,
                                         10.4,
                                         18.4,
                                         9.9,
                                         15.5,
                                         14.1,
                                         17.4),
                            stringsAsFactors = FALSE)
#####################################################################
# GTDB data for 16S analysis
#####################################################################
anno_gtdb_16s <- data.frame(species = c("Pseudomonas_aeruginosa",
                                        "Escherichia_flexneri",
                                        "Salmonella_enterica",
                                        "Lactobacillus_H_fermentum",
                                        "Enterococcus_faecalis",
                                        "Staphylococcus_aureus",
                                        "Listeria_monocytogenes_B",
                                        "Bacillus_marinus"),
                            genus = c("Pseudomonas",
                                      "Escherichia",
                                      "Salmonella",
                                      "Lactobacillus_H",
                                      "Enterococcus",
                                      "Staphylococcus",
                                      "Listeria",
                                      "Bacillus"),
                            family = c("Pseudomonadaceae",
                                       "Enterobacteriaceae",
                                       "Enterobacteriaceae",
                                       "Lactobacillaceae",
                                       "Enterococcaceae",
                                       "Staphylococcaceae",
                                       "Listeriaceae",
                                       "Bacillaceae"),
                            stringsAsFactors = FALSE)
anno_gtdb_16s$Expected <- anno_ncbi_16s$Expected

#####################################################################
# NCBI data for shotgun
#####################################################################
anno_ncbi_shotgun <- anno_ncbi_16s
anno_ncbi_shotgun$Expected <- c(12, 12, 12, 12, 12, 12, 12, 12)
# add fungi
fungi <- data.frame(species = c("Saccharomyces cerevisiae", "Cryptococcus neoformans"),
                    genus = c("Saccharomyces", "Cryptococcus"),
                    family = c("Saccharomycetaceae","Tremellaceae"),
                    Expected = c(2, 2),
                    stringsAsFactors = FALSE)
anno_ncbi_shotgun <- bind_rows(anno_ncbi_shotgun, fungi)

#####################################################################
# GTDB data for shotgun
#####################################################################
# NB - The fungi aren't present in the GTDB so have to
# sort the abundances accordingly
anno_gtdb_shotgun <- anno_gtdb_16s
gtdb_abundance <- anno_ncbi_shotgun$Expected + 0.5
anno_gtdb_shotgun$Expected <- gtdb_abundance[1:8]

zymobiomics <- list("anno_gtdb_16s" = anno_gtdb_16s,
                    "anno_ncbi_16s" = anno_ncbi_16s,
                    "anno_gtdb_shotgun" = anno_gtdb_shotgun,
                    "anno_ncbi_shotgun" = anno_ncbi_shotgun)
#
# # write RData file
# save(anno_gtdb_shotgun, anno_ncbi_shotgun, anno_gtdb_16s, anno_ncbi_16s,
#      file="./data/zymobiomics.RData")

usethis::use_data(zymobiomics, overwrite=TRUE)
