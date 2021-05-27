## code to prepare `dss_example` dataset goes here

# save_dir <- "/gfs/devel/syen/OCMSlooksy/data-raw"
save_dir <- "./data-raw"

# database from copied from /gfs/work/syen/test_dss/2_dada2/dss_db
# metadata copied from /gfs/devel/syen/OCMSlooksy/data-raw/dss_metadata.tsv

# read database
con <- dbConnect(RSQLite::SQLite(), file.path(save_dir, "dss_db"))

# writing database table into list to save as RData file
# extract data tables
table_ls <- RSQLite::dbListTables(con)

dss_example <- list()
for(i in 1:length(table_ls)) {
  query <- sprintf("SELECT * FROM %s", table_ls[i])
  entry <- RSQLite::dbGetQuery(con, query)

  dss_example[[table_ls[i]]] <- entry
}
# close connection
dbDisconnect(con)

dss_met <- read.csv(file.path(save_dir, "dss_metadata.tsv"), sep="\t")
# keep relevant columns
keep <- apply(dss_met, 2, function(x) !all(is.na(x)))

dss_met <- dss_met[colnames(dss_example$merged_abundance_id),which(keep)]

# add metadata to list
dss_example$metadata <- dss_met

save(dss_example, 'dss_example', file = './data/dss_example.RData')
