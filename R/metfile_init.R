#' function to initiate a metadata file from dada2 db file
#'
#' metfile_init

#' @param db_file rsqlite database file
#' @param out_dir output directory
#' @param ref_table = name of table in database from which sampleID are generated.
#'           default NULL. when NULL assumes table is merged_abundance_id
#'           which is the count table produced from the ocms_16s dada2 pipeline
#' @param id_orient = indicates orientation of sample IDs. default 'col' indicates
#'           samples are in columns. 'row' indicates samples are in rows and
#'           sample IDs are the rownames
#' @param dummy = default NULL. when string supplied, adds a dummy column of NAs with the string as the column name.
#' @param export = boolean. default TRUE. writes metadata file as tsv file.
#' @returns dataframe of metadata


metfile_init <- function(db_file, out_dir, ref_table=NULL, id_orient = 'col',
                         dummy = NULL, export=TRUE) {

  # load libraries--------------------------------------------------------------
  require(RSQLite)

  # check inputs----------------------------------------------------------------
  if(!file.exists(db_file)) {
    stop("db_file not found")
  }

  if(!file.exists(out_dir)) {
    stop("out_dir not found")
  }

  if(!id_orient %in% c('col','row')) {
    stop("id_orient should be either 'col' or 'row'")
  }

  # connect to database
  con <- dbConnect(SQLite(), db_file)

  # table in database
  table_ls <- dbListTables(con)

  if(is.null(ref_table)) {
    ref_table <- "merged_abundance_id"
  } else {

    if(!ref_table %in% table_ls) {
      stop("ref_table supplied not found in database file")
    }
  }

  # read in data----------------------------------------------------------------
  query <- sprintf("SELECT * FROM %s;", ref_table)

  # read in count data
  ref_df <- dbGetQuery(con, query)

  # disconenct dabase
  dbDisconnect(con)

  # get sample IDs--------------------------------------------------------------
  if(id_orient == 'col') {
    ids <- colnames(ref_df)
  } else {
    ids <- rownames(ref_df)
  }

  # initiate metadata table
  met <- data.frame(sampleID = ids)

  if(!is.null(dummy)) {
    met[[dummy]] <- NA
  }

  if(export) {
    write.table(met, file.path(out_dir, "metadata.tsv"),
                sep = "\t", row.names=FALSE, quote=FALSE)
  }

  return(met)
}