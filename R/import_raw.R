#' import_raw
#' Import raw data from MaxQuant into a SQLite database
#' 
#' @param peptides_file A character string of the path to the peptides file
#' @param proteins_file A character string of the path to the proteins file
#' @param db_file A character string of the path to the SQLite database file
#' 
#' @return NULL
#' @export
#' @importFrom DBI dbConnect
#' @importFrom DBI dbDisconnect
#' @importFrom DBI dbWriteTable
#' @importFrom RSQLite SQLite
#' @importFrom utils read.table
import_raw<- function(peptides_file, proteins_file, db_file) {
  peptides <- read.table (peptides_file, header = TRUE, sep = "\t")
  proteins <- read.table(proteins_file, header = TRUE, sep = "\t")
  
  #new SQLites database
  con <- dbConnect(SQLite(), dbname = db_file)
  
  #write the data to the database 
  dbWriteTable(con, 'peptides', peptides, overwrite = TRUE)
  dbWriteTable(con, 'proteins', peptides, overwrite = TRUE)
  
  #close the connection
  dbDisconnect(con)
}

