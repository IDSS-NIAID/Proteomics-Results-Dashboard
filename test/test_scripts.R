library(ProtResDash)

# move default data to sqlite file
extdata <- system.file("extdata", package = 'ProtResDash')

import_raw ( file.path(extdata, "peptides.txt"), 
             file.path(extdata, "proteinGroups.txt"), 
             "data/protein_peptidedb.sqlite")
