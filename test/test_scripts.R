library(ProtResDash)
library(RSQLite)
library(dplyr)
library(dbplyr)

# move default data to sqlite file
extdata <- system.file("extdata", package = 'ProtResDash')

import_raw ( file.path(extdata, "peptides.txt"), 
             file.path(extdata, "proteinGroups.txt"), 
             "data/protein_peptidedb.sqlite")

# load sqlite database
con <- dbConnect(SQLite(), dbname = "data/protein_peptidedb.sqlite")
peptides <- tbl(con, "peptides")
proteins <- tbl(con, "proteins")

# plot PCA of proteins
# It's a TMT experiment. 
# It looks like LM_control_1 is tagged with 126 reporter (which is most likely
# reporter 1 in the table); LM_control_2: 127N; LM_control_3: 127C; LM_control_4: 128N;
# LM_control_5: 128C; LM_IgM_1: 129N; LM_IgM_2:129C; LM_IgM_3: 130N; LM_IgM_4: 130C;
# LM_IgM_5: 131. These are then fractionated for a total of 12 acquisitions
# (files PI22-27419.raw to PI22-27430.raw).

data <- proteins |>
  select(starts_with('Reporter.intensity.corrected')) |>
  select(1:10) |>
  collect()

library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

##### Distribution of intensity scores #####
pivot_longer(data, everything(), names_to = 'sample', values_to = 'intensity') |>
  mutate(sample = str_replace(sample, 'Reporter.intensity.corrected.', ''),
         intensity = intensity + 1) |>
  
  ggplot(aes(x = sample, y = intensity, group = sample)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal()


##### Missing Values by Sample #####

# count missing values
data |>
  pivot_longer(everything(), names_to = 'sample', values_to = 'intensity') |>
  mutate(sample = str_replace(sample, 'Reporter.intensity.corrected.', ''),
         missing = is.na(intensity)) |>
  
  ggplot(aes(sample, fill = missing)) +
  geom_bar() +
  theme_minimal() +
  scale_fill_brewer(palette = 'Paired')




