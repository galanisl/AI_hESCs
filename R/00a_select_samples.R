
# This script selects the Epi and PE samples from the Petropoulos et al. study
# using sample annotations from Stirparo et al. (2018) Development

library(dplyr)

# Read the sample meta-data
smdata <- readr::read_tsv("PRJEB11202.txt")

# Read the sample labels
slabels <- readxl::read_excel("Stirparo_TableS4.xlsx")

# Focus on EPI and PE samples from Petropoulos et al. 
slabels <- slabels[grep("Petropoulos", slabels$Study), ] %>% 
  filter(`Revised lineage (this study)` %in% 
           c("epiblast", "primitive_endoderm")) %>% 
  select(Cell, `Revised lineage (this study)`) %>% 
  rename(sample_title = Cell, lineage = `Revised lineage (this study)`) %>% 
  mutate(sample_title = stringr::str_replace_all(sample_title, "_", ".")) %>% 
  arrange(lineage)

# Get the FTP link to the FASTQ files for the samples of interest
sdata <- left_join(slabels, smdata, by = "sample_title")

readr::write_tsv(sdata, "sample_data.tsv")
