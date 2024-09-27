#####################################################
## Sequencing analysis of nanopore 16S sequencing  ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

{
  library(readxl)
  library(ggplot2)
  library(tidyverse)
  library(ggfortify)
  library(reshape2)
  library(remotes)
  library(DEGreport)
  library(pheatmap)
  library(vsn)
  library(ggforce)
  library(ggpubr)
  library(emmeans)
  library(ggforce)
  library(ggrepel)
  library(RColorBrewer)
  library(dendextend)
  library(openxlsx)
}

#Dataimport and metadata extraction
######################################################################################################################
meta <- read.csv("tables/DB044_meta.csv")
meta <- meta %>%
  mutate(barcode = str_c("barcode", str_pad(barcode, width = 2, pad = "0")))

#Define run number
run_nr <- "3"

#Read in run data
if(run_nr == "1"){
  run <- read.csv("tables/combined_estimated_counts_run1.csv")
}else if(run_nr == "2"){
  run <- read.csv("tables/combined_estimated_counts_run2.csv")
}else if(run_nr =="3"){
  run <- read.csv("tables/combined_estimated_counts_run3.csv")
}else if(run_nr == "4"){
  run <- read.csv("tables/combined_estimated_counts_run4.csv")
} else{
  print("Not applicable")}


run <- run %>% 
  filter(species != "")
rownames(run) <- run$species

#Check which rows have undefined bacteria to remove

if(run_nr == "1"){
  run <- run[-c(3,11),-c(1:12)]
}else if(run_nr == "2"){
  run <- run[-11,-c(1:12)]
}else if(run_nr =="3"){
  run <- run[-c(2,10),-c(1:12)]
}else if(run_nr == "4"){
  run <- run[-9,-c(1:12)]
} else{
  print("Not applicable")
}

run_t <- as.data.frame(t(run))


#select metadata for run 1
meta_run <- meta %>% 
  filter(PCR.run == run_nr)

#match metadata with run based on barcode column
run_t$barcode <- rownames(run_t)
run_combined <- left_join(run_t, meta_run, by = "barcode")
run_combined <- run_combined %>% 
  select(-c(barcode, sample_nr, PCR.run))

#Make long df
run_long <- reshape2::melt(run_combined, id.vars = c("day", "hour", "replicate", "unique_identification"))
colnames(run_long)[5:6] <- c("bacteria", "counts")

write_csv(run_long, file = paste0("tables/run", run_nr, "_combined_counts_meta.csv"))

run1 <- read.csv("tables/run1_combined_counts_meta.csv")
run2 <- read.csv("tables/run2_combined_counts_meta.csv")
run3 <- read.csv("tables/run3_combined_counts_meta.csv")
run4 <- read.csv("tables/run4_combined_counts_meta.csv")


run_long_combined <- rbind(run1, run2, run3, run4)
write.csv(run_long_combined, "tables/run_combined_counts_meta.csv")

