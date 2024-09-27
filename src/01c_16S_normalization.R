#####################################################
## Sequencing analysis of nanopore 16S sequencing  ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

## IMPORTANT: Run 01b instead of 01 for counts!!

#Normalize counts with 16S copies per bug
######################################################################################################################

copies_16S <- read.csv("tables/16S_copies.csv")
counts <- read.csv("tables/run_combined_counts_meta.csv")
colnames(copies_16S) [1:2] <- c("bacteria", "16S_copies")

#Rename bacteria
counts <- counts %>%
  mutate(bacteria = recode(bacteria,
                           "Enterococcus faecalis" = "E. faecalis KB1",
                           "Lactobacillus reuteri" = "L. reuteri I49",
                           "Bacteroides caecimuris" = "B. caecimuris I48",
                           "Blautia sp. YL58" = "B. coccoides YL58",
                           "Muribaculum intestinale" = "M. intestinale YL27",
                           "Lachnoclostridium sp. YL32" = "C. clostridium YL32",
                           "Burkholderiales bacterium YL45" = "T. muris YL45",
                           "Erysipelotrichaceae bacterium I46" = "C. innocuum I46",
                           "Akkermansia muciniphila" = "A. muciniphila YL44",
                           "Flavonifractor plautii" = "F. plautii YL31"
  ))

copies_16S <- copies_16S %>%
  mutate(bacteria = recode(bacteria,
                           "KB1" = "E. faecalis KB1",
                           "I49" = "L. reuteri I49",
                           "I48" = "B. caecimuris I48",
                           "YL58" = "B. coccoides YL58",
                           "YL27" = "M. intestinale YL27",
                           "YL32" = "C. clostridium YL32",
                           "YL45" = "T. muris YL45",
                           "I46" = "C. innocuum I46",
                           "YL44" = "A. muciniphila YL44",
                           "YL31" = "F. plautii YL31"
  ))

# Pivot counts wide, so that bacteria are rows and samples are columns
counts <- counts %>%
  select(-X, -day, -hour, -replicate)

counts_wide <- counts %>%
  pivot_wider(names_from = unique_identification, values_from = counts)

counts_wide <- inner_join(copies_16S, counts_wide, by = "bacteria")

# Normalize each column by 16S copies
counts_wide_normalized <- counts_wide %>%
  mutate(across(3:90, ~ . / `16S_copies`))

# Add row with colsums
col_sums <- colSums(counts_wide_normalized[3:90], na.rm = TRUE)

# Create a new row with column sums
new_row <- c("Total_counts", NA, col_sums)

# Add the new row to the data frame
counts_wide_normalized <- rbind(counts_wide_normalized, new_row)

write.csv(counts_wide_normalized, "tables/normalised_counts.csv")

# Calculate relative abundances
# For each column, divide the counts per bacteria by the total_counts & *100 --> done in excel

relative_abundances <- read.csv("tables/relative_abundances.csv")

# Replace all instances of "#VALUE!" with NA
relative_abundances[relative_abundances == "#VALUE!"] <- NA
species <- relative_abundances$bacteria

# Convert dataframe to numeric (if necessary)
relative_abundances <- apply(relative_abundances, 2, function(x) as.numeric(as.character(x)))
relative_abundances <- relative_abundances[,-1]
rownames(relative_abundances) <- species
relative_abundances <- as.data.frame(relative_abundances)

# Convert rownames to a column and reset index
relative_abundances <- tibble::rownames_to_column(relative_abundances, var = "bacteria")

# Convert to long format
abundance_long <- pivot_longer(relative_abundances, 
                        cols = -bacteria, 
                        names_to = "unique_identification", 
                        values_to = "abundance")

abundance_long <- abundance_long %>% 
  filter(!is.na(abundance))


abundance_long <- abundance_long %>%
  separate(unique_identification, into = c("day", "hour", "replicate"), sep = "_", remove = F) %>% 
  relocate(unique_identification, day, hour, replicate)

write.csv(abundance_long, "tables/relative_abundances_long.csv")
