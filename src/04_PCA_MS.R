#####################################################
## Metabolomic analysis of OMM12 timecourse.       ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################


# Import of cleaned skyline results
######################################################################################################################
set.seed(5)
ms_results <- read.csv("results/molecule_quantification_filtered_DB044.csv")

ms_results <- ms_results[,c(1:2, 5:6)]
colnames(ms_results) [1:4] <- c("file_name", "Molecule", "Area", "RT")

#remove all qc, water, lead & tail samples
ms_results <- ms_results[!grepl("water|qc|lead|tail|spacer", ms_results$file_name, ignore.case = TRUE), ]

# split up file_name into metadata columns
# First create string column containing all the info
ms_results <- ms_results %>%
  mutate(string = sub(".*__S\\d+\\s(.*)\\.d", "\\1", file_name))

# Second: split day, hour & replicate

ms_results <- ms_results %>%
  separate(string, into = c("day", "hour", "replicate"), sep = "_", remove = F)

ms_results <- ms_results %>% 
  mutate(hour = case_when(
    hour == "48h" ~ "24h",
    hour == "72h" ~ "24h",
    hour == "96h" ~ "24h",
    hour == "120h" ~ "24h",
    hour == "144h" ~ "24h",
    hour == "168h" ~ "24h",
    TRUE ~ hour),
    day = case_when(
      day == "media" ~ "blank",
      TRUE ~ day
    )
  )
write.csv(ms_results, "results/ms_results.csv")


ms_results_summary <- ms_results %>%
  group_by(string, Molecule) %>%
  summarise(Area = mean(Area), .groups = 'drop')

ms_results_matrix <- ms_results_summary %>%
  select(Molecule, string, Area) %>%
  pivot_wider(names_from = Molecule, values_from = Area) %>%
  column_to_rownames(var = "string") %>%
  as.matrix()

pca_data <- ms_results_matrix

ms_results_matrix <- as.data.frame(ms_results_matrix)
ms_results_matrix$string <- rownames(ms_results_matrix)

ms_results_matrix <- ms_results_matrix %>%
  separate(string, into = c("day", "hour", "replicate"), sep = "_", remove = F) %>% 
  relocate(string, day, hour, replicate)

ms_results_matrix <- ms_results_matrix %>% 
  mutate(hour = case_when(
    hour == "48h" ~ "24h",
    hour == "72h" ~ "24h",
    hour == "96h" ~ "24h",
    hour == "120h" ~ "24h",
    hour == "144h" ~ "24h",
    hour == "168h" ~ "24h",
    TRUE ~ hour),
    day = case_when(
      day == "media" ~ "blank",
      TRUE ~ day
    )
)

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)
screeplot(pca_result)

# Extracting scores
pca_scores <- as.data.frame(pca_result$x)

# Define custom shapes (you can add more shapes if needed)
custom_shapes <- c(15, 16, 1, 4, 3, 17, 18, 25)

# Add metadata
pca_scores$day <- as.factor(ms_results_matrix$day)
pca_scores$hour <- as.factor(ms_results_matrix$hour)
pca_scores$replicate <- as.factor(ms_results_matrix$replicate)
pca_scores$day <- factor(pca_scores$day, levels = c("d0", "d1", "d2", "d3", "d4", "d5", "d6", "blank"))


autoplot(pca_result, data = pca_scores, colour = "hour", shape = "day",  size = 6, loadings = F) + 
  theme_pubr() +
  scale_shape_manual(values = custom_shapes)
ggsave("plots/PCA_metabolites_all.png")


