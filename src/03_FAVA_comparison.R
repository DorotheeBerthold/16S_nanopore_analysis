#####################################################
## Sequencing analysis of nanopore 16S sequencing  ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

library(FAVA)
library(patchwork)

# Compare conditions with in vivo using FAVA
######################################################################################################################

#Metadata needs to be on left-handside of matrix
# All relative abundances on right-handside of matrix
#each row = 1 sample, summing up to 1

#in vivo##############################################################################################################
vivo2023 <- read.csv("tables/Weiss2023_invivo.csv")

#Extract origin
vivo2023$origin <- str_extract(vivo2023$relab_16Scorr, "TV\\d+_D\\d+_([A-Za-z]+)")
vivo2023$origin <- str_replace(vivo2023$origin, "TV\\d+_D\\d+_", "")

vivo2023 <- vivo2023 %>% 
  mutate(origin = recode(origin,
                         "Je" = "Jejunum",
                         "Il" = "Ileum",
                         "Ce" = "Cecum",
                         "Co" = "Colon",
                         "F" = "Feces"))

#Extract community
vivo2023$community <- str_extract(vivo2023$relab_16Scorr, "(OMM[\\w-]+)$")

# Extract replicate
vivo2023$replicate <- str_extract(vivo2023$relab_16Scorr, "(?<=_)[0-9]+(?=_)")

vivo2023 <- vivo2023[,-1]

vivo2023$origin <- factor(vivo2023$origin, levels = c("Jejunum", "Ileum", "Cecum", "Colon", "Feces"))

vivo2023 <- vivo2023 %>% 
  filter(origin == "Cecum",
         community == "OMM12")

combined_Weiss <- vivo2023

# Make sure all rowSums are 1 (by summing up & dividing each column by rowsums)
combined_Weiss$rowsums <- rowSums(combined_Weiss[, 1:12])
combined_Weiss[, 1:12] <- combined_Weiss[, 1:12] / combined_Weiss$rowsums
rowSums(combined_Weiss[, 1:12])

# Relocate metadata columns so they are on left-handside of matrix
combined_Weiss <- combined_Weiss %>% 
  select(-community, -rowsums, -origin) %>% 
  relocate(replicate) 

#Column names are hardcoded for the plot
ordered_bacteria_short <- c("YL44", "YL2", "I48", "YL58", "YL32", "I46", 
                            "KB1", "YL31", "I49", "YL27", "YL45", "KB18")

# Define a palette of colors with enough colors for all unique bacteria
palette2<- paletteer_d("MetBrewer::Signac", length(ordered_bacteria_short))

# Create a named vector to map each bacteria to a color
color_mapping_short <- setNames(palette2, ordered_bacteria_short)


# Calculate FAVA for cecum
fava_result <- fava(relab_matrix = combined_Weiss, K = 12)


# Compare cecum average to each condition of screen
##############################################################################################################
avg_cecum <- colMeans(combined_Weiss[2:13], na.rm = TRUE)
avg_cecum <- data.frame(t(avg_cecum))
avg_cecum$group <- "cecum"
avg_cecum <- avg_cecum %>% 
  relocate(group)

# Prepare dataframe accordingly: select desired endpoint
endpoint_fava <- run_combined %>% 
  filter(day == "d0" & time == "24h") %>% 
  filter(!media %in% c("carrier", "carrier_2"))
endpoint_fava <- endpoint_fava[,c(1:11,17)]

endpoint_fava <- endpoint_fava %>% 
  relocate(group)

#######
# Doublecheck colnames before and adjust following vector accordingly!!
new_colnames <- c("group", "YL44", "I48", "YL2", "YL58", "YL45", "KB1", "I46", "YL31", "I49", "YL27", "YL32")

colnames(endpoint_fava) <- new_colnames

# bring in same order

avg_cecum <- avg_cecum[, new_colnames]

colnames(avg_cecum) == colnames(endpoint_fava)

# Join the two df
cecum_endpoint <- rbind(endpoint_fava, avg_cecum)

# Make sure all rowSums are 1 (by summing up & dividing each column by rowsums)
cecum_endpoint[is.na(cecum_endpoint)] <- 0
cecum_endpoint$rowsums <- rowSums(cecum_endpoint[, 2:12])
cecum_endpoint[, 2:12] <- cecum_endpoint[, 2:12] / cecum_endpoint$rowsums
rowSums(cecum_endpoint[, 2:12])

cecum_endpoint <- cecum_endpoint[,-13]

# Loop over each group combined with cecum:
##############################################################################################################
# Extract cecum row
cecum_row <- cecum_endpoint[cecum_endpoint$group == "cecum", ]

# Initialize an empty data frame to store results
fava_results <- data.frame(group = integer(), fava_value = numeric())


# Loop through groups 1 to 7
for (grp in 1:7) {
  # Extract group row
  group_row <- cecum_endpoint[cecum_endpoint$group == as.character(grp), ]
  
  # Combine cecum row with group row
  combined_matrix <- rbind(cecum_row, group_row)
  
  # Remove columns with all NA values
  combined_matrix <- combined_matrix[, colSums(!is.na(combined_matrix)) > 0]
  
  # Calculate fava
  fava_value <- fava(relab_matrix = combined_matrix, K = 12)
  
  # Store group and fava value in results
  fava_results <- rbind(fava_results, data.frame(group = grp, fava_value = fava_value))
}

fava_results$media <- c("APF", "APF-xy-in", "APF", "APF-xy-in", "APF", "APF-xy-in", "APF-xy-in")

# Plot fava_results
ggplot(fava_results, aes(x = factor(group), y = fava_value, fill = media)) +
  geom_bar(stat = "identity", colour = "black", alpha = 0.7) + 
  geom_hline(yintercept = 0.0645, color = "red", linetype = "dashed") + # Bar plot
  theme_pubclean() +                                 # Minimal theme
  labs(
    title = "FAVA 24h in vitro compared to average cecum",
    x = "Group",
    y = "FAVA"
  )
ggsave("plots/FAVA_24h_cecum.png")
