#####################################################
## Sequencing analysis of nanopore 16S sequencing  ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

library(FAVA)
library(patchwork)


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


#in vitro##############################################################################################################
vitro2021 <- read.csv("tables/Weiss2021_invitro.csv")

# Extract community, media, and replicate
vitro2021$community <- ifelse(grepl("OMM-1", vitro2021$sample), "OMM-1", "OMM")
vitro2021$origin <- sub(".*_(M[0-9]+|Inok)_.*", "\\1", vitro2021$sample)
vitro2021$replicate <- sub(".*_(W[0-9]+|Inok).*", "\\1", vitro2021$sample)
vitro2021$experiment <- sub("^(E[0-9]+)_.*", "\\1", vitro2021$sample)

vitro2021 <- vitro2021 %>% 
  filter(origin != "Inok",
         origin != "M3")

vitro2021$replicate <- paste0(vitro2021$experiment, "_", vitro2021$replicate)

vitro2021 <- vitro2021 %>% 
  select(-X, -experiment, -sample)

vitro2021 <- vitro2021 %>% 
  mutate(community = recode(community,
                            "OMM" = "OMM12"))


# Combine in vivo & in vitro df
combined_Weiss <- rbind(vitro2021, vivo2023)

# Make sure all rowSums are 1 (by summing up & dividing each column by rowsums)
combined_Weiss$rowsums <- rowSums(combined_Weiss[, 1:12])
combined_Weiss[, 1:12] <- combined_Weiss[, 1:12] / combined_Weiss$rowsums
rowSums(combined_Weiss[, 1:12])

# Relocate metadata columns so they are on left-handside of matrix
combined_Weiss <- combined_Weiss %>% 
  filter(community == "OMM12") %>% 
  select(-community, -rowsums) %>% 
  relocate(origin, replicate) 

#Column names are hardcoded for the plot
ordered_bacteria_short <- c("YL44", "YL2", "I48", "YL58", "YL32", "I46", 
                      "KB1", "YL31", "I49", "YL27", "YL45", "KB18")

# Define a palette of colors with enough colors for all unique bacteria
palette2<- paletteer_d("MetBrewer::Signac", length(ordered_bacteria_short))

# Create a named vector to map each bacteria to a color
color_mapping_short <- setNames(palette2, ordered_bacteria_short)


# Check that rowSums are 1
rowSums(combined_Weiss[, 3:14], na.rm = TRUE)


# Plot relative abundances
###################################################################################################################
plot_relabund(combined_Weiss,
              group = "origin",
              arrange = "vertical",
              K = 12) + # specify amount of numerical columns
  ggplot2::scale_color_manual(values = color_mapping_short) +
  ggplot2::scale_fill_manual(values = color_mapping_short)

ggsave("plots/FAVA_relabundance_Weiss_combined.png")


# Compute FAVA (unweighted)
###################################################################################################################
# First, split into cecum & each media separately to compare to cecum


# Filter to only include "M" values and Cecum
unique_origins <- combined_Weiss %>%
  filter(grepl("^M", origin)) %>%
  pull(origin) %>%
  unique()


# Function to calculate fava for each M
compute_fava <- function(m_origin) {
  Weiss_fava <- combined_Weiss %>%
    filter(origin %in% c(m_origin, "Cecum"))
  
  # Calculate the fava result
  fava_result <- fava(relab_matrix = Weiss_fava, K = 12)
  
  # Return as a tibble with the media (origin) and fava result
  tibble(media = m_origin, fava = fava_result)
}

# Iterate over all M origins and compute fava, combine into a table
fava_results <- map_df(unique_origins, compute_fava)

fava_results$media <- factor(fava_results$media, levels = c("M1", "M2", "M4", "M8", "M10", "M11"))

#{fava_results <- fava_results %>% 
  mutate(media2 = recode(media,
                        "M1" = "AF + glucose",
                        "M2" = "AF - glucose",
                        "M4" = "AF - glucose + mucin",
                        "M8" = "AF - glucose + mucin + C5/C6 sugars",
                        "M10" = "AF - glucose + mucin + C5/C6 sugars + xylan/inulin",
                        "M11" = "AF - glucose + starch"))
#}

ggplot(fava_results, aes(media, fava)) +
  geom_point(aes(fill = media), color = "black", shape = 21, size = 8, stroke = 1) +
  labs(title = "FAVA comparison to Cecum", y = "FAVA", x = "Media") +
  theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x = element_blank())
ggsave("plots/fava_cecum.png")

Weiss_fava <- combined_Weiss %>%
  filter(origin %in% c(unique_origins, "Cecum"))

plot_relabund(Weiss_fava,
              group = "origin",
              arrange = "vertical",
              K = 12) + # specify amount of numerical columns
  ggplot2::scale_color_manual(values = color_mapping_short) +
  ggplot2::scale_fill_manual(values = color_mapping_short)
ggsave("plots/relabund_cecum_media.png")


fava(relab_matrix = Weiss_fava,
     group = "origin",
     K = 12)

rowSums(Weiss_fava[, 3:14], na.rm = TRUE)

# Bootstrapping to get an idea of the error
###################################################################################################################
library(boot)

combine_boot <- function(m_origin) {
  Weiss_boot <- Weiss_fava %>%
    filter(origin %in% c(m_origin, "Cecum"))
  
  # Add column based on each media-cecum combination
  Weiss_boot$media <- m_origin
  
  return(Weiss_boot)  # Ensure the function returns a data frame
}

# Iterate over all M origins and compute fava, combine into a table
boot_results <- map_df(unique_origins, combine_boot)

boot_results[, 3:14] <- boot_results[, 3:14] / rowSums(boot_results[, 3:14])

boot_results <- boot_results %>% 
  relocate(media)


# Modify the function to return data instead of the plot
bootstrap_and_return_data <- function(m_origin, data) {
  # Filter data and perform bootstrapping as before
  boot_results <- data %>% 
    filter(media == m_origin)
  
  statistic_function <- function(data, indices) {
    sample_data <- data[indices, ]
    result <- fava(relab_matrix = sample_data, K = 12)
    return(result)
  }
  
  results <- boot(data = boot_results, statistic = statistic_function, R = 1000)
  results_ci <- boot.ci(results, 0.95, type = "norm")
  
  # Extract values into a data frame
  plot_data <- data.frame(
    t0 = results_ci$t0,
    ci_lower = results_ci$normal[2],
    ci_upper = results_ci$normal[3],
    media = m_origin
  )
  
  return(plot_data)
}

# Generate a data frame for all media origins
combined_data <- map_df(unique_origins, ~bootstrap_and_return_data(.x, boot_results))

combined_data$media <- factor(combined_data$media, levels = c("M1", "M2", "M4", "M8", "M10", "M11"))

# Create a combined plot for all media
ggplot(combined_data, aes(x = media, y = t0, colour = media)) +
  geom_point(aes(fill = media), size = 8, shape = 21, color = "black", stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(x = "Media", y = "FAVA", title = "Bootstrapped not stratified for Cecum & all media") +
  theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x = element_blank())
ggsave("plots/bootstrapped_FAVA_media.png")

# Stratify by media or Cecum for bootstrapping
###################################################################################################################

# Define the statistic function for bootstrapping, with indices --> here you can't use boot() as the function itself resamples

# Define the function for bootstrapping
bootstrap_fava <- function(data, m_origin, n_bootstraps = 1000) {
  # Initialize a vector to store the results
  boot_results <- numeric(n_bootstraps)
  
  # Number of samples in the data
  n <- nrow(data)
  
  # Perform bootstrapping
  for (i in 1:n_bootstraps) {
    # Sample with replacement
    sampled_indices <- sample(1:n, n, replace = TRUE)
    sample_data <- data[sampled_indices, ]
    
    # Split the resampled data into "Cecum" and the media group
    media_data <- sample_data %>% filter(origin == m_origin)
    cecum_data <- sample_data %>% filter(origin == "Cecum")
    
    # Determine the number of samples to draw from each group
    n_samples <- min(nrow(media_data), nrow(cecum_data))
    
    # Sample equal numbers from both "Cecum" and media
    media_sample <- media_data[sample(1:nrow(media_data), n_samples, replace = TRUE), ]
    cecum_sample <- cecum_data[sample(1:nrow(cecum_data), n_samples, replace = TRUE), ]

    # Combine the sampled data
    combined_sample <- rbind(media_sample, cecum_sample)
    
    # Calculate the statistic (using your fava function)
    result <- fava(relab_matrix = combined_sample, K = 12)
    
    # Store the result
    boot_results[i] <- result
  }
  
  # Calculate confidence intervals
  ci <- quantile(boot_results, probs = c(0.025, 0.975))
  
  return(list(
    mean = mean(boot_results),
    ci_lower = ci[1],
    ci_upper = ci[2],
    media = m_origin
  ))
}

# Loop through each media type and apply the bootstrap function
unique_origins <- unique(boot_results$media) # Extract unique media types
bootstrap_results <- lapply(unique_origins, function(m_origin) {
  # Filter data for current media type
  filtered_data <- boot_results %>% filter(media == m_origin)
  
  # Perform bootstrapping
  bootstrap_fava(filtered_data, m_origin)
})

# Convert results to a data frame for easier plotting
results_df <- do.call(rbind, lapply(bootstrap_results, function(result) {
  data.frame(
    media = result$media,
    mean = result$mean,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper
  )
}))

results_df$media <- factor(results_df$media, levels = c("M1", "M2", "M4", "M8", "M10", "M11"))

ggplot(results_df, aes(x = media, y = mean, color = media)) +
  geom_point(aes(fill = media), size = 8, shape = 21, color = "black", stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(x = "Media", y = "FAVA", title = "Bootstrapped stratified for Cecum & all media") +
  theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x = element_blank()) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2)
  
ggsave("plots/bootstrapped_FAVA_media_stratified.png")
