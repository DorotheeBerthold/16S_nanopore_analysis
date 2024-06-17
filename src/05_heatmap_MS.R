#####################################################
## Metabolomic analysis of OMM12 timecourse.       ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#Matrix for heatmap
######################################################################################################################

# Use object ms_results from 04_PCA
ms_results <- read.csv("results/ms_results.csv")
ms_results_mean <- ms_results %>% 
  group_by(Molecule, day, hour) %>% 
  mutate(mean_Area = mean(Area)) %>% 
  select(-replicate)

ms_results_mean$annotation <- paste0(ms_results_mean$day, "_", ms_results_mean$hour)
ms_results_mean2 <- ms_results_mean %>% 
  group_by(Molecule, mean_Area) %>% 
  distinct(annotation)

ms_wide <- reshape2::acast(ms_results_mean2, formula = Molecule ~ annotation,
                           fill = 0, value.var = "mean_Area",
                           fun.aggregate = sum, drop = F)
ms_scaled <- t(scale(t(ms_wide)))


#create annotation df

ann_df <- ms_results_mean %>% 
  group_by(day, hour) %>% 
  distinct(annotation)
ann_df <- as.data.frame(ann_df)
rownames(ann_df) <- ann_df$annotation
ann_df <- ann_df[,-3]
ann_df$day <- factor(ann_df$day, levels = c("d0", "d1", "d2", "d3", "d4", "d5", "d6", "blank"))
ann_df$hour <- factor(ann_df$hour, levels = c("0h", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "24h", "blank"))

#put in same order as scaled matrix
ms_scaled <- ms_scaled[, match(rownames(ann_df), colnames(ms_scaled))]

#create heatmap
library(RColorBrewer)


annotation_colors_list <- list(
  day = setNames(colorRampPalette(c("black", "white"))(length(levels(ann_df$day))), levels(ann_df$day)),
  hour = setNames(brewer.pal(length(levels(ann_df$hour)), "PiYG"), levels(ann_df$hour))
)


# Extract day and hour information from the column names
day_hour <- strsplit(colnames(ms_scaled), "_")
day <- sapply(day_hour, `[`, 1)
hour <- sapply(day_hour, `[`, 2)

# Convert day and hour to numeric for sorting
day_numeric <- as.numeric(sub("d", "", day))
hour_numeric <- as.numeric(sub("h", "", hour))

# Combine day and hour into a data frame
day_hour_df <- data.frame(day = day_numeric, hour = hour_numeric, colname = colnames(ms_scaled))

# Sort the data frame based on day and hour
sorted_day_hour_df <- day_hour_df[order(day_hour_df$day, day_hour_df$hour), ]

# Reorder the columns in ms_scaled
ms_scaled_sorted <- ms_scaled[, sorted_day_hour_df$colname]

# Reorder the columns in ann_df to match the sorted order
ann_df_sorted <- ann_df[sorted_day_hour_df$colname, ]

# Create the heatmap with the sorted columns
heatmap <- pheatmap(ms_scaled_sorted, show_rownames = F, show_colnames = F, cluster_cols = F, annotation_col = ann_df_sorted, annotation_colors = annotation_colors_list, cutree_rows = 4, cutree_cols = 5)
save_pheatmap_pdf(heatmap, filename = "plots/heatmap_metabolites_sorted.pdf", width = 7, height = 10)
