#####################################################
## Sequencing analysis of nanopore 16S sequencing  ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

library(viridis)
library(pals)


#create stacked barplot faceted by day & hour
######################################################################################################################

#Read in combined df generated in 01_dataimport.R
run_long_combined <- read.csv("tables/run_combined_abundances_meta.csv")

#Factor hours so they are in right order on the plot
run_long_combined$hour <- factor(run_long_combined$hour, levels = c("inoculum", "0h", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "24h"))

#Rename bacteria
run_long_combined <- run_long_combined %>%
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

# Create column with only day & hour
run_long_combined$time <- paste0(run_long_combined$day, "_", run_long_combined$hour)
run_long_combined$time <- factor(run_long_combined$time, levels = unique(run_long_combined$time))

#set distinct colours for the bacteria
unique_bacteria <- unique(run_long_combined$bacteria)

# Define a palette of colors with enough colors for all unique bacteria
palette <- viridis(length(unique_bacteria))
palette2 <- cols25(length(unique_bacteria))

# Create a named vector to map each bacteria to a color
color_mapping <- setNames(palette2, unique_bacteria)

run_24h <- run_long_combined %>% 
  filter(hour == "24h")

run_days_with24h <- run_long_combined %>%
  filter(!(hour %in% c("0h", "inoculum")),
         day %in% c("d0", "d1", "d6"))

run_days_without24h <- run_long_combined %>%
  filter(!(hour %in% c("24h", "0h", "inoculum")))

inoculum <- run_long_combined %>% 
  filter(hour == "0h") %>%
  mutate(hour = recode(hour, "0h" = "inoculum"))

# Plot inoculum (aka 0h in the replicates)
ggplot(inoculum, aes(x = factor(replicate), y = abundance, fill = bacteria)) +
  geom_bar(stat = "identity") +
  labs(x = "Replicate", y = "Abundance", fill = "Bacteria", title = "Inoculum") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_mapping) +
  theme_pubclean() +
  theme(legend.position = "right")
ggsave("plots/relative_abundance_inoculum_DB044.png")

# Plot changes over d0, d1, d6
ggplot(run_days_without24h, aes(x = factor(replicate), y = abundance, fill = bacteria)) +
  geom_bar(stat = "identity") +
  labs(x = "Replicate", y = "Abundance", fill = "Bacteria", title = "Dynamic changes over feeding cycle") +
  facet_grid(day ~ hour, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_mapping) +
  theme_pubclean() +
  theme(legend.position = "right")
ggsave("plots/relative_abundance_DB044.png")

# Plot 24h endpoints
ggplot(run_24h, aes(x = factor(replicate), y = abundance, fill = bacteria)) +
  geom_bar(stat = "identity") +
  labs(x = "Replicate", y = "Abundance", fill = "Bacteria", title = "Stability of 24h timepoint after each day") +
  facet_wrap(~ day, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_mapping) +
  theme_pubclean() +
  theme(legend.position = "right")
ggsave("plots/relative_abundance_DB044_24h.png")



# Plot lineplots to visualize changes
ggplot(run_days_with24h, aes(x = hour, y = abundance, color = bacteria, group = bacteria)) +
  geom_point() +
  geom_smooth(aes(fill = bacteria), se = TRUE, method = "loess") +
  labs(title = "Abundance of Bacteria Over Time",
       x = "Time",
       y = "Abundance") +
  theme_pubclean() +
  facet_wrap(~ day, scales = "free_x") +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/changes_during_day_lineplot.png") 

ggplot(run_24h, aes(x = time, y = abundance, group = bacteria)) +
  geom_point() +
  geom_smooth(aes(fill = bacteria, color = bacteria), se = TRUE, method = "loess") +
  labs(title = "Stability of community at 24h timepoints",
       x = "Time",
       y = "Abundance") +
  theme_pubclean() +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/changes_24hlineplot.png") 

