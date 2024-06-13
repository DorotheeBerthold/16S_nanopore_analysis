#####################################################
## Metabolomic analysis of OMM12 timecourse.       ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

library(readxl)

# Import of skyline results list
######################################################################################################################

skyline <- read.csv("tables/Molecule_Quantification_DB044_updated.csv")
skyline_molecules <- skyline %>% 
  distinct(Molecule)
# 303 metabolites

# Import negative retention time library
library <- read_xlsx("tables/library_MSMM.xlsx", sheet = 1, col_names = TRUE)
library_molecules <- library %>% 
  distinct(Molecule)
# 287 metabolites
matching_molecules <- inner_join(skyline_molecules, library_molecules, by = join_by(Molecule))
#177 that match

skyline_unique <- anti_join(skyline_molecules, matching_molecules, by = join_by(Molecule)) # 126
library_unique <- anti_join(library_molecules, matching_molecules, by = join_by(Molecule)) #110

skyline_library <- inner_join(skyline, library, by = "Molecule", relationship = "many-to-many")

# Filter out molecules where area is below 1000
# check relative retention times
# Divide RT from analyte by RT from standard library *100
# if higher than 10-15% difference, investigate peak

skyline_RRT <- skyline_library %>% 
  mutate(RRT = Best.Retention.Time / RT) %>%
  group_by(Molecule) %>% 
  mutate(mean_RRT = mean(RRT),
         sd_rt = sd(Best.Retention.Time),
         mean_area = mean(Total.Area)) %>%
  select(Molecule, RT, mean_RRT, sd_rt, mean_area)

skyline_RRT <- skyline_RRT %>% 
  group_by(mean_RRT, sd_rt, mean_area, RT) %>% 
  distinct(Molecule)

area_filter <- skyline_RRT %>% 
  filter(mean_area<1000) %>% 
  select(Molecule)

skyline_RRT <- skyline_RRT %>% 
  filter(mean_area > 1000)


#list of metabolites to check on skyline

skyline_check <- skyline_RRT %>% 
  group_by(Molecule) %>% 
  filter(mean_RRT > 10 | sd_rt > 1)

write_csv(skyline_check, "results/skyline_check_B044.csv")
write.xlsx(skyline_check, "results/skyline_check_B044.xlsx")

# Plot RRT against sd of RRT --> higher RRT has also higher sd

ggplot(skyline_check, aes(sd_rt, mean_RRT)) +
  geom_point()

# Remove molecules that were checked in skyline but still have negative peak
molecules <- skyline_check$Molecule
area <- area_filter$Molecule

skyline_library <- skyline_library[!skyline_library$Molecule %in% area,]
skyline_library <- skyline_library[!skyline_library$Molecule %in% molecules,]
write_csv(skyline_library, "results/molecule_quantification_filtered_DB044.csv")


# make long df
skyline_long <- skyline_library[,c(1:2, 5:6)]

#plot RT per molecule
ggplot(skyline_long, aes(Molecule, Best.Retention.Time)) + 
  geom_boxplot(outliers = F) +
  theme_classic()

ggsave("plots/retention_time_shifts_DB044.png")
