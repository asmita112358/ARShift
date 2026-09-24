# -----------------------------------------------------------------------------
# Prepare the Schurch et al. mIF data in CRC for analysis
# -----------------------------------------------------------------------------

# Packages
library(readxl)
library(tidyverse)
library(data.table)

# Cell-level data downloaded from: https://data.mendeley.com/datasets/mpjzbtfgfr/1

# Load in the cell-level data
raw_data <- data.table::fread("/Users/ssamorod/Documents/Schurch_CRC_Data/CRC_clusters_neighborhoods_markers.csv")

# Load in the clinical data
clinical <- read_xlsx("/Users/ssamorod/Documents/Schurch_CRC_Data/CRC_TMAs_patient_annotations.xlsx")

# Remove the bottom rows
clinical <- clinical[1:35,]

# Load in tumor composition data for the multi-tumor TMA (don't use this?)
# composition <- read_xlsx("/Users/ssamorod/Documents/Schurch_CRC_Data/Multi-tumor_TMA_composition.xlsx")

# Check the columns
identical(raw_data$`X:X`, raw_data$`X_withinTile:X_withinTile`)
identical(raw_data$`Y:Y`, raw_data$`Y_withinTile:Y_withinTile`)

# Select columns of interest
cells <- raw_data %>%
  select(patients, groups, `File Name`, `X:X`, `Y:Y`, ClusterName, `neighborhood name`) %>%
  rename(
    "PID" = "patients",
    "x" = "X:X",
    "y" = "Y:Y",
    "Neighborhood" = `neighborhood name`,
    "group" = "groups",
    "Region" = `File Name`
  )

# How many images per patient? How many cells per image?
cells %>%
  group_by(PID, Region) %>%
  summarise(
    n_cells = n()
  ) %>%
  View()
  
# Plot some example images
cells %>%
  filter(Region == "reg001_A") %>%
  ggplot(aes(x = x, y = y, color = ClusterName)) + 
  geom_point() + 
  theme_bw() 

# What were the detected cell types?
unique(cells$ClusterName)

# What were the neighborhoods?
unique(cells$Neighborhood)

# Plot just the T cells
cells %>%
  filter(Region == "reg001_A") %>%
  filter(grepl("T cells", ClusterName)) %>%
  ggplot(aes(x = x, y = y, color = ClusterName)) + 
  geom_point() + 
  theme_bw() 

cells %>%
  filter(Region == "reg002_A") %>%
  filter(grepl("T cells", ClusterName)) %>%
  ggplot(aes(x = x, y = y, color = ClusterName)) + 
  geom_point() + 
  theme_bw() 

# Plot the tertiary lymphoid structures
tls <- c("B cells", "CD11c+ DCs", "CD4+ T cells", "CD4+ T cells CD45RO+", "CD8+ T cells")

cells %>%
  filter(Region == "reg003_B") %>%
  filter(ClusterName %in% tls) %>%
  ggplot(aes(x = x, y = y, color = ClusterName)) + 
  geom_point() + 
  theme_bw() 

# ---
# Merge in the clinical data
# ---

cells <- left_join(
  cells,
  clinical,
  by = join_by("PID" == "Patient")
)

# cells <- left_join(
#   cells,
#   composition,
#   by = join_by("Region" == "region") 
# )

# Save the data
save(
  cells,
  file = "/Users/ssamorod/Documents/Schurch_CRC_Data/Schurch_CRC_Data_Cleaned.rda"
)
