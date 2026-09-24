load("Schurch_CRC_Data_Cleaned.rda")
library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)
library(alphahull)

# Plot the tertiary lymphoid structures
bm <- c("B cells", "CD11b+CD68+ macrophages", "CD68+CD163+ macrophages", "CD68+ macrophages", "CD68+ macrophages GzmB+")

cells %>%
  filter(Region == "reg001_A") %>%
  #filter(ClusterName %in% bm) %>%
  ggplot(aes(x = x, y = y, color = ClusterName)) + 
  geom_point() + 
  theme_bw() 


regions_with_both <- cells %>%
  # 1. Create the celltype column
  mutate(
    celltype = case_when(
      ClusterName == "B cells" ~ "B cells",
      stringr::str_detect(ClusterName, "macrophages") ~ "macrophages",
      TRUE ~ "Other"
    )
  ) %>%
  # 2. Filter for only the two cell types of interest
  filter(celltype %in% c("B cells", "macrophages")) %>%
  # 3. Count cells per Region and celltype
  group_by(PID, Group, Region, celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  # 4. Pivot wider so B cells and macrophages are separate columns
  pivot_wider(names_from = celltype, values_from = n_cells, values_fill = 0) %>%
  # 5. Filter for regions where both cell types strictly exceed 50
  filter(`B cells` > 50 & macrophages > 50) 
  # 6. Extract just the unique Region names
  #pull(Region)

# View the result
print(regions_with_both, n = 33)

slide_list <- regions_with_both %>% pull(Region)


##7. Extract the cells in the short listed region list
cells_shortlist <- cells %>% filter(Region %in% slide_list) %>% mutate(celltype = case_when(
  ClusterName == "B cells" ~ "B cells",
  stringr::str_detect(ClusterName, "macrophages") ~ "macrophages",
  TRUE ~ "Other"
)) %>% filter(celltype %in% c("B cells", "macrophages")) %>% select(PID, Region, Group, x,y, celltype)


# 8. Extract windows for shortlisted cells

windows <- list()
for(i in seq_along(slide_list)){
  region <- slide_list[i]
  all_cells <- cells %>%
    filter(Region == region) %>%
   select(x,y)
  ppp_obj <- ppp(all_cells$x, all_cells$y, window = owin(range(all_cells$x), range(all_cells$y)))
  lambda <- density(ppp_obj, bw = "bw.diggle", at = "points")
  
}

