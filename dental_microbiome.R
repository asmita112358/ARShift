library(spatstat)
library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome3")

muco.list <- c("2023_02_08_hsdm_group_2_sample_06_fov_01_centroid_sciname.csv",            
               "2023_02_08_hsdm_group_2_sample_06_fov_02_centroid_sciname.csv" ,           
               "2023_02_18_hsdm_group_II_patient_6_fov_01_centroid_sciname.csv"  ,         
               "2023_10_16_hsdm_slide_IIB_fov_01_centroid_sciname.csv",                    
               "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv" ,                   
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_01_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_02_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_03_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_04_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_05_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_06_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_14_aspect_MB_fov_01_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_14_aspect_MB_fov_02_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_15_aspect_MB_fov_01_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_15_aspect_MB_fov_02_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_01_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_02_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_03_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_04_centroid_sciname.csv")

healthy.list <- c( "2023_02_08_hsdm_group_1_sample_06_fov_01_centroid_sciname.csv",           
                   "2023_02_08_hsdm_group_1_sample_11_fov_01_centroid_sciname.csv",           
                   "2023_02_08_hsdm_group_1_sample_12_fov_01_centroid_sciname.csv",           
                   "2023_02_18_hsdm_group_I_patient_11_fov_01_centroid_sciname.csv",          
                   "2023_02_18_hsdm_group_I_patient_11_fov_02_centroid_sciname.csv",          
                   "2023_02_18_hsdm_group_I_patient_13_fov_01_centroid_sciname.csv",          
                   "2023_02_18_hsdm_group_I_patient_6_fov_01_centroid_sciname.csv",           
                   "2023_10_16_hsdm_slide_IL_fov_01_centroid_sciname.csv",                    
                   "2023_10_16_hsdm_slide_IL_fov_02_centroid_sciname.csv",                    
                   "2023_10_16_hsdm_slide_IL_fov_03_centroid_sciname.csv",                    
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_01_centroid_sciname.csv",
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_02_centroid_sciname.csv",
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_03_centroid_sciname.csv",
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_04_centroid_sciname.csv")

process_df <- function(df, tile_num=NULL){
  if (is.null(tile_num)) tile <- df
  else
    tile <- df[df$tile==tile_num,]
  
  #coords <- gsub("\\[|\\]", "", tile$coord)  # Remove brackets
  coords <- gsub("\\[|\\]|\\(|\\)", "", tile$coord)
  split_coords <- strsplit(coords, ",")  # Split into x and y
  
  x <- as.numeric(sapply(split_coords, `[`, 1))  # Extract x values
  y <- as.numeric(sapply(split_coords, `[`, 2))  # Extract y values
  
  # adjust x,y coordinates
  x <- x-min(x)
  y <- y-min(y)
  
  tile.df <- data.frame(x=x,y=y,sciname=tile$sciname, tile = tile$tile)
  return(tile.df)
}

##Read all data into list
df_muco <- list()
for(i in seq_along(muco.list)){
  path <- muco.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=as.factor(img.df$sciname))
  ppp_df <- data.frame(ppp_img)
  df_muco[[i]] <- ppp_df
}
df_healthy <- list()
for(i in seq_along(healthy.list)){
  path <- healthy.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=as.factor(img.df$sciname))
  ppp_df <- data.frame(ppp_img)
  df_healthy[[i]] <- ppp_df
}
