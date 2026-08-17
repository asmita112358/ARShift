library(spatstat)
library(ggplot2)
library(here)
library(stringr)
library(mcprogress)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/main_fun.R", echo = FALSE)

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
ppp_muco <- list()
for(i in seq_along(muco.list)){
  path <- muco.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=as.factor(img.df$sciname))
  ppp_muco[[i]] <- ppp_img
}
ppp_healthy <- list()
for(i in seq_along(healthy.list)){
  path <- healthy.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=as.factor(img.df$sciname))
 
  ppp_healthy[[i]] <- ppp_img
}

##Save all real data results in dir
dir <- "~/Library/CloudStorage/OneDrive-JohnsHopkins/spatial_assoc_test/ARShift/real_data_results/"
#saveRDS(W_new, paste0(dir,"test.RDS"))

muco_taxa <- lapply(ppp_muco, function(x)levels(x$marks))
muco_taxa <- Reduce(intersect, muco_taxa)

healthy_taxa <- lapply(ppp_healthy, function(x)levels(x$marks))
healthy_taxa <- Reduce(intersect, healthy_taxa)

setdiff(healthy_taxa, muco_taxa)   ##Report this in the manuscript

all_taxa <- intersect(muco_taxa, healthy_taxa)
d <- length(all_taxa)

analyse_slide <- function(slide = ppp_healthy[[5]], lite = FALSE, taxa_subset = all_taxa) {
  
  
  d <- length(taxa_subset)
  ARshift <- VCshift <- Torshift <- Rshift <- matrix(NA_real_, nrow = d, ncol = d)
  rownames(ARshift) <- rownames(VCshift) <- rownames(Torshift) <- rownames(Rshift) <- taxa_subset
  colnames(ARshift) <- colnames(VCshift) <- colnames(Torshift) <- colnames(Rshift) <- taxa_subset
  
  for(i in 1:(d-1)) {
    for(j in (i + 1):d ){
      
      rmax <- 0.3 * incircle(slide$window)$r
      r <- seq(0, rmax, length.out = 50)
      
      obj <- tryCatch(
        {
          test_spatial_association(
            slide,
            base_taxa = taxa_subset[i],
            shift_taxa = taxa_subset[j],
            r = r,
            type = "hom",
            lite = lite
          )
        },
        error = function(e) {
          message("Failed for i = ", i, ", j = ", j, ": ", e$message)
          return(NULL)
        }
      )
      
      if (is.null(obj)) {
        ARshift[i, j] <- NA_real_
        Rshift[i, j]  <- NA_real_
        VCshift[i, j] <- NA_real_
        Torshift[i, j] <- NA_real_
      } else {
        ARshift[i, j] <- obj$pval_Kcross_area
        Rshift[i, j]  <- obj$pval_Kcross_uncorr
        VCshift[i, j] <- obj$pval_Kcross_vc
        Torshift[i, j] <- obj$pval_Kcross_tor
      }
      
     
      cat(i, j, "\n")
    }
  }
  
  return(list(
    ARshift = ARshift,
    Rshift = Rshift,
    VCshift = VCshift,
    Torshift = Torshift
  ))
}
library(parallel)


##Analysis for slides with mucositis

results_muco <- pmclapply(ppp_muco, analyse_slide, mc.cores = 7L, lite = FALSE)
saveRDS(results_muco, file = paste0(dir,"resultsmuco.rds"))

sapply(results_muco, function(x)x$ARshift[1,2])
sapply(results_muco, function(x)x$VCshift[1,2])
sapply(results_muco, function(x)x$Rshift[1,2])
sapply(results_muco, function(x)x$Torshift[1,2])


results_muco_qvalue <- list()
d <- nrow(results_muco[[1]]$ARshift)
taxa_names <- rownames(results_muco[[1]]$ARshift)
for(i in 1:length(results_muco)){
 
  #qARshift <- qVCshift <- qTorshift <- qRshift <- matrix(NA_real_, nrow = d, ncol = d)
 
  
  qARshift <- matrix(p.adjust(results_muco[[i]]$ARshift, "BH"), nrow = d)
  qVCshift <- matrix(p.adjust(results_muco[[i]]$VCshift, "BH"), nrow = d)
  qTorshift <- matrix(p.adjust(results_muco[[i]]$Torshift, "BH"), nrow = d)
  qRshift <- matrix(p.adjust(results_muco[[i]]$Rshift, "BH"), nrow = d)
  rownames(qARshift) <- rownames(qVCshift) <- rownames(qTorshift) <- rownames(qRshift) <- taxa_names
  colnames(qARshift) <- colnames(qVCshift) <- colnames(qTorshift) <- colnames(qRshift) <- taxa_names
  
  results_muco_qvalue[[i]] <- list(
    ARshift = qARshift,
    Rshift = qRshift,
    VCshift = qVCshift,
    Torshift = qTorshift
  )
  
}

sapply(results_muco_qvalue, function(x)x$ARshift[1,2])
sapply(results_muco_qvalue, function(x)x$VCshift[1,2])
sapply(results_muco_qvalue, function(x)x$Rshift[1,2])
sapply(results_muco_qvalue, function(x)x$Torshift[1,2])


saveRDS(results_muco_qvalue, file = paste0(dir,"Qresultsmuco.rds"))


##Analysis for healthy slides

results_healthy <- pmclapply(ppp_healthy, analyse_slide, mc.cores = 7L, lite = FALSE)
saveRDS(results_healthy, file = paste0(dir,"resultshealthy.rds"))

sapply(results_healthy, function(x)x$ARshift[1,2])
sapply(results_healthy, function(x)x$VCshift[1,2])
sapply(results_healthy, function(x)x$Rshift[1,2])
sapply(results_healthy, function(x)x$Torshift[1,2])

results_healthy_qvalue <- list()
d <- nrow(results_healthy[[1]]$ARshift)

taxa_names <- rownames(results_healthy[[1]]$ARshift)

for(i in 1:length(results_healthy)){
  
  #qARshift <- qVCshift <- qTorshift <- qRshift <- matrix(NA_real_, nrow = d, ncol = d)
  
  
  qARshift <- matrix(p.adjust(results_healthy[[i]]$ARshift, "BH"), nrow = d)
  qVCshift <- matrix(p.adjust(results_healthy[[i]]$VCshift, "BH"), nrow = d)
  qTorshift <- matrix(p.adjust(results_healthy[[i]]$Torshift, "BH"), nrow = d)
  qRshift <- matrix(p.adjust(results_healthy[[i]]$Rshift, "BH"), nrow = d)
  rownames(qARshift) <- rownames(qVCshift) <- rownames(qTorshift) <- rownames(qRshift) <- taxa_names
  colnames(qARshift) <- colnames(qVCshift) <- colnames(qTorshift) <- colnames(qRshift) <- taxa_names
  
  results_healthy_qvalue[[i]] <- list(
    ARshift = qARshift,
    Rshift = qRshift,
    VCshift = qVCshift,
    Torshift = qTorshift
  )
  
}

sapply(results_healthy_qvalue, function(x)x$ARshift[1,2])
sapply(results_healthy_qvalue, function(x)x$VCshift[1,2])
sapply(results_healthy_qvalue, function(x)x$Rshift[1,2])
sapply(results_healthy_qvalue, function(x)x$Torshift[1,2])


saveRDS(results_healthy_qvalue, file = paste0(dir,"Qresultshealthy.rds"))


##plot

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

extract_pairwise_results <- function(results_list,
                                     group_label,
                                     method = "ARshift",
                                     value_name = "pvalue") {
  
  out <- lapply(seq_along(results_list), function(s) {
    
    mat <- results_list[[s]][[method]]
    
    taxa <- rownames(mat)
    d <- length(taxa)
    
    # Extract upper triangle only
    inds <- which(upper.tri(mat), arr.ind = TRUE)
    
    data.frame(
      slide = s,
      group = group_label,
      taxon1 = taxa[inds[, 1]],
      taxon2 = taxa[inds[, 2]],
      pair = paste(taxa[inds[, 1]], taxa[inds[, 2]], sep = " vs "),
      value = mat[inds],
      method = method,
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(out)
  names(out)[names(out) == "value"] <- value_name
  
  return(out)
}

df_healthy_p <- extract_pairwise_results(
  results_list = results_healthy,
  group_label = "Healthy",
  method = "ARshift",
  value_name = "pvalue"
)

df_muco_p <- extract_pairwise_results(
  results_list = results_muco,
  group_label = "Mucositis",
  method = "ARshift",
  value_name = "pvalue"
)

df_p <- bind_rows(df_healthy_p, df_muco_p)
df_healthy_q <- extract_pairwise_results(
  results_list = results_healthy_qvalue,
  group_label = "Healthy",
  method = "ARshift",
  value_name = "qvalue"
)

df_muco_q <- extract_pairwise_results(
  results_list = results_muco_qvalue,
  group_label = "Mucositis",
  method = "ARshift",
  value_name = "qvalue"
)

df_q <- bind_rows(df_healthy_q, df_muco_q)
plot_one_pair_log <- function(df,
                              taxon1,
                              taxon2,
                              value_col = "qvalue",
                              alpha_line = 0.05) {
  
  pair_name_1 <- paste(taxon1, taxon2, sep = " vs ")
  pair_name_2 <- paste(taxon2, taxon1, sep = " vs ")
  
  df_pair <- df %>%
    filter(pair %in% c(pair_name_1, pair_name_2)) %>%
    filter(!is.na(.data[[value_col]])) %>%
    mutate(log_value = -log10(pmax(.data[[value_col]], 1e-6)))
  
  ggplot(df_pair, aes(x = group, y = log_value, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
    geom_jitter(aes(color = group),
                width = 0.12,
                size = 2.2,
                alpha = 0.85,
                show.legend = FALSE) +
    geom_hline(yintercept = -log10(alpha_line),
               linetype = "dashed",
               color = "gray30") +
    scale_fill_manual(values = c("Healthy" = "#4DAF4A",
                                 "Mucositis" = "#E41A1C")) +
    scale_color_manual(values = c("Healthy" = "#4DAF4A",
                                  "Mucositis" = "#E41A1C")) +
    labs(
      title = paste0(taxon1, " vs ", taxon2),
      x = NULL,
      y = ifelse(value_col == "qvalue",
                 expression(-log[10]("BH q-value")),
                 expression(-log[10]("p-value")))
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}

plot_one_pair_log(df_q, taxon1 = all_taxa[1], taxon2 = all_taxa[2], value_col = "qvalue")
plot_dir <- paste0(dir, "pairwise_plots_ARshift_qvalue/")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

taxa_pairs <- combn(all_taxa, 2, simplify = FALSE)


for (pair in taxa_pairs) {
  
  taxon1 <- pair[1]
  taxon2 <- pair[2]
  
  p <- plot_one_pair(
    df = df_q,
    taxon1 = taxon1,
    taxon2 = taxon2,
    value_col = "qvalue",
    alpha_line = 0.05
  )
  
  file_name <- paste0(
    str_replace_all(taxon1, "[^A-Za-z0-9]+", "_"),
    "_vs_",
    str_replace_all(taxon2, "[^A-Za-z0-9]+", "_"),
    "_ARshift_qvalue.png"
  )
  
  ggsave(
    filename = file.path(plot_dir, file_name),
    plot = p,
    width = 5,
    height = 4
  )
}

plot_dir_p <- paste0(dir, "pairwise_plots_ARshift_pvalue/")
dir.create(plot_dir_p, showWarnings = FALSE, recursive = TRUE)

for (pair in taxa_pairs) {
  
  taxon1 <- pair[1]
  taxon2 <- pair[2]
  
  p <- plot_one_pair(
    df = df_p,
    taxon1 = taxon1,
    taxon2 = taxon2,
    value_col = "pvalue",
    alpha_line = 0.05
  )
  
  file_name <- paste0(
    str_replace_all(taxon1, "[^A-Za-z0-9]+", "_"),
    "_vs_",
    str_replace_all(taxon2, "[^A-Za-z0-9]+", "_"),
    "_ARshift_pvalue.png"
  )
  
  ggsave(
    filename = file.path(plot_dir_p, file_name),
    plot = p,
    width = 5,
    height = 4
  )
}

##

results_healthy_qvalue <- readRDS("real_data_results/Qresultshealthy.rds")
results_muco_qvalue <- readRDS("real_data_results/Qresultsmuco.rds")

sapply(results_healthy_qvalue, function(x)x$ARshift[1,2])
sapply(results_healthy_qvalue, function(x)x$VCshift[1,2])
sapply(results_healthy_qvalue, function(x)x$Rshift[1,2])
sapply(results_healthy_qvalue, function(x)x$Torshift[1,2])

sapply(results_muco_qvalue, function(x)x$ARshift[1,2])
sapply(results_muco_qvalue, function(x)x$VCshift[1,2])
sapply(results_muco_qvalue, function(x)x$Rshift[1,2])
sapply(results_muco_qvalue, function(x)x$Torshift[1,2])

