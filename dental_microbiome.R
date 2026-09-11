library(spatstat)
library(ggplot2)
library(here)
library(stringr)
library(mcprogress)
library(ggpp)
library(phyloseq)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/main_fun_fast.R", echo = FALSE)


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
list.files("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/windows")

##Read all data into list
ppp_muco <- list()
for(i in seq_along(muco.list)){
  win_path <- paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/windows/", muco.list[i], ".RData")
  load(win_path)
  path <- muco.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, window = W,
              marks=as.factor(img.df$sciname))
  ppp_muco[[i]] <- ppp_img
}
ppp_healthy <- list()
for(i in seq_along(healthy.list)){
  win_path <- paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/windows/", healthy.list[i], ".RData")
  load(win_path)
  path <- healthy.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, window = W,
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

two_taxa <- setdiff(healthy_taxa, muco_taxa)   ##Report this in the manuscript

#how many muco slides do not contain two_taxa?
which(sapply(ppp_muco, function(x) all(two_taxa %notin% levels(x$marks))))  #slide 6
table(ppp_muco[[4]]$marks)


all_taxa <- intersect(muco_taxa, healthy_taxa)
d <- length(all_taxa)

#table of all_taxa frequency in healthy and mucositis slides
healthy_taxa_table <- lapply(ppp_healthy, function(x)table(factor(x$marks, levels = all_taxa)))
healthy_taxa_table <- do.call(rbind, healthy_taxa_table)

muco_taxa_table <- lapply(ppp_muco, function(x)table(factor(x$marks, levels = all_taxa)))
muco_taxa_table <- do.call(rbind, muco_taxa_table)

all_taxa_table <- rbind(healthy_taxa_table, muco_taxa_table)
all_taxa_prop <- all_taxa_table / rowSums(all_taxa_table)
rownames(all_taxa_prop) <- c(paste0("Healthy_", seq_along(ppp_healthy)), paste0("Mucositis_", seq_along(ppp_muco)))

heatmap(all_taxa_prop, Rowv = NA, margins = c(10, 10))
OTU  <- otu_table(as.matrix(all_taxa_table), taxa_are_rows = FALSE) 
META <- sample_data(data.frame(group_label = as.factor(group_label)))

# Combine into a single phyloseq object
pseq <- phyloseq(OTU, META)

# 3. Run ANCOM-BC (or ancombc2)
# Replace "GroupVariable" with the column name in your metadata you want to test
out <- ancombc2(
  data = pseq, 
  fix_formula = "group_label", 
  p_adj_method = "fdr"
)
results_df <- out$res
head(results_df)
results_df %>% select(taxon, diff_group_labelmuco, diff_robust_group_labelmuco)
#wrapper for all functions

analyse_slide <- function(slide = ppp_healthy[[5]], lite = FALSE, taxa_subset = all_taxa) {
  rmax <- 0.3 * incircle(slide$window)$r
  r <- seq(0, rmax, length.out = 50)
  
  d <- length(taxa_subset)
  ARshift <- VCshift <- Torshift <- Rshift <- matrix(NA_real_, nrow = d, ncol = d)
  rownames(ARshift) <- rownames(VCshift) <- rownames(Torshift) <- rownames(Rshift) <- taxa_subset
  colnames(ARshift) <- colnames(VCshift) <- colnames(Torshift) <- colnames(Rshift) <- taxa_subset
  
  for(i in 1:(d-1)) {
    for(j in (i + 1):d ){
     
        
       obj <- tryCatch(
        {
          test_spatial_association_pmclapply(
            slide,
            base_taxa = taxa_subset[i],
            shift_taxa = taxa_subset[j],
            r = r,
            type = "hom",
            lite = lite,
            correction = "border"
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
library(mcprogress)

##Analysis for slides with mucositis
#taxa_subset <- c("Actinomyces", "Campylobacter", "Veillonella", "Porphyromonas", "Streptococcus", "Fusobacterium", "Gemella")

analyse_slide(slide = ppp_muco[[9]], lite = FALSE, taxa_subset = taxa_subset)




for(i in 5:19){
  results_muco <- analyse_slide(slide = ppp_muco[[i]], lite = FALSE, taxa_subset = all_taxa)
  cat("Finished slide ", i, "\n")
  saveRDS(results_muco, file = paste0(dir,"resultsmuco_slide_",i,".rds"))
}

beepr::beep(4)

##lite analysis for slide 4. 

results_muco <- analyse_slide(ppp_muco[[4]], lite = TRUE, taxa_subset = all_taxa)
saveRDS(results_muco, file = paste0(dir,"resultsmuco_slide_4.rds"))


##get qvalues
results_muco_qvalue <- list()
d <- length(all_taxa)
taxa_names <- all_taxa
for(i in 1:length(ppp_muco)){

  #qARshift <- qVCshift <- qTorshift <- qRshift <- matrix(NA_real_, nrow = d, ncol = d)

  #read data
  results_muco <- readRDS(file = paste0(dir,"resultsmuco_slide_",i,".rds"))
  qARshift <- matrix(p.adjust(results_muco$ARshift, "BH"), nrow = d)
  qVCshift <- matrix(p.adjust(results_muco$VCshift, "BH"), nrow = d)
  qTorshift <- matrix(p.adjust(results_muco$Torshift, "BH"), nrow = d)
  qRshift <- matrix(p.adjust(results_muco$Rshift, "BH"), nrow = d)
  rownames(qARshift) <- rownames(qVCshift) <- rownames(qTorshift) <- rownames(qRshift) <- taxa_names
  colnames(qARshift) <- colnames(qVCshift) <- colnames(qTorshift) <- colnames(qRshift) <- taxa_names

  results_muco_qvalue[[i]] <- list(
    ARshift = qARshift,
    Rshift = qRshift,
    VCshift = qVCshift,
    Torshift = qTorshift
  )

}



saveRDS(results_muco_qvalue, file = paste0(dir,"Qresultsmuco.rds"))


##Analysis for healthy slides

for(i in seq_along(ppp_healthy)){
  results_healthy <- analyse_slide(slide = ppp_healthy[[i]], lite = FALSE, taxa_subset = all_taxa)
  cat("Finished slide ", i, "\n")
  saveRDS(results_healthy, file = paste0(dir,"resultshealthy_slide_",i,".rds"))
  
}
saveRDS(results_healthy, file = paste0(dir,"resultshealthy_subset.rds"))
beepr::beep(4)

#qvalue for healthy slides
results_healthy_qvalue <- list()

for(i in 1:length(ppp_healthy)){

  #qARshift <- qVCshift <- qTorshift <- qRshift <- matrix(NA_real_, nrow = d, ncol = d)

  results_healthy <- readRDS(file = paste0(dir,"resultshealthy_slide_",i,".rds"))
  qARshift <- matrix(p.adjust(results_healthy$ARshift, "BH"), nrow = d)
  qVCshift <- matrix(p.adjust(results_healthy$VCshift, "BH"), nrow = d)
  qTorshift <- matrix(p.adjust(results_healthy$Torshift, "BH"), nrow = d)
  qRshift <- matrix(p.adjust(results_healthy$Rshift, "BH"), nrow = d)
  rownames(qARshift) <- rownames(qVCshift) <- rownames(qTorshift) <- rownames(qRshift) <- taxa_names
  colnames(qARshift) <- colnames(qVCshift) <- colnames(qTorshift) <- colnames(qRshift) <- taxa_names

  results_healthy_qvalue[[i]] <- list(
    ARshift = qARshift,
    Rshift = qRshift,
    VCshift = qVCshift,
    Torshift = qTorshift
  )

}

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


# df_healthy_q <- extract_pairwise_results(
#   results_list = results_healthy_qvalue,
#   group_label = "Healthy",
#   method = "ARshift",
#   value_name = "qvalue"
# )
# 
# df_muco_q <- extract_pairwise_results(
#   results_list = results_muco_qvalue,
#   group_label = "Mucositis",
#   method = "ARshift",
#   value_name = "qvalue"
# )

#df_q <- bind_rows(df_healthy_q, df_muco_q)
library(ggrepel)

plot_one_pair_log <- function(df,
                              taxon1,
                              taxon2,
                              value_col = "pvalue",
                              alpha_line = 0.05,
                              slide_col = "slide") {
  
  pair_name_1 <- paste(taxon1, taxon2, sep = " vs ")
  pair_name_2 <- paste(taxon2, taxon1, sep = " vs ")
  
  df_pair <- df %>%
    filter(pair %in% c(pair_name_1, pair_name_2)) %>%
    filter(!is.na(.data[[value_col]])) %>%
    mutate(
      log_value = -log10(pmax(.data[[value_col]], 1e-6)),
      significant = .data[[value_col]] < alpha_line
    )
  jitter_pos <- position_jitter(width = 0.1, height = 0, seed = 123)
  
  ggplot(df_pair, aes(x = group, y = log_value, fill = group)) +
    geom_violin(alpha = 0.5, color = "black") +
    geom_point(
      aes(color = group),
      position = jitter_pos,
      size = 2.2,
      alpha = 1,
      show.legend = FALSE
    ) +
    geom_text_s(
      data = df_pair %>% filter(significant),
      aes(label = .data[[slide_col]]),
      color = "black",
      size = 3,
      nudge_y = 0.5
    ) +
    geom_hline(
      yintercept = -log10(alpha_line),
      linetype = "dashed",
      color = "gray30"
    ) +
    scale_fill_manual(values = c("Healthy" = "navy",
                                 "Mucositis" = "darkred")) +
    scale_color_manual(values = c("Healthy" = "navy",
                                  "Mucositis" = "darkred")) +
    labs(
      title = paste0(taxon1, " vs ", taxon2),
      x = NULL,
      y = ifelse(
        value_col == "qvalue",
        expression(-log[10]("BH q-value")),
        expression(-log[10]("p-value"))
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}

plot_one_pair_log <- function(df,
                              taxon1,
                              taxon2,
                              value_col = "pvalue",
                              alpha_line = 0.05,
                              slide_col = "slide",
                              jitter_width = 0.1) {
  
  pair_name_1 <- paste(taxon1, taxon2, sep = " vs ")
  pair_name_2 <- paste(taxon2, taxon1, sep = " vs ")
  
  set.seed(123)
  
  df_pair <- df %>%
    filter(pair %in% c(pair_name_1, pair_name_2)) %>%
    filter(!is.na(.data[[value_col]])) %>%
    mutate(
      log_value = -log10(pmax(.data[[value_col]], 1e-6)),
      significant = .data[[value_col]] < alpha_line,
      
      # numeric x-position for each group
      group_num = as.numeric(factor(group, levels = c("Healthy", "Mucositis"))),
      
      # pre-jittered x-position
      xjit = jitter(group_num, amount = jitter_width)
    )
  
  ggplot(df_pair, aes(y = log_value)) +
    
    geom_violin(
      aes(x = group_num, group = group, fill = group),
      alpha = 0.5,
      color = "black"
    ) +
    
    geom_point(
      aes(x = xjit, color = group),
      size = 2.2,
      alpha = 1,
      show.legend = FALSE
    ) +
    
    ggrepel::geom_text_repel(
      data = df_pair %>% filter(significant),
      aes(
        x = xjit,
        label = .data[[slide_col]]
      ),
      color = "black",
      size = 3,
      show.legend = FALSE,
      nudge_x = 0.08,
      nudge_y = 0.05,
      force = 0.3,
      box.padding = 0.1,
      point.padding = 0.05,
      
      min.segment.length = 0,
      segment.color = "black",
      segment.size = 0.4,
      
      max.overlaps = Inf
    ) +
    
    geom_hline(
      yintercept = -log10(alpha_line),
      linetype = "dashed",
      color = "gray30"
    ) +
    
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Healthy", "Mucositis")
    ) +
    
    scale_fill_manual(values = c("Healthy" = "#26818e",
                                 "Mucositis" = "#c33b4f")) +
    scale_color_manual(values = c("Healthy" = "#26818e",
                                  "Mucositis" = "#c33b4f")) +
    
    labs(
      title = paste0(taxon1, " vs ", taxon2),
      x = NULL,
      y = ifelse(
        value_col == "qvalue",
        expression(-log[10]("BH q-value")),
        expression(-log[10]("p-value"))
      )
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}
plot_one_pair_log(df_p, taxon1 = taxa_subset[1], taxon2 = taxa_subset[2], value_col = "pvalue")

plot_Kcross <- function(ppp_healthy, ppp_muco, taxon1, taxon2){
  Kcross_healthy <- lapply(ppp_healthy, function(x) Kcross(x, i = taxon1, j = taxon2, correction = "border", rmax = r_max)$bor)
  Kcross_muco <- lapply(ppp_muco, function(x) Kcross(x, i = taxon1, j = taxon2, correction = "border", rmax = r_max)$bor)
  r <- Kcross(ppp_healthy[[1]], i = taxon1, j = taxon2, correction = "border", rmax = r_max)$r
  
  # ##get L(r) - r for each slide
   Lcross_healthy <- lapply(Kcross_healthy, function(K) sqrt(K/pi) - r)
  Lcross_muco <- lapply(Kcross_muco, function(K) sqrt(K/pi) - r)
  
  ##Plot muco and healthy Lcross curve
  d<- do.call(rbind, c(Lcross_muco, Lcross_healthy))
  d_df <- as.data.frame(d)
  colnames(d_df) <- paste0("r_", seq_along(r))
  
  d_df$sample_id <- seq_len(nrow(d_df))
  d_df$group_label <- group_label
  
  d_long <- d_df %>%
    pivot_longer(
      cols = starts_with("r_"),
      names_to = "r_index",
      values_to = "d"
    ) %>%
    mutate(
      r_index = as.integer(sub("r_", "", r_index)),
      r = r[r_index]
    )
  
  ggplot(d_long, aes(x = r, y = d, color = group_label, group = sample_id)) +
    geom_line(alpha = 0.5) +
    scale_color_manual(values = c("#26818e",
                                  "#c33b4f")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6)+
    labs(y = expression(L(r) - r))+
  theme_minimal() +
    labs(
      color = "Group"
    )
    
}

plot_dir_p <- paste0(dir, "violin_plots_ARshift_pvalue/")
dir.create(plot_dir_p, showWarnings = FALSE, recursive = TRUE)

taxa_pairs_subset <- taxa_pairs[c(1,5, 7, 11, 14, 21)]
plot_list <- list()
Kcross_plots <- list()
k = 1
group_label <- rep(c("muco","healthy"), c(length(muco.list), length(healthy.list)))
n <- length(muco.list) + length(healthy.list)
incircle_radius <- c()
ppp_data <- c(ppp_healthy, ppp_muco)
for(i in seq_along(ppp_data)){
  incircle_radius[i] <- incircle(ppp_data[[i]]$window)$r
}

r_max <- 0.5*min(incircle_radius)

for (pair in taxa_pairs_subset) {
  
  taxon1 <- pair[1]
  taxon2 <- pair[2]
  
  plot_list[[k]] <- plot_one_pair_log(
    df = df_p,
    taxon1 = taxon1,
    taxon2 = taxon2,
    value_col = "pvalue",
    alpha_line = 0.05
  )
  
    Kcross_plots[[k]] <- plot_Kcross(ppp_healthy, ppp_muco, taxon1, taxon2)
  
  # file_name <- paste0(
  #   str_replace_all(taxon1, "[^A-Za-z0-9]+", "_"),
  #   "_vs_",
  #   str_replace_all(taxon2, "[^A-Za-z0-9]+", "_"),
  #   "_ARshift_pvalue.png"
  # )
  # 
  # ggsave(
  #   filename = file.path(plot_dir_p, file_name),
  #   plot = p,
  #   width = 5,
  #   height = 4
  # )
  k = k + 1
}

ggpubr::ggarrange(plotlist = c(plot_list, Kcross_plots),  nrow = length(taxa_pairs_subset),ncol = 2, common.legend = TRUE, legend = "top")
ggsave("real_data_results2/combined_violin_Kcross_plots.png", width = 8, height = 16, units = "in")



##Plot ppp healthy 3 and ppp_muco 18 for taxa_subset
taxa_subset <- c("Actinomyces", "Campylobacter", "Veillonella", "Fusobacterium", "Gemella")

df <- data.frame(ppp_healthy[[8]])
df <- df[df$marks %in% taxa_subset, ]
p1 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 0.5) +
  theme_minimal()+
  labs(title = "Healthy Slide 8") 
p1

built_plot <- ggplot_build(p1)
extracted_colors <- unique(built_plot$data[[1]]$colour)
extracted_colors
df <- data.frame(ppp_muco[[18]])
df <- df[df$marks %in% taxa_subset, ]
p2 <- ggplot(df, aes(x = x, y = y, color = marks))+
  geom_point(size = 0.5) +
  theme_minimal()+
  labs(title = "Mucositis Slide 18")
p2
ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave("real_data_results2/exampleslides.png", width = 8, height = 4, units = "in")
