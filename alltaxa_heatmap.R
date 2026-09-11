##heatmap for all taxa
#run untill line 91 in dental microbiome.R
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
dir <- "~/Library/CloudStorage/OneDrive-JohnsHopkins/spatial_assoc_test/ARShift/real_data_results/"
results_muco <- list()
for(i in 1:length(ppp_muco)){
  
  #qARshift <- qVCshift <- qTorshift <- qRshift <- matrix(NA_real_, nrow = d, ncol = d)
  
  #read data
  data <- readRDS(file = paste0(dir,"resultsmuco_slide_",i,".rds"))
  results_muco[[i]] <- data$ARshift
  
}

results_healthy <- list()
for(i in 1:length(ppp_healthy)){
  
  data <- readRDS(file = paste0(dir, "resultshealthy_slide_",i,".rds"))
  results_healthy[[i]] <- data$ARshift
}


m <- results_healthy[[1]]
idx <- upper.tri(m)
taxa_pair <- paste(rownames(m)[row(m)[idx]], colnames(m)[col(m)[idx]], sep = "_")

out <- matrix(nrow = 91, ncol = 14 + 19)
rownames(out) <- taxa_pair
for(i in seq_along(ppp_muco)){
  idx <- upper.tri(results_muco[[i]])
  pvals <- results_muco[[i]][idx]
  out[,i] <- as.numeric(pvals)
}
k <- length(ppp_muco)
for(i in seq_along(ppp_healthy)){
  idx <- upper.tri(results_healthy[[i]])
  pvals <- results_healthy[[i]][idx]
  out[,i+k] <- as.numeric(pvals)
}

colnames(out) <- c(paste0("Mucositis_", seq_along(ppp_muco)), paste0("Healthy_", seq_along(ppp_healthy)))
out <- -log10(out)

mid <- -log10(0.05)
n <- 100

# Same style as pheatmap default: blue -> white -> red
my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(n)

# Split colors equally around the midpoint
minv <- min(out, na.rm = TRUE)
maxv <- max(out, na.rm = TRUE)

low_n  <- floor(n / 2)
high_n <- n - low_n

# Breaks so that mid is exactly the center
my_breaks <- c(
  seq(minv, mid, length.out = low_n + 1),
  seq(mid, maxv, length.out = high_n + 1)[-1]
)

# identify columns by group
muco_cols <- grep("^Mucositis_", colnames(out))
healthy_cols <- grep("^Healthy_", colnames(out))

# cluster within each group
hc_muco <- hclust(dist(t(out[, muco_cols, drop = FALSE])))
hc_healthy <- hclust(dist(t(out[, healthy_cols, drop = FALSE])))

# reorder columns within each group
out2 <- cbind(
  out[, muco_cols[hc_muco$order], drop = FALSE],
  out[, healthy_cols[hc_healthy$order], drop = FALSE]
)

# column annotations
ann <- data.frame(
  Group = c(rep("Mucositis", length(muco_cols)),
            rep("Healthy", length(healthy_cols)))
)
ann_colors <- list(
  Group = c(
    Mucositis = "darkred",
    Healthy = "navy"
  )
)
rownames(ann) <- colnames(out2)

# plot
ph1 <- pheatmap(out2,
         cluster_cols = FALSE,
         cluster_rows = TRUE,   # or TRUE if you want row clustering
         annotation_col = ann,
         annotation_colors = ann_colors,
         annotation_legend = FALSE,
         gaps_col = length(muco_cols),
         treeheight_row = 0,
         treeheight_col = 0,
         margins = c(10, 10),
         color = my_colors,
         breaks = my_breaks,
         show_colnames = FALSE)

ph1
ggsave("plots/heatmap_alltaxa.png", width = 6, height = 10,unit = "in", dpi = 300)



results_healthy_qvalue <- readRDS("real_data_results/Qresultshealthy.rds")
results_muco_qvalue <- readRDS("real_data_results/Qresultsmuco.rds")


sum(sapply(results_healthy_qvalue, function(x) x$ARshift["Actinomyces", "Rothia"]) <= 0.05)
sum(sapply(results_muco_qvalue, function(x) x$ARshift["Actinomyces", "Rothia"]) <= 0.05)


###Heatmap for q value

outm <- matrix(nrow = 91, ncol = 19)
outh <- matrix(nrow = 91, ncol = 14)

qval_muco_ARshift <- lapply(results_muco_qvalue, function(x) x$ARshift)
qval_healthy_ARshift <- lapply(results_healthy_qvalue, function(x) x$ARshift)
rownames(outm) <- rownames(outh) <- taxa_pair
for(i in seq_along(ppp_muco)){
  idx <- upper.tri(qval_muco_ARshift[[i]])
  pvals <- qval_muco_ARshift[[i]][idx]
  outm[,i] <- as.numeric(pvals)
}
k <- length(ppp_muco)
for(i in seq_along(ppp_healthy)){
  idx <- upper.tri(qval_healthy_ARshift[[i]])
  pvals <- qval_healthy_ARshift[[i]][idx]
  outh[,i] <- as.numeric(pvals)
}

##Fisher's test between the number of significant pairs in mucositis and healthy groups
sig_m <- rowSums(outm <= 0.05)
sig_h <- rowSums(outh <= 0.05)
p_diffs <- c()
for(i in 1:91){
  TAB <- rbind(c(sig_m[i], 19 - sig_m[i]), c(sig_h[i], 14 - sig_h[i]))
  obj <- fisher.test(TAB)
  p_diffs[i] <- obj$p.value
}
qval_diffs <- data.frame(taxa_pair, sig_m, sig_h, p_diffs)

taxa_pairs_considered <- qval_diffs %>% filter(p_diffs <= 0.1)
saveRDS(qval_diffs, "real_data_results/qval_diffs.rds")


##t test between the qvalues from mucositis and healthy groups
p_t <- c()
for(i in 1:91){
  obj <- t.test(outm[i,], outh[i,])
  p_t[i] <- obj$p.value
}
qval_diffs_t <- data.frame(taxa_pair, sig_m, sig_h, p_t)
taxa_pairs_considered <- qval_diffs_t %>% filter(p_t <= 0.05)

##Extract the rows that have Fusobacterium in taxa pair

Fuso_taxa_pairs <- qval_diffs_t %>% filter(grepl("Fusobacterium", taxa_pair))
saveRDS(Fuso_taxa_pairs, "real_data_results/Fuso_taxa_pairs.rds")
out <- cbind(outm, outh)
colnames(out) <- c(paste0("Mucositis_", seq_along(ppp_muco)), paste0("Healthy_", seq_along(ppp_healthy)))
out <- -log10(out)
mid <- -log10(0.05)
n <- 100

# Same style as pheatmap default: blue -> white -> red
my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(n)

# Split colors equally around the midpoint
minv <- min(out, na.rm = TRUE)
maxv <- max(out, na.rm = TRUE)

low_n  <- floor(n / 2)
high_n <- n - low_n

# Breaks so that mid is exactly the center
my_breaks <- c(
  seq(minv, mid, length.out = low_n + 1),
  seq(mid, maxv, length.out = high_n + 1)[-1]
)

# identify columns by group
muco_cols <- grep("^Mucositis_", colnames(out))
healthy_cols <- grep("^Healthy_", colnames(out))


# reorder columns within each group by the pvalue clustering
out2 <- cbind(
  out[, muco_cols[hc_muco$order], drop = FALSE],
  out[, healthy_cols[hc_healthy$order], drop = FALSE]
)

# column annotations
ann <- data.frame(
  Group = c(rep("Mucositis", length(muco_cols)),
            rep("Healthy", length(healthy_cols)))
)
ann_colors <- list(
  Group = c(
    Mucositis = "darkred",
    Healthy = "navy"
  )
)
rownames(ann) <- colnames(out2)

# plot
ph2 <- pheatmap(
  out2,
  cluster_cols = FALSE,
  cluster_rows = ph1$tree_row,  # reuse existing hclust object
  annotation_col = ann,
  annotation_colors = ann_colors,
  annotation_legend = FALSE,
  gaps_col = length(muco_cols),
  treeheight_row = 0,
  treeheight_col = 0,
  margins = c(10, 10),
  color = my_colors,
  breaks = my_breaks,
  show_colnames = FALSE
)
ph2


library(patchwork)
library(ggplotify)
ph1 <- as.ggplot(ph1)
ph2 <- as.ggplot(ph2)
ph1 + ph2 +
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave("plots/pq_heatmap_alltaxa.png", width = 12, height = 10, unit = "in", dpi = 1200)

