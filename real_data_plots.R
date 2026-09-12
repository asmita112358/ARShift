##Plots for real data analysis
#Run until line 85 in dental_microbiome.R
library(pals)
library(ggvenn)
library(ggVennDiagram)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyr)
library(dplyr)
library(grid)
#Example plot, muco i
i = 17
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome_twosample")
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))

path <- muco.list[i]
img <- read.csv(file=paste0("centroid_sciname_tables/",path))
img.df <- process_df(img)
ppp_img=ppp_img=rthin(ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
                    marks=as.factor(img.df$sciname)), 0.2)
#ppp_img <- rthin(ppp_img, 0.2)
ppp_df <- data.frame(ppp_img) %>% filter(marks %in% all_taxa)
colnames(ppp_df)[3] <- "Taxa"
taxon_colors <- c(
  "Actinomyces"    = "#F8766D",  # salmon (specified)
  "Campylobacter"  = "#D99000",  # warm ochre
  "Capnocytophaga" = "#A68B00",  # muted gold
  "Fusobacterium"  = "#7A9A01",  # olive green
  "Gemella"        = "#53B400",  # green (specified)
  "Lautropia"      = "#001F3F",  # deep navy (specified)
  "Leptotrichia"   = "#00A98F",  # sea green
  "Porphyromonas"  = "#00A6A6",  # teal
  "Prevotella"     = "#1B9ECA",  # soft cyan-blue
  "Rothia"         = "#4C78A8",  # muted blue
  "Selenomonas"    = "#8E7CC3",  # soft purple
  "Streptococcus"  = "#C77CFF",  # lavender
  "Treponema"      = "#E76BF3",  # orchid
  "Veillonella"    = "#FB61D7"   # pink (specified)
)
muco1 <- ggplot(ppp_df, aes(x = x, y = y, color = Taxa)) +
  geom_point(size = 0.5) +
  ggtitle("Peri Implant Mucositis, Image 17")+
  theme_void() +
  scale_color_manual(values = taxon_colors) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, hjust = 0.5))+
  theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13)) 

muco1
dir <- "~/Library/CloudStorage/OneDrive-JohnsHopkins/spatial_assoc_test/ARShift/real_data_results/"
#ggsave(paste0(dir,"muco19.png"), muco1, width = 5, height = 6, units = "in")
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift")
results_healthy_qvalue <- readRDS("real_data_results/Qresultshealthy.rds")
results_muco_qvalue <- readRDS("real_data_results/Qresultsmuco.rds")

##Venn diagram for slide 17, mucositis

results_healthy_qvalue <- readRDS(paste0(dir,"Qresultshealthy.rds"))
results_muco_qvalue <- readRDS(paste0(dir,"Qresultsmuco.rds"))


ARShift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$ARshift))
VC_shift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$VCshift))
Rshift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$Rshift))
Torshift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$Torshift))

df <- data.frame(ARShift_q, VC_shift_q, Rshift_q, Torshift_q)
colSums(df <= 0.05)
taxa_names <- rownames(results_muco_qvalue[[i]]$ARshift)
taxa_pair <- combn(taxa_names, 2)
taxa_pair_list <- paste0(taxa_pair[1,], "_", taxa_pair[2,])

x = list(ARShift = taxa_pair_list[ARShift_q <= 0.05], VCShift = taxa_pair_list[VC_shift_q <= 0.05], 
         RShift = taxa_pair_list[Rshift_q <= 0.05], Torshift = taxa_pair_list[Torshift_q <= 0.05])
vplot <- ggVennDiagram(x, label = "count", label_size = 7) + #scale_fill_gradient(low = "#D73027", high = "gray90") + 
  ggtitle("Significantly co-clustered taxa pairs")+theme(plot.title = element_text(hjust = 0.5, size = 16), legend.position = "bottom")
vplot
muco1 + vplot +
plot_layout(ncol = 2, widths = c(3,3)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.99)
  )
ggsave("plots/slide17_results.png", width = 14.5, height = 8, units = "in")


##Complete table of detections

rejections_muco <- matrix(nrow = 19, ncol = 4)
for(i in 1:length(results_muco_qvalue)){
  ARShift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$ARshift))
  VC_shift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$VCshift))
  Rshift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$Rshift))
  Torshift_q <- na.omit(as.vector(results_muco_qvalue[[i]]$Torshift))
  
  df <- data.frame(ARShift_q, VC_shift_q, Rshift_q, Torshift_q)
  rejections_muco[i,] = colSums(df <= 0.05)
  
}
colnames(rejections_muco) <- c("ARShift", "VCShift", "Rshift", "Torshift")
rejections_muco
df <- as.data.frame(rejections_muco)
df$Slide <- 1:nrow(df)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = c(ARShift, VCShift, Rshift, Torshift),
    names_to = "Method",
    values_to = "Detections"
  )
pastel_cols <- c(
  "ARShift" = "#F4A6A6",
  "VCShift"  = "#B7E4C7",
  "Rshift"   = "#CDB4DB",
  "Torshift" = "#FFE5A5"
)
# Grouped barplot
pm <- ggplot(df_long, aes(x = factor(Slide), y = Detections, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = pastel_cols) +
  coord_cartesian(ylim = c(0, 55)) +
  labs(x = "Slide", y = "Number of detections", fill = "Method") +
  ggtitle("Donors with mucositis")
pm
rejections_healthy <- matrix(nrow = length(results_healthy_qvalue), ncol = 4)
for(i in 1:length(results_healthy_qvalue)){
  ARShift_q <- na.omit(as.vector(results_healthy_qvalue[[i]]$ARshift))
  VC_shift_q <- na.omit(as.vector(results_healthy_qvalue[[i]]$VCshift))
  Rshift_q <- na.omit(as.vector(results_healthy_qvalue[[i]]$Rshift))
  Torshift_q <- na.omit(as.vector(results_healthy_qvalue[[i]]$Torshift))
  
  df <- data.frame(ARShift_q, VC_shift_q, Rshift_q, Torshift_q)
  rejections_healthy[i,] = colSums(df <= 0.05)
}
colnames(rejections_healthy) <- c("ARShift", "VCShift", "Rshift", "Torshift")
rejections_healthy

df <- as.data.frame(rejections_healthy)
df$Slide <- 1:nrow(df)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = c(ARShift, VCShift, Rshift, Torshift),
    names_to = "Method",
    values_to = "Detections"
  )
pastel_cols <- c(
  "ARShift" = "#F4A6A6",
  "VCShift"  = "#B7E4C7",
  "Rshift"   = "#CDB4DB",
  "Torshift" = "#FFE5A5"
)
# Grouped barplot
ph <- ggplot(df_long, aes(x = factor(Slide), y = Detections, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = pastel_cols) +
  coord_cartesian(ylim = c(0, 55)) + 
  labs(x = "Slide", y = "Number of detections", fill = "Method") +
  ggtitle("Healthy donors")
ph 


# pm + ph + 
# plot_layout(ncol = 2, widths = c(4,3), guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(
#     plot.tag = element_text(size = 16, face = "bold"),
#     plot.tag.position = c(0.02, 0.99),
#     legend.position = "bottom"
#   )

library(ggpubr)
all_rej <- ggarrange(pm + rremove("ylab"), ph + rremove("ylab"), widths = c(4,3), common.legend = TRUE, legend = "bottom")
annotate_figure(all_rej, left = textGrob("Number of rejections", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
ggsave("plots/all_rejections.png", width = 14, height = 4, units = "in", dpi = 1200)



##setworking directory in ARShift for these plots
results_healthy <- readRDS("real_data_results/resultshealthy.rds")
results_muco <- readRDS("real_data_results/resultsmuco.rds")
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
df_p <- df_p %>% select(slide, group, pair, pvalue)
