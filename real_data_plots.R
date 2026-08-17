##Plots for real data analysis
#Run until line 85 in dental_microbiome.R
library(pals)
library(ggvenn)
library(ggVennDiagram)
library(ggplot2)
library(tidyr)
library(dplyr)
#Example plot, muco i
i = 9
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome_twosample")
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))

path <- muco.list[i]
img <- read.csv(file=paste0("centroid_sciname_tables/",path))
img.df <- process_df(img)
ppp_img=ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
                    marks=as.factor(img.df$sciname))
#ppp_img <- rthin(ppp_img, 0.2)
ppp_df <- data.frame(ppp_img)
colnames(ppp_df)[3] <- "Taxa"
# Custom 18-color distinct pastel palette
my_18_colors <- c(
  "#FFB4B4", "#FFF5C2", "#C0EEF2", "#B5DEFF", "#D5B4B4", "#E1BEE7",
  "#F9DBBD", "#FFE0AC", "#E8F9FD", "#E3BEC6", "#E4DCCF", "#F3E8EE",
  "#F67280", "#F9B208", "#97DECE", "#6CAEED", "#A7BBC7", "#BB9CC0"
)
muco1 <- ggplot(ppp_df, aes(x = x, y = y, color = Taxa)) +
  geom_point(size = 1) +
  scale_color_manual(values = my_18_colors) + 
  coord_equal() + # Important for spatial data
  ggtitle("Slide 9, mucositis")+
  theme_minimal() +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))

muco1
dir <- "~/Library/CloudStorage/OneDrive-JohnsHopkins/spatial_assoc_test/ARShift/real_data_results/"
ggsave(paste0(dir,"muco9.png"), muco5, width = 5, height = 6, units = "in")
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift")
results_healthy_qvalue <- readRDS("real_data_results/Qresultshealthy.rds")
results_muco_qvalue <- readRDS("real_data_results/Qresultsmuco.rds")

##Venn diagram for slide 9, mucositis

results_healthy_qvalue <- readRDS(paste0(dir,"Qresultshealthy.rds"))
results_muco_qvalue <- readRDS(paste0(dir,"Qresultsmuco.rds"))

i = 9
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
vplot <- ggVennDiagram(x) + scale_fill_gradient(low = "lightblue", high = "darkblue") + 
  ggtitle("Taxa pairs detected")+theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
vplot
muco1 + vplot +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.99)
  )
ggsave("real_data_results/slide9_results.png", width = 11, height = 6.5, units = "in")


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
  labs(x = "Slide", y = "Number of detections", fill = "Method") +
  ggtitle("Healthy donors")
ph 
pm + ph + 
plot_layout(ncol = 2, widths = c(4,3), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.99),
    legend.position = "bottom"
  )

ggsave("real_data_results/all_tests.png", width = 14, height = 4, units = "in", dpi = 1200)

