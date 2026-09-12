##Post hoc analysis of between group differences based on p-values
library(spatstat)
library(ggplot2)
library(here)
library(stringr)
library(mcprogress)
library(ggpp)
library(phyloseq)
library(dplyr)
library(qvalue)
library(tidyr)


#qvalue for healthy slides

#read pvalue data
results_healthy_pvalue = list()
for(i in 1:14){
  obj <- readRDS(file = paste0("real_data_results/resultshealthy_slide_",i,".rds"))
  results_healthy_pvalue[[i]] = obj$ARshift
}
taxa_pair <- c()

qvalues_healthy <- matrix(nrow = 91, ncol = 14)
counter = 1
for(i in 1:13){
  for(j in (i+1):14){
   #print(c(i,j))
    taxa_pair <- c(taxa_pair, paste(rownames(results_healthy_pvalue[[1]])[i], colnames(results_healthy_pvalue[[1]])[j], sep = "_"))
    pvalue_healthy = c()
    for(k in 1:14){
      pvalue_healthy = c(pvalue_healthy, results_healthy_pvalue[[k]][i,j])
    }
    qvalues_healthy[counter, ] <- p.adjust(pvalue_healthy, method = "BH")
    counter = counter + 1
    }
}
qvalues_healthy_df <- data.frame(taxa_pair, qvalues_healthy)

taxa_pair = c()
qvalues_muco <- matrix(nrow = 91, ncol = 19)
counter = 1
for(i in 1:13){
  for(j in (i+1):14){
   
    taxa_pair <- c(taxa_pair, paste(rownames(results_healthy_pvalue[[1]])[i], colnames(results_healthy_pvalue[[1]])[j], sep = "_"))
    pvalue_muco = c()
    for(k in 1:19){
      obj <- readRDS(file = paste0("real_data_results/resultsmuco_slide_",k,".rds"))
      pvalue_muco = c(pvalue_muco, obj$ARshift[i,j])
    }
    qvalues_muco[counter, ] <- p.adjust(pvalue_muco, method = "BH")
    counter = counter + 1
  }
  
  
}

qvalues_muco_df <- data.frame(taxa_pair, qvalues_muco)

##Post hoc 1: Fisher's test between the number of significant pairs in mucositis and healthy groups
#get rejections based on qvalues

rej_muco <- rowSums(qvalues_muco <= 0.05)
names(rej_muco) <- taxa_pair

rej_healthy <- rowSums(qvalues_healthy <= 0.05)
names(rej_healthy) <- taxa_pair


##Fisher's test between the number of significant pairs in mucositis and healthy groups

p_diffs <- c()
for(i in 1:91){
  TAB <- rbind(c(rej_muco[i], 19 - rej_muco[i]), c(rej_healthy[i], 14 - rej_healthy[i]))
  obj <- fisher.test(TAB)
  p_diffs[i] <- obj$p.value
}
names(p_diffs) <- taxa_pair
qval_diffs <- p.adjust(p_diffs, method = "BH")
names(qval_diffs) <- taxa_pair

df <- data.frame(taxa_pair, rej_muco, rej_healthy, p_diffs, qval_diffs)
saveRDS(df, "real_data_results/fisher_posthoc_results.rds")


##post hoc 2: ks test between the qvalues from mucositis and healthy groups

# ks.p <- c()
# for(i in 1:91){
#   obj <- ks.test(qvalues_muco[i,], qvalues_healthy[i,])
#   ks.p[i] <- obj$p.value
# }
# names(ks.p) <- taxa_pair
# qval_ks <- p.adjust(ks.p, method = "BH")
# names(qval_ks) <- taxa_pair
# 
# df <- data.frame(taxa_pair, rej_muco, rej_healthy, ks.p, qval_ks, p_diffs, qval_diffs)
# df %>% filter(qval_ks <= 0.05) %>% filter(p_diffs <= 0.05) %>% arrange(qval_ks)


##violin plots for the two significant pairs

df <- df %>% filter(qval_diffs <= 0.1) %>% arrange(qval_diffs)

sig_pairs <- df$taxa_pair
df_q <- qvalues_muco_df %>% filter(taxa_pair %in% sig_pairs) %>% pivot_longer(cols = -taxa_pair, names_to = "slide", values_to = "qvalue") %>%
  mutate(group = "Mucositis")
df_q$slide <- as.numeric(gsub("X", "", df_q$slide))

df_h<- qvalues_healthy_df %>% filter(taxa_pair %in% sig_pairs) %>% pivot_longer(cols = -taxa_pair, names_to = "slide", values_to = "qvalue") %>%
  mutate(group = "Healthy")

df_h$slide <- as.numeric(gsub("X", "", df_h$slide))

df <- rbind(df_q, df_h)

df$qvalue <- -log10(df$qvalue)
df$significant <- df$qvalue > -log10(0.05)
df$xjit = jitter(as.numeric(factor(df$group)), amount = 0.1)

p1 <- ggplot(df, aes(x = group, y = qvalue, fill = group)) +
  geom_violin(alpha = 0.5, color = "black") +
  geom_point(
    aes(x = xjit, color = group),
    size = 2.2,
    alpha = 1,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = df %>% filter(significant),
    aes(
      x = xjit,
      label = .data[["slide"]]
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
  )+
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray30"
  )+
  scale_fill_manual(values = c("Healthy" = "navy",
                               "Mucositis" = "darkred")) +
  scale_color_manual(values = c("Healthy" = "navy",
                                "Mucositis" = "darkred")) +
  
  facet_wrap(~ taxa_pair, nrow = 2, labeller = as_labeller(
    c("Actinomyces_Gemella" = " Actinomyces - Gemella",
      "Lautropia_Veillonella" = " Lautropia - Veillonella")
  )) +
  theme_minimal() +
  theme(strip.text = element_text(size = 14), axis.title = element_text(size = 13), axis.text = element_text(size = 12)) +
  labs(y = "-log10(BH q-value)", x = "Group") +
  theme(legend.position = "none")
p1
ggsave(p1, file = "plots/violin_plot_posthoc.png", width = 4, height = 5, units = "in", dpi = 1200)


##Plot ppp healthy 8,13 and ppp_muco 12,18 for taxa_subset
taxa1 <- c("Actinomyces",  "Lautropia")
taxa2 <- c("Gemella", "Veillonella")
taxa_subset <- unique(c(taxa1, taxa2))


# Match each original taxon label to the color used in the first geom_point layer

taxon_colors <- c(
  "Actinomyces"   = "#F8766D",  # salmon
  "Gemella"       = "#53B400",  # green
  "Lautropia"     = "#001F3F",  # deep navy -- changed
  "Veillonella"   = "#FB61D7"   # pink
)

df <- data.frame(ppp_healthy[[8]])
df <- df[df$marks %in% taxa_subset, ]
p1 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = taxon_colors, drop = FALSE) +
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))+
  labs(title = "Peri Implant Health, 8") 
p1

df <- data.frame(ppp_healthy[[13]])
df <- df[df$marks %in% taxa_subset, ]
p2 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = taxon_colors, drop = FALSE) +
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))+
  labs(title = "Peri Implant Health, 13")

p2

df <- data.frame(ppp_muco[[12]])
df <- df[df$marks %in% taxa_subset, ]
p3 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = taxon_colors, drop = FALSE) +
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))+
  labs(title = "Peri Implant Mucositis, 12")
p3

df <- data.frame(ppp_muco[[5]])
df <- df[df$marks %in% taxa_subset, ]
p4 <- ggplot(df, aes(x = x, y = y, color = marks))+
  geom_point(size = 0.5) +
  scale_color_manual(values = taxon_colors, drop = FALSE) +
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))+
  labs(title = "Peri Implant Mucositis, 5")
p4


library(patchwork)

p1 + p2 + p3 + p4 + plot_layout(ncol = 2, guides = "collect") +plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 13, face = "bold"),
                                                                                                         plot.tag.position = c(0.02, 0.97), legend.position = "bottom", 
                                                                                                         legend.title = element_blank(), legend.text = element_text(size = 12))

ggsave("plots/exampleslides.png", width = 5.5, height = 6, units = "in", dpi = 1200)
#plot Kcross
r_max <- 370
plot_Kcross <- function(ppp_healthy, ppp_muco, taxon1, taxon2){
  Kcross_healthy <- lapply(ppp_healthy, function(x) Kcross(x, i = taxon1, j = taxon2, correction = "isotropic", rmax = r_max)$iso)
  Kcross_muco <- lapply(ppp_muco, function(x) Kcross(x, i = taxon1, j = taxon2, correction = "isotropic", rmax = r_max)$iso)
  r <- Kcross(ppp_healthy[[1]], i = taxon1, j = taxon2, correction = "border", rmax = r_max)$r
  group_label <- c(rep("Mucositis", length(ppp_muco)), rep("Healthy", length(ppp_healthy)))
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
  
  d_summary <- d_long %>%
    group_by(group_label, r) %>%
    summarise(
      mean_d = mean(d, na.rm = TRUE),
      lower  = quantile(d, probs = 0.025, na.rm = TRUE),
      upper  = quantile(d, probs = 0.975, na.rm = TRUE),
      n      = sum(!is.na(d)),
      .groups = "drop"
    )
  
  
  ggplot() +
    # Quantile bands: same group colors, transparent fill
    geom_ribbon(
      data = d_summary,
      aes(
        x = r,
        ymin = lower,
        ymax = upper,
        fill = group_label,
        group = group_label
      ),
      alpha = 0.15,
      colour = NA
    ) +
    
    # Individual slide curves
    geom_line(
      data = d_long,
      aes(
        x = r,
        y = d,
        color = group_label,
        linetype = group_label,
        group = sample_id
      ),
      alpha = 1,
      linewidth = 0.5
    ) +
    scale_color_manual(values = c("navy",
                                  "darkred")) +
    scale_fill_manual(values = c("navy",
                                 "darkred")) +
    scale_linetype_manual(
      values = c(
        "Healthy" = "solid",
        "Mucositis" = "dashed"
      )
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6)+
    theme_minimal() +
    labs(
      y = expression(L(r) - r),
    )
  
}
p3 <- plot_Kcross(ppp_healthy, ppp_muco, "Actinomyces", "Gemella") + 
  ggtitle("Actinomyces - Gemella") +
  ylim(-200,500)+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), axis.title = element_text(size = 13), axis.text = element_text(size = 12)) 
p4 <- plot_Kcross(ppp_healthy, ppp_muco, "Lautropia", "Veillonella") + 
  ggtitle("Lautropia - Veillonella") +
  ylim(-200,500)+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), axis.title = element_text(size = 13), axis.text = element_text(size = 12))
library(patchwork)
pK <- p3 + p4 + plot_layout(ncol = 2, widths = c(1,1), axes = "collect") 
pK


win_area <- c()
for(i in 1:19){
  win_area[i] <- area.owin(ppp_muco[[i]]$window)
}
win_area_healthy <- c()
for(i in 1:14){
  win_area_healthy[i] <- area.owin(ppp_healthy[[i]]$window)
}
##overlapped histogram, color coded by healthy vs muco
hist_df <- data.frame(area = c(win_area, win_area_healthy), group = c(rep("Mucositis", 19), rep("Healthy", 14)))

hist(win_area, breaks = 20, col = "darkred",   ylim = c(0,9e-08),freq = FALSE)
hist(win_area_healthy, breaks = 20, col = "navy", ylim = c(0,9e-08),add = TRUE, freq = FALSE)



