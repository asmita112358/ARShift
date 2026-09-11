##violin and Kcross plots for differentially co-clustering taxa

library(ggrepel)
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
plot_Kcross <- function(ppp_healthy, ppp_muco, taxon1, taxon2){
  Kcross_healthy <- lapply(ppp_healthy, function(x) Kcross(x, i = taxon1, j = taxon2, correction = "isotropic", rmax = r_max)$iso)
  Kcross_muco <- lapply(ppp_muco, function(x) Kcross(x, i = taxon1, j = taxon2, correction = "isotropic", rmax = r_max)$iso)
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


df_muco <- extract_pairwise_results(results_muco_qvalue, group_label = "Mucositis", value_name = "qvalue")
df_healthy <- extract_pairwise_results(results_healthy_qvalue, group_label = "Healthy", value_name = "qvalue")
df <- bind_rows(df_healthy, df_muco)
taxa1 <- c("Actinomyces", "Campylobacter", "Porphyromonas", "Lautropia")
taxa2 <- c("Gemella", "Gemella", "Treponema", "Veillonella")

df %>% filter(pair %in% paste0(taxa1, " vs ", taxa2)) 

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

for (k in 1:4) {
  
  taxon1 <- taxa1[k]
  taxon2 <- taxa2[k]
  
  plot_list[[k]] <- plot_one_pair_log(
    df = df,
    taxon1 = taxon1,
    taxon2 = taxon2,
    value_col = "qvalue",
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
 
}
ggpubr::ggarrange(plotlist = plot_list,  nrow = 4,ncol = 1, common.legend = TRUE, legend = "top", byrow = FALSE)

ggsave("real_data_results2/post-hoc-violins.png", width = 4, height = 10, units = "in", dpi = 1200)


obj <- Kcross(ppp_muco[[6]], i = taxon1, j = taxon2)
