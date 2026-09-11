##violin and Kcross plots for differentially co-clustering taxa

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

