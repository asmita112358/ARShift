##Powerplots
library(tidyverse)

# Data
df <- tribble(
  ~window, ~P,   ~RShift, ~VCShift, ~ARShift, ~TorShift,
  "a", "P1", 0.57, 0.26, 0.94, 0.96,
  "b", "P1", 0.345, 0.077, 0.716, 0.149,
  "a", "P2", 0.714, 0.272, 0.974, 0.968,
  "b", "P2", 0.466, 0.174, 0.836, 0.212,
  "a", "P3", 0.612, 0.272, 0.94, 0.944,
  "b", "P3", 0.38, 0.158, 0.724, 0.086,
  "a", "P4", 0.265, 0.089, 0.613, 0.552,
  "b", "P4", 0.216, 0.092, 0.482, 0.104,
  "a", "P5", 0.454, 0.206, 0.796, 0.752,
  "b", "P5", 0.272, 0.142, 0.592, 0.164,
  "a", "P6", 0.307, 0.166, 0.632, 0.673,
  "b", "P6", 0.182, 0.098, 0.42,  0.108,
  "a", "P7", 0.655, 0.094, 0.952, 0.961,
  "b", "P7", 0.362, 0.065, 0.731, 0.083,
  "a", "P8", 0.621, 0.105, 0.928, 0.949,
  "b", "P8", 0.348, 0.069, 0.705, 0.042,
  "a", "P9", 0.254, 0.062, 0.577, 0.722,
  "b", "P9", 0.156, 0.045, 0.392, 0.016
)

# Put P1-P9 into 3 row groups
df2 <- df %>%
  mutate(
    P = factor(P, levels = paste0("P", 1:9)),
    group = case_when(
      P %in% c("P1", "P2", "P3") ~ "Thomas",
      P %in% c("P4", "P5", "P6") ~ "LGCP",
      TRUE ~ "Neyman Scott"
    ),
    group = factor(group, levels = c("Thomas", "LGCP", "Neyman Scott")),
    window = factor(window, levels = c("a", "b"),
                    labels = c("a (Square windows)", "b (Irregular windows)"))
  ) %>%
  pivot_longer(
    cols = c(RShift, VCShift, ARShift, TorShift),
    names_to = "method",
    values_to = "power"
  ) %>%
  mutate(
    method = factor(method, levels = c("RShift", "VCShift", "ARShift", "TorShift"))
  )

pastel_cols <- c(
  "ARShift" = "#F4A6A6",
  "VCShift"  = "#B7E4C7",
  "RShift"   = "#CDB4DB",
  "TorShift" = "#FFE5A5"
)

# Plot
ggplot(df2, aes(x = P, y = power, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
  facet_grid(window ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pastel_cols) +
  labs(
    x = NULL,
    y = "Empirical power",
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
ggsave("plots/powerplots.png", height = 7, width = 7, units = "in")
