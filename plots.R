#Plots for simulation settings
library(ggplot2)
library(ggpubr)
data1 <-  generate_Thomas(win, kappa = c(12,10), mu = c(30, 20), scale = c(0.02, 0.01), dependant = FALSE)
df <- data.frame(data1)
p1 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S1") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))

p1
data2 <-  generate_Thomas(win, kappa = c(12,10), mu = c(30, 20), scale = c(0.05, 0.04), dependant = FALSE)
df <- data.frame(data2)
p2 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S2") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p2

data3 <- generate_Thomas(win, kappa = c(12,10), mu = c(30, 20), scale = c(0.15, 0.1), dependant = FALSE)
df <- data.frame(data3)
p3 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S3") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p3

data4 <- generate_LGCP(win, corrfun = "gauss", mu = 4.5, scale = 0.1, dependant = FALSE)
df <- data.frame(data4)
p4 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S4") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p4


data5 <- generate_LGCP(win, corrfun = "gauss", mu = 4.5, scale = 0.3, dependant = FALSE)
df <- data.frame(data5)
p5 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S5") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p5

data6 <- generate_LGCP(win, corrfun = "gauss", mu = 4.5, scale = 0.5, dependant = FALSE)
df <- data.frame(data6)
p6 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S6") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p6

p_all <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3,nrow = 2, legend = "bottom",common.legend = TRUE)
ggsave("plots/sim_examples.png", width = 6, height = 4,units = "in", dpi = 1200)
