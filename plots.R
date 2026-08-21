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

data7 <- generate_rCauchy(kappa = c(20,15), mu = c(30, 20), scale = c(0.02, 0.01))
df <- data.frame(data7)
p7 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S7") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p7

data8 <- generate_rCauchy(kappa = c(20,15), mu = c(30, 20), scale = c(0.05, 0.04))
df <- data.frame(data8)
p8 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S8") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p8


data9 <- generate_rCauchy(kappa = c(20,15), mu = c(30, 20), scale = c(0.15, 0.1))
df <- data.frame(data9)
p9 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S9") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p9


p_all <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3,nrow = 3, legend = "bottom",common.legend = TRUE)
p_all
ggsave("plots/sim_examples.png", width = 6, height = 6,units = "in", dpi = 1200)



###Plots for inhomogeneous simulations

data10 <- generate_inhom_thomas(win,  scale = c(0.05, 0.04), dependant = FALSE)
df <- data.frame(data10)
p10 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S10") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p10


data11 <- generate_inhom_LGCP(win, mu = 6, corrfun = "matern", dependant = FALSE)
df <- data.frame(data11)
p11 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S11") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p11

data12 <- generate_inhom_cauchy(
     win,
       scale = c(0.05, 0.04),
       dependant = FALSE
   )
df <- data.frame(data12)
p12 <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_minimal() +
  labs(title = "S12") +
  scale_color_manual(values = c("coral3", "darkolivegreen"))
p12

p_supp <- ggarrange(p10,p11,p12, ncol = 3,nrow = 1, legend = "bottom",common.legend = TRUE)
p_supp

ggsave("plots/inhom_sim_examples.png", width = 6, height = 2.5,units = "in", dpi = 1200)
