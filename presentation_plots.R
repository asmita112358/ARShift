##Plots for presentation
library(ggplot2)
library(ggpubr)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/generate_data.R", echo = FALSE)
#Load Irregular window if needed
path <- "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv"
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))

W_new <- W
W_new$bdry[[1]]$x <- scale(W$bdry[[1]]$x, center = 0, scale = 6000)
W_new$bdry[[1]]$y <- scale(W$bdry[[1]]$y, center = 0, scale = 6000)
W_new$xrange <- range(W_new$bdry[[1]]$x)
W_new$yrange <- range(W_new$bdry[[1]]$y)


data <- generate_LGCP(win = W_new, corrfun = "gauss", mu = 6.5, scale = 0.1, dependant = FALSE)
df <- data.frame(data)
p <- ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point() +
  theme_void() + 
  scale_color_manual(
    name = "taxon",
    values = c("1" = "coral3", "2" = "darkolivegreen"),
    labels = c("1" = "Taxa i", "2" = "Taxa j")
  )
p


jump_radius <- incircle(data$window)$r

# Random shift vector from disc centered at origin
shift_pp <- runifdisc(1, radius = jump_radius, centre = c(0, 0))
shift_vec <- c(shift_pp$x, shift_pp$y)

shift_vec
taxa_i <- data[marks(data) == "1"]
taxa_j <- data[marks(data) == "2"]
taxa_j_shifted <- shift(taxa_j, vec = shift_vec)
W_original <- Window(data)
W_shifted  <- shift(W_original, vec = shift_vec)

W_overlap <- intersect.owin(W_original, W_shifted)
owin_to_df <- function(W) {
  Wp <- as.polygonal(W)
  bd <- Wp$bdry
  
  do.call(rbind, lapply(seq_along(bd), function(i) {
    data.frame(
      x = bd[[i]]$x,
      y = bd[[i]]$y,
      id = i
    )
  }))
}

df_i <- as.data.frame(taxa_i)
df_i$taxon <- "Taxa i"

df_j_shifted <- as.data.frame(taxa_j_shifted)
df_j_shifted$taxon <- "Taxa j shifted"

df_shifted <- rbind(df_i, df_j_shifted)

overlap_df <- owin_to_df(W_overlap)
p1 <- ggplot(df_shifted, aes(x = x, y = y, color = taxon)) +
  geom_point() +
  geom_path(
    data = overlap_df,
    aes(x = x, y = y, group = id),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Taxa i" = "coral3",
      "Taxa j shifted" = "darkolivegreen"
    )
  ) +
  coord_equal() +
  theme_void()

p1

taxa_i_overlap <- taxa_i[W_overlap]
taxa_j_shifted_overlap <- taxa_j_shifted[W_overlap]
df_i_overlap <- as.data.frame(taxa_i_overlap)
df_i_overlap$taxon <- "Taxa i"

df_j_overlap <- as.data.frame(taxa_j_shifted_overlap)
df_j_overlap$taxon <- "Taxa j shifted"

df_overlap_points <- rbind(df_i_overlap, df_j_overlap)
p2 <- ggplot(df_overlap_points, aes(x = x, y = y, color = taxon)) +
  geom_point() +
  geom_path(
    data = overlap_df,
    aes(x = x, y = y, group = id),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Taxa i" = "coral3",
      "Taxa j shifted" = "darkolivegreen"
    )
  ) +
  coord_equal() +
  theme_void()

p2
