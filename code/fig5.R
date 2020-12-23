# Code for figure 5: relative values of eigenvalues and distances to border/vertex
# for the 3-species experimental systems of Friedman et al (2017) Nat Ecol Evol

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)}
if(!require(ppcor)) {install.packages("ppcor"); library(ppcor)}

# load results ------------------------------
# whether to save plots
save_plots <- FALSE
# number of species
S <- 3
# number of fixed points to sample
fixed_points <- 2000
# which norm to use
norm <- "l2"
# load R data
load(file = paste("results/friedman2017_data/eigenvalues_distances_", fixed_points, "_fixed_points_", 
                  norm, "_norm.RData", sep = ""))
# check if all eigenvalues are negative
all(df$large_eigen < 0)
all(df$large_eigen_emp < 0)

# compute correlations and partial correlations ------------------------------
eigen1_border <- c()
eigen2_vertex <- c()
eigen1_eigen2 <- c()
border_vertex <- c()
for (i in unique(df$matrix)) {
  # data frame for current matrix
  sub_df <- subset(df, matrix == i)
  sub_df$matrix <- factor(sub_df$matrix)
  # compute most abundant species
  sub_df$most_abundant <- apply(sub_df[ , 1:S], 1, function(x) which.max(x))
  # compute correlations
  eigen1_border[i] <- cor(sub_df$large_eigen, sub_df$min_distance2border_r)
  eigen2_vertex[i] <- cor(sub_df$second_large_eigen, sub_df$min_distance2vertex_r)
  eigen1_eigen2[i] <- pcor.test(sub_df$large_eigen, sub_df$second_large_eigen,
                                rank(sub_df$large_eigen))$estimate
  border_vertex[i] <- pcor.test(sub_df$min_distance2border_r, sub_df$min_distance2vertex_r,
                                rank(sub_df$min_distance2vertex_r))$estimate
}
# mean and sd of correlations for Tables S6 and S7
mean(eigen1_border)
sd(eigen1_border)
mean(eigen2_vertex)
sd(eigen2_vertex)
mean(eigen1_eigen2)
sd(eigen1_eigen2)
mean(border_vertex)
sd(border_vertex)
# empirical values of largest and second smallest eigenvalue
large_eigen_emp <- unique(df$large_eigen_emp)
second_large_eigen_emp <- unique(df$second_large_eigen_emp)
# empirical values of distances to border and vertex
dist2border_emp <- unique(df$min_distance2border_r_emp)
dist2vertex_emp <- unique(df$min_distance2vertex_r_emp)
# compute relative indicators
min_large_eigen <- tapply(df$large_eigen, df$matrix, min)
ratio_large_eigen <- large_eigen_emp / min_large_eigen
min_second_large_eigen <- tapply(df$second_large_eigen, df$matrix, min)
ratio_second_large_eigen <- second_large_eigen_emp / min_second_large_eigen
max_dist2border <- tapply(df$min_distance2border_r, df$matrix, max)
ratio_dist2border <- dist2border_emp / max_dist2border
max_dist2vertex <- tapply(df$min_distance2vertex_r, df$matrix, max)
ratio_dist2vertex <- dist2vertex_emp / max_dist2vertex

# plot relationship between eigenvalues for a single example ------------------------------
i = 6
sub_df <- subset(df, matrix == i)[ , c("large_eigen", "second_large_eigen")]
sub_df$ratio_large_eigen <- sub_df$large_eigen / min(sub_df$large_eigen)
sub_df$ratio_second_large_eigen <- sub_df$second_large_eigen / min(sub_df$second_large_eigen)
sub_df <- rbind(sub_df, c(large_eigen_emp[i], second_large_eigen_emp[i],
                          ratio_large_eigen[i], ratio_second_large_eigen[i]))
sub_df$type = c(rep("Random", nrow(sub_df)-1), "Experimental")
# plot
fig_5A <- ggplot(data = sub_df, aes(x = ratio_large_eigen, y = ratio_second_large_eigen,
                                    shape = type, color = type, fill = type, alpha = type, size = type)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = c(c(pal_material("orange")(9)[8], "gray"))) +
  scale_alpha_manual(values = c(1, 0.7)) +
  scale_size_manual(values = c(6, 3)) +
  xlab(expression(paste("Relative full recovery ", (tilde(lambda)[1])))) +
  ylab(expression(paste("Relative partial recovery ", (tilde(lambda)[2])))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.7, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave("figs/fig5/fig5_A.pdf", fig_5A, width = 16, height = 16, units = "cm")
}

# plot relationship between distances for a single example ------------------------------
i = 6
sub_df <- subset(df, matrix == i)[ , c("min_distance2border_r", "min_distance2vertex_r")]
sub_df$ratio_dist2border <- sub_df$min_distance2border_r / max(sub_df$min_distance2border_r)
sub_df$ratio_dist2vertex <- sub_df$min_distance2vertex_r / max(sub_df$min_distance2vertex_r)
sub_df <- rbind(sub_df, c(dist2border_emp[i], dist2vertex_emp[i],
                          ratio_dist2border[i], ratio_dist2vertex[i])) 
sub_df$type = c(rep("Random", nrow(sub_df)-1), "Experimental")
# plot
fig_5B <- ggplot(data = sub_df, aes(x = ratio_dist2border, y = ratio_dist2vertex,
                                    shape = type, color = type, fill = type, alpha = type, size = type)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = c(pal_material("orange")(9)[8], "gray")) +
  scale_alpha_manual(values = c(1, 0.7)) +
  scale_size_manual(values = c(6, 3)) +
  xlab(expression(paste("Relative full resistance ", (min(group("{", tilde(d)[b], "}")))))) +
  ylab(expression(paste("Relative partial resistance ", (min(group("{", tilde(d)[v], "}")))))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.7, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave("figs/fig5/fig5_B.pdf", fig_5B, width = 16, height = 16, units = "cm")
}

# plot relative value of eigenvalues for all 17 systems ------------------------------
plot_df <- data.frame(ratio_large_eigen, ratio_second_large_eigen,
                      ratio_dist2border, ratio_dist2vertex)
# reorder systems
plot_df <- plot_df[order(plot_df$ratio_large_eigen, decreasing = TRUE), ]
rownames(plot_df) <- 1:nrow(plot_df)
# plot
fig_5C <- ggplot(data = plot_df, aes(x = ratio_large_eigen, y = ratio_second_large_eigen,
                                     label = rownames(plot_df))) +
  geom_point(size = 6, shape = 21, fill = pal_material("orange")(9)[8]) +
  geom_text_repel(size = 5, segment.color = NA, point.padding = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab(expression(paste("Relative full recovery ", (tilde(lambda)[1])))) +
  ylab(expression(paste("Relative partial recovery ", (tilde(lambda)[2])))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.7, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave("figs/fig5/fig5_C.pdf", fig_5C, width = 16, height = 16, units = "cm")
}

# plot relative value of distances for all 17 systems ------------------------------
fig_5D <- ggplot(data = plot_df, aes(x = ratio_dist2border, y = ratio_dist2vertex,
                                     label = rownames(plot_df))) +
  geom_point(size = 6, shape = 21, fill = pal_material("orange")(9)[8]) +
  geom_text_repel(size = 5, segment.color = NA, point.padding = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab(expression(paste("Relative full resistance ", (min(group("{", tilde(d)[b], "}")))))) +
  ylab(expression(paste("Relative partial resistance ", (min(group("{", tilde(d)[v], "}")))))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.7, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave("figs/fig5/fig5_D.pdf", fig_5D, width = 16, height = 16, units = "cm")
}
