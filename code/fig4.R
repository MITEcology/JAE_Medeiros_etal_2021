# Code for figure 4: correlations between eigenvalues and distances to border/vertex 
# for the 3-species experimental systems of Friedman et al (2017) Nat Ecol Evol

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)}
if(!require(ppcor)) {install.packages("ppcor"); library(ppcor)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

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
# empirical largest and second smallest eigenvalue
large_eigen_emp <- unique(df$large_eigen_emp)
second_large_eigen_emp <- unique(df$second_large_eigen_emp)
# empirical distances to border and vertex
dist2border_emp <- unique(df$min_distance2border_r_emp)
dist2vertex_emp <- unique(df$min_distance2vertex_r_emp)

# plot full recovery vs full resistance for a single example ------------------------------
i = 6
sub_df <- subset(df, matrix == i)[ , c("large_eigen", "min_distance2border_r")]
print(cor(sub_df$large_eigen, sub_df$min_distance2border_r))
sub_df <- rbind(sub_df, c(large_eigen_emp[i], dist2border_emp[i]))
sub_df$type = c(rep("Random", nrow(sub_df)-1), "Experimental")
# plot
fig_4A <- ggplot(data = sub_df, aes(x = min_distance2border_r, y = large_eigen,
                                    shape = type, color = type, fill = type, alpha = type, size = type)) +
  geom_point() +
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = c(c(pal_material("orange")(9)[8], "gray"))) +
  scale_alpha_manual(values = c(1, 0.7)) +
  scale_size_manual(values = c(6, 3)) +
  ylab(expression(paste("Full recovery ", (lambda[1])))) +
  xlab(expression(paste("Full resistance ", (min(group("{", d[b], "}")))))) +
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
  ggsave("figs/fig4/fig4_A.pdf", fig_4A, width = 16, height = 16, units = "cm")
}

# plot partial recovery vs partial resistance for a single example ------------------------------
i = 6
sub_df <- subset(df, matrix == i)[ , c("second_large_eigen", "min_distance2vertex_r")]
print(cor(sub_df$second_large_eigen, sub_df$min_distance2vertex_r))
sub_df <- rbind(sub_df, c(second_large_eigen_emp[i], dist2vertex_emp[i]))
sub_df$type = c(rep("Random", nrow(sub_df)-1), "Experimental")
# plot
fig_4B <- ggplot(data = sub_df, aes(x = min_distance2vertex_r, y = second_large_eigen,
                                    shape = type, color = type, fill = type, alpha = type, size = type)) +
  geom_point() +
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = c(pal_material("orange")(9)[8], "gray")) +
  scale_alpha_manual(values = c(1, 0.7)) +
  scale_size_manual(values = c(6, 3)) +
  scale_x_continuous(breaks = c(0, 0.025, 0.05, 0.075)) +
  ylab(expression(paste("Partial recovery ", (lambda[2])))) +
  xlab(expression(paste("Partial resistance ", (min(group("{", d[v], "}")))))) +
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
  ggsave("figs/fig4/fig4_B.pdf", fig_4B, width = 16, height = 16, units = "cm")
}

# compute correlations and partial correlations for all experimental systems ------------------------------
# loading data
A <- as.matrix(read.table("data/friedman2017_data/friedman2017_matrix.txt"))
r <- as.matrix(read.table("data/friedman2017_data/friedman2017_r.txt"))
K <- as.matrix(read.table("data/friedman2017_data/friedman2017_k.txt"))
# transform matrix to r formalism
A <- A * as.numeric(r / K)
# number of species
S <- 3
# index combinations
q <- combn(1:nrow(A), S)
A_list <- list()
r_list <- list()
# extract subsets of matrices and parameters
for (j in 1:ncol(q)) { 
  index <- q[ , j]
  A_list[[j]] <- -A[index, index]
  r_list[[j]] <- r[index ,]
  names(r_list[[j]]) <- rownames(r)[index]
}
# systems to use (only locally stable empirical systems)
matrix_ids <- c(2, 3, 12, 13, 14, 20, 27, 29, 30, 32, 47, 48, 49, 50, 51, 52, 54)
A_list <- A_list[matrix_ids]
r_list <- r_list[matrix_ids]
# compute most abundant species and its total interaction effect 
df$matrix <- factor(df$matrix)
df$most_abundant <- apply(df[ , 1:S], 1, function(x) which.max(x))
mean_int_strength <- lapply(A_list, function(mat) (apply(mat, 2, sum) - diag(mat)) / (nrow(mat) - 1))
df_mean_int_strength <- as.data.frame(matrix(unlist(mean_int_strength), nrow = length(mean_int_strength), 
                                             byrow = TRUE))
names(df_mean_int_strength) <- as.character(1:S)
df_mean_int_strength <- df_mean_int_strength %>% slice(rep(1:n(), each = fixed_points))
df$int_strength_most_abund <- df_mean_int_strength[cbind(1:nrow(df_mean_int_strength), 
                                                         df$most_abundant)]
# compute correlations
plot_df <- data.frame()
eigen1_min_dist2border_r <- c()
eigen2_min_dist2vertex_r <- c()
eigen1_eigen2 <- c()
dist2border_dist2vertex <- c()
for (l in 1:length(A_list)) {
  eigen1_min_dist2border_r[l] <- cor(subset(df, matrix == l)$large_eigen,
                                     subset(df, matrix == l)$min_distance2border_r)
  eigen2_min_dist2vertex_r[l] <- cor(subset(df, matrix == l)$second_large_eigen,
                                     subset(df, matrix == l)$min_distance2vertex_r)
  eigen1_eigen2[l] <- pcor.test(subset(df, matrix == l)$large_eigen,
                                subset(df, matrix == l)$second_large_eigen,
                                rank(subset(df, matrix == l)$large_eigen))$estimate
  dist2border_dist2vertex[l] <- pcor.test(subset(df, matrix == l)$min_distance2border_r,
                                          subset(df, matrix == l)$min_distance2vertex_r,
                                          rank(subset(df, matrix == l)$min_distance2vertex_r))$estimate
}
plot_df <- data.frame(eigen1_min_dist2border_r, eigen2_min_dist2vertex_r,
                      eigen1_eigen2, dist2border_dist2vertex)

# boxplots with correlations for all experimental systems ------------------------------
# plot
plot_df$system <- 1:nrow(plot_df)
plot_df <- gather(plot_df, "correlation", "value", -system)
plot_df$correlation <- factor(plot_df$correlation, 
                              levels = c("eigen1_min_dist2border_r", "eigen2_min_dist2vertex_r",
                                         "eigen1_eigen2", "dist2border_dist2vertex"))
fig_4C <- ggplot(data = plot_df, aes(x = as.factor(correlation), y = value)) +
  geom_boxplot(size = 1.5, outlier.shape = NA, color = "black") +
  geom_jitter(size = 3.6, alpha = 0.6, position = position_jitter(0.2)) +
  xlab("") +
  ylab("Correlation") +
  scale_y_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        legend.position = "none")
if (save_plots) {
  ggsave("figs/fig4/fig4_C.pdf", fig_4C, width = 38, height = 14, units = "cm")
}
