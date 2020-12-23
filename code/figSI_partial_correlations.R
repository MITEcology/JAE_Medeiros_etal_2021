# Code for figures 2b and 3b with partial correlations between eigenvalues and distances 
# to border/vertex for multiple matrices, number of species and interaction types

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(ppcor)) {install.packages("ppcor"); library(ppcor)}
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

# settings ------------------------------
# whether to save plots
save_plots <- FALSE
# number of matrices
n <- 100
# interaction strength
strength <- "low"
# connectance
connectance <- 1
# which norm to use
norm <- "l2"
# color palette
pal <- c(pal_material("red")(9)[8], pal_material("blue")(9)[8], pal_material("yellow")(9)[8])

# compute partial correlations for all numbers of species and interaction types ------------------------------
# number of species
S_all <- c(3, 4, 5)
# number of fixed points to sample
fixed_points_all <- 100 * 2^(S_all-2)
# type of interaction
interaction_all <- c("negative", "both", "positive")
# compute partial correlations
plot_df <- data.frame()
for (i in 1:length(S_all)) {
  for (j in 1:length(interaction_all)) {
    # load eigenvalues and distances for fixed points
    load(file = paste("results/random_matrices/", S_all[i], "sp/eigenvalues_evenness_distances_", n, "_random_matrices_", 
                      fixed_points_all[i], "_fixed_points_", norm, "_norm_", strength, "_strength_", 
                      interaction_all[j], "_sign", ".RData", sep = ""))
    # load matrices
    load(paste("data/random_matrices/", S_all[i], "sp/", n, "_random_matrices_",
               strength, "_strength_", connectance, "_connectance_", 
               interaction_all[j], "_sign.RData", sep = ""))
    # removing equilibrium states with positive largest eigenvalue
    df <- df[df$lambda1 < 0, ]
    # obtaining indexes of remaining matrices
    remaining_mat <- unique(df$matrix)
    matrix_list <- matrix_list[remaining_mat]
    # compute partial recovery (second smallest eigenvalue)
    if (S_all[i] == 3) {
      df$lambda_penultimate <- df$lambda2
    }
    if (S_all[i] == 4) {
      df$lambda_penultimate <- df$lambda3
    }
    if (S_all[i] == 5) {
      df$lambda_penultimate <- df$lambda4
    }
    # compute most abundant species and its total interaction effect 
    df$matrix <- factor(df$matrix)
    df$most_abundant <- apply(df[ , 1:S_all[i]], 1, function(x) which.max(x))
    mean_int_strength <- lapply(matrix_list, function(mat) (apply(mat, 2, sum) - diag(mat)) / (nrow(mat) - 1))
    df_mean_int_strength <- as.data.frame(matrix(unlist(mean_int_strength), nrow = length(mean_int_strength), 
                                                 byrow = TRUE))
    names(df_mean_int_strength) <- as.character(1:S_all[i])
    df_mean_int_strength <- df_mean_int_strength[rep(1:length(matrix_list), 
                                                     tapply(df$N1, df$matrix, length)), ]
    df$int_strength_most_abund <- df_mean_int_strength[cbind(1:nrow(df_mean_int_strength), 
                                                             df$most_abundant)]
    # compute partial correlations
    lambda1_min_dist2border_r <- c()
    lambda_penultimate_min_dist2vertex_r <- c()
    lambda1_lambda_penultimate <- c()
    dist2border_dist2vertex <- c()
    for (l in 1:length(remaining_mat)) {
      lambda1_min_dist2border_r[l] <- pcor.test(subset(df, matrix == remaining_mat[l])$lambda1,
                                                subset(df, matrix == remaining_mat[l])$min_distance2border_r,
                                                subset(df, matrix == remaining_mat[l])$most_abundant)$estimate
      lambda_penultimate_min_dist2vertex_r[l] <- pcor.test(subset(df, matrix == remaining_mat[l])$lambda_penultimate,
                                                           subset(df, matrix == remaining_mat[l])$min_distance2vertex_r,
                                                           subset(df, matrix == remaining_mat[l])$most_abundant)$estimate
      lambda1_lambda_penultimate[l] <- pcor.test(subset(df, matrix == remaining_mat[l])$lambda1,
                                                 subset(df, matrix == remaining_mat[l])$lambda_penultimate,
                                                 rank(subset(df, matrix == remaining_mat[l])$lambda1))$estimate
      dist2border_dist2vertex[l] <- pcor.test(subset(df, matrix == remaining_mat[l])$min_distance2border_r,
                                              subset(df, matrix == remaining_mat[l])$min_distance2vertex_r,
                                              rank(subset(df, matrix == remaining_mat[l])$min_distance2vertex_r))$estimate
    }
    # add current data frame to full data frame
    curr_df <- data.frame(S = rep(S_all[i], length(remaining_mat)),
                          interaction = rep(interaction_all[j], length(remaining_mat)),
                          strength = rep(strength, length(remaining_mat)),
                          lambda1_min_dist2border_r, lambda_penultimate_min_dist2vertex_r,
                          lambda1_lambda_penultimate, dist2border_dist2vertex)
    plot_df <- rbind(plot_df, curr_df)
  }
}

# plot of largest eigenvalue and minimum distance to border ------------------------------
# change factor level order
plot_df$interaction <- as.character(plot_df$interaction)
plot_df$interaction[plot_df$interaction == "negative"] <- "Competition system"
plot_df$interaction[plot_df$interaction == "both"] <- "Antagonistic system"
plot_df$interaction[plot_df$interaction == "positive"] <- "Mutualistic system"
plot_df$interaction <- factor(plot_df$interaction, levels = c("Competition system", "Mutualistic system",
                                                              "Antagonistic system"))
# partial correlation between full and partial resistance
fig <- ggplot(data = plot_df, aes(x = as.factor(S), y = dist2border_dist2vertex, color = interaction)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_color_manual(values = pal) +
  facet_wrap(~interaction, ncol = length(interaction_all)) +
  xlab(expression(paste("Number of species ", (S)))) +
  ylab(expression(atop("Partial correlation between full and partial resistance",
                       ~ paste("(", rho, "(", min(group("{", d[b], "}")), ", ", min(group("{", d[v], "}")), 
                               "| rank of ", min(group("{", d[v], "}")), "))", sep = "")))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/figSI/partial_correlations/partial_cor_border_vertex.pdf", sep = ""),
         fig, width = 38, height = 17, units = "cm")
}
# partial correlation between full and partial recovery
fig <- ggplot(data = plot_df, aes(x = as.factor(S), y = lambda1_lambda_penultimate, color = interaction)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_color_manual(values = pal) +
  facet_wrap(~interaction, ncol = length(interaction_all)) +
  xlab(expression(paste("Number of species ", (S)))) +
  ylab(expression(atop("Partial correlation between full and partial recovery",
                       ~ paste("(", rho, "(", lambda[1], ", ", lambda[S-1], 
                               "| rank of ", lambda[1], "))", sep = "")))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/figSI/partial_correlations/partial_cor_eigen1_eigen_penultimate.pdf", sep = ""),
         fig, width = 38, height = 17, units = "cm")
}
# partial correlation between full recovery and full resistance
fig <- ggplot(data = plot_df, aes(x = as.factor(S), y = lambda1_min_dist2border_r, color = interaction)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_color_manual(values = pal) +
  facet_wrap(~interaction, ncol = length(interaction_all)) +
  xlab(expression(paste("Number of species ", (S)))) +
  ylab(expression(atop("Partial correlation between full recovery and full resistance",
                       ~ paste("(", rho, "(", lambda[1], ", ", min(group("{", d[b], "}")), 
                               "| most abundant species))", sep = "")))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/figSI/partial_correlations/partial_cor_eigen1_border.pdf", sep = ""),
         fig, width = 38, height = 17, units = "cm")
}
# partial correlation between partial recovery and partial resistance
fig <- ggplot(data = plot_df, aes(x = as.factor(S), y = lambda_penultimate_min_dist2vertex_r, color = interaction)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_color_manual(values = pal) +
  facet_wrap(~interaction, ncol = length(interaction_all)) +
  xlab(expression(paste("Number of species ", (S)))) +
  ylab(expression(atop("Partial correlation between partial recovery and partial resistance",
                       ~ paste("(", rho, "(", lambda[S-1], ", ", min(group("{", d[v], "}")), 
                               "| most abundant species))", sep = "")))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/figSI/partial_correlations/partial_cor_eigen_penultimate_vertex.pdf", sep = ""),
         fig, width = 38, height = 17, units = "cm")
}
