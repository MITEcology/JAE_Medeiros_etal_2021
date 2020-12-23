# Code for figures 2 and 3: correlations between eigenvalues and distances 
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
# interaction strength ("low" or "high")
strength <- "low"
# connectance
connectance <- 1
# which norm to use
norm <- "l2"
# color palette
pal <- c(pal_material("red")(9)[8], pal_material("blue")(9)[8], pal_material("yellow")(9)[8])

# example scatterplot for 3 species ------------------------------
# number of species
S <- 3
# number of fixed points to sample
fixed_points <- 100 * 2^(S-2)
# types of interaction
interaction_all <- c("negative", "both", "positive")
# example matrix for low strength
if (strength == "low")
  matrix_id = 30
# example matrix for high strength
if (strength == "high")
  matrix_id = 9
# compute variables for multiple interaction types
plot_df <- data.frame()
for (i in 1:length(interaction_all)) {
  # load eigenvalues and distances for fixed points
  load(file = paste("results/random_matrices/", S, "sp/eigenvalues_evenness_distances_", n, "_random_matrices_", 
                    fixed_points, "_fixed_points_", norm, "_norm_", strength, "_strength_", 
                    interaction_all[i], "_sign", ".RData", sep = ""))
  # load matrices
  load(paste("data/random_matrices/", S, "sp/", n, "_random_matrices_",
             strength, "_strength_", connectance, "_connectance_", 
             interaction_all[i], "_sign.RData", sep = ""))
  # removing equilibrium states with positive largest eigenvalue
  df <- df[df$lambda1 < 0, ]
  # obtaining indexes of remaining matrices
  remaining_mat <- unique(df$matrix)
  matrix_list <- matrix_list[remaining_mat]
  # use one matrix as example
  A <- matrix_list[[matrix_id]]
  # compute most abundant species and its total interaction effect 
  df$matrix <- factor(df$matrix)
  df$most_abundant <- apply(df[ , 1:S], 1, function(x) which.max(x))
  mean_int_strength <- lapply(matrix_list, function(mat) (apply(mat, 2, sum) - diag(mat)) / (nrow(mat) - 1))
  df_mean_int_strength <- as.data.frame(matrix(unlist(mean_int_strength), nrow = length(mean_int_strength), 
                                               byrow = TRUE))
  names(df_mean_int_strength) <- as.character(1:S)
  df_mean_int_strength <- df_mean_int_strength[rep(1:length(matrix_list), 
                                                   tapply(df$N1, df$matrix, length)), ]
  df$int_strength_most_abund <- df_mean_int_strength[cbind(1:nrow(df_mean_int_strength), 
                                                           df$most_abundant)]
  # print correlations for this matrix
  print(interaction_all[i])
  print("Correlation between full recovery and full resistance:")
  print(cor(subset(df, matrix == matrix_id)$lambda1, 
            subset(df, matrix == matrix_id)$min_distance2border_r))
  print("Correlation between partial recovery and partial resistance:")
  print(cor(subset(df, matrix == matrix_id)$lambda2, 
            subset(df, matrix == matrix_id)$min_distance2vertex_r))
  # add current data frame to full data frame
  curr_df <- data.frame(S = rep(S, n),
                        interaction = rep(interaction_all[i], n),
                        strength = rep(strength, n),
                        lambda1 = subset(df, matrix == matrix_id)$lambda1,
                        lambda2 = subset(df, matrix == matrix_id)$lambda2,
                        min_distance2border = subset(df, matrix == matrix_id)$min_distance2border_r,
                        min_distance2vertex = subset(df, matrix == matrix_id)$min_distance2vertex_r,
                        int_strength_most_abund = subset(df, matrix == matrix_id)$int_strength_most_abund,
                        most_abundant = subset(df, matrix == matrix_id)$most_abundant)
  plot_df <- rbind(plot_df, curr_df)
}

# plots of largest eigenvalue and minimum distance to border ------------------------------
# change factor level order
plot_df$interaction <- as.character(plot_df$interaction)
plot_df$interaction[plot_df$interaction == "negative"] <- "Competition system"
plot_df$interaction[plot_df$interaction == "both"] <- "Antagonistic system"
plot_df$interaction[plot_df$interaction == "positive"] <- "Mutualistic system"
plot_df$interaction <- factor(plot_df$interaction, levels = c("Competition system", "Mutualistic system",
                                                              "Antagonistic system"))
# axis limits
y_lim <- c(min(plot_df$lambda1), max(plot_df$lambda1))
# plot for competition
fig_2A_comp <- ggplot(data = subset(plot_df, interaction == "Competition system"), 
                      aes(x = min_distance2border, y = lambda1,
                          shape = as.factor(most_abundant))) +
  geom_point(size = 4, fill = pal[1]) +
  scale_shape_manual(values = c(21, 22, 24), name = "Most abundant species") +
  facet_wrap( ~ interaction) +
  xlab(expression(paste("Full resistance ", (min(group("{", d[b], "}")))))) +
  ylab(expression(paste("Full recovery ", (lambda[1])))) +
  scale_y_continuous(limits = y_lim) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_text(size = 18),
        legend.text = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_A_comp_", strength, ".pdf", sep = ""), 
         fig_2A_comp, width = 16, height = 16, units = "cm")
}
# plot for antagonism
fig_2A_ant <- ggplot(data = subset(plot_df, interaction == "Antagonistic system"), 
                      aes(x = min_distance2border, y = lambda1,
                          shape = as.factor(most_abundant))) +
  geom_point(size = 4, fill = pal[3]) +
  scale_shape_manual(values = c(21, 22, 24), name = "Most abundant species") +
  facet_wrap( ~ interaction) +
  xlab(expression(paste("Full resistance ", (min(group("{", d[b], "}")))))) +
  ylab("") +
  scale_y_continuous(limits = y_lim) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_text(size = 18),
        legend.text = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_A_ant_", strength, ".pdf", sep = ""), 
         fig_2A_ant, width = 16, height = 16, units = "cm")
}
# plot for mutualism
fig_2A_mut <- ggplot(data = subset(plot_df, interaction == "Mutualistic system"), 
                      aes(x = min_distance2border, y = lambda1,
                          shape = as.factor(most_abundant))) +
  geom_point(size = 4, fill = pal[2]) +
  scale_shape_manual(values = c(21, 22, 24), name = "Most abundant species") +
  facet_wrap( ~ interaction) +
  xlab(expression(paste("Full resistance ", (min(group("{", d[b], "}")))))) +
  ylab("") +
  scale_y_continuous(limits = y_lim) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_text(size = 18),
        legend.text = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_A_mut_", strength, ".pdf", sep = ""), 
         fig_2A_mut, width = 16, height = 16, units = "cm")
}

# plots of second smallest eigenvalue and minimum distance to vertex ------------------------------
# axis limits
y_lim <- c(min(plot_df$lambda2), max(plot_df$lambda2))
# plot for competition
fig_3A_comp <- ggplot(data = subset(plot_df, interaction == "Competition system"), 
                      aes(x = min_distance2vertex, y = lambda2,
                          shape = as.factor(most_abundant))) +
  geom_point(size = 4, fill = pal[1]) +
  scale_shape_manual(values = c(21, 22, 24), name = "Most abundant species") +
  facet_wrap( ~ interaction) +
  xlab(expression(paste("Partial resistance ", (min(group("{", d[v], "}")))))) +
  ylab(expression(paste("Partial recovery ", (lambda[2])))) +
  scale_y_continuous(limits = y_lim) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_text(size = 18),
        legend.text = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig3/fig3_A_comp_", strength, ".pdf", sep = ""), 
         fig_3A_comp, width = 16, height = 16, units = "cm")
}
# plot for antagonism
fig_3A_ant <- ggplot(data = subset(plot_df, interaction == "Antagonistic system"), 
                      aes(x = min_distance2vertex, y = lambda2,
                          shape = as.factor(most_abundant))) +
  geom_point(size = 4, fill = pal[3]) +
  scale_shape_manual(values = c(21, 22, 24), name = "Most abundant species") +
  facet_wrap( ~ interaction) +
  xlab(expression(paste("Partial resistance ", (min(group("{", d[v], "}")))))) +
  ylab("") +
  scale_y_continuous(limits = y_lim) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_text(size = 18),
        legend.text = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig3/fig3_A_ant_", strength, ".pdf", sep = ""), 
         fig_3A_ant, width = 16, height = 16, units = "cm")
}
# plot for mutualism
fig_3A_mut <- ggplot(data = subset(plot_df, interaction == "Mutualistic system"), 
                     aes(x = min_distance2vertex, y = lambda2,
                         shape = as.factor(most_abundant))) +
  geom_point(size = 4, fill = pal[2]) +
  scale_shape_manual(values = c(21, 22, 24), name = "Most abundant species") +
  facet_wrap( ~ interaction) +
  xlab(expression(paste("Partial resistance ", (min(group("{", d[v], "}")))))) +
  ylab("") +
  scale_y_continuous(limits = y_lim) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.title = element_text(size = 18),
        legend.text = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig3/fig3_A_mut_", strength, ".pdf", sep = ""), 
         fig_3A_mut, width = 16, height = 16, units = "cm")
}

# compute correlations for all numbers of species and interaction types ------------------------------
# number of species
S_all <- c(3, 4, 5)
# number of fixed points to sample
fixed_points_all <- 100 * 2^(S_all-2)
# type of interaction
interaction_all <- c("negative", "both", "positive")
# compute correlations
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
    # compute correlations and partial correlations
    lambda1_min_dist2border_r <- c()
    lambda_penultimate_min_dist2vertex_r <- c()
    lambda1_lambda_penultimate <- c()
    dist2border_dist2vertex <- c()
    for (l in 1:length(remaining_mat)) {
      lambda1_min_dist2border_r[l] <- cor(subset(df, matrix == remaining_mat[l])$lambda1,
                                          subset(df, matrix == remaining_mat[l])$min_distance2border_r)
      lambda_penultimate_min_dist2vertex_r[l] <- cor(subset(df, matrix == remaining_mat[l])$lambda_penultimate,
                                             subset(df, matrix == remaining_mat[l])$min_distance2vertex_r)
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
# plot
fig_2B <- ggplot(data = plot_df, aes(x = as.factor(S), y = lambda1_min_dist2border_r, color = interaction)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_color_manual(values = pal) +
  facet_wrap(~interaction, ncol = length(interaction_all)) +
  xlab(expression(paste("Number of species ", (S)))) +
  ylab(expression(atop("Correlation between full recovery",
                       "and full resistance" ~ (rho(lambda[1], min(group("{", d[b], "}"))))))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_B_", strength, ".pdf", sep = ""),
         fig_2B, width = 38, height = 14, units = "cm")
}

# plot of second smallest eigenvalue and minimum distance to vertex ------------------------------
fig_3B <- ggplot(data = plot_df, aes(x = as.factor(S), y = lambda_penultimate_min_dist2vertex_r, color = interaction)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_color_manual(values = pal) +
  facet_wrap(~interaction, ncol = length(interaction_all)) +
  xlab(expression(paste("Number of species ", (S)))) +
  ylab(expression(atop("Correlation between partial recovery",
                       "and partial resistance" ~ (rho(lambda[S-1], min(group("{", d[v], "}"))))))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        strip.background = element_rect(size = 1, fill = "white"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig3/fig3_B_", strength, ".pdf", sep = ""), 
         fig_3B, width = 38, height = 14, units = "cm")
}
