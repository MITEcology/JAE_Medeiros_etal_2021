# Code for figure 1: distribution of eigenvalues and distances to border/vertex
# across the feasibility domain of a 3-species experimental system (Ea-Pp-Sm) 
# from Friedman et al (2017) Nat Ecol Evol

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(corrplot)) {install.packages("corrplot"); library(corrplot)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
source("code/functions/simplex_sampling.R")

# settings ------------------------------
# whether to save plots
save_plots <- FALSE
# number of species
S <- 3
# number of fixed points
fixed_points <- 20000 * 2^(S-2)
# which norm to use
norm <- "l2"
# matrix to use (this is Ea-Pp-Sm system)
matrix_id <- 6
# load R data
load(file = paste("results/friedman2017_data/eigenvalues_distances_", fixed_points, "_fixed_points_", 
                  norm, "_norm.RData", sep = ""))
# check if all eigenvalues are negative
all(df$large_eigen < 0)
all(df$large_eigen_emp < 0)

# loading data ------------------------------
A <- as.matrix(read.table("data/friedman2017_data/friedman2017_matrix.txt"))
r <- as.matrix(read.table("data/friedman2017_data/friedman2017_r.txt"))
K <- as.matrix(read.table("data/friedman2017_data/friedman2017_k.txt"))
# transform matrix to r formalism
A <- A * as.numeric(r / K)

# extract matrices ------------------------------
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

# Fig A: fixed points (N space) colored by largest eigenvalue ------------------------------
# create data frame with borders of feasibility domain
A <- matrix(c(1, 0, 0,
              0, 1, 0, 
              0, 0, 1), nrow = S, ncol = S)
names(A) <- paste("N", 1:S, sep = "")
t <- seq(0, 1, 0.05)
v1 <- list()
v2 <- list()
v3 <- list()
for (i in 1:length(t)) {
  x <- (A[1, ] * t[i]) + (A[2, ] * (1 - t[i]))
  v1[[i]] <- x / sqrt(sum(x^2))
  x <- (A[2, ] * t[i]) + (A[3, ] * (1 - t[i]))
  v2[[i]] <- x / sqrt(sum(x^2))
  x <- (A[3, ] * t[i]) + (A[1, ] * (1 - t[i]))
  v3[[i]] <- x / sqrt(sum(x^2))
}
A1 <- data.frame(matrix(unlist(v1), nrow = length(v1), byrow = TRUE))
names(A1) <- paste("N", 1:S, sep = "") 
A1$large_eigen <- NA
A2 <- data.frame(matrix(unlist(v2), nrow = length(v2), byrow = TRUE))
names(A2) <- paste("N", 1:S, sep = "") 
A2$large_eigen <- NA
A3 <- data.frame(matrix(unlist(v3), nrow = length(v3), byrow = TRUE))
names(A3) <- paste("N", 1:S, sep = "")
A3$large_eigen <- NA
# make plot
fig_A <- plot_ly(x = ~N2, y = ~N1, z = ~N3, color = ~large_eigen,
                 colors = viridis_pal(option = "B", direction = -1)(fixed_points),
                 width = 700, height = 450) %>% 
  add_markers(data = subset(df, matrix == matrix_id), 
              marker = list(size = 2.7)) %>% 
  add_lines(data = A1, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 17)) %>% 
  add_lines(data = A2, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 17)) %>% 
  add_lines(data = A3, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 17)) %>% 
  layout(scene = list(xaxis = list(title = "N<sub>2</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"), 
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0, 1),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      yaxis = list(title = "N<sub>1</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0, 1),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      zaxis = list(title = "N<sub>3</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0, 1),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
         margin = list(l = 0,
                       r = 0,
                       b = 0,
                       t = 0),
         showlegend = FALSE) %>%
  colorbar(len = 0.38,
           y = 0.8,
           x = 0.81,
           xanchor = "center",
           yanchor = "top",
           thickness = 18,
           outlinecolor = "black",
           bordercolor = "black",
           tickfont = list(size = 14, 
                           family = "Arial, sans-serif",
                           color = "black"),
           ticklen = 4.5,
           tickcolor = "black",
           title = list(text = "Largest\neigenvalue",
                        font = list(size = 18, 
                                    family = "Arial, sans-serif",
                                    color = "black"))) 
# save plot
if (save_plots) {
  orca(fig_A, "figs/figSI/fig1/fig1_A_EaPpSm.pdf", format = "pdf",
       width = 700, height = 450)
}

# Fig B: fixed points (N space) colored by second smallest eigenvalue ------------------------------
A1$second_large_eigen <- NA
A2$second_large_eigen <- NA
A3$second_large_eigen <- NA
# make plot
fig_B <- plot_ly(x = ~N2, y = ~N1, z = ~N3, color = ~second_large_eigen,
                 colors = viridis_pal(option = "B", direction = -1)(fixed_points),
                 width = 700, height = 450) %>% 
  add_markers(data = subset(df, matrix == matrix_id), 
              marker = list(size = 2.7)) %>% 
  add_lines(data = A1, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 17)) %>% 
  add_lines(data = A2, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 17)) %>% 
  add_lines(data = A3, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 17)) %>% 
  layout(scene = list(xaxis = list(title = "N<sub>2</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"), 
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0, 1),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      yaxis = list(title = "N<sub>1</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0, 1),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      zaxis = list(title = "N<sub>3</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0, 1),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
         margin = list(l = 0,
                       r = 0,
                       b = 0,
                       t = 0),
         showlegend = FALSE) %>%
  colorbar(len = 0.38,
           y = 0.8,
           x = 0.83,
           xanchor = "center",
           yanchor = "top",
           thickness = 18,
           outlinecolor = "black",
           bordercolor = "black",
           tickfont = list(size = 14, 
                           family = "Arial, sans-serif",
                           color = "black"),
           ticklen = 4.5,
           tickcolor = "black",
           title = list(text = "Second smallest\neigenvalue",
                        font = list(size = 18, 
                                    family = "Arial, sans-serif",
                                    color = "black"))) 
# save plot
if (save_plots) {
  orca(fig_B, "figs/figSI/fig1/fig1_B_EaPpSm.pdf", format = "pdf",
       width = 700, height = 450)
}

# Fig C: fixed points (r space) colored by minimum distance to border ------------------------------
# create data frame with borders of feasibility domain
A <- A_list[[matrix_id]]
if (norm == "l1") 
  A <- as.data.frame(apply(A, 2, function(x) x / abs(sum(x))))
if (norm == "l2") 
  A <- as.data.frame(apply(A, 2, function(x) x / sqrt(sum(x^2))))
A <- as.data.frame(-t(A))
names(A) <- paste("r", 1:S, sep = "")
t <- seq(0, 1, 0.5)
v1 <- list()
v2 <- list()
v3 <- list()
for (i in 1:length(t)) {
  x <- (A[1, ] * t[i]) + (A[2, ] * (1 - t[i]))
  v1[[i]] <- x / sqrt(sum(x^2))
  x <- (A[2, ] * t[i]) + (A[3, ] * (1 - t[i]))
  v2[[i]] <- x / sqrt(sum(x^2))
  x <- (A[3, ] * t[i]) + (A[1, ] * (1 - t[i]))
  v3[[i]] <- x / sqrt(sum(x^2))
}
A1 <- data.frame(matrix(unlist(v1), nrow = length(v1), byrow = TRUE))
names(A1) <- paste("r", 1:S, sep = "") 
A1$min_distance2border_r <- NA
A2 <- data.frame(matrix(unlist(v2), nrow = length(v2), byrow = TRUE))
names(A2) <- paste("r", 1:S, sep = "") 
A2$min_distance2border_r <- NA
A3 <- data.frame(matrix(unlist(v3), nrow = length(v3), byrow = TRUE))
names(A3) <- paste("r", 1:S, sep = "")
A3$min_distance2border_r <- NA
# make plot
fig_C <- plot_ly(x = ~r2, y = ~r1, z = ~r3, color = ~min_distance2border_r,
                 colors = viridis_pal(option = "B", direction = -1)(fixed_points),
                 width = 700, height = 450) %>% 
  add_markers(data = subset(df, matrix == matrix_id), 
              marker = list(size = 2.7)) %>% 
  add_lines(data = A1, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 14)) %>% 
  add_lines(data = A2, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 14)) %>% 
  add_lines(data = A3, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 14)) %>% 
  layout(scene = list(xaxis = list(title = "r<sub>2</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"), 
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0.65, 0.8),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      yaxis = list(title = "r<sub>1</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0.45, 0.65),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      zaxis = list(title = "r<sub>3</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0.3, 0.45),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
         margin = list(l = 0,
                       r = 0,
                       b = 0,
                       t = 0),
         showlegend = FALSE) %>% 
  colorbar(len = 0.38,
           y = 0.8,
           x = 0.83,
           xanchor = "center",
           yanchor = "top",
           thickness = 18,
           outlinecolor = "black",
           bordercolor = "black",
           tickfont = list(size = 14, 
                           family = "Arial, sans-serif",
                           color = "black"),
           ticklen = 4.5,
           tickcolor = "black",
           title = list(text = "Distance to\nclosest border",
                        font = list(size = 18, 
                                    family = "Arial, sans-serif",
                                    color = "black")))
# save plot
if (save_plots) {
  orca(fig_C, "figs/figSI/fig1/fig1_C_EaPpSm.pdf", format = "pdf",
       width = 700, height = 450)
}

# Fig D: fixed points (r space) colored by minimum distance to vertex ------------------------------
A1$min_distance2vertex_r <- NA
A2$min_distance2vertex_r <- NA
A3$min_distance2vertex_r <- NA
# make plot
fig_D <- plot_ly(x = ~r2, y = ~r1, z = ~r3, color = ~min_distance2vertex_r,
                 colors = viridis_pal(option = "B", direction = -1)(fixed_points),
                 width = 700, height = 450) %>% 
  add_markers(data = subset(df, matrix == matrix_id), 
              marker = list(size = 2.7)) %>% 
  add_lines(data = A1, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 14)) %>% 
  add_lines(data = A2, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 14)) %>% 
  add_lines(data = A3, type = 'scatter3d', mode = 'lines',
            line = list(color = "black", width = 14)) %>% 
  layout(scene = list(xaxis = list(title = "r<sub>2</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"), 
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0.65, 0.8),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      yaxis = list(title = "r<sub>1</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0.45, 0.65),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      zaxis = list(title = "r<sub>3</sub>",
                                   titlefont = list(size = 24, 
                                                    family = "Arial, sans-serif",
                                                    color = "black"),
                                   tickfont = list(size = 14,
                                                   family = "Arial, sans-serif",
                                                   color = "black"),
                                   range = c(0.3, 0.45),
                                   ticklen = 4.5,
                                   tickcolor = "#7F7F7F",
                                   gridcolor = "#7F7F7F",
                                   gridwidth = 0.6,
                                   zerolinewidth = 2.4),
                      aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
         margin = list(l = 0,
                       r = 0,
                       b = 0,
                       t = 0),
         showlegend = FALSE) %>% 
  colorbar(len = 0.38,
           y = 0.8,
           x = 0.83,
           xanchor = "center",
           yanchor = "top",
           thickness = 18,
           outlinecolor = "black",
           bordercolor = "black",
           tickfont = list(size = 14, 
                           family = "Arial, sans-serif",
                           color = "black"),
           ticklen = 4.5,
           tickcolor = "black",
           title = list(text = "Distance to\nclosest vertex",
                        font = list(size = 18, 
                                    family = "Arial, sans-serif",
                                    color = "black")))
# save plot
if (save_plots) {
  orca(fig_D, "figs/figSI/fig1/fig1_D_EaPpSm.pdf", format = "pdf",
       width = 700, height = 450)
}
