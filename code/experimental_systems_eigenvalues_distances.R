# Computes full and partial recovery (largest and second smallest eigenvalue) and
# full and partial resistance (minimum distance to border and vertex) for 17 experimental
# 3-species microbial systems (data from Friedman et al 2017 Nat Ecol Evol) as well as
# for 2000 random samples for each experimental system

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(xtable)) {install.packages("xtable"); library(xtable)}
source("code/functions/hypersphere_sampling.R")
source("code/functions/eigenvalues_jacobian.R")
source("code/functions/lotka_volterra.R")
source("code/functions/distance2vertex.R")
source("code/functions/distance2border.R")
source("code/functions/arc_length.R")
source("code/functions/simplex_sampling.R")
set.seed(4)

# settings ------------------------------
# number of species
S <- 3
# number of fixed points to sample
fixed_points <- 2000
# defining model
func <- lotka_volterra
# which norm to use
norm <- "l2"
# whether to save latex tables
save_tex <- FALSE

# loading data ------------------------------
A <- as.matrix(read.table("data/friedman2017_data/friedman2017_matrix.txt"))
r <- as.matrix(read.table("data/friedman2017_data/friedman2017_r.txt"))
K <- as.matrix(read.table("data/friedman2017_data/friedman2017_k.txt"))
# transform matrix to r formalism
A <- A * as.numeric(r / K)

# extract matrices ------------------------------
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
matrix_id <- c(2, 3, 12, 13, 14, 20, 27, 29, 30, 32, 47, 48, 49, 50, 51, 52, 54)
A_list <- A_list[matrix_id]
r_list <- r_list[matrix_id]

# save latex tables ------------------------------
if (save_tex) {
  # A matrix
  A_tex = xtable(-A)
  # save latex table
  print.xtable(A_tex, type = "latex", file = "results/friedman2017_data/friedman2017_matrix.tex")
  # r vector
  colnames(r) <- "r"
  r_tex = xtable(r)
  # save latex table
  print.xtable(r_tex, type = "latex", file = "results/friedman2017_data/friedman2017_r.tex")
  # list of species names
  A_list_names <- lapply(A_list, rownames)
  df_A_list_names <- data.frame(matrix(unlist(A_list_names), nrow = length(A_list_names), 
                                       byrow = TRUE))
  # rename systems according to fig 4
  df_A_list_names <- df_A_list_names[c(6, 12, 2, 1, 11, 13, 7, 16, 8, 10, 9, 3, 4, 5, 15, 17, 14), ]
  rownames(df_A_list_names) <- 1:17 
  names(df_A_list_names) <- c("species 1", "species 2", "species 3")
  A_list_names_tex = xtable(df_A_list_names)
  # save latex table
  print.xtable(A_list_names_tex, type = "latex", file = "results/friedman2017_data/friedman2017_sp_names.tex")
}

# sample points on unit hypersphere ------------------------------
if (norm == "l1") 
  N <- simplex_sampling(fixed_points, S)
if (norm == "l2")
  N <- replicate(fixed_points, hypersphere_sampling(S, positive = TRUE, within = FALSE), simplify = FALSE)

# run analysis for each interaction matrix ------------------------------
df <- data.frame()
for (i in 1:length(A_list)) {
  print(i)
  curr_df <- data.frame()
  # read interaction matrix
  A_curr <- A_list[[i]]
  # read growth rate vector
  r_curr <- r_list[[i]]
  # compute abundance values using equilibrium equation
  N_curr <- -c(solve(A_curr) %*% r_curr)
  # normalize abundance vector to unit norm
  if (norm == "l1") 
    N_curr <- N_curr / sum(N_curr)
  if (norm == "l2")
    N_curr <- N_curr / sqrt(sum(N_curr^2))
  # recompute r vector
  r_curr <- -c(A_curr %*% N_curr)
  # LV parameters
  parms <- as.numeric(c(r_curr, A_curr))
  names(parms) <- c(paste("r", 1:S, sep = ""),
                    paste(paste("a", 1:S, sep = ""), 
                          rep(1:S, each = S), sep = ""))
  # compute empirical eigenvalues (diag(N[[j]]) %*% A gives exactly same result)
  out <- eigenvalues_jacobian(y = as.numeric(N_curr), func = func, parms = parms)
  large_eigen_emp <- sort(out[[1]])[S]
  second_large_eigen_emp <- sort(out[[1]])[S-1]
  # compute eigenvalues for randomly sampled abundance vectors
  large_eigen <- c()
  second_large_eigen <- c()
  r_curr_list <- list()
  for (j in 1:length(N)) {
    # compute r vector using equilibrium equation
    r_curr_list[[j]] <- -c(A_curr %*% N[[j]])
    # LV parameters
    parms <- as.numeric(c(r_curr_list[[j]], A_curr))
    names(parms) <- c(paste("r", 1:S, sep = ""),
                      paste(paste("a", 1:S, sep = ""), 
                            rep(1:S, each = S), sep = ""))
    # compute largest and second smallest Jacobian eigenvalues
    out <- eigenvalues_jacobian(y = as.numeric(N[[j]]), func = func, parms = parms)
    large_eigen[j] <- sort(out[[1]])[S]
    second_large_eigen[j] <- sort(out[[1]])[S-1]
  }
  # normalize r vector to unit norm
  if (norm == "l1") 
    r_curr_dist <- r_curr / sum(r_curr)
  if (norm == "l2")
    r_curr_dist <- r_curr / sqrt(sum(r_curr^2))
  # normalize column vectors of A to unit norm
  if (norm == "l1") 
    A_curr_dist <- apply(A_curr, 2, function(x) x / abs(sum(x)))
  if (norm == "l2") 
    A_curr_dist <- apply(A_curr, 2, function(x) x / sqrt(sum(x^2)))
  # compute minimum distance to vertex and border
  min_distance2vertex_r_emp <- min(distance2vertex(r_curr_dist, -A_curr_dist, norm))
  min_distance2border_r_emp <- min(distance2border(r_curr_dist, -A_curr_dist, n = 100 * 2^(S-2),
                                                   norm, all = FALSE))
  # compute distances for randomly sampled abundance vectors
  min_distance2vertex_r <- c()
  min_distance2border_r <- c()
  for (j in 1:length(N)) {
    # normalize r vectors to unit norm
    if (norm == "l1") 
      r_curr_list[[j]] <- r_curr_list[[j]] / sum(r_curr_list[[j]])
    if (norm == "l2")
      r_curr_list[[j]] <- r_curr_list[[j]] / sqrt(sum(r_curr_list[[j]]^2))
    # compute minimum distance to vertex and border
    min_distance2vertex_r[j] <- min(distance2vertex(r_curr_list[[j]], -A_curr_dist, norm))
    min_distance2border_r[j] <- min(distance2border(r_curr_list[[j]], -A_curr_dist, n = 100 * 2^(S-2),
                                                    norm, all = FALSE))
  }
  # build results data frame
  curr_df <- data.frame(matrix = rep(i, fixed_points),
                        large_eigen_emp = rep(large_eigen_emp, fixed_points),
                        large_eigen = large_eigen,
                        second_large_eigen_emp = rep(second_large_eigen_emp, fixed_points),
                        second_large_eigen = second_large_eigen,
                        min_distance2border_r_emp = rep(min_distance2border_r_emp, fixed_points),
                        min_distance2border_r = min_distance2border_r,
                        min_distance2vertex_r_emp = rep(min_distance2vertex_r_emp, fixed_points),
                        min_distance2vertex_r = min_distance2vertex_r)
  df_N <- data.frame(matrix(unlist(N), nrow = length(N), byrow = TRUE))
  names(df_N) <- paste("N", 1:S, sep = "")
  df_r <- data.frame(matrix(unlist(r_curr_list), nrow = length(r_curr_list), byrow = TRUE))
  names(df_r) <- paste("r", 1:S, sep = "")
  curr_df <- cbind(df_N, df_r, curr_df)
  df <- rbind(df, curr_df)
}
# save results data frame
save(df, file = paste("results/friedman2017_data/eigenvalues_distances_", fixed_points, "_fixed_points_", 
                      norm, "_norm.RData", sep = ""))
