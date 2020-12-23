# Computes full and partial recovery (largest and second smallest eigenvalue) and
# full and partial resistance (minimum distance to border and vertex) for multiple
# interaction matrices with a given number of species, interaction type and strength

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
source("code/functions/lotka_volterra.R")
source("code/functions/distance2border.R")
source("code/functions/distance2vertex.R")
source("code/functions/arc_length.R")
source("code/functions/euclidean_distance.R")
source("code/functions/simplex_sampling.R")
source("code/functions/hypersphere_sampling.R")
source("code/functions/eigenvalues_jacobian.R")

# settings ------------------------------
# number of matrices
n <- 100
# number of species (3, 4, or 5)
S <- 3
# interaction strength ("low" or "high")
strength <- "low"
if (strength == "low")
  stren <- 1/S
if (strength == "high")
  stren <- 1/sqrt(S)
# connectance
connectance <- 1
# interaction sign ("negative" for competition, "both" for antagonism, or "positive" for mutualism)
sign <- "negative"
# number of fixed points to sample
fixed_points <- 100 * 2^(S-2)
# number of points on each border to sample
n_points_border <- 100 * 2^(S-2)
# defining model
func <- lotka_volterra
# which norm to use
norm <- "l2"
# load random matrices 
load(paste("data/random_matrices/", S, "sp/", n, "_random_matrices_",
           strength, "_strength_", connectance, "_connectance_", 
           sign, "_sign.RData", sep = ""))

# perform analyzes for each interaction matrix ------------------------------
df <- data.frame()
for (i in 1:length(matrix_list)) {
  print(i)
  A <- matrix_list[[i]]
  # sampling feasible species abundances
  if (norm == "l1") 
    N <- simplex_sampling(fixed_points, S)
  if (norm == "l2")
    N <- replicate(fixed_points, hypersphere_sampling(S, positive = TRUE, within = FALSE), simplify = FALSE)
  # obtain r vectors using equilibrium equation
  r <- lapply(N, function(mat, x) -c(mat %*% x), mat = A)
  # compute eigenvalues, trace and determinant of Jacobian for each equilibrium
  trace_J <- c()
  det_J <- c()
  eigen_J <- list()
  for (j in 1:fixed_points) {
    # LV parameters
    parms <- as.numeric(c(r[[j]], A))
    names(parms) <- c(paste("r", 1:S, sep = ""),
                      paste(paste("a", 1:S, sep = ""), 
                            rep(1:S, each = S), sep = ""))
    # computing eigenvalues, trace, and determinant Jacobian (diag(N[[j]]) %*% A gives exactly same result)
    out <- eigenvalues_jacobian(y = as.numeric(N[[j]]), func = func, parms = parms)
    eigen_J[[j]] <- sort(out[[1]], decreasing = TRUE)
    trace_J[j] <- out[[2]]
    det_J[j] <- out[[3]]
  }
  # normalize r vectors to unit norm
  r <- lapply(r, function(x) {
    if (norm == "l1") 
      x <- x / sum(x)
    if (norm == "l2")
      x <- x / sqrt(sum(x^2))
    return(x) })
  # normalize column vectors of A to unit norm
  if (norm == "l1") 
    A <- apply(A, 2, function(x) x / abs(sum(x)))
  if (norm == "l2") 
    A <- apply(A, 2, function(x) x / sqrt(sum(x^2)))
  # compute evenness of each N and r
  evenness_N <- sapply(N, function(x) -sum(x * log(x)) / log(S))
  evenness_r <- sapply(r, function(x) -sum(x * log(x)) / log(S))
  # compute asymmetry of feasibility domain
  dist_cols <- as.numeric(dist(t(-A)))
  asymmetry <- sd(dist_cols)
  # compute distances to vertices of the feasibility domain
  distance2vertex_r <- lapply(r, distance2vertex, A = -A,
                              norm = norm)
  min_distance2vertex_r <- sapply(distance2vertex_r, min, na.rm = TRUE)
  max_distance2vertex_r <- sapply(distance2vertex_r, max, na.rm = TRUE)
  avg_distance2vertex_r <- sapply(distance2vertex_r, mean, na.rm = TRUE)
  # compute distances to borders of the feasibility domain
  distance2border_r <- lapply(r, distance2border, A = -A, 
                              n = n_points_border, norm = norm, all = FALSE)
  min_distance2border_r <- sapply(distance2border_r, min, na.rm = TRUE)
  max_distance2border_r <- sapply(distance2border_r, max, na.rm = TRUE)
  avg_distance2border_r <- sapply(distance2border_r, mean, na.rm = TRUE)
  # add results to full data frame
  df_N <- data.frame(matrix(unlist(N), nrow = length(N), byrow = TRUE))
  names(df_N) <- paste("N", 1:S, sep = "")
  df_r <- data.frame(matrix(unlist(r), nrow = length(r), byrow = TRUE))
  names(df_r) <- paste("r", 1:S, sep = "")
  df_eigen <- data.frame(matrix(unlist(eigen_J), nrow = length(eigen_J), byrow = TRUE))
  names(df_eigen) <- paste("lambda", 1:S, sep = "")
  curr_df <- cbind(df_N, df_r, df_eigen)
  curr_df$trace_J <- trace_J
  curr_df$det_J <- det_J
  curr_df$evenness_N <- evenness_N
  curr_df$evenness_r <- evenness_r
  curr_df$asymmetry <- asymmetry
  curr_df$min_distance2vertex_r <- min_distance2vertex_r
  curr_df$max_distance2vertex_r <- max_distance2vertex_r
  curr_df$avg_distance2vertex_r <- avg_distance2vertex_r
  curr_df$min_distance2border_r <- min_distance2border_r
  curr_df$max_distance2border_r <- max_distance2border_r
  curr_df$avg_distance2border_r <- avg_distance2border_r
  curr_df$matrix <- i
  df <- rbind(df, curr_df)
}
# save results data frame
save(df, file = paste("results/random_matrices/", S, "sp/eigenvalues_evenness_distances_", n, "_random_matrices_", 
                      fixed_points, "_fixed_points_", norm, "_norm_", strength, "_strength_", 
                      sign, "_sign", ".RData", sep = ""))
