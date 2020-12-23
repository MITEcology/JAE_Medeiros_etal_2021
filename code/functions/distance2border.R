# Function that calculates n distances from a vector (r vector) 
# to the convex hull formed by the column vectors of matrix A (faces of 
# the feasibility domain)
# Inputs: r = vector of intrinsic growth rates; A = interaction matrix; 
#         n = number of samples for each convex hull; 
#         norm = which norm to use (l1 or l2);
#         all = whether to output all distances or only minimum distance
# Output: list of distances for each convex hull
distance2border <- function(r, A, n, norm, all) {
  vertices <- combn(1:nrow(A), nrow(A)-1, simplify = FALSE)
  distances_list <- list()
  for (i in 1:length(vertices)) {
    vertex <- vertices[[i]]
    if (norm == "l1") {
      distances_list[[i]] <- 1:n %>% 
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(A)-1))
          border_point <- c(A[, vertex] %*% t)
          euclidean_distance(r, border_point)
        })
    }
    if (norm == "l2") {
      distances_list[[i]] <- 1:n %>% 
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(A)-1))
          border_point <- c(A[, vertex] %*% t)
          border_point_norm <- border_point / sqrt(sum(border_point^2))
          arc_length(r, border_point_norm)
        })
    }
  } 
  names(distances_list) <- unlist(lapply(vertices, paste, collapse = "_"))
  if (all) {
    return(distances_list)
  } else {
    return(sapply(distances_list, min))
  }
}