# Function that calculates the distances from a vector (r vector) 
# to the column vectors of matrix A (vertices of the feasibility domain)
# Inputs: r = vector of intrinsic growth rates; A = interaction matrix;
#         norm = which norm to use (l1 or l2)
# Output: vector of distances
distance2vertex <- function(r, A, norm) {
  if (norm == "l1") {
    distances <- 1:nrow(A) %>% 
      map_dbl(~euclidean_distance(r, A[, .x]))
  }
  if (norm == "l2") {
    distances <- 1:nrow(A) %>% 
      map_dbl(~arc_length(r, A[, .x]))
  }
  names(distances) <- 1:nrow(A)
  return(distances)
}