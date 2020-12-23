# Function that samples m vectors uniformly on a (n-1)-dimensional simplex
# Inputs: m = number of samples; n = dimension of the space
# Output: m vectors of length n that sum up to 1
simplex_sampling <- function(m, n) {
  r <- list()
  for (j in 1:m) {
    dist <- c(sort(runif(n-1, 0, 1)), 1)
    r[[j]] <- c(dist[1], diff(dist))
  }
  return(r)
}