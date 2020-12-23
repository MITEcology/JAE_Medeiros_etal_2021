# Function that samples a d-dimensional vector uniformly 
# within or on the d-dimensional hypersphere
# Inputs: d = dimension of the space to sample from;
#         positive = whether only vectors in the positive orthant are needed;
#         within = whether the points should also be sampled inside the hypersphere
# Output: the sampled vector
hypersphere_sampling <- function(d, positive, within) {
  x <- rnorm(d, 0, 1)
  if (positive)
    x <- abs(x)
  norm_x <- sqrt(sum(x^2))
  x_scaled <- x / norm_x
  if (within) {
    u <- runif(1, 0, 1)^(1/d)
    x_scaled <- x_scaled * u
  }
  return(x_scaled)
}