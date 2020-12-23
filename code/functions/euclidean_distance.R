# Function that calculates the euclidean distance between two vectors
# Input: a, b = the two vectors
# Output: the euclidean distance
euclidean_distance <- function(a, b) {
  sqrt(sum((a - b)^2))
}