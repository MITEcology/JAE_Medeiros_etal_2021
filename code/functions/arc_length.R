# Function that calculates the arc length between two vectors
# Input: a, b = two vectors with unit l2 norm
# Output: the arc length in radians
arc_length <- function(a, b) {
  acos(sum(a * b))
}