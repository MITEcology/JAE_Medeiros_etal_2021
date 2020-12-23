# Function that generates a random interaction matrix
# Inputs: num = community size; stren = interaction strength; 
#         conne = connectance; sign = sign of interactions;
#         diag = what to put in diagonal elements
# Output: random interaction matrix
generate_interaction_matrix <- function(num, stren, conne, sign, diag) {
  n <- num * (num-1)
  Inte <- rnorm(n, mean = 0, sd = stren)
  zeroes <- sample(c(rep.int(1, floor(n * conne)), 
                     rep.int(0,(n - floor(n * conne)))))
  Inte[which(zeroes == 0)] <- 0
  mat <- matrix(NA, nrow = num, ncol = num)
  if (sign == "positive") {
    Inte <- abs(Inte)
    mat[upper.tri(mat) | lower.tri(mat)] <- Inte
  }
  if (sign == "negative") {
    Inte <- -abs(Inte)
    mat[upper.tri(mat) | lower.tri(mat)] <- Inte
  }
  if (sign == "both") {
    Inte_pos <- abs(Inte[1:(length(Inte)/2)])
    mat[upper.tri(mat)] <- Inte_pos
    Inte_neg <- -abs(Inte[((length(Inte)/2) + 1):length(Inte)])
    mat[lower.tri(mat)] <- Inte_neg
  }
  if (diag == "one")
    diag(mat) <- -1
  return(mat)
}
