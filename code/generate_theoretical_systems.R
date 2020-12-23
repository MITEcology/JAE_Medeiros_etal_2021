# Generates n random interaction matrices

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}
source("code/functions/generate_interaction_matrix.R")

# generate and save matrices ------------------------------
# number of species (3, 4, or 5)
S <- 3
# interaction strength ("low" or "high")
strength <- "low"
if (strength == "low")
  stren <- 1/S
if (strength == "high")
  stren <- 1/sqrt(S)
# connectance
conne <- 1
# interaction signs ("negative" for competition, "both" for antagonism, or "positive" for mutualism)
sign <- "negative"
# number of matrices to generate
n <- 100
# generate matrices
matrix_list <- list()
for (i in 1:n) {
  print(i)
  # sample matrix
  matrix_list[[i]] <- generate_interaction_matrix(num = S, stren = stren, 
                                                  conne = conne, sign = sign, 
                                                  diag = "one")
  # print message if matrix is not positive definite (i.e. globally stable)
  if(!is.positive.definite(-matrix_list[[i]] + t(-matrix_list[[i]]))) {
    print("not positive definite")
  }
}
# save matrices
save(matrix_list, file = paste("data/random_matrices/", S, "sp/", n, "_random_matrices_",
                               strength, "_strength_", conne, "_connectance_", 
                               sign, "_sign.RData", sep = ""))
