# Function that computes the eigenvalues, trace, and determinant of the Jacobian 
# matrix of a dynamical system for given state and parameter values
# Inputs: y = vector of state variables; func = function containing the ODE system; 
#         parms = vector of parameters;
# Output: list containing vector of eigenvalues, trace (sum of eigenvalues), and
#         determinant (product of eigenvalues)
eigenvalues_jacobian <- function(y, func, parms) {
  J <- jacobian.full(y = y, func = func, parms = parms)
  eigen_J <- sort(Re(eigen(J)$values), decreasing = TRUE)
  trace_J <- sum(diag(J))
  det_J <- det(J)
  return(list(eigen_J, trace_J, det_J))
}