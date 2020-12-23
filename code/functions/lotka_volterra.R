# Function containing the system of ordinary differential equations
# given by the generalized Lotka-Volterra dynamics 
# Inputs: t = vector of time points;
#         state = vector of state variable; 
#         parameters = vector of parameters (first r vector than vectorized A matrix);
# Output: list containing vector of rates of change
lotka_volterra <- function(t, state, parameters) {
  r <- parameters[1:length(state)]
  A <- matrix(parameters[(length(state) + 1):(length(state) + length(state) * length(state))],
              nrow = length(state), ncol = length(state))
  dx_dt <- state * (r + c(A %*% state))
  list(dx_dt)
}
