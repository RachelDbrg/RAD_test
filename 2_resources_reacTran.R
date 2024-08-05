# Install and load the ReacTran package
library(ReacTran)

# Define the model parameters
params <- list(
  resource1_growth_rate = 0.06,  # Growth rate of resource 1
  resource2_growth_rate = 0.05,  # Growth rate of resource 2
  resource1_diffusion_rate = 0.05, # Diffusion rate of resource 1
  resource2_diffusion_rate = 0.04, # Diffusion rate of resource 2
  prey_growth_rate = 10,        # Prey growth rate
  predation_rate = 0.01,        # Rate at which predators eat prey
  browsing_rate = 0.5,          # Rate at which prey consume resources
  predator_death_rate = 0.001,  # Predator death rate
  predator_efficiency = 0.01,   # Efficiency of converting prey to predators
  diffusion_rate_prey = 10,     # Diffusion rate of prey
  diffusion_rate_predator = 0.1, # Diffusion rate of predator
  advection_rate_prey = 0,      # Advection rate of prey
  advection_rate_predator = 0.1 # Advection rate of predator
)

# Set up the numerical grid
grid_size <- 50 # Size of the grid
time_steps <- 5000 # Number of time steps
dx <- 1.0 # Spatial step size
dy <- 1.0 # Spatial step size
dt <- 0.01 # Time step size

# Initialize population densities
prey <- matrix(0, nrow = grid_size, ncol = grid_size)
predator <- matrix(0, nrow = grid_size, ncol = grid_size)
resource1 <- matrix(0, nrow = grid_size, ncol = grid_size)
resource2 <- matrix(0, nrow = grid_size, ncol = grid_size)

# Define initial positions randomly
set.seed(123) # For reproducibility
num_prey <- 0
num_predators <- 0
num_resources <- 1000
prey_positions <- sample(1:(grid_size * grid_size), num_prey)
predator_positions <- sample(1:(grid_size * grid_size), num_predators)
resource_positions <- sample(1:(grid_size * grid_size), num_resources)

prey[prey_positions] <- 1.0
predator[predator_positions] <- 1.0
resource1[resource_positions] <- 1.0
resource2[resource_positions] <- 0.5 # Different initial density for resource 2

# Function to define the reaction-advection-diffusion model
rad_model <- function(time, state, params) {
  # Extract population densities
  prey <- state[1:grid_size^2]
  predator <- state[(grid_size^2 + 1):(2 * grid_size^2)]
  resource1 <- state[(2 * grid_size^2 + 1):(3 * grid_size^2)]
  resource2 <- state[(3 * grid_size^2 + 1):(4 * grid_size^2)]
  
  # Reshape to 2D matrices
  prey <- matrix(prey, nrow = grid_size, ncol = grid_size)
  predator <- matrix(predator, nrow = grid_size, ncol = grid_size)
  resource1 <- matrix(resource1, nrow = grid_size, ncol = grid_size)
  resource2 <- matrix(resource2, nrow = grid_size, ncol = grid_size)
  
  # Reaction terms
  resource1_growth <- params$resource1_growth_rate * (1 - resource1 - 0.5 * resource2)  # Assuming logistic growth with competition
  resource2_growth <- params$resource2_growth_rate * (1 - resource2 - 0.5 * resource1)  # Assuming logistic growth with competition
  prey_growth <- params$prey_growth_rate * prey * (resource1 + resource2)
  resource1_consumption <- params$browsing_rate * prey * resource1
  resource2_consumption <- params$browsing_rate * prey * resource2
  predation <- params$predation_rate * prey * predator
  predator_growth <- params$predator_efficiency * predation
  predator_death <- params$predator_death_rate * predator
  
  # Diffusion terms
  prey_diffusion <- tran.2D(C = prey, D.x = params$diffusion_rate_prey, D.y = params$diffusion_rate_prey, dx = dx, dy = dy)$dC
  predator_diffusion <- tran.2D(C = predator, D.x = params$diffusion_rate_predator, D.y = params$diffusion_rate_predator, dx = dx, dy = dy)$dC
  resource1_diffusion <- tran.2D(C = resource1, D.x = params$resource1_diffusion_rate, D.y = params$resource1_diffusion_rate, dx = dx, dy = dy)$dC
  resource2_diffusion <- tran.2D(C = resource2, D.x = params$resource2_diffusion_rate, D.y = params$resource2_diffusion_rate, dx = dx, dy = dy)$dC
  
  # Advection terms
  prey_advection <- tran.2D(C = prey, v.x = params$advection_rate_prey, v.y = params$advection_rate_prey, D.x = 0.1, D.y = 0.1, dx = dx, dy = dy)$dC
  predator_advection <- tran.2D(C = predator, v.x = params$advection_rate_predator, v.y = params$advection_rate_predator, D.x = 0.1, D.y = 0.1, dx = dx, dy = dy)$dC
  
  # Update population densities
  dprey <- prey_growth - predation + prey_diffusion + prey_advection
  dpredator <- predator_growth - predator_death + predator_diffusion + predator_advection
  dresource1 <- resource1_growth - resource1_consumption + resource1_diffusion
  dresource2 <- resource2_growth - resource2_consumption + resource2_diffusion
  
  # Return the derivatives as a vector
  list(c(as.vector(dprey), as.vector(dpredator), as.vector(dresource1), as.vector(dresource2)))
}

# Initial state vector
state <- c(as.vector(prey), as.vector(predator), as.vector(resource1), as.vector(resource2))

# Time points for the output
times <- seq(0, dt * time_steps, by = dt)

# Set larger workspace sizes
lrw <- 10000000 # Real workspace
liw <- 100000  # Integer workspace

# Run the simulation
out <- ode.2D(y = state, times = times, func = rad_model, parms = params, 
              nspec = 4, dimens = c(grid_size, grid_size),
              lrw = lrw, liw = liw)

# Extract prey, predator, and resource populations at different times
prey_initial <- matrix(out[1, 2:(grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
prey_intermediate1 <- matrix(out[200, 2:(grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
prey_intermediate2 <- matrix(out[300, 2:(grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
prey_intermediate3 <- matrix(out[400, 2:(grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
prey_final <- matrix(out[length(times), 2:(grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)

predator_initial <- matrix(out[1, (grid_size^2 + 2):(2 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
predator_final <- matrix(out[length(times), (grid_size^2 + 2):(2 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)




resource1_initial <- matrix(out[1, (2 * grid_size^2 + 2):(3 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource1_int1 <- matrix(out[100, (2 * grid_size^2 + 2):(3 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource1_int2 <- matrix(out[200, (2 * grid_size^2 + 2):(3 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource1_int3 <- matrix(out[300, (2 * grid_size^2 + 2):(3 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource1_final <- matrix(out[length(times), (2 * grid_size^2 + 2):(3 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)


# Plot resource populations at initial and final times
image(resource1_initial, main = "Resource Population at Time = 0", col = heat.colors(256))
image(resource1_int1, main = "Resource Population at Time = 100", col = heat.colors(256))
image(resource1_int2, main = "Resource Population at Time = 200", col = heat.colors(256))
image(resource1_int3, main = "Resource Population at Time = 300", col = heat.colors(256))
image(resource1_final, main = "Resource Population at Time = 500", col = heat.colors(256))




resource2_initial <- matrix(out[1, (3 * grid_size^2 + 2):(4 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource2_int1 <- matrix(out[100, (3 * grid_size^2 + 2):(4 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource2_int2 <- matrix(out[200, (3 * grid_size^2 + 2):(4 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource2_int3 <- matrix(out[300, (3 * grid_size^2 + 2):(4 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)
resource2_final <- matrix(out[length(times), (3 * grid_size^2 + 2):(4 * grid_size^2 + 1)], nrow = grid_size, ncol = grid_size)


# Plot resource populations at initial and final times
image(resource2_initial, main = "Resource Population at Time = 0", col = heat.colors(256))
image(resource2_int1, main = "Resource Population at Time = 100", col = heat.colors(256))
image(resource2_int2, main = "Resource Population at Time = 200", col = heat.colors(256))
image(resource2_int3, main = "Resource Population at Time = 300", col = heat.colors(256))
image(resource2_final, main = "Resource Population at Time = 500", col = heat.colors(256))


par(mfrow = c(2, 2))
