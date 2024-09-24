source(here::here("environmental-module/load_libraries.R"))
source(here::here("foodborne-module/utilities/estimate_variables.R"))
source(here::here("farm-module/run_farm_module_parallel.R"))
source(here::here("environmental-module/utilities/decay_model.R"))
source(here::here("environmental-module/utilities/visualization.R"))

## Initialization

Runs <- 1000 # number of simulation to be performed

# Create dataframe with specified number of rows and column names
data <- data.frame(matrix(1:Runs, nrow = Runs, ncol = 1))
colnames(data) <- "Runs"

# Load estimated variables metadata 
# TODO: separate constants and variables
input <- read_excel(here("environmental-module/data-input/estimated_variables.xlsx"))

# Simulation of estimated variables
data <- estimate_variables(data, input, N=Runs)

# Run farm module to get initial load and prevalence
parallel_output <- batch_simulator_parallel(n_sim = Runs)

# Load farm module inputs
input_list_farm = load_inputs()

data$C_barn <- parallel_output[36, 1,]/(parallel_output[36, 4,] + input_list_farm$farm_size * input_list_farm$litter_mass)

#computing decay rate in watershed
data$k <- get_decay_rate(data)

# Risk calculation for each module
source(here::here("environmental-module/swat_model.R"))
source(here::here("environmental-module/submodule_drinking.R"))
source(here::here("environmental-module/submodule_swimming.R"))
source(here::here("environmental-module/submodule_lettuce.R"))

# Plot ECDFs
plot_ecdfs(data)



