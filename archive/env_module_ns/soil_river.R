# SOIL MODEL #
library(tidyverse)
library(nls2)
library(readxl)
library(triangle)
set.seed(123)  # For reproducibility

# Define number of simulations
n_simulations <- 100

## Estimation of decay rate

# Load the data to estimate the decay rate
# Sharma et al. 2019 - supplementary material

data <- read_excel("sharma_2.xlsx")

# Define the exponential decay function
exp_decay <- function(concentration, decay_rate, days) {
  concentration * exp(-decay_rate * days)
}

# Fit the exponential decay model
# Here, 'y' is the log-transformed E. coli count, and 'dpi' is days post-inoculation.
fit <- nls(y ~ exp_decay(concentration, decay_rate, dpi), 
           data = data, 
           start = list(concentration = max(data$y), decay_rate = 0.001),
           control = nls.control(maxiter = 100))

# Print model summary
summary(fit)

# Extract the estimated decay rate
decay_rate_est <- coef(fit)[["decay_rate"]]

#####

# Input parameters from farm module
e_coli_concentration <- 1.6 * 10^4  # CFU/g, mean E. coli concentration in manure
std_dev <- 0.02 * 10^4  # Standard deviation

# Reduction parameters
stocking_reduction <- 1  # 100%
composting_reduction <- 1  # 100%
anaerobic_reduction_30 <- 0.184272811  # 18%
anaerobic_reduction_37 <- 0.554397820  # 55%

# Calculate initial concentrations and standard deviations for all scenarios
e_coli_concentration_composted <- e_coli_concentration * (1-composting_reduction)
std_dev_composted <- std_dev * (1-composting_reduction)

e_coli_concentration_stocked <- e_coli_concentration * (1-stocking_reduction)
std_dev_stocked <- std_dev * (1-stocking_reduction)

e_coli_concentration_anaerobic_30 <- e_coli_concentration * (1-anaerobic_reduction_30)
std_dev_anaerobic_30 <- std_dev * (1-anaerobic_reduction_30)

e_coli_concentration_anaerobic_37 <- e_coli_concentration * (1-anaerobic_reduction_37)
std_dev_anaerobic_37 <- std_dev * (1-anaerobic_reduction_37)

# Assumed
manure_application_rate <- 2  # kg/m2 - 2 kg per square meter

# Time period (days)
days <- 0:200  # From day 0 to day 100

# Total field area (m²)
field_area <- 4000  # Example field area

# Soil properties
soil_depth_m <- 0.1  # Depth of the topsoil layer where manure is mixed (in meters) - https://www.mdpi.com/2073-4395/14/2/278
soil_bulk_density <- 1.47  # g/cm³ -  https://extension.okstate.edu/fact-sheets/basics-of-soil-bulk-density.html
soil_bulk_density_kg_per_m3 <- soil_bulk_density * 1000  # Convert to kg/m³

# Calculate soil volume and mass per square meter
soil_volume_per_m2 <- soil_depth_m * 1  # Volume of soil per m² (in cubic meters)
soil_mass_per_m2 <- soil_volume_per_m2 * soil_bulk_density_kg_per_m3  # Mass of soil per m² (in kg)
soil_mass_per_m2_grams <- soil_mass_per_m2 * 1000  # Convert soil mass to grams

# Calculate total manure applied per square meter (g/m2)
manure_applied_per_m2_grams <- manure_application_rate * 1000  # kg/m² to g/m²

# Simulate random E. coli concentrations per square meter based on mean and std deviation
n_simulations <- 1000  # Number of simulations to run

# Scenario 1: Fresh manure
e_coli_concentrations_per_gram_fresh <- rnorm(n_simulations, e_coli_concentration, std_dev) # cfu/g manure
initial_e_coli_per_m2_fresh <- e_coli_concentrations_per_gram_fresh * manure_applied_per_m2_grams
# CFU per g/m² (fresh manure)
e_coli_cfu_g_m2_fresh <- e_coli_concentrations_per_gram_fresh * manure_applied_per_m2_grams / (soil_mass_per_m2_grams + manure_applied_per_m2_grams)

# Scenario 2: Composted manure
e_coli_concentrations_per_gram_composted <- rnorm(n_simulations, e_coli_concentration_composted, std_dev_composted) # cfu/g composted manure
initial_e_coli_per_m2_composted <- e_coli_concentrations_per_gram_composted * manure_applied_per_m2_grams
# CFU per g/m² (composted manure)
e_coli_cfu_g_m2_composted <- e_coli_concentrations_per_gram_composted * manure_applied_per_m2_grams / (soil_mass_per_m2_grams + manure_applied_per_m2_grams)

# Scenario 3: Stocked manure
e_coli_concentrations_per_gram_stocked <- rnorm(n_simulations, e_coli_concentration_stocked, std_dev_stocked) # cfu/g stocked manure
initial_e_coli_per_m2_stocked <- e_coli_concentrations_per_gram_stocked * manure_applied_per_m2_grams
# CFU per g/m² (stocked manure)
e_coli_cfu_g_m2_stocked <- e_coli_concentrations_per_gram_stocked * manure_applied_per_m2_grams / (soil_mass_per_m2_grams + manure_applied_per_m2_grams)

# Scenario 4: Anaerobic digestion at 30°C
e_coli_concentrations_per_gram_anaerobic_30 <- rnorm(n_simulations, e_coli_concentration_anaerobic_30, std_dev_anaerobic_30) # cfu/g anaerobic manure 30°C
initial_e_coli_per_m2_anaerobic_30 <- e_coli_concentrations_per_gram_anaerobic_30 * manure_applied_per_m2_grams
# CFU per g/m² (anaerobic manure 30°C)
e_coli_cfu_g_m2_anaerobic_30 <- e_coli_concentrations_per_gram_anaerobic_30 * manure_applied_per_m2_grams / (soil_mass_per_m2_grams + manure_applied_per_m2_grams)

# Scenario 5: Anaerobic digestion at 37°C
e_coli_concentrations_per_gram_anaerobic_37 <- rnorm(n_simulations, e_coli_concentration_anaerobic_37, std_dev_anaerobic_37) # cfu/g anaerobic manure 37°C
initial_e_coli_per_m2_anaerobic_37 <- e_coli_concentrations_per_gram_anaerobic_37 * manure_applied_per_m2_grams
# CFU per g/m² (anaerobic manure 37°C)
e_coli_cfu_g_m2_anaerobic_37 <- e_coli_concentrations_per_gram_anaerobic_37 * manure_applied_per_m2_grams / (soil_mass_per_m2_grams + manure_applied_per_m2_grams)

# Dataframe to store the decay data over time using the mean decay rate for multiple simulations

create_decay_data <- function(initial_e_coli_per_m2, e_coli_cfu_g_m2) {
  tibble(
    simulation = rep(1:n_simulations, each = length(days)),
    day = rep(days, times = n_simulations),
    e_coli_concentration_per_m2 = unlist(lapply(initial_e_coli_per_m2, function(init_conc) exp_decay(init_conc, decay_rate_est, days))),
    total_e_coli_on_field = e_coli_concentration_per_m2 * field_area,  # Total E. coli on the field
    # E. coli concentration per g·m² (adjusted with the soil dilution)
    e_coli_cfu_per_g_m2 = unlist(lapply(e_coli_cfu_g_m2, function(init_conc) exp_decay(init_conc, decay_rate_est, days)))
  )
}

# Create decay data for each scenario
decay_data_fresh <- create_decay_data(initial_e_coli_per_m2_fresh, e_coli_cfu_g_m2_fresh)
decay_data_composted <- create_decay_data(initial_e_coli_per_m2_composted, e_coli_cfu_g_m2_composted)
decay_data_stocked <- create_decay_data(initial_e_coli_per_m2_stocked, e_coli_cfu_g_m2_stocked)
decay_data_anaerobic_30 <- create_decay_data(initial_e_coli_per_m2_anaerobic_30, e_coli_cfu_g_m2_anaerobic_30)
decay_data_anaerobic_37 <- create_decay_data(initial_e_coli_per_m2_anaerobic_37, e_coli_cfu_g_m2_anaerobic_37)

# Function to summarize decay results
summarize_decay_results <- function(decay_data) {
  decay_data %>%
    group_by(day) %>%
    summarise(
      # CFU/m²
      mean_e_coli_concentration_per_m2 = mean(e_coli_concentration_per_m2),
      std_e_coli_concentration_per_m2 = sd(e_coli_concentration_per_m2),
      # Total CFU on field
      mean_total_e_coli_on_field = mean(total_e_coli_on_field),
      std_total_e_coli_on_field = sd(total_e_coli_on_field),
      # CFU/g·m²
      mean_e_coli_cfu_per_g_m2 = mean(e_coli_cfu_per_g_m2),
      std_e_coli_cfu_per_g_m2 = sd(e_coli_cfu_per_g_m2)
    )
}

# Summarize results for each scenario
soil_results_summary_fresh <- summarize_decay_results(decay_data_fresh)
soil_results_summary_composted <- summarize_decay_results(decay_data_composted)
soil_results_summary_stocked <- summarize_decay_results(decay_data_stocked)
soil_results_summary_anaerobic_30 <- summarize_decay_results(decay_data_anaerobic_30)
soil_results_summary_anaerobic_37 <- summarize_decay_results(decay_data_anaerobic_37)

# Plot results for all scenarios, m2
ggplot() +
  geom_line(data = soil_results_summary_fresh, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Fresh Manure")) +
  geom_line(data = soil_results_summary_composted, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Composted Manure")) +
  geom_line(data = soil_results_summary_stocked, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Stocked Manure")) +
  geom_line(data = soil_results_summary_anaerobic_30, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Anaerobic Digestion 30°C")) +
  geom_line(data = soil_results_summary_anaerobic_37, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Anaerobic Digestion 37°C")) +
  labs(title = "E. coli Decay per Square Meter in Soil",
       x = "Day",
       y = "CFU/m²",
       color = "Intervention") +
  theme_minimal()

# log scale
ggplot() +
  geom_line(data = soil_results_summary_fresh, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Fresh Manure"), size = 1) +
  geom_line(data = soil_results_summary_composted, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Composted Manure"), size = 1) +
  geom_line(data = soil_results_summary_stocked, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Stocked Manure"), size = 1) +
  geom_line(data = soil_results_summary_anaerobic_30, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Anaerobic Digestion 30°C"), size = 1) +
  geom_line(data = soil_results_summary_anaerobic_37, aes(x = day, y = mean_e_coli_concentration_per_m2, color = "Anaerobic Digestion 37°C"), size = 1) +
  scale_y_log10() +  # Apply log scale to better highlight smaller differences
  labs(title = "E. coli Decay per Square Meter in Soil Over 200 Days",
       x = "Days",
       y = "E. coli Concentration (CFU/m², log scale)",
       color = "Scenario") +
  theme_minimal()


# RIVER MODEL #

# Fixed parameters, from Sowah 2020 (SWAT model)
bacteria_partition_coeff <- 0.95  # Fraction of bacteria bound to manure
bacteria_active_fraction <- 0.75  # Fraction of bacteria that are active
wash_off_fraction <- 0.5  # Fraction of bacteria washed off during runoff
field_area <- 8093.71  # Field area (m²)
total_ecoli_bathing_site_after_decay <- 0

# Function to run simulations for different soil scenarios
run_river_simulation <- function(soil_results_summary) {
  map_dfr(1:n_simulations, function(sim) {
    # Environmental parameters for decay rate (Mancini's Equation)
    T <- runif(1, 21, 28)  # River water temperature (°C)
    seawater <- runif(1, 0.035, 0.075)  # Salinity (%)
    IA <- rtriangle(1, 17.32, 25.38, 22.73)  # Solar radiation levels (ly/hr)
    et <- runif(1, 0.26, 0.31)  # Light extinction coefficient (m⁻¹)
    H <- runif(1, 0.5, 6)  # Water depth (m)
    
    bathing_site_volumes <- runif(1, min = 67500, max = 82500) * 1000  # Convert to liters
    
    # Parameters for the SWAT-like model
    runoff_coefficients <- runif(1, 0.05, 0.2)  # Random runoff coefficient for this simulation
    
    # Mancini's equation for decay rate
    k <- 0.8 + 0.006 * seawater * 1.07^(T - 20) + (IA * (1 - exp(-et * H))) / (et * H)
    
    # River flow and velocity parameters
    #Dist_x <- runif(1, 800, 1000)  # Distance from source to river mouth (m)
    #F_River_x <- rtriangle(1, 67500, 82500, 75000)  # Flow rate of river (m³/day)
    #Width_x <- rtriangle(1, 4.33, 21.71, 14.00)  # Width of river (m)
    #Vel_x <- F_River_x / (Width_x * H)  # River velocity (m/s)
    #t_x <- Dist_x / Vel_x  # Exposure time in river (s)
    #Fr_Day_x <- rtriangle(1, 67500, 82500)
    
    # Merging with rainfall events & soil model 
    # Calculating concentration of E. coli in bathing site
    
    all_days <- tibble(day = 1:200) %>%
      # Merge with results_summary to bring in mean_total_e_coli_on_field
      left_join(soil_results_summary %>% select(day, mean_total_e_coli_on_field), by = "day") %>%
      mutate(
        e_coli_free_for_transport = ifelse(is.na(mean_total_e_coli_on_field), 0, mean_total_e_coli_on_field * (1 - bacteria_partition_coeff)),  # Bacteria free for transport
        e_coli_washoff = e_coli_free_for_transport * wash_off_fraction,  # Bacteria washed off
        # Cumulative E. coli in the bathing site (accounting for decay)
        total_ecoli_bathing_site = lag(total_ecoli_bathing_site_after_decay, default = 0) + e_coli_washoff,  # Total E. coli in the bathing site
        
        # Apply Mancini decay to the total E. coli in the bathing site
        total_ecoli_bathing_site_after_decay = total_ecoli_bathing_site * exp(-k), #removed t_x
        
        # E. coli concentration in the bathing site
        e_coli_concentration_bathing_site = total_ecoli_bathing_site_after_decay  / mean(bathing_site_volumes)  # E. coli concentration in bathing site (CFU/L)
      )
    
    # Return relevant results for this simulation
    all_days %>%
      select(day, e_coli_concentration_bathing_site) %>%
      mutate(simulation = sim)  # Tag each result with the simulation number
  })
}

# Run river simulations for each soil scenario
simulation_results_fresh <- run_river_simulation(soil_results_summary_fresh)
simulation_results_composted <- run_river_simulation(soil_results_summary_composted)
simulation_results_stocked <- run_river_simulation(soil_results_summary_stocked)
simulation_results_anaerobic_30 <- run_river_simulation(soil_results_summary_anaerobic_30)
simulation_results_anaerobic_37 <- run_river_simulation(soil_results_summary_anaerobic_37)

# Summarize river simulation results for each scenario
summarize_river_results <- function(simulation_results) {
  simulation_results %>%
    group_by(day) %>%
    summarize(
      mean_concentration = mean(e_coli_concentration_bathing_site, na.rm = TRUE) / 1000,
      sd_concentration = sd(e_coli_concentration_bathing_site, na.rm = TRUE) / 1000
    )
}

river_summary_results_fresh <- summarize_river_results(simulation_results_fresh)
river_summary_results_composted <- summarize_river_results(simulation_results_composted)
river_summary_results_stocked <- summarize_river_results(simulation_results_stocked)
river_summary_results_anaerobic_30 <- summarize_river_results(simulation_results_anaerobic_30)
river_summary_results_anaerobic_37 <- summarize_river_results(simulation_results_anaerobic_37)

# Plot results for river concentration comparison
ggplot() +
  geom_line(data = river_summary_results_fresh, aes(x = day, y = mean_concentration , color = "Fresh Manure")) +
  geom_ribbon(data = river_summary_results_fresh, aes(x = day, ymin = (mean_concentration - sd_concentration) , ymax = (mean_concentration + sd_concentration)), alpha = 0.2, fill = "red") +
  geom_line(data = river_summary_results_composted, aes(x = day, y = mean_concentration , color = "Composted Manure")) +
  geom_ribbon(data = river_summary_results_composted, aes(x = day, ymin = (mean_concentration - sd_concentration) , ymax = (mean_concentration + sd_concentration) ), alpha = 0.2, fill = "green") +
  geom_line(data = river_summary_results_stocked, aes(x = day, y = mean_concentration , color = "Stocked Manure")) +
  geom_line(data = river_summary_results_anaerobic_30, aes(x = day, y = mean_concentration , color = "Anaerobic Digestion 30°C")) +
  geom_line(data = river_summary_results_anaerobic_37, aes(x = day, y = mean_concentration , color = "Anaerobic Digestion 37°C")) +
  labs(title = "E. coli Concentration in Watershed",
       x = "Day",
       y = "CFU/mL",
       color = "Intervention") +
  theme_minimal()

# log scale
ggplot() +
  geom_line(data = river_summary_results_fresh, aes(x = day, y = mean_concentration , color = "Fresh Manure")) +
  geom_line(data = river_summary_results_composted, aes(x = day, y = mean_concentration , color = "Composted Manure")) +
  geom_line(data = river_summary_results_stocked, aes(x = day, y = mean_concentration , color = "Stocked Manure")) +
  geom_line(data = river_summary_results_anaerobic_30, aes(x = day, y = mean_concentration , color = "Anaerobic Digestion 30°C")) +
  geom_line(data = river_summary_results_anaerobic_37, aes(x = day, y = mean_concentration , color = "Anaerobic Digestion 37°C")) +
  labs(title = "E. coli Concentration in Watershed",
       x = "Day",
       y = "CFU/mL",
       color = "Scenario") +
  theme_minimal() +
  scale_y_log10() +  # Apply log scale to better highlight smaller differences
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
