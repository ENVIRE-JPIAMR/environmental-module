library(tidyverse)

# Decay rates from studies (k values)
results <- c(
  0.083117706,  # Habteselassie et al.
  0.074276938,  # Islam et al.
  0.067798339,  # Merchant et al.
  0.063084523,  # Pang et al.
  0.009767825,  # Overbeek M0611H
  0.011912521,  # Overbeek M0611L
  0.011537529   # Sharma et al.
)

# Calculate the mean decay rate from the 'results' dataframe
# We can find a better way than the mean maybe?
decay_mean <- mean(results)

# Input parameters from farm module
mean_e_coli_concentration <- 1.6 * 10^4  # CFU/g, mean E. coli concentration in manure
std_dev <- 0.02 * 10^4  # Standard deviation

# Assumed
manure_application_rate <- 2  # kg/m2 - 2 kg per square meter

##  Other
# Time period (days)
days <- 0:200  # From day 0 to day 200

# Total field area (m²)
field_area <- 4000  # Example field area



# Soil properties
soil_depth_m <- 0.1  # Depth of the topsoil layer where manure is mixed (in meters) , https://www.mdpi.com/2073-4395/14/2/278
soil_bulk_density <- 1.47  # g/cm³, https://extension.okstate.edu/fact-sheets/basics-of-soil-bulk-density.html
soil_bulk_density_kg_per_m3 <- soil_bulk_density * 1000  # Convert to kg/m³




# Calculate soil volume and mass per square meter
soil_volume_per_m2 <- soil_depth_m * 1  # Volume of soil per m² (in cubic meters)
soil_mass_per_m2 <- soil_volume_per_m2 * soil_bulk_density_kg_per_m3  # Mass of soil per m² (in kg)

# Calculate total manure applied per square meter (g/m2)
manure_applied_per_m2_grams <- manure_application_rate * 1000  # kg/m² to g/m²

# Simulate random E. coli concentrations per square meter based on mean and std deviation
n_simulations <- 1000  # Number of simulations to run
e_coli_concentrations_per_m2 <- rnorm(n_simulations, mean_e_coli_concentration, std_dev)

# Adjust initial E. coli load based on dilution with soil, per square meter
initial_e_coli_per_m2 <- e_coli_concentrations_per_m2 * manure_applied_per_m2_grams / (soil_mass_per_m2 + manure_applied_per_m2_grams)




# Function to calculate E. coli concentration over time per m² using the mean decay rate
calculate_decay <- function(initial_concentration, decay_rate, days) {
  initial_concentration * exp(-decay_rate * days)
}

# Dataframe to store the decay data over time using the mean decay rate for multiple simulations
decay_data <- tibble(
  simulation = rep(1:n_simulations, each = length(days)),
  day = rep(days, times = n_simulations),
  e_coli_concentration_per_m2 = unlist(lapply(initial_e_coli_per_m2, function(init_conc) calculate_decay(init_conc, decay_mean, days))),
  total_e_coli_on_field = e_coli_concentration_per_m2 * field_area  # Total E. coli on the field
)

# Results decay model
results_summary <- decay_data %>%
  group_by(day) %>%
  summarise(
    mean_e_coli_concentration_per_m2 = mean(e_coli_concentration_per_m2),
    std_e_coli_concentration_per_m2 = sd(e_coli_concentration_per_m2),
    mean_total_e_coli_on_field = mean(total_e_coli_on_field),
    std_total_e_coli_on_field = sd(total_e_coli_on_field)
  )

# Extract mean E. coli on the field for each day from the decay model
mean_total_e_coli_on_field_by_day <- results_summary %>%
  select(day, mean_total_e_coli_on_field)

## Input parameters for the  kind-of-SWAT model

# Random sample of runoff coefficients
runoff_coefficients <- runif(n_simulations, min = 0.05, max = 0.2)

# Bacteria partition coefficient for manure (BACTKDDB = 0.95, from Sowah)
bacteria_partition_coeff <- 0.95

# Apply bacterial parameters (BACT_SWF, WOF_P, also from Sowah paper)
bacteria_active_fraction <- 0.75  # BACT_SWF
wash_off_fraction <- 0.5  # WOF_P

# Apply rainfall data
total_rainfall_mm <- 100  # example rainfall, Berlin July 2022
total_rainfall_m <- total_rainfall_mm / 1000  # Convert to meters

# Bathing site volume (from O'Flaherty 2019)
bathing_site_volumes <- runif(n_simulations, min = 67500, max = 82500) * 1000  # Convert to liters

# Results using the daily mean from the decay model
results_day_by_day <- mean_total_e_coli_on_field_by_day %>%
  mutate(
    e_coli_free_for_transport = mean_total_e_coli_on_field * (1 - bacteria_partition_coeff),  # Bacteria free for transport
    e_coli_washoff = e_coli_free_for_transport * wash_off_fraction,  # Bacteria washed off
    runoff_volume = total_rainfall_m * field_area * mean(runoff_coefficients),  # Mean runoff volume in cubic meters
    e_coli_concentration_in_runoff = e_coli_washoff / runoff_volume,  # E. coli concentration (CFU/m3) in runoff water
    e_coli_concentration_bathing_site = e_coli_concentration_in_runoff / mean(bathing_site_volumes)  # E. coli concentration in bathing site (CFU/L)
  )


# Plot
ggplot(results_day_by_day, aes(x = day, y = e_coli_concentration_bathing_site)) +
  geom_line(color = "blue") +
  labs(title = "E. coli Concentration in Bathing Site Over 200 Days",
       x = "Day",
       y = "E. coli Concentration in Bathing Site (CFU/L)") +
  theme_minimal()



