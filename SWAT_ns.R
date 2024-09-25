library(tidyverse)

# Input parameters
mean_e_coli_concentration <- 1.6 * 10^4  # CFU/g, mean E. coli concentration from farm module
std_dev <- 0.02 * 10^4  # Standard deviation
manure_application_rate <- 2  # kg/m2 - example
field_area <- 4000  # m2 - example

# Runoff coefficient distribution for the "Marins" catchment (https://www.sciencedirect.com/science/article/pii/S2095633921000885)
mean_runoff_coefficient <- 0.12
sd_runoff_coefficient <- 0.1

# Generate a random sample of runoff coefficients based on a normal distribution
n_simulations <- 1000
runoff_coefficients <- rnorm(n_simulations, mean_runoff_coefficient, sd_runoff_coefficient)

# Total manure applied in kg and convert to grams
total_manure_applied_grams <- manure_application_rate * field_area * 1000

# Generate random E. coli concentrations based on mean and standard deviation
e_coli_concentrations <- rnorm(n_simulations, mean_e_coli_concentration, std_dev)

# Bacteria partition coefficient for manure (BACTKDDB = 0.95, from Sowah)
bacteria_partition_coeff <- 0.95

# Apply bacterial parameters (BACT_SWF, WOF_P, also from Sowah paper)
bacteria_active_fraction <- 0.75  # BACT_SWF
wash_off_fraction <- 0.5  # WOF_P

# Apply rainfall data
total_rainfall_mm <- 100  # example rainfall Berlin July 2022
total_rainfall_m <- total_rainfall_mm / 1000

# Bathing site volume ( from O'Flaherty 2019)
bathing_site_volumes <- runif(n_simulations, min = 67500, max = 82500) * 1000  # Convert to liters

# Estimate E. coli transport for each simulation
results <- data.frame(runoff_coefficient = runoff_coefficients, e_coli_concentration = e_coli_concentrations, bathing_site_volume = bathing_site_volumes) %>%
  mutate(
    total_e_coli_load = e_coli_concentration * total_manure_applied_grams * bacteria_active_fraction,  # CFU in manure
    e_coli_free_for_transport = total_e_coli_load * (1 - bacteria_partition_coeff),  # Bacteria free for transport
    e_coli_washoff = e_coli_free_for_transport * wash_off_fraction,  # Bacteria washed off
    runoff_volume = total_rainfall_m * field_area * runoff_coefficient,  # Runoff volume in cubic meters
    e_coli_in_runoff = e_coli_washoff * runoff_coefficient,  # E. coli transported in runoff
    e_coli_concentration_in_runoff = e_coli_in_runoff / runoff_volume,  # E. coli concentration in runoff water (CFU/m3)
    e_coli_concentration_bathing_site = e_coli_in_runoff / bathing_site_volume  # E. coli concentration in bathing site (CFU/L)
  )

# Summarize results
summary_results <- results %>%
  summarise(
    mean_runoff_coefficient = mean(runoff_coefficient),
    sd_runoff_coefficient = sd(runoff_coefficient),
    mean_e_coli_concentration_in_runoff = mean(e_coli_concentration_in_runoff),
    sd_e_coli_concentration_in_runoff = sd(e_coli_concentration_in_runoff),
    mean_e_coli_concentration_bathing_site = mean(e_coli_concentration_bathing_site),
    sd_e_coli_concentration_bathing_site = sd(e_coli_concentration_bathing_site)
  )

# Display results
print(summary_results)
