library(tidyverse)
library(triangle)
set.seed(123)  # For reproducibility

# Define number of simulations
n_simulations <- 100

# Rainfall events - if it rains every day
rainfall_events <- tibble(
  day = 1:200,  # Days from 1 to 200
  rainfall_mm = rep(50, 200)  # Use the value 50 mm for all 200 days
)

# Fixed parameters
bacteria_partition_coeff <- 0.95  # Fraction of bacteria bound to manure
bacteria_active_fraction <- 0.75  # Fraction of bacteria that are active
wash_off_fraction <- 0.5  # Fraction of bacteria washed off during runoff
field_area <- 4000
total_ecoli_bathing_site_after_decay <- 0
# Simulate 100 runs using stochastic parameters
simulation_results <- map_dfr(1:n_simulations, function(sim) {
  
  # Environmental parameters for decay rate (Mancini's Equation)
  T <- runif(1, 21, 28)  # River water temperature (°C)
  seawater <- runif(1, 0.035, 0.075)  # Salinity (%)
  IA <- rtriangle(1, 17.32, 25.38, 22.73)  # Solar radiation levels (ly/hr)
  et <- runif(1, 0.26, 0.31)  # Light extinction coefficient (m⁻¹)
  H <- runif(1, 0.5, 6)  # Water depth (m)
  
  ## From O'Flaherty 2020 (Swimming model)
  bathing_site_volumes <- runif(1, min = 67500, max = 82500) * 1000  # Convert to liters
  
  # Parameters for the SWAT-like model
  runoff_coefficients <- runif(1, 0.05, 0.2)  # Random runoff coefficient for this simulation
  
  # Mancini's equation for decay rate
  k <- 0.8 + 0.006 * seawater * 1.07^(T - 20) + (IA * (1 - exp(-et * H))) / (et * H)
  
  # River flow and velocity parameters
  Dist_x <- runif(1, 800, 1000)  # Distance from source to river mouth (m)
  F_River_x <- rtriangle(1, 67500, 82500, 75000)  # Flow rate of river (m³/day)
  Width_x <- rtriangle(1, 4.33, 21.71, 14.00)  # Width of river (m)
  Vel_x <- F_River_x / (Width_x * H)  # River velocity (m/s)
  t_x <- Dist_x / Vel_x  # Exposure time in river (s)
  Fr_Day_x <- rtriangle(1, 67500, 82500)
  
  # Merging with rainfall events & soil model 
  # Calculating concentration of E. coli in bathing site

  all_days <- tibble(day = 1:200) %>%
    # Merge with results_summary to bring in mean_total_e_coli_on_field
    left_join(results_summary %>% select(day, mean_total_e_coli_on_field), by = "day") %>%
    left_join(rainfall_events, by = "day") %>%
    mutate(
      # Set rainfall to 0 for non-rainfall days
      rainfall_mm = replace_na(rainfall_mm, 0),  # Set rainfall to 0 for non-rainfall days
      total_rainfall_m = rainfall_mm / 1000,  # Convert rainfall to meters
      e_coli_free_for_transport = ifelse(is.na(mean_total_e_coli_on_field), 0, mean_total_e_coli_on_field * (1 - bacteria_partition_coeff)),  # Bacteria free for transport
      e_coli_washoff = e_coli_free_for_transport * wash_off_fraction,  # Bacteria washed off
      runoff_volume = total_rainfall_m * field_area * runoff_coefficients,  # Mean runoff volume in cubic meters
      e_coli_in_runoff = ifelse(rainfall_mm > 0, e_coli_washoff, 0),  # Only add E. coli to runoff on rainfall days
      # Cumulative E. coli in the bathing site (accounting for decay)
      total_ecoli_bathing_site = lag(total_ecoli_bathing_site_after_decay, default = 0) + e_coli_in_runoff,  # Total E. coli in the bathing site
      
      # Apply Mancini decay to the total E. coli in the bathing site
      total_ecoli_bathing_site_after_decay = total_ecoli_bathing_site * exp(-k * t_x), 
      
      # E. coli concentration in the bathing site
      e_coli_concentration_bathing_site = total_ecoli_bathing_site_after_decay  / mean(bathing_site_volumes)  # E. coli concentration in bathing site (CFU/L)
    )
  
    
  # Return relevant results for this simulation
  all_days %>%
    select(day, e_coli_concentration_bathing_site) %>%
    mutate(simulation = sim)  # Tag each result with the simulation number
})

summary_results <- simulation_results %>%
  group_by(day) %>%
  summarize(
    mean_concentration = mean(e_coli_concentration_bathing_site, na.rm = TRUE),
    sd_concentration = sd(e_coli_concentration_bathing_site, na.rm = TRUE)
  )

# Plot results: Mean and variability across simulations
ggplot(simulation_results, aes(x = day, y = e_coli_concentration_bathing_site * 1000, color = as.factor(simulation))) +
  geom_line() +
  labs(title = "E. coli Concentration in Bathing Site (Mean ± SD over 100 simulations)",
       x = "Day",
       y = "CFU/L") +
  theme_minimal() +
  guides(color = "none")



  
# Plot results: Mean and variability across simulations
ggplot(summary_results, aes(x = day)) +
  geom_line(aes(y = mean_concentration * 1000), color = "blue") +
  geom_ribbon(aes(ymin = (mean_concentration - sd_concentration) * 1000, 
                  ymax = (mean_concentration + sd_concentration) * 1000), 
              alpha = 0.2) +
  labs(title = "E. coli Concentration in Bathing Site (Mean ± SD over 1000 simulations)",
       x = "Day",
       y = "CFU/L") +
  theme_minimal()


ggplot() +
  geom_line(data = simulation_results, 
            aes(x = day, y = e_coli_concentration_bathing_site * 1000, group = simulation), 
            color = "lightgray", alpha = 0.2) +
  geom_line(data = summary_results, 
            aes(x = day, y = mean_concentration * 1000), 
            color = "blue", size = 1) +
  # Overlay ribbon for mean ± SD
  geom_ribbon(data = summary_results, 
              aes(x = day, ymin = (mean_concentration - sd_concentration) * 1000, 
                  ymax = (mean_concentration + sd_concentration) * 1000), 
              fill = "blue", alpha = 0.2) +
  labs(title = "E. coli Concentration in Bathing Site (Mean ± SD over 100 simulations)",
       x = "Day", 
       y = "CFU/L") +
  theme_minimal()

