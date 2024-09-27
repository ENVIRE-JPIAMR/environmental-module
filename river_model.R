library(tidyverse)
library(triangle)
### RIVER MODEL ###

# Rainfall events - if it rains 5 time
rainfall_events <- tibble(
  day = c(10, 30, 50, 80, 150),  # Example days when rainfall occurs
  rainfall_mm = c(200, 200, 200, 200, 200)  # Example rainfall amounts in mm
)

# Rainfall events - if it rains every day
rainfall_events <- tibble(
  day = 1:200,  # Days from 1 to 200
  rainfall_mm = rep(200, 200)  # Use the value 200 mm for all 200 days
)

# Input parameters for the SWAT-like model
runoff_coefficients <- runif(n_simulations, min = 0.05, max = 0.2)

## From Sowah 2020
bacteria_partition_coeff <- 0.95 # BACTKDDB - Fraction of bacteria bound to manure
bacteria_active_fraction <- 0.75  # BACT_SWF - Fraction of bacteria that are active
wash_off_fraction <- 0.5  # WOF_P - Fraction of bacteria washed off during runoff

## From O'Flaherty 2020 (Swimming model)
bathing_site_volumes <- runif(n_simulations, min = 67500, max = 82500) * 1000  # Convert to liters

 # Environmental parameters for decay rate (Mancini's Equation)
T <- runif(1, 21, 28)  # Temperature of river water (°C)
seawater <- runif(1, 0.035, 0.075)  # Salinity (%)
IA <- rtriangle(1, 17.32, 25.38, 22.73)  # Solar radiation levels (ly/hr)
et <- runif(1, 0.26, 0.31)  # Light extinction coefficient (m⁻¹)
H <- runif(1, 0.5, 6)  # Depth of water (m)

# Mancini's equation for decay rate
k <- 0.8 + 0.006 * seawater * 1.07^(T - 20) + (IA * (1 - exp(-et * H))) / (et * H)

# River flow and velocity parameters
Dist_x <- runif(1, 800, 1000)  # Distance from source to river mouth (m)
F_River_x <- rtriangle(1, 67500, 82500, 75000)  # Flow rate of river (m³/day)
Width_x <- rtriangle(1, 4.33, 21.71, 14.00)  # Width of river (m)
Vel_x <- F_River_x / (Width_x * H)  # Velocity of river (m/s)

# Exposure time in river
t_x <- Dist_x / Vel_x  # Exposure time (s)

# Create a full dataset for all 200 days
all_days <- tibble(day = 0:200)

# Initialize total E. coli in the bathing site and track bacteria added each day
total_ecoli_bathing_site <- 0

# Concentration in the bathing site with decay applied using Mancini's equation
results_with_rainfall <- all_days %>%
  left_join(results_summary, by = "day") %>%
  left_join(rainfall_events, by = "day") %>%
  mutate(
    rainfall_mm = replace_na(rainfall_mm, 0),  # Set rainfall to 0 for non-rainfall days
    total_rainfall_m = rainfall_mm / 1000,  # Convert rainfall to meters
    e_coli_free_for_transport = ifelse(is.na(mean_total_e_coli_on_field), 0, mean_total_e_coli_on_field * (1 - bacteria_partition_coeff)),  # Bacteria free for transport
    e_coli_washoff = e_coli_free_for_transport * wash_off_fraction,  # Bacteria washed off
    runoff_volume = total_rainfall_m * field_area * mean(runoff_coefficients),  # Mean runoff volume in cubic meters
    e_coli_in_runoff = ifelse(rainfall_mm > 0, e_coli_washoff, 0),  # Only add E. coli to runoff on rainfall days
    
   EC_Env_x = e_coli_in_runoff * exp(-k * t_x),  # E. coli concentration after decay in river
    
    
    bacteria_added_bathing_site = ifelse(rainfall_mm > 0, EC_Env_x / mean(bathing_site_volumes), 0),  # Bacteria added to the bathing site
    
    # Cumulative E. coli in the bathing site (accounting for decay)
    total_ecoli_bathing_site = lag(total_ecoli_bathing_site, default = 0) + bacteria_added_bathing_site,  # Total E. coli in the bathing site
    
    # Apply Mancini decay to the total E. coli in the bathing site
    total_ecoli_bathing_site_after_decay = total_ecoli_bathing_site * exp(-k * t_x), 
    
    # E. coli concentration in the bathing site
    e_coli_concentration_bathing_site = total_ecoli_bathing_site_after_decay / mean(bathing_site_volumes)  # E. coli concentration in bathing site (CFU/L)
  )

# Plot 
ggplot(results_with_rainfall, aes(x = day, y = e_coli_concentration_bathing_site)) +
  geom_line(color = "blue") +
  labs(title = "Cumulative E. coli Concentration in Bathing Site with Decay Over 200 Days",
       x = "Day",
       y = "Cumulative E. coli Concentration in Bathing Site (CFU/L)") +
  theme_minimal()

