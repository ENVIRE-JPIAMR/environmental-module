script_start_time <- Sys.time()

# lettuce model
# Load necessary libraries
library(tidyverse)
library(triangle)
library(mc2d)  # For rtri function

# O'Flaherty 2019
# Parameters for the calculations
f <- 0.00007642         # Unitless
k1 <- 4.45              # d^-1
k2 <- 0.06981           # d^-1
LC100g <- 100           # g

# Supplementary material table 3 - O'Flaherty 2019
Diffxy <-  0.826  # Post-harvesting washing with water only - Probability before expiry date

# Set number of simulations
n_simulations <- 100  # Adjust as needed

# Number of days
n_days <- nrow(river_summary_results_fresh)

# Define the days of lettuce planting and the duration of lettuce growth
planting_days <- c(1, 20, 50, 100, 150)  # Different planting days to simulate
lettuce_growth_duration <- 35  # Lettuce stays in the field for 35 days

set.seed(123)  # For reproducibility

# Function to simulate results for each day
simulate_day <- function(day, planting_day, irrigation_end_day, river_summary_results, sim_results) {
  if (day >= planting_day && day <= irrigation_end_day) {
    # Sample ECwater using mean_concentration and sd_concentration
    ECwater <- rnorm(1, mean = river_summary_results$mean_concentration[day], sd = river_summary_results$sd_concentration[day])
    # Ensure ECwater is non-negative
    ECwater <- max(ECwater, 0)
    
    # Sample Waterattach from Lognormal3 distribution (threshold at 0.006 ml/g)
    Waterattach <- 0.006 + rlnorm(1, meanlog = -4.75, sdlog = 0.50)  # ml/g
    
    # Calculate ECadhesion (CFU/g)
    ECadhesion <- ECwater * Waterattach  # CFU/g
    
    # Time since irrigation started
    t <- day - planting_day + 1  # Days since lettuce was planted
    
    # Calculate log R (Log CFU)
    log_R <- f * k1 * t + (1 - f) * k2 * t
    
    # Calculate ECdecay (CFU/g)
    ECdecay <- ECadhesion * 10^(-log_R)
    
    # Store ECadhesion and ECdecay before accumulation
    sim_results$ECadhesion[day] <- ECadhesion
    sim_results$ECdecay[day] <- ECdecay
    
    # Accumulate ECadhesion and ECdecay over the entire production period
    ECadhesion_accum <- sum(sim_results$ECadhesion[planting_day:day], na.rm = TRUE)
    ECdecay_accum <- sum(sim_results$ECdecay[planting_day:day], na.rm = TRUE)
    
    # Calculate PredECHarvest (CFU/g)
    PredECHarvest <- max(ECadhesion_accum - ECdecay_accum, 0)
    
    # Store the calculated value in the dataframe
    sim_results$PredECHarvest[day] <- PredECHarvest
  }
  return(sim_results)
}

# Function to simulate each scenario
simulate_scenario <- function(planting_day, n_days, n_simulations, river_summary_results) {
  irrigation_end_day <- planting_day + lettuce_growth_duration - 1
  exposure_results_list <- map(1:n_simulations, function(sim) {
    sim_results <- data.frame(
      day = 1:n_days,
      ECadhesion = NA_real_,
      ECdecay = NA_real_,
      PredECHarvest = NA_real_,
      PredECxy = NA_real_,
      PredECxyz = NA_real_,
      HExy = NA_real_,
      HExyz = NA_real_
    )
    for (day in 1:n_days) {
      sim_results <- simulate_day(day, planting_day, irrigation_end_day, river_summary_results, sim_results)
    }
    sim_results$simulation <- sim
    
    # Apply post-harvest treatment after harvest
    PredECHarvest <- sim_results$PredECHarvest[irrigation_end_day]
    PredECxy <- PredECHarvest * Diffxy  # CFU/g
    
    # Sample CWashing from triangular distribution (min 0.65, mode 0.99, max 0.99)
    CWashing <- rtriang(1, min = 0.65, mode = 0.99, max = 0.99)
    
    # Calculate PredECxyz (CFU/g)
    PredECxyz <- PredECxy * (1 - CWashing)
    
    # Calculate Human Exposure without washing (HExy)
    HExy <- LC100g * PredECxy  # CFU/100g
    
    # Calculate Human Exposure with washing (HExyz)
    HExyz <- LC100g * PredECxyz  # CFU/100g
    
    # Store post-harvest results in the dataframe
    sim_results$PredECxy[irrigation_end_day] <- PredECxy
    sim_results$PredECxyz[irrigation_end_day] <- PredECxyz
    sim_results$HExy[irrigation_end_day] <- HExy
    sim_results$HExyz[irrigation_end_day] <- HExyz
    
    return(sim_results)
  })
  scenario_results <- bind_rows(exposure_results_list)
  scenario_results$planting_day <- planting_day
  return(scenario_results)
}

# Simulate both scenarios: fresh and composted manure
scenarios <- list(
  fresh = river_summary_results_fresh,
  composted = river_summary_results_composted,
  stocked = river_summary_results_stocked,
  anaerobic_30 = river_summary_results_anaerobic_30,
  anaerobic_37 = river_summary_results_anaerobic_37
  
)

all_scenarios_results <- map_df(names(scenarios), function(scenario_name) {
  river_summary_results <- scenarios[[scenario_name]]
  scenario_results <- map(planting_days, simulate_scenario, n_days = n_days, n_simulations = n_simulations, river_summary_results = river_summary_results)
  scenario_results <- bind_rows(scenario_results)
  scenario_results$scenario <- scenario_name
  return(scenario_results)
})

# Create a data frame for exposure values at harvest (end of growth period) for all scenarios
exposure_at_harvest <- all_scenarios_results %>%
  filter(day %in% (planting_days + lettuce_growth_duration - 1)) %>%
  select(simulation, planting_day, scenario, HExy, HExyz) %>%
  pivot_longer(cols = c(HExy, HExyz), names_to = "ExposureType", values_to = "Exposure") %>%
  filter(!is.na(Exposure))  # Remove rows with NA values

# Plotting Human Exposure with and without washing as a boxplot for different planting days and scenarios
exposure_at_harvest %>%
  ggplot(aes(x = factor(planting_day), y = Exposure, fill = ExposureType)) +
  geom_boxplot() +
  facet_wrap(~scenario) +
  labs(
    title = "Human Exposure from Consuming 100 g of Lettuce",
    x = "Planting Day",
    y = "CFU/100g",
    fill = "Consumption Method"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c("HExy" = "red", "HExyz" = "blue"),
    labels = c("Not washing", "Washing")
  )


exposure_at_harvest %>%
  ggplot(aes(x = factor(planting_day), y = Exposure, fill = ExposureType)) +
  geom_boxplot(outlier.size = 1.5) +  # Slightly bigger outliers for visibility
  facet_wrap(~scenario, scales = "free_y", ncol = 2, strip.position = "top") +  # Two facets in one row
  labs(
    title = "Human Exposure from Consuming 100 g of Lettuce",
    x = "Planting Day",
    y = "CFU/100g",
    fill = "Consumption Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase title size
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14),  # Increase y-axis label size
    strip.text = element_text(size = 14),  # Increase facet title size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12),  # Increase legend text size
    panel.spacing = unit(1, "lines")  # Add more space between facets
  ) +
  scale_fill_manual(
    values = c("HExy" = "red", "HExyz" = "blue"),
    labels = c("Not washing", "Washing")
  )




exposure_at_harvest %>%
  ggplot(aes(x = factor(planting_day), y = Exposure, fill = ExposureType)) +
  geom_boxplot(outlier.size = 1.5) +  # Slightly bigger outliers for visibility
  facet_wrap(~scenario, ncol = 2, strip.position = "top", labeller = as_labeller(
    c(
      anaerobic_30 = "Anaerobic 30°C",
      anaerobic_37 = "Anaerobic 37°C",
      composted = "Composted",
      fresh = "Fresh",
      stocked = "Stocked"
    )
  )) +  # Update facet labels
  labs(
    title = "Human Exposure from Consuming 100 g of Lettuce",
    x = "Planting Day",
    y = "CFU/100g",
    fill = "Consumption Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase title size
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14),  # Increase y-axis label size
    strip.text = element_text(size = 14),  # Increase facet title size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12),  # Increase legend text size
    panel.spacing = unit(1, "lines")  # Add more space between facets
  ) +
  scale_fill_manual(
    values = c("HExy" = "red", "HExyz" = "blue"),
    labels = c("Not washing", "Washing")
  )


# Calculate reduction in exposure for each intervention compared to fresh manure
exposure_reduction <- all_scenarios_results %>%
  filter(day %in% (planting_days + lettuce_growth_duration - 1)) %>%
  group_by(simulation, planting_day, scenario) %>%
  summarise(HExy = mean(HExy, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = scenario, values_from = HExy) %>%
  mutate(
    composted_reduction = 100 * (fresh - composted) / fresh,
    stocked_reduction = 100 * (fresh - stocked) / fresh,
    anaerobic_30_reduction = 100 * (fresh - anaerobic_30) / fresh,
    anaerobic_37_reduction = 100 * (fresh - anaerobic_37) / fresh
  ) %>%
  select(simulation, planting_day, composted_reduction, stocked_reduction, anaerobic_30_reduction, anaerobic_37_reduction) %>%
  pivot_longer(cols = c(composted_reduction, stocked_reduction, anaerobic_30_reduction, anaerobic_37_reduction),
               names_to = "Condition", values_to = "Reduction")

# Plot reduction in human exposure compared to fresh manure using bar plot with error bars
exposure_summary <- exposure_reduction %>%
  group_by(planting_day, Condition) %>%
  summarise(
    mean_reduction = mean(Reduction, na.rm = TRUE),
    sd_reduction = sd(Reduction, na.rm = TRUE),
    .groups = 'drop'
  )

# plot for presentation

exposure_summary %>%
  ggplot(aes(x = factor(planting_day), y = mean_reduction, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_reduction - sd_reduction, ymax = mean_reduction + sd_reduction),
                position = position_dodge(width = 0.8), width = 0.3) +
  labs(
    x = "Planting Day",
    y = "Reduction (%)",
    fill = "Intervention"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_fill_manual(
    values = c("composted_reduction" = "blue",
               "stocked_reduction" = "green",
               "anaerobic_30_reduction" = "purple",
               "anaerobic_37_reduction" = "orange"),
    labels = c("composted_reduction" = "Composted",
               "stocked_reduction" = "Stocked",
               "anaerobic_30_reduction" = "Anaerobic 30°C",
               "anaerobic_37_reduction" = "Anaerobic 37°C")
  )

exposure_summary %>%
  filter(planting_day == 50) %>%  # Filter for planting day 50
  ggplot(aes(x = "", y = mean_reduction, fill = Condition)) +  # Set x-axis to an empty string
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_reduction - sd_reduction, ymax = mean_reduction + sd_reduction),
                position = position_dodge(width = 0.8), width = 0.3) +
  labs(
    x = "",  # Remove x-axis label
    y = "Reduction (%)",
    fill = "Intervention"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_blank(),  # Hide x-axis title
    axis.text.x = element_blank(),  # Hide x-axis text labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  scale_fill_manual(
    values = c("composted_reduction" = "blue",
               "stocked_reduction" = "green",
               "anaerobic_30_reduction" = "purple",
               "anaerobic_37_reduction" = "orange"),
    labels = c("composted_reduction" = "Composted",
               "stocked_reduction" = "Stocked",
               "anaerobic_30_reduction" = "Anaerobic 30°C",
               "anaerobic_37_reduction" = "Anaerobic 37°C"))
