# Load the required package for the triangular distribution
library(mc2d)

# Define two scenarios: fresh and composted manure
scenarios <- list(
  fresh = river_summary_results_fresh,
  composted = river_summary_results_composted,
  stocked = river_summary_results_stocked,
  anaerobic_30 = river_summary_results_anaerobic_30,
  anaerobic_37 = river_summary_results_anaerobic_37
)

# Add human exposure calculations to each scenario with 100 simulations
swimming_simulation_results <- map_df(names(scenarios), function(scenario_name) {
  simulation_results <- scenarios[[scenario_name]]
  scenario_results <- map_df(1:1000, function(sim) {
    simulation_results %>%
      rowwise() %>%
      mutate(
        # Generate HC_y for each exposure scenario
        # Uniform distribution
        HC_SwimAdult = runif(1, min = 0, max = 70.67),         # Uniform distribution
        HC_SwimChild = runif(1, min = 0, max = 205.33),        # Uniform distribution
        
        # Convert E. coli concentration to CFU/mL, ensuring non-negative values
        EC_Bathing_mL = max(0, rnorm(1, mean = mean_concentration, sd = sd_concentration)),
        
        # Compute HE_xy for each scenario (HE_xy = EC_Bathing_x * HC_y)
        
        HE_SwimAdult = EC_Bathing_mL * HC_SwimAdult,
        HE_SwimChild = EC_Bathing_mL * HC_SwimChild
      ) %>%
      ungroup()
  })
  scenario_results$scenario <- scenario_name
  return(scenario_results)
})

# Summarize HE_xy for each day and exposure scenario
swimming_exposure_summary <- swimming_simulation_results %>%
  group_by(day, scenario) %>%
  summarize(
    
    mean_HE_SwimAdult = mean(HE_SwimAdult, na.rm = TRUE),
    sd_HE_SwimAdult = sd(HE_SwimAdult, na.rm = TRUE),
    mean_HE_SwimChild = mean(HE_SwimChild, na.rm = TRUE),
    sd_HE_SwimChild = sd(HE_SwimChild, na.rm = TRUE)
  )

# View the summarized exposure results
print(swimming_exposure_summary)

# Reshape the exposure_summary dataframe to long format
swimming_exposure_long <- swimming_exposure_summary %>%
  pivot_longer(
    cols = -c(day, scenario),
    names_to = c(".value", "scenario_type"),
    names_pattern = "(mean_HE|sd_HE)_(.*)"
  )

# View the reshaped data
head(swimming_exposure_long)

# Filter data for specific days
swim <- swimming_exposure_summary %>% filter(day %in% c(1, 100))
print(swim)

# Create the plot
swimming_exposure_long %>%
  filter(day <= 100) %>%
  ggplot(aes(x = day, y = mean_HE, color = scenario_type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_HE - sd_HE, ymax = mean_HE + sd_HE, fill = scenario_type), alpha = 0.5) +  # Increase alpha for better visibility
  facet_wrap(~scenario, scales = "free_y") +
  labs(
    title = "Human Exposure from Recreational Swimming",
    x = "Day",
    y = "CFU/event",
    color = "Exposure Scenario", 
    fill = "Exposure Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis")

library(scales)

swimming_exposure_long %>%
  filter(day <= 100) %>%
  ggplot(aes(x = day, y = mean_HE, color = scenario_type)) +
  geom_line(size = 1) +  # Plot only the mean lines
  facet_wrap(~scenario, scales = "free_y") +
  labs(
    title = "Human Exposure from Recreational Swimming",
    x = "Day",
    y = "CFU/event",
    color = "Exposure Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_viridis_d(option = "viridis") +
  scale_y_continuous(labels = scales::scientific_format())  # Correct scientific notation


swimming_exposure_long %>%
  filter(day <= 100) %>%
  filter(scenario_type == "SwimAdult") %>%
  ggplot(aes(x = day, y = mean_HE, color = scenario)) +
  geom_line(size = 1) +  # Plot only the mean lines
  labs(
    title = "Human Exposure from Recreational Swimming",
    x = "Day",
    y = "CFU/event",
    color = "Exposure Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(
    values = c(
      "anaerobic_30" = "purple",
      "anaerobic_37" = "blue",
      "composted" = "green",
      "fresh" = "darkgreen",
      "stocked" = "gold"
    ),
    labels = c(
      "anaerobic_30" = "Anaerobic 30°C",
      "anaerobic_37" = "Anaerobic 37°C",
      "composted" = "Composted",
      "fresh" = "Fresh",
      "stocked" = "Stocked"
    )
  ) +
  scale_y_continuous(labels = scales::scientific_format())  # Correct scientific notation



swimming_exposure_long %>%
  filter(day <= 100) %>%
  filter(scenario_type == "SwimAdult") %>%
  ggplot(aes(x = day, y = mean_HE, color = scenario, linetype = scenario)) +
  geom_line(size = 1) +  # Use different linetypes to distinguish scenarios
  labs(
    title = "Human Exposure from Recreational Swimming",
    x = "Day",
    y = "CFU/event",
    color = "Exposure Scenario",
    linetype = "Exposure Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(
    values = c(
      "anaerobic_30" = "purple",
      "anaerobic_37" = "blue",
      "composted" = "green",
      "fresh" = "darkgreen",
      "stocked" = "gold"
    ),
    labels = c(
      "anaerobic_30" = "Anaerobic 30°C",
      "anaerobic_37" = "Anaerobic 37°C",
      "composted" = "Composted",
      "fresh" = "Fresh",
      "stocked" = "Stocked"
    )
  ) +
  scale_y_continuous(labels = scales::scientific_format())  # Correct scientific notation

swimming_exposure_long %>%
  filter(day <= 100) %>%
  filter(scenario_type == "SwimAdult") %>%
  ggplot(aes(x = day, y = mean_HE, color = scenario, linetype = scenario)) +
  geom_line(size = 1) +  # Use different linetypes to distinguish scenarios
  labs(
    title = "Human Exposure from Recreational Swimming",
    x = "Day",
    y = "CFU/event",
    color = "Exposure Scenario",
    linetype = "Exposure Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(
    values = c(
      "anaerobic_30" = "purple",
      "anaerobic_37" = "blue",
      "composted" = "green",
      "fresh" = "darkgreen",
      "stocked" = "gold"
    ),
    labels = c(
      "anaerobic_30" = "Anaerobic 30°C",
      "anaerobic_37" = "Anaerobic 37°C",
      "composted" = "Composted",
      "fresh" = "Fresh",
      "stocked" = "Stocked"
    )
  ) +
  scale_y_continuous(labels = scales::scientific_format()) +  # Correct scientific notation
  guides(linetype = "none")  # Remove the lower right legend


### Calculate reduction in exposure for SwimAdult compared to the fresh manure baseline
swimming_exposure_long_adult <- swimming_exposure_long %>%
  filter(scenario_type == "SwimAdult") %>%
  group_by(day) %>%
  mutate(
    baseline_HE = mean_HE[scenario == "fresh"],
    reduction = if_else(scenario == "fresh", 0, (baseline_HE - mean_HE) / baseline_HE * 100)
  ) %>%
  ungroup()

# Create the reduction plot
swimming_exposure_long_adult %>%
  filter(day <= 100) %>%
  ggplot(aes(x = day, y = reduction, color = scenario)) +
  geom_line(size = 1) +
  labs(
    title = "Reduction in Human Exposure Compared to Fresh Manure",
    x = "Day",
    y = "Reduction (%)",
    color = "Intervention"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(
    values = c("composted" = "blue",
               "stocked" = "red",
               "anaerobic_30" = "purple",
               "anaerobic_37" = "orange"),
    labels = c("composted" = "Composted",
               "stocked" = "Stocked",
               "anaerobic_30" = "Anaerobic 30°C",
               "anaerobic_37" = "Anaerobic 37°C")) + 
            theme(
        legend.position = "right",
        plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      )
 


