#####################################################################
# Integrated Model for Farm, Soil, River, Swimming, & Lettuce Modules #
# Authors: Nunzio Sarnino and Subhasish Basask                        #
# Date: November 2025                                                 #
# Purpose: To estimate human exposure to ESBL-producing E. coli       #
#          originating from broiler production via two exposure paths:#
#          - Consumption of lettuce                                  #
#          - Recreational swimming in contaminated water              #
#####################################################################

#############################
# Setup and Libraries  #
#############################
# Load required libraries; note that many tidyverse functions are used.
library(tidyverse)    # Data manipulation and plotting
library(mc2d)         # For PERT and other distributions
library(future)       # Parallelization
library(furrr)        # Parallel map using future
library(scales)       # For extra ggplot scaling functions
library(here)         # For relative file paths
library(purrr)        # Functional programming tools
library(nls2)         # Non-linear least squares
library(readxl)       # Reading Excel files
library(triangle)     # Triangular distributions
library(gridExtra)    # For grid graphics
library(writexl)      # To export the results in excel
library(lhs)          # Latin Hypercube Sampling
library(sensitivity)  # for pcc()
library(progressr)    # for progress bar





###################################
# FARM MODEL AND SIMULATIONS  #
###################################
# This section sets up the farm model that simulates the dynamics of ESBL E. coli
# in a broiler flock. It reads model parameters from a CSV input file and creates
# functions to initialize the flock, update bacterial growth,
# simulate transmission, feces production, ingestion, excretion, and environmental decay.

### Load Model Inputs ###
load_inputs <- function(input_manual = list()) {
  # Construct file path relative to the working directory.
  input_file <- here("inputs.csv")
  
  # Read the CSV which must have columns: Variable, Value, and Type.
  df_read <- read.csv(input_file, header = TRUE, sep = ';')
  
  # Parse "OBJECT" type entries (more complex structures)
  input_objects <- list(
    water_consum.min  = eval(parse(text = df_read$Value[df_read$Variable == "water_consum.min"])),
    water_consum.max  = eval(parse(text = df_read$Value[df_read$Variable == "water_consum.max"])),
    water_consum.mean = eval(parse(text = df_read$Value[df_read$Variable == "water_consum.mean"])),
    weight            = eval(parse(text = df_read$Value[df_read$Variable == "weight"])),
    daily_gain        = eval(parse(text = df_read$Value[df_read$Variable == "daily_gain"])),
    daily_intake      = eval(parse(text = df_read$Value[df_read$Variable == "daily_intake"])),
    phages_reduction  = eval(parse(text = df_read$Value[df_read$Variable == "phages_reduction"]))
  )
  
  # Parse numeric entries from the CSV that are not "OBJECT"
  idx_double <- df_read$Type != "OBJECT"
  df_double  <- data.frame(
    id = df_read$Variable[idx_double],
    val = unlist(lapply(df_read$Value[idx_double], function(x) eval(parse(text = x))))
  )
  named_vector <- with(df_double, setNames(val, id))
  input_doubles <- lapply(split(named_vector, names(named_vector)), unname)
  
  # Combine and allow manual overrides if provided
  input_list <- modifyList(c(input_objects, input_doubles), input_manual)
  
  return(input_list)
}
input_list <- load_inputs()


### Define the Farm Module Environment ###
# The function new.farm_module() creates a dedicated environment that stores model parameters
# and a set of functions to simulate daily dynamics in the flock.
new.farm_module <- function(input_list = load_inputs()){
  fm <- new.env()
  fm$params <- input_list
  
  # Initialize the flock data frame with healthy and infected birds.
  fm$initialize_df <- function() {
    healthy <- tribble(~infection_duration, ~age, -1, 1)
    sick    <- tribble(~infection_duration, ~age, 1, 1)
    
    # Calculate the number of animals based on density, target weight, and farm size.
    n_animals <- round(fm$params$farm_density / fm$params$target_weight * fm$params$farm_size)
    n_animals_infected <- rbinom(1, n_animals, fm$params$prevalence)
    n_animals_healthy  <- n_animals - n_animals_infected
    
    animals <- rbind(
      sick[rep(1, n_animals_infected), ],
      healthy[rep(1, n_animals_healthy), ]
    ) %>%
      mutate(
        feces_gut         = 0,
        sum_feces_gut     = 0,
        sum_feces_env     = 0,
        sum_feces_cont_env= 0,
        C_esbl_gut        = ifelse(infection_duration == -1, fm$params$esbl.min, fm$params$esbl.max),
        C_sum_esbl_env    = 0,
        B_infection_status= infection_duration != -1,
        ingested_feces    = 0,
        C_esbl_excreted   = 0
      )
    return(animals)
  }
  
  # Logistic growth model for ESBL E. coli in infected birds.
  fm$logistic_growth <- function(animals) {
    K <- fm$params$K * animals$feces_gut
    r <- 10 ^ runif(1, fm$params$r.min, fm$params$r.max)
    animals <- animals %>%
      mutate(
        C_esbl_gut = ifelse(
          infection_duration != -1,
          K / (1 + ((K - C_esbl_gut) / C_esbl_gut) * exp(-r * infection_duration)),
          C_esbl_gut
        )
      )
    return(animals)
  }
  
  # Transmission of ESBL from litter to broiler
  fm$transmission <- function(animals) {
    esbl_conc_feces <- ifelse(
      sum(animals$sum_feces_cont_env) > 0,
      sum(animals$C_sum_esbl_env) / sum(animals$sum_feces_cont_env),
      0
    )
    foi <- ifelse(esbl_conc_feces == 0, 0, log10(esbl_conc_feces) * fm$params$beta.mean)
    
    N_susceptible  <- sum(animals$infection_duration == -1)
    N_new_infected <- round(N_susceptible * (1 - exp(-foi * fm$params$Dt)))
    N_new_infected <- min(max(N_new_infected, 0), N_susceptible)
    
    if (N_new_infected > 0) {
      susceptible_indices <- which(animals$infection_duration == -1)
      newly_infected_indices <- sample(susceptible_indices, N_new_infected, replace = FALSE)
      animals$infection_duration[newly_infected_indices] <- 0
    }
    
    total_esbl <- sum(esbl_conc_feces * animals$ingested_feces[animals$infection_duration == 0])
    total_doner_feces <- sum(animals$ingested_feces[animals$infection_duration > 0])
    
    animals <- animals %>%
      mutate(
        B_infection_status = infection_duration != -1,
        C_esbl_gut = ifelse(
          infection_duration == 0,
          C_esbl_gut + esbl_conc_feces * ingested_feces,
          C_esbl_gut
        ),
        C_sum_esbl_env = ifelse(
          infection_duration > 0,
          C_sum_esbl_env - (total_esbl * ingested_feces) / total_doner_feces,
          C_sum_esbl_env
        )
      )
    return(animals)
  }
  
  # Calculate daily feces production
  fm$feces_production <- function(day, animals) {
    feces_amount <- runif(
      nrow(animals),
      min = fm$params$water_consum.min[day],
      max = fm$params$water_consum.max[day]
    ) * fm$params$water_reduction +
      fm$params$daily_intake[day] - fm$params$daily_gain[day]
    
    animals$feces_gut     <- feces_amount
    animals$sum_feces_gut <- animals$sum_feces_gut + feces_amount
    return(animals)
  }
  
  # Determine daily feces ingestion
  fm$feces_ingestion <- function(day, animals) {
    animals$ingested_feces <- rep(fm$params$daily_intake[day] * fm$params$ingestion_rate, nrow(animals))
    return(animals)
  }
  
  # Excrete ESBL from gut to the environment
  fm$excretion <- function(animals) {
    animals <- animals %>%
      mutate(
        C_esbl_excreted = ifelse(infection_duration != -1, C_esbl_gut * fm$params$e_rate, 0),
        C_sum_esbl_env  = C_sum_esbl_env + C_esbl_excreted,
        C_esbl_gut      = C_esbl_gut - C_esbl_excreted,
        sum_feces_env   = sum_feces_gut,
        sum_feces_cont_env = ifelse(infection_duration != -1, sum_feces_cont_env + feces_gut, sum_feces_cont_env)
      )
    return(animals)
  }
  
  # Apply environmental decay to the ESBL in the environment
  fm$environmental_decay <- function(animals) {
    animals <- animals %>%
      mutate(
        C_sum_esbl_env = ifelse(infection_duration != -1, C_sum_esbl_env * (1 - fm$params$ed_rate), C_sum_esbl_env)
      )
    return(animals)
  }
  
  # (Optional) Phages intervention function on litter
  fm$phages_intervention_litter <- function(day, animals) {
    animals <- animals %>%
      mutate(
        C_sum_esbl_env = C_sum_esbl_env * (1 - fm$params$phages_reduction[day])
      )
    return(animals)
  }
  
  # Run the simulation over the production period
  fm$run <- function(animals, total_feces_input = NA, total_C_esbl_input = NA, 
                     phages = FALSE, thinning = FALSE, day, until) {
    day_idx <- fm$params$day.min
    animals_full <- list()
    animals$day <- rep(day_idx, nrow(animals))
    animals_full[[day_idx]] <- animals
    
    while (day_idx < fm$params$day.max) {
      animals <- fm$feces_production(day_idx, animals)
      animals <- fm$feces_ingestion(day_idx, animals)
      animals <- fm$excretion(animals)
      animals <- fm$logistic_growth(animals)
      animals <- fm$transmission(animals)
      animals <- fm$environmental_decay(animals)
      
      # Update age and infection duration.
      animals$age <- animals$age + 1
      animals$infection_duration <- ifelse(animals$infection_duration != -1,
                                           animals$infection_duration + 1,
                                           -1)
      
      if (phages) {
        animals <- fm$phages_intervention_litter(day_idx, animals)
      }
      
      if (thinning) {
        # Thinning intervention logic would be added here.
      }
      
      day_idx <- day_idx + 1
      animals$day <- rep(day_idx, nrow(animals))
      animals_full[[day_idx]] <- animals
    }
    return(list(animals_full))
  }
  
  return(fm)
}

### Running the Farm Model ###
run_single_simulation <- function(phages = FALSE, 
                                  total_feces_input = NA, 
                                  total_C_esbl_input = NA, 
                                  thinning = FALSE) {
  fm <- new.farm_module(input_list)
  simulation_result <- fm$run(
    animals = fm$initialize_df(), 
    total_feces_input, total_C_esbl_input, 
    phages, thinning, day = NA, until = NA
  )
  return(simulation_result)
}

# Run multiple simulations and combine results.
n_simulations <- 1000
multiple_simulations <- map(1:n_simulations, ~ run_single_simulation(phages = FALSE))
all_results <- map_dfr(multiple_simulations, ~ bind_rows(.x[[1]], .id = "day"), .id = "simulation") %>%
  mutate(day = as.numeric(day))

# Summarize key metrics (e.g., cumulative ESBL in environment).
results <- all_results %>%
  group_by(day, simulation) %>%
  summarise(
    C_sum_esbl_env = sum(C_sum_esbl_env, na.rm = TRUE),
    sum_feces_env  = sum(sum_feces_env),
    C_esbl_gut     = sum(C_esbl_gut),
    mean_gut       = mean(C_esbl_gut),
    .groups = "drop"
  )

# Calculate manure concentration (CFU per gram) using litter mass and farm size.
litter <- input_list$litter_mass * input_list$farm_size
results <- results %>% mutate(cfu_g = C_sum_esbl_env / (sum_feces_env + litter))
farm_cfu <- results %>% 
  group_by(day) %>% 
  summarize(
    mean_cfu_g = mean(cfu_g, na.rm = TRUE),
    sd_cfu_g   = sd(cfu_g, na.rm = TRUE),
    lower_cfu_g   = quantile(cfu_g, probs = 0.025, na.rm = TRUE),
    upper_cfu_g   = quantile(cfu_g, probs = 0.975, na.rm = TRUE),
    .groups    = "drop"
  )
write

# For subsequent modules, extract the final day statistics
mean_last_day <- farm_cfu %>% filter(day == 36) %>% pull(mean_cfu_g)
sd_last_day   <- farm_cfu %>% filter(day == 36) %>% pull(sd_cfu_g)

# Additional output: Infection prevalence over time.
all_animals <- all_results
infected <- all_animals %>% group_by(day, simulation) %>% select(infection_duration)
daily_infected <- infected %>%
  filter(infection_duration != -1) %>%
  group_by(day, simulation) %>%
  summarise(num_infected = n(), .groups = "drop")
total_animals <- all_animals %>% group_by(day, simulation) %>% summarise(total = n(), .groups = "drop")
daily_infected <- left_join(daily_infected, total_animals, by = c("day", "simulation")) %>%
  mutate(prop_infected = num_infected / total)
mean_prop_infected <- daily_infected %>% group_by(day) %>% 
  summarise(mean_prop_infected = mean(prop_infected), 
            sd_prop        = sd  (prop_infected),
            lower_prop    = quantile(prop_infected, 0.025, na.rm = TRUE),
            upper_prop    = quantile(prop_infected, 0.975, na.rm = TRUE),
            .groups = "drop")


# Optional - Delete all_animals and all_results to free memory

rm(all_animals, all_results)
gc()  # Optional: triggers garbage collection to free up memory immediately



#########################
# SOIL MODEL        #
#########################
# This module estimates the decay of ESBL E. coli in soil following manure application.
# It uses experimental data with an exponential decay model fitted via non-linear least squares.

# From Sharma et al. (2020) supplementary material
data <- read_excel(here("sharma_2.xlsx"))

# Define an exponential decay function
exp_decay <- function(concentration, decay_rate, days) {
  concentration * exp(-decay_rate * days)
}

# Fit the decay model to the experimental data
fit <- nls(y ~ exp_decay(concentration, decay_rate, dpi), 
           data = data, 
           start = list(concentration = max(data$y), decay_rate = 0.001),
           control = nls.control(maxiter = 100))
decay_rate_est <- coef(fit)[["decay_rate"]]
decay_rate_est
data


# Define manure treatment scenarios (reductions based on intervention)
e_coli_concentration <- mean_last_day   # from farm model
std_dev             <- sd_last_day


# Define input parameters
manure_application_rate <- 2        # kg/m2
days <- 0:200                       # simulation time horizon (days)
field_area <- 4046.86                  # m² # 1 Ha
soil_depth_m <- 0.1                 # topsoil depth (m)
soil_bulk_density <- 1.47           # g/cm³ -> convert to kg/m³ below
soil_bulk_density_kg_per_m3 <- soil_bulk_density * 1000
soil_volume_per_m2 <- soil_depth_m * 1  
soil_mass_per_m2 <- soil_volume_per_m2 * soil_bulk_density_kg_per_m3  
soil_mass_per_m2_grams <- soil_mass_per_m2 * 1000
manure_applied_per_m2_grams <- manure_application_rate * 1000

# Monte Carlo simulations for initial soil E. coli per scenario
n_simulations <- 1000

# Function for  initial soil E. coli
simulate_soil <- function(e_coli_concentration, std_dev) {
  conc_per_gram <- rnorm(n_simulations, e_coli_concentration, std_dev)
  initial_e_coli_per_m2 <- conc_per_gram * manure_applied_per_m2_grams
  e_coli_cfu_g_m2 <- initial_e_coli_per_m2 / (soil_mass_per_m2_grams + manure_applied_per_m2_grams)
  list(initial = initial_e_coli_per_m2, cfu = e_coli_cfu_g_m2)
}


# Function to create decay data given initial concentrations
create_decay_data <- function(initial_e_coli_per_m2, e_coli_cfu_g_m2) {
  tibble(
    simulation = rep(1:n_simulations, each = length(days)),
    day = rep(days, times = n_simulations),
    e_coli_concentration_per_m2 = unlist(lapply(initial_e_coli_per_m2, 
                                                function(init_conc) exp_decay(init_conc, decay_rate_est, days))),
    total_e_coli_on_field = e_coli_concentration_per_m2 * field_area, 
    e_coli_cfu_per_g_m2   = unlist(lapply(e_coli_cfu_g_m2, 
                                          function(init_conc) exp_decay(init_conc, decay_rate_est, days)))
  )
}


# Summarize decay results by day
summarize_decay_results <- function(decay_data) {
  decay_data %>%
    group_by(day) %>%
    summarise(
      mean_e_coli_concentration_per_m2 = mean(e_coli_concentration_per_m2),
      std_e_coli_concentration_per_m2  = sd(e_coli_concentration_per_m2),
      lower_cfu_m2 = quantile(e_coli_concentration_per_m2, 0.025, na.rm = TRUE),
      upper_cfu_m2 = quantile(e_coli_concentration_per_m2, 0.975, na.rm = TRUE),
      mean_total_e_coli_on_field       = mean(total_e_coli_on_field),
      std_total_e_coli_on_field        = sd(total_e_coli_on_field),
      mean_e_coli_cfu_per_g_m2         = mean(e_coli_cfu_per_g_m2),
      std_e_coli_cfu_per_g_m2          = sd(e_coli_cfu_per_g_m2),
      .groups = "drop"
    )
}

# Run baseline soil simulation
soil_baseline <- simulate_soil(e_coli_concentration, std_dev)
decay_data_baseline <- create_decay_data(soil_baseline$initial, soil_baseline$cfu)
soil_results_summary <- summarize_decay_results(decay_data_baseline)




#########################
# RIVER MODEL       #
#########################
# This model uses a simplified SWAT-like approach to simulate the transport of ESBL E. coli
# from the field to a river (or bathing site). It accumulates daily load and applies a decay rate
# estimated via Mancini's equation.

# Global river parameters
bacteria_partition_coeff <- 0.95  
wash_off_fraction        <- 0.5    
n_simulations            <- 1000    

run_river_simulation <- function(soil_results_summary) {
  # Ensure days start from 1
  daily_soil <- soil_results_summary %>% filter(day >= 1) %>% arrange(day)
  
  sim_results <- map_dfr(1:n_simulations, function(sim) {
    # Sample environmental parameters
    T <- runif(1, 21, 28)                    # water temperature (°C)
    seawater <- runif(1, 0.035, 0.075)         # salinity
    IA <- rtriangle(1, 17.32, 25.38, 22.73)    # solar radiation
    et <- runif(1, 0.26, 0.31)                # light extinction coefficient
    H <- runif(1, 0.5, 6)                     # water depth (m)
    # Bathing site volume in liters
    bathing_site_volume <- runif(1, 67500, 82500) * 1000
    
    # Decay rate calculated using Mancini's equation
    k <- 0.8 + 0.006 * seawater * 1.07^(T - 20) +
      (IA * (1 - exp(-et * H))) / (et * H)
    
    daily_data <- daily_soil %>%
      select(day, mean_total_e_coli_on_field) %>%
      mutate(
        e_coli_free_for_transport = ifelse(is.na(mean_total_e_coli_on_field),
                                           0,
                                           mean_total_e_coli_on_field * (1 - bacteria_partition_coeff)),
        e_coli_washoff = e_coli_free_for_transport * wash_off_fraction
      ) %>% arrange(day)
    
    n_days <- nrow(daily_data)
    cumulative_load <- numeric(n_days)
    cumulative_load[1] <- daily_data$e_coli_washoff[1]
    if(n_days > 1) {
      for (i in 2:n_days) {
        cumulative_load[i] <- cumulative_load[i - 1] * exp(-k) + daily_data$e_coli_washoff[i]
      }
    }
    
    daily_data <- daily_data %>%
      mutate(
        total_ecoli_bathing_site = cumulative_load,
        e_coli_concentration_bathing_site = total_ecoli_bathing_site / bathing_site_volume,
        simulation = sim
      ) %>%
      select(simulation, day, e_coli_concentration_bathing_site)
    
    daily_data
  })
  
  # Summarize simulation results across Monte Carlo iterations
  river_summary <- sim_results %>%
    group_by(day) %>%
    summarize(
      mean_concentration = mean(e_coli_concentration_bathing_site, na.rm = TRUE) / 1000,
      sd_concentration   = sd(e_coli_concentration_bathing_site, na.rm = TRUE) / 1000,
      lower_CFU_mL  = quantile(e_coli_concentration_bathing_site / 1000, 0.025, na.rm = TRUE),
      upper_CFU_mL  = quantile(e_coli_concentration_bathing_site / 1000, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  return(river_summary)
}

# Run river simulations for each soil scenario.
simulation_results <- run_river_simulation(soil_results_summary)



###############################
# SWIMMING EXPOSURE MODEL #
###############################
# This module estimates human exposure via recreational swimming in the contaminated water.
# It performs Monte Carlo simulations that account for variability in water flow, volume,
# and ingestion rates.




# Define scenario names and associate with river simulation results.
swimming_simulation_results <- map_df(1:1000, function(sim) {
  # 1) draw the one‐time downstream flow & volume for this simulation
  Fr_Day_sim <- rtriangle(1, a = 67500, b = 82500, c = (67500 + 82500)/2)  # L/day
  V_sim      <- runif(1, 67500, 82500)      # L
  

  # 2) compute a single dilution factor
  dil_factor <- Fr_Day_sim / V_sim
  
  # 3) now loop over days *without* rowwise()
  simulation_results %>%
    mutate(
      simulation    = sim,
      EC_river_mL   = pmax(0, rnorm(n(), mean_concentration, sd_concentration)),
      EC_swim_mL    = (EC_river_mL * 1000 * dil_factor) / 1000,  # CFU/mL
      HC_SwimAdult  = runif(n(), 0, 70.67),
      HC_SwimChild  = runif(n(), 0, 205.33),
      HE_SwimAdult  = EC_swim_mL * HC_SwimAdult,
      HE_SwimChild  = EC_swim_mL * HC_SwimChild
    )
})


# Summarize daily exposures.
swimming_exposure_summary <- swimming_simulation_results %>%
  group_by(day) %>%
  summarize(
    mean_HE_SwimAdult = mean(HE_SwimAdult, na.rm = TRUE),
    sd_HE_SwimAdult   = sd(HE_SwimAdult, na.rm = TRUE),
    lower_swim_adult  = quantile(HE_SwimAdult, 0.025, na.rm = TRUE),
    upper_swim_adult  = quantile(HE_SwimAdult, 0.975, na.rm = TRUE),
    mean_HE_SwimChild = mean(HE_SwimChild, na.rm = TRUE),
    sd_HE_SwimChild   = sd(HE_SwimChild, na.rm = TRUE),
    lower_swim_child  = quantile(HE_SwimChild, 0.025, na.rm = TRUE),
    upper_swim_child  = quantile(HE_SwimChild, 0.975, na.rm = TRUE),
    .groups = "drop"
  )



##################################
# LETTUCE EXPOSURE MODULE    #
##################################
# This module estimates human exposure from consuming 100g of lettuce irrigated with
# contaminated river water. The model simulates daily adhesion of E. coli to lettuce, its
# subsequent decay via a biphasic model, and the final exposure after post-harvest washing.

# Key parameters from O’Flaherty et al., 2019.
f      <- 0.00007642  # fraction for biphasic decay
k1     <- 4.45        # fast decay rate (d^-1)
k2     <- 0.06981     # slow decay rate (d^-1)
LC100g <- 100         # ingestion amount (100g lettuce)
Diffxy <- 0.826       #  post-harvest washing reduction
# User defined parameters
n_days <- nrow(simulation_results)  # number of simulation days
planting_days <- c(1, 20, 50, 100, 150)     # available planting days
lettuce_growth_duration <- 35              # days from planting to harvest


# Use the river summary results from the baseline simulation only
river_summary_results <- simulation_results

# Function: simulate_day
simulate_day <- function(day, planting_day, irrigation_end_day, river_summary_results, sim_results) {
  if (day >= planting_day && day <= irrigation_end_day) {
    ECwater <- rnorm(1, mean = river_summary_results$mean_concentration[day], 
                     sd   = river_summary_results$sd_concentration[day])
    ECwater <- max(ECwater, 0)
    
    # Volume of water attachment (ml/g) from a lognormal offset.
    Waterattach <- 0.006 + rlnorm(1, meanlog = -4.75, sdlog = 0.50)
    ECadhesion  <- ECwater * Waterattach  # CFU per gram
    
    t <- day - planting_day + 1
    log_R   <- f * k1 * t + (1 - f) * k2 * t
    ECdecay <- ECadhesion * 10^(-log_R)
    
    sim_results$ECadhesion[day] <- ECadhesion
    sim_results$ECdecay[day]    <- ECdecay
    
    ECadhesion_accum <- sum(sim_results$ECadhesion[planting_day:day], na.rm = TRUE)
    ECdecay_accum    <- sum(sim_results$ECdecay[planting_day:day], na.rm = TRUE)
    PredECHarvest    <- max(ECadhesion_accum - ECdecay_accum, 0)
    
    sim_results$PredECHarvest[day] <- PredECHarvest
  }
  return(sim_results)
}

# Function: simulate_scenario
simulate_scenario <- function(planting_day, n_days, n_simulations, river_summary_results) {
  irrigation_end_day <- planting_day + lettuce_growth_duration - 1
  
  exposure_results_list <- map(1:n_simulations, function(sim) {
    sim_results <- data.frame(
      day           = 1:n_days,
      ECadhesion    = NA_real_,
      ECdecay       = NA_real_,
      PredECHarvest = NA_real_,
      PredECxy      = NA_real_,
      PredECxyz     = NA_real_,
      HExy          = NA_real_,
      HExyz         = NA_real_
    )
    for (day in 1:n_days) {
      sim_results <- simulate_day(day, planting_day, irrigation_end_day, 
                                  river_summary_results, sim_results)
    }
    sim_results$simulation <- sim
    # Post-harvest washing treatment
    PredECHarvest <- sim_results$PredECHarvest[irrigation_end_day]
    PredECxy      <- PredECHarvest * Diffxy
    CWashing      <- rtriang(1, min = 0.65, mode = 0.99, max = 0.99)
    PredECxyz     <- PredECxy * (1 - CWashing)
   
    
    sim_results$PredECxy[irrigation_end_day]  <- PredECxy
    sim_results$PredECxyz[irrigation_end_day] <- PredECxyz
    sim_results$HExy[irrigation_end_day]  <- LC100g * PredECxy
    sim_results$HExyz[irrigation_end_day] <- LC100g * PredECxyz
    sim_results
  })
  scenario_results <- bind_rows(exposure_results_list)
  scenario_results$planting_day <- planting_day
  scenario_results
}

# Run lettuce model for all planting days
all_scenarios_results <- map_df(planting_days, function(planting_day) {
  simulate_scenario(planting_day, n_days, 100, river_summary_results)
})

# Extract final exposures at harvest
exposure_at_harvest <- all_scenarios_results %>%
  filter(day %in% (planting_days + lettuce_growth_duration - 1)) %>%
  select(simulation, planting_day, HExy, HExyz) %>%
  pivot_longer(cols = c(HExy, HExyz), names_to = "ExposureType", values_to = "Exposure") %>%
  filter(!is.na(Exposure))
exposure_at_harvest


lettuce_summary <- exposure_at_harvest %>%
  group_by(planting_day, ExposureType) %>%
  summarise(mean_exp = mean(Exposure), sd_exp = sd(Exposure), lower_CFU   = quantile(Exposure, 0.025, na.rm = TRUE),
            upper_CFU   = quantile(Exposure, 0.975, na.rm = TRUE), .groups="drop") %>%
  mutate(Treatment = recode(ExposureType,
                            HExy  = "Not Washing",
                            HExyz = "Washing"))
lettuce_summary



##### UNIFORM PLOT THEME #####


my_theme <- theme_minimal() +
  theme(
    plot.title      = element_text(size = 14, face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10),
    panel.spacing   = unit(1, "lines")
  )

##### 1. Farm: CFU per g litter ####
p1 <- ggplot(farm_cfu, aes(x = day, y = mean_cfu_g)) +
  geom_line(size = 1, color = "darkred") +
  geom_ribbon(aes(ymin = lower_cfu_g, ymax = upper_cfu_g), alpha = 0.2) +
    labs(
        x = "Day",
    y = "CFU/g"
  ) +
  my_theme



##### 2. Farm: Proportion Infected #####
p2 <- ggplot(mean_prop_infected, aes(x = day, y = mean_prop_infected)) +
  geom_line(size = 1, color = "steelblue") +
  geom_ribbon(aes(ymin = lower_prop, ymax = upper_prop), alpha = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Day",
    y = "Prevalence"
  ) +
  my_theme

##### 3. Soil: E. coli per g soil #####
p3 <- ggplot(soil_results_summary, aes(x = day, y = mean_e_coli_concentration_per_m2)) +
  geom_line(size = 1, color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_cfu_m2, ymax = upper_cfu_m2), alpha = 0.2) +
    labs(
    
    x = "Day",
    y = expression(CFU / m^2)
  ) +
  my_theme



##### 4. River: Concentration #####
p4 <- ggplot(simulation_results, aes(x = day, y = mean_concentration)) +
  geom_line(size = 1, color = "darkblue") +
  geom_ribbon(aes(ymin = lower_CFU_mL, ymax = upper_CFU_mL), alpha = 0.2) +
  
  labs(
    x = "Day",
    y = "CFU/mL"
  ) +
  my_theme


##### 5. Swimming Exposure #####
swim_long <- swimming_exposure_summary %>%
  pivot_longer(c(mean_HE_SwimAdult, mean_HE_SwimChild),
               names_to = "group", values_to = "exposure") %>%
  mutate(group = recode(group,
                        mean_HE_SwimAdult = "Adult",
                        mean_HE_SwimChild = "Child"))
p5 <- ggplot(swim_long, aes(x = day, y = exposure, color = group)) +
  geom_line(size = 1) +
  labs(
    title = "Swimming Exposure Over Time",
    x = "Day",
    y = "CFU per Event",
    color = "Age Group"
  ) +
  my_theme



##### 6. Lettuce Exposure #####

p9 <- ggplot(lettuce_summary, aes(x = planting_day, y = mean_exp, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_CFU, ymax = upper_CFU), width = 5)  +
  labs(
       x = "Planting Day",
    y = "CFU per 100 g",
    color = "Treatment"
  ) +
  my_theme

#########################
# DOSE–RESPONSE & DALY MODULE
#########################

##  Parameters  ---------------------------------------------
# From Heida et al. 2025, paper and Table S2 
r_dose_resp    <- 1.95e-6                     # dose-response 
f_upec      <- 0.10                           # fraction of E. coli that are uPEC
# Conditional probabilities (
p_uti_given_col<- 0.067                       # P(UTI | urinary colon)
# extra parameters in the loop (swimming_risk and lettuce_risk)

## Helper functions  ---------------------------------------
dose_to_risk <- function(dose){
  1 - exp(-r_dose_resp * dose * f_upec)
}


##  Swimming pathway  ---------------------------------------
swimming_risk <- swimming_simulation_results %>% 
  mutate(
    p_colon_urine  = runif(n(), 0.35, 0.46) ,                       # P(urinary colon | gut colon) 
    daly_per_case  = runif(n(), 3.7, 12.84)     ,                # DALY x case⁻¹ 
    
    Risk_GI_Adult  = dose_to_risk(HE_SwimAdult),
    Risk_GI_Child  = dose_to_risk(HE_SwimChild),
    Risk_UTI_Adult = dose_to_risk(HE_SwimAdult) * p_colon_urine * p_uti_given_col,
    Risk_UTI_Child = dose_to_risk(HE_SwimChild) * p_colon_urine * p_uti_given_col,
    DALY_Adult     = Risk_UTI_Adult * daly_per_case,
    DALY_Child     = Risk_UTI_Child * daly_per_case
  )

# mean (± sd if needed) per day  
swimming_risk_summary <- swimming_risk %>% 
  group_by(day) %>% 
  summarise(across(
    c(Risk_GI_Adult:DALY_Child),
    list(mean = mean, sd = sd),
    .names = "{.fn}_{.col}"
  ),
  .groups = "drop"
  )

##  Lettuce consumption pathway  -----------------------------
lettuce_risk <- exposure_at_harvest %>% 
  mutate(
    p_colon_urine  = runif(n(), 0.35, 0.46) ,                       # P(urinary colon | gut colon) 
    daly_per_case  = runif(n(), 3.7, 12.84)     ,                # DALY x case⁻¹ 
    Risk_GI  = dose_to_risk(Exposure),
    Risk_UTI = dose_to_risk(Exposure) * p_colon_urine * p_uti_given_col,
    DALY     = Risk_UTI * daly_per_case
  )

lettuce_risk_summary <- lettuce_risk %>% 
  group_by(planting_day, ExposureType) %>% 
  summarise(
    mean_Risk_GI  = mean(Risk_GI),
    sd_Risk_GI = sd(Risk_GI),
    mean_Risk_UTI = mean(Risk_UTI),
    sd_Risk_UTI = sd(Risk_UTI),
    mean_DALY     = mean(DALY),
    sd_DALY = sd(DALY),
    .groups = "drop"
  )



lettuce_risk_summary <- lettuce_risk %>% 
  group_by(planting_day, ExposureType) %>% 
  summarise(
    mean_Risk_GI   = mean(Risk_GI),
    sd_Risk_GI     = sd(Risk_GI),
    lower_Risk_GI  = quantile(Risk_GI, 0.025),
    upper_Risk_GI  = quantile(Risk_GI, 0.975),
    
    mean_Risk_UTI  = mean(Risk_UTI),
    sd_Risk_UTI    = sd(Risk_UTI),
    lower_Risk_UTI = quantile(Risk_UTI, 0.025),
    upper_Risk_UTI = quantile(Risk_UTI, 0.975),
    
    mean_DALY      = mean(DALY),
    sd_DALY        = sd(DALY),
    lower_DALY     = quantile(DALY, 0.025),
    upper_DALY     = quantile(DALY, 0.975),
    
    .groups = "drop"
  )



### PLOTS ###

#  Risk + DALY long format ---------------------------------------------
risk_long <- swimming_risk_summary %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to    = c("metric", "group"),
    names_pattern = "mean_([A-Za-z_]+)_(Adult|Child)",
    values_to   = "value"
  ) %>%
  mutate(
    group = recode(group, Adult = "Adult", Child = "Child"),
    metric = recode(metric,
                    "Risk_GI"  = "Gut Colonization",
                    "Risk_UTI" = "UTI",
                    "DALY"     = "DALY")
  )

#  Gut-colonization risk -----------------------------------------------
p6 <- risk_long %>%
  filter(metric == "Gut Colonization") %>%
  ggplot(aes(day, value, color = group)) +
  geom_line(size = 1) +
  labs(
    title = "Risk of Gut Colonization",
    x     = "Day",
    y     = "Probability",
    color = "Age Group"
  ) +
  my_theme

# UTI risk ------------------------------------------------------------
p7 <- risk_long %>%
  filter(metric == "UTI") %>%
  ggplot(aes(day, value, color = group)) +
  geom_line(size = 1) +
  labs(
    title = "Risk of UTI",
    x     = "Day",
    y     = "Probability",
    color = "Age Group"
  ) +
  my_theme

# DALY per swim -------------------------------------------------------
p8 <- risk_long %>%
  filter(metric == "DALY") %>%
  ggplot(aes(day, value, color = group)) +
  geom_line(size = 1) +
  labs(
    title = "DALY Lost per Swim",
    x     = "Day",
    y     = "DALY per Event",
    color = "Age Group"
  ) +
  my_theme

lettuce_long <- lettuce_risk_summary %>%
  pivot_longer(
    cols = c(mean_Risk_GI, mean_Risk_UTI, mean_DALY),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    Treatment = recode(ExposureType,
                       HExy  = "Not Washing",
                       HExyz = "Washing"),
    metric = recode(metric,
                    mean_Risk_GI  = "Gut Colonization",
                    mean_Risk_UTI = "UTI",
                    mean_DALY     = "DALY")
  )

#  Gut‐colonization risk plot ------------------------------------------
p10 <- lettuce_long %>%
  filter(metric == "Gut Colonization") %>%
  ggplot(aes(x = planting_day, y = value, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Risk of Gut Colonization",
    x     = "Planting Day",
    y     = "Probability",
    color = "Treatment"
  ) +
  my_theme

#  UTI risk plot -------------------------------------------------------
p11 <- lettuce_long %>%
  filter(metric == "UTI") %>%
  ggplot(aes(x = planting_day, y = value, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Risk of UTI",
    x     = "Planting Day",
    y     = "Probability",
    color = "Treatment"
  ) +
  my_theme

# DALY per serving plot ----------------------------------------------
p12 <- lettuce_long %>%
  filter(metric == "DALY") %>%
  ggplot(aes(x = planting_day, y = value, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "DALY Lost per 100 g Serving",
    x     = "Planting Day",
    y     = "DALY per Serving",
    color = "Treatment"
  ) +
  my_theme

# Write results
write_xlsx(farm_cfu,                  "farm_cfu.xlsx")
write_xlsx(mean_prop_infected,        "mean_prop_infected.xlsx")
write_xlsx(soil_results_summary,      "soil_results_summary.xlsx")
write_xlsx(simulation_results,        "river_results_summary.xlsx")
write_xlsx(lettuce_summary,    "lettuce_summary.xlsx")
write_xlsx(swimming_exposure_summary, "swimming_exposure_summary.xlsx")
write_xlsx(lettuce_risk_summary,      "lettuce_risk_summary.xlsx")
write_xlsx(swimming_risk_summary,      "swimming_risk_summary.xlsx")
write_xlsx(lettuce_risk,                "lettuce_risk.xlsx")

########################################
#  Sensitivity Analysis (Farm Module)  #
########################################

# set up progressr
handlers("txtprogressbar")


#  Define parameters to vary and their ranges -------------
param_defs <- list(
  farm_density          = c(min = 0.5 * input_list$farm_density,      max = 1.5 * input_list$farm_density),
  farm_size             = c(min = 0.5 * input_list$farm_size,         max = 1.5 * input_list$farm_size),
  target_weight         = c(min = 0.5 * input_list$target_weight,      max = 1.5 * input_list$target_weight),
  prevalence            = c(min = 0.5 * input_list$prevalence,         max = 1.5 * input_list$prevalence),
  K                     = c(min = 0.5 * input_list$K,                 max = 1.5 * input_list$K),
  beta.mean             = c(min = 0.5 * input_list$beta.mean,         max = 1.5 * input_list$beta.mean),
  ed_rate               = c(min = 0.5 * input_list$ed_rate,           max = 1.5 * input_list$ed_rate),
  r.max               = c(min = 0.5 * input_list$r.max,             max = 1.5 * input_list$r.max),
   water_reduction = c(min = 0.5 * input_list$water_reduction,  max = 1.5 * input_list$water_reduction),
  ingestion_rate        = c(min = 0.5 * input_list$ingestion_rate,    max = 1.5 * input_list$ingestion_rate),
  e_rate                = c(min = 0.5 * input_list$e_rate,           max = 1.5 * input_list$e_rate),
  litter_mass             = c(min=0.5*input_list$litter_mass,       max=1.5*input_list$litter_mass)
)

#   Latin Hypercube Sampling ------------- 
n_sens <- 1000
lhs_mat <- randomLHS(n_sens, length(param_defs))
colnames(lhs_mat) <- names(param_defs)
param_samples <- as_tibble(lhs_mat)
for(p in names(param_defs)) {        
  param_samples[[p]] <- param_defs[[p]]["min"] +
    lhs_mat[,p] * (param_defs[[p]]["max"] - param_defs[[p]]["min"])
}

# 8.4  Plan parallel execution -------------
plan(multisession, workers = parallel::detectCores() - 1)


output_metric <- NULL

with_progress({
  p <- progressor(along = seq_len(n_sens))
#  Parallel sensitivity loop -------------
# Metric: CFU/g on day 36
output_metric <- future_map_dbl(seq_len(n_sens), function(i) {
  p()  # advance bar
  # 1) Override farm inputs
  mods <- list(
    farm_density    = param_samples$farm_density[i],
    farm_size       = param_samples$farm_size[i],
    target_weight   = param_samples$target_weight[i],
    prevalence      = param_samples$prevalence[i],
    K               = param_samples$K[i],
   r.max           = param_samples$`r.max`[i],
    beta.mean       = param_samples$beta.mean[i],
    ed_rate         = param_samples$ed_rate[i],
      water_reduction= param_samples$water_reduction[i],
    ingestion_rate  = param_samples$ingestion_rate[i],
    e_rate          = param_samples$e_rate[i]
  )
  input_mod <- modifyList(input_list, mods)
  
  #  Run simulation
  fm <- new.farm_module(input_mod)
  sim_list <- fm$run(animals = fm$initialize_df(), phages=FALSE, thinning=FALSE)[[1]]
  
  #  Compute CFU/g by day
  df <- bind_rows(sim_list, .id = "day") %>%
    mutate(day = as.integer(day)) %>%
    group_by(day) %>%
    summarise(
      total_ESBL  = sum(C_sum_esbl_env, na.rm = TRUE),
      total_feces = sum(sum_feces_env, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cfu_g = total_ESBL / (total_feces + input_mod$litter_mass * input_mod$farm_size)
    )
  
  #  Extract day 36 CFU/g
  df %>% filter(day == 36) %>% pull(cfu_g)
}, .options = furrr_options(seed = TRUE)
)
})

# return to sequential plan
plan(sequential)

#  Compute PRCC and plot results -------------
prcc_res_fm <- pcc(
  X     = param_samples,
  y     = output_metric,
  rank  = TRUE,
  nboot = 100


prcc_df_fm <- as.data.frame(prcc_res_fm$PRCC)
colnames(prcc_df_fm) <- c("prcc", "bias", "std.error", "min.ci", "max.ci")
prcc_df_fm$parameter <- rownames(prcc_df_fm)
rownames(prcc_df_fm) <- NULL

# Customize label for plotting
label_map <- c(
  farm_density   = "Farm stocking density",
  farm_size      = "Farm size",
  target_weight  = "Target bird weight",
  prevalence     = "Initial flock prevalence",
  K              = "Gut carrying capacity",
  beta.mean      = "Transmission coefficient",
  ed_rate        = "On-farm decay rate",
  r.max          = "Maximum growth rate",
  water_reduction= "Water loss rate",
  ingestion_rate = "Bird ingestion rate",
  e_rate         = "Shedding rate",
  litter_mass    = "Litter mass"
)


# Matching equation symbols (plotmath expressions)
symbol_map  <- c(
  farm_density   = expression(rho[farm]),
  farm_size      = expression(A),
  target_weight  = expression(w[tar]),
  prevalence     = expression(p[init]),
  K              = expression(K),
  beta.mean      = expression(beta),
  ed_rate        = expression(phi),
  r.max          = expression(r[max]),
  water_reduction= expression(alpha),
  ingestion_rate = expression(rho[ingest]),
  e_rate         = expression(epsilon),
  litter_mass    = expression(L)
)

#Combine “text (symbol)”; each element becomes an expression
nice_labels <- map2(label_map, symbol_map,
                    ~ bquote(.(.x)~"("~.(.y)~")")) |>
  unlist()         # keep a named vector

# Plot
prcc_df_fm$parameter <- factor(
  prcc_df_fm$parameter,
  levels = prcc_df_fm$parameter[order(prcc_df_fm$prcc)]
)

ggplot(prcc_df_fm, aes(parameter, prcc)) +
  geom_col() +
  geom_errorbar(aes(ymin = min.ci, ymax = max.ci), width = .2) +
  coord_flip() +
  scale_x_discrete(labels = nice_labels) +   # <- new labels here
  labs(
    x = NULL,
    y = "Partial Rank Correlation\nCoefficient (PRCC)"
  ) +
  theme_minimal(base_size = 12)

########################################
# Sensitivity Analysis (Lettuce Exposure)
########################################
# start progress bar
handlers("txtprogressbar")

#Define which parameters to vary and their ranges 
param_defs <- list(
  # Farm‐module parameters (±50%)
  farm_density            = c(min = 0.5 * input_list$farm_density,           max = 1.5 * input_list$farm_density),
  target_weight           = c(min = 0.5 * input_list$target_weight,          max = 1.5 * input_list$target_weight),
  prevalence              = c(min = 0.5 * input_list$prevalence,             max = 1.5 * input_list$prevalence),
  K                       = c(min = 0.5 * input_list$K,                      max = 1.5 * input_list$K),
  beta.mean               = c(min = 0.5 * input_list$beta.mean,              max = 1.5 * input_list$beta.mean),
  ed_rate                 = c(min = 0.5 * input_list$ed_rate,                max = 1.5 * input_list$ed_rate),
  ingestion_rate          = c(min = 0.5 * input_list$ingestion_rate,         max = 1.5 * input_list$ingestion_rate),
  e_rate                  = c(min = 0.5 * input_list$e_rate,                 max = 1.5 * input_list$e_rate),
  litter_mass             = c(min = 0.5 * input_list$litter_mass,            max = 1.5 * input_list$litter_mass),
  # Soil‐module parameters
  decay_rate_est          = c(min = 0.5 * decay_rate_est,                    max = 1.5 * decay_rate_est),
  manure_application_rate = c(min = 0.5 * manure_application_rate,           max = 1.5 * manure_application_rate),
  soil_bulk_density       = c(min = 0.5 * soil_bulk_density,                 max = 1.5 * soil_bulk_density),
  # River‐module parameters
  bacteria_partition_coeff= c(min = 0.5 * bacteria_partition_coeff,           max = 1.5 * bacteria_partition_coeff),
  wash_off_fraction       = c(min = 0.5 * wash_off_fraction,                 max = 1.5 * wash_off_fraction),
  # Lettuce‐module parameters
  f                       = c(min = 0.5 * f,                                max = 1.5 * f),
  k1                      = c(min = 0.5 * k1,                               max = 1.5 * k1),
  k2                      = c(min = 0.5 * k2,                               max = 1.5 * k2),
  Diffxy                  = c(min = 0.5 * Diffxy,                           max = 1.5 * Diffxy)
)

#  Generate Latin Hypercube Samples -------------
n_sens   <- 1000
lhs_mat  <- randomLHS(n_sens, length(param_defs))
colnames(lhs_mat) <- names(param_defs)

param_samples <- as_tibble(lhs_mat)
for(p in names(param_defs)) {
  param_samples[[p]] <- param_defs[[p]]["min"] +
    lhs_mat[, p] * (param_defs[[p]]["max"] - param_defs[[p]]["min"])
}

#  Run full model chain in parallel and extract final lettuce exposure -------------
plant_day <- 50  # fixing planting day for sensitivity metric

# set up parallel plan
plan(multisession, workers = parallel::detectCores() - 1)

output_metric <- NULL

with_progress({
  p <- progressor(along = seq_len(n_sens))
  output_metric <- future_map_dbl(seq_len(n_sens), function(i) {
    p()
  # 1) Override farm inputs
  mods_farm <- list(
    farm_density   = param_samples$farm_density[i],
    target_weight  = param_samples$target_weight[i],
    prevalence     = param_samples$prevalence[i],
    K              = param_samples$K[i],
    beta.mean      = param_samples$beta.mean[i],
    ed_rate        = param_samples$ed_rate[i],
    ingestion_rate = param_samples$ingestion_rate[i],
    e_rate         = param_samples$e_rate[i]
  )
  input_mod <- modifyList(input_list, mods_farm)

  # 2) Assign module globals
  assign("litter_mass",               param_samples$litter_mass[i],            envir = .GlobalEnv)
  assign("decay_rate_est",            param_samples$decay_rate_est[i],         envir = .GlobalEnv)
  assign("manure_application_rate",   param_samples$manure_application_rate[i],envir = .GlobalEnv)
  assign("soil_bulk_density",         param_samples$soil_bulk_density[i],      envir = .GlobalEnv)
  assign("bacteria_partition_coeff",  param_samples$bacteria_partition_coeff[i],envir = .GlobalEnv)
  assign("wash_off_fraction",         param_samples$wash_off_fraction[i],      envir = .GlobalEnv)
  assign("f",                         param_samples$f[i],                      envir = .GlobalEnv)
  assign("k1",                        param_samples$k1[i],                     envir = .GlobalEnv)
  assign("k2",                        param_samples$k2[i],                     envir = .GlobalEnv)
  assign("Diffxy",                    param_samples$Diffxy[i],                 envir = .GlobalEnv)

  # -- Farm module
  fm      <- new.farm_module(input_mod)
  farm_out<- fm$run(animals = fm$initialize_df())[[1]]
  farm_df <- bind_rows(farm_out, .id = "day") %>%
               mutate(day = as.numeric(day))
  last_day<- max(farm_df$day)
  mean_last<- farm_df %>%
                filter(day == last_day) %>%
                summarise(mean_env = mean(C_sum_esbl_env /
                                          (sum_feces_env + litter_mass * input_mod$farm_size))) %>%
                pull(mean_env)
  sd_last  <- farm_df %>%
                filter(day == last_day) %>%
                summarise(sd_env = sd(C_sum_esbl_env /
                                       (sum_feces_env + litter_mass * input_mod$farm_size))) %>%
                pull(sd_env)

  # -- Soil module
  soil_base  <- simulate_soil(mean_last, sd_last)
  decay_data <- create_decay_data(soil_base$initial, soil_base$cfu)
  soil_sum   <- summarize_decay_results(decay_data)

  # -- River module
  river_sum <- run_river_simulation(soil_sum)

  # -- Lettuce module for chosen planting day
  scen      <- simulate_scenario(plant_day, n_days, 1000, river_sum)
  final_val <- scen %>%
                 filter(day == (plant_day + lettuce_growth_duration - 1)) %>%
                 summarise(mean_HExyz = mean(HExyz, na.rm = TRUE)) %>%
                 pull(mean_HExyz)

  final_val
} , .options = furrr_options(seed = TRUE)
)
})
# return to sequential plan
plan(sequential)

#  Compute PRCC and plot results -------------
prcc_res_lettuce <- pcc(
  X     = param_samples,
  y     = output_metric,
  rank  = TRUE,
  nboot = 100
)

print(prcc_res)
plot(prcc_res)

# turn into data.frame and label plot nicely
prcc_res_lettuce <- as.data.frame(prcc_res$PRCC)
colnames(prcc_res_lettuce) <- c("prcc", "bias", "std.error", "min.ci", "max.ci")
prcc_res_lettuce$parameter <- rownames(prcc_res_lettuce)
rownames(prcc_res_lettuce) <- NULL

# define labels
label_map <- c(
  farm_density            = "Farm stocking density",
  target_weight           = "Target weight",
  prevalence              = "Initial flock prevalence",
  K                       = "Gut carrying capacity",
  beta.mean               = "Transmission coefficient",
  ed_rate                 = "Environmental decay rate",
  ingestion_rate          = "Bird ingestion rate",
  e_rate                  = "Shedding rate",
  litter_mass             = "Litter mass",
  decay_rate_est          = "Soil decay rate",
  manure_application_rate = "Manure application",
  soil_bulk_density       = "Soil bulk density",
  bacteria_partition_coeff= "Partition coefficient",
  wash_off_fraction       = "Wash-off fraction",
  f                       = "Biphasic decay fraction",
  k1                      = "Fast decay rate",
  k2                      = "Slow decay rate",
  Diffxy                  = "Post-harvest wash"
)

## Plot-math symbols that will go in parentheses 
symbol_map <- c(
  farm_density             = expression(rho[farm]),
  target_weight            = expression(w[tar]),
  prevalence               = expression(p[init]),
  K                        = expression(K),
  beta.mean                = expression(beta),
  ed_rate                  = expression(phi),          # φ  (environmental decay)
  ingestion_rate           = expression(rho[ingest]),
  e_rate                   = expression(epsilon),
  litter_mass              = expression(L),
  decay_rate_est           = expression(lambda[s]),    # soil decay – change if needed
  manure_application_rate  = expression(M[app]),
  soil_bulk_density        = expression(rho[soil]),
  bacteria_partition_coeff = expression(K[d]),         # partition coefficient
  wash_off_fraction        = expression(omega),        # pick the letter you prefer
  f                        = expression(f),
  k1                       = expression(k[1]),
  k2                       = expression(k[2]),
  Diffxy                   = expression(D[wash])
)

##  Combine text + symbol into a single expression for each tick -----------
nice_labels <- map2(label_map, symbol_map,
                    ~ bquote(.(.x)~"("~.(.y)~")")) |>
  unlist()


prcc_res_lettuce$parameter <- factor(
  prcc_res_lettuce$parameter,
  levels = prcc_res_lettuce$parameter[order(prcc_res_lettuce$prcc)]
)

## Plot -
ggplot(prcc_res_lettuce, aes(parameter, prcc)) +
  geom_col() +
  geom_errorbar(aes(ymin = min.ci, ymax = max.ci), width = .2) +
  coord_flip() +
  scale_x_discrete(labels = nice_labels) +   # << the new labels
  labs(
    x = NULL,
    y = "Partial Rank Correlation\nCoefficient (PRCC)"
  ) +
  theme_minimal(base_size = 12)

########################################
# Sensitivity Analysis (Swimming Adult)
########################################


# set up progress bar
handlers("txtprogressbar")


#Define parameters and ranges 
param_defs <- list(
  # Farm‐module (±50%)
  farm_density            = c(min=0.5*input_list$farm_density,      max=1.5*input_list$farm_density),
  target_weight           = c(min=0.5*input_list$target_weight,     max=1.5*input_list$target_weight),
  prevalence              = c(min=0.5*input_list$prevalence,        max=1.5*input_list$prevalence),
  K                       = c(min=0.5*input_list$K,                 max=1.5*input_list$K),
  beta.mean               = c(min=0.5*input_list$beta.mean,         max=1.5*input_list$beta.mean),
  ed_rate                 = c(min=0.5*input_list$ed_rate,           max=1.5*input_list$ed_rate),
  ingestion_rate          = c(min=0.5*input_list$ingestion_rate,    max=1.5*input_list$ingestion_rate),
  e_rate                  = c(min=0.5*input_list$e_rate,            max=1.5*input_list$e_rate),
  litter_mass             = c(min=0.5*input_list$litter_mass,       max=1.5*input_list$litter_mass),
  
  # Soil‐module
  decay_rate_est          = c(min=0.5*decay_rate_est,               max=1.5*decay_rate_est),
  manure_application_rate = c(min=0.5*manure_application_rate,      max=1.5*manure_application_rate),
  soil_bulk_density       = c(min=0.5*soil_bulk_density,            max=1.5*soil_bulk_density),
  
  # River‐module
  bacteria_partition_coeff= c(min=0.5*bacteria_partition_coeff,      max=1.5*bacteria_partition_coeff),
  wash_off_fraction       = c(min=0.5*wash_off_fraction,            max=1.5*wash_off_fraction),
  
  # Swimming‐module
  flow_rate               = c(min = 0.5 * 67500,   # L/day
                              max = 1.5 * 82500),  # L/day
  bathing_site_volume     = c(min = 0.5 * 67500,   # L
                              max = 1.5 * 82500),  # L
  HC_Adult_max             = c(min = 0.5 * 70.67,                   max = 1.5 * 70.67)
)

# Latin Hypercube Sampling
n_sens  <- 1000
lhs_mat <- randomLHS(n_sens, length(param_defs))
colnames(lhs_mat) <- names(param_defs)

param_samples <- as_tibble(lhs_mat)
for (p in names(param_defs)) {
  param_samples[[p]] <- param_defs[[p]]["min"] +
    lhs_mat[,p] * (param_defs[[p]]["max"] - param_defs[[p]]["min"])
}

# 8.4  Parallel model runs with progress bar -------------
plan(multisession, workers = parallel::detectCores() - 1)

output_metric <- NULL
with_progress({
  p <- progressor(along = seq_len(n_sens))
  
  output_metric <- future_map_dbl(
    seq_len(n_sens),
    function(i) {
      p()  # advance bar
      
      # -- override farm & soil & river globals --
      mods_farm <- list(
        farm_density   = param_samples$farm_density[i],
        target_weight  = param_samples$target_weight[i],
        prevalence     = param_samples$prevalence[i],
        K              = param_samples$K[i],
        beta.mean      = param_samples$beta.mean[i],
        ed_rate        = param_samples$ed_rate[i],
        ingestion_rate = param_samples$ingestion_rate[i],
        e_rate         = param_samples$e_rate[i]
      )
      input_mod <- modifyList(input_list, mods_farm)
      assign("litter_mass",               param_samples$litter_mass[i],            envir = .GlobalEnv)
      assign("decay_rate_est",            param_samples$decay_rate_est[i],         envir = .GlobalEnv)
      assign("manure_application_rate",   param_samples$manure_application_rate[i],envir = .GlobalEnv)
      assign("soil_bulk_density",         param_samples$soil_bulk_density[i],      envir = .GlobalEnv)
      assign("bacteria_partition_coeff",  param_samples$bacteria_partition_coeff[i],envir = .GlobalEnv)
      assign("wash_off_fraction",         param_samples$wash_off_fraction[i],      envir = .GlobalEnv)
      
      # -- override swim params locally --
      Fr_i           <- param_samples$flow_rate[i]            # [L/day]
      V_i            <- param_samples$bathing_site_volume[i]  # [L]
      HC_adult_max_i  <- param_samples$HC_Adult_max[i]
      
      # 1) Farm module
      fm       <- new.farm_module(input_mod)
      farm_out <- fm$run(fm$initialize_df())[[1]]
      farm_df  <- bind_rows(farm_out, .id="day") %>% mutate(day = as.numeric(day))
      last_day <- max(farm_df$day)
      mean_last<- farm_df %>%
        filter(day == last_day) %>%
        summarise(mean_env = mean(C_sum_esbl_env /
                                    (sum_feces_env + litter_mass * input_mod$farm_size))) %>%
        pull(mean_env)
      sd_last  <- farm_df %>%
        filter(day == last_day) %>%
        summarise(sd_env = sd(C_sum_esbl_env /
                                (sum_feces_env + litter_mass * input_mod$farm_size))) %>%
        pull(sd_env)
      
      # 2) Soil module
      soil_base  <- simulate_soil(mean_last, sd_last)
      soil_sum   <- summarize_decay_results(create_decay_data(soil_base$initial, soil_base$cfu))
      
      # 3) River module
      river_sum  <- run_river_simulation(soil_sum)
      n_days_run <- nrow(river_sum)
      
      # 4) Swimming module (100 replicates)
      swim_df <- map_df(1:1000, function(sim) {
      Fr_sim <- Fr_i      # [L/day]
        V_sim  <- V_i       # [L]
        dilute <- Fr_sim / V_sim
        
        river_sum %>%
          mutate(simulation    = sim,
                 EC_river_mL   = pmax(0, rnorm(n(), mean_concentration, sd_concentration)),
                 EC_swim_mL    = EC_river_mL * dilute,
                 HC_SwimAdult  = runif(n(), 0, HC_adult_max_i),
                 HE_SwimAdult  = EC_swim_mL * HC_SwimAdult)
      })
      
      # average exposure over all days & sims
      mean(swim_df$HE_SwimAdult, na.rm = TRUE)
    },
    .options = furrr_options(seed = TRUE)
  )
})
plan(sequential)

# Compute PRCC and plot results -------------
prcc_res_swimming <- pcc(
  X     = param_samples,
  y     = output_metric,
  rank  = TRUE,
  nboot = 100
)


# turn into data.frame and label/plot nicely
prcc_res_swimming <- as.data.frame(prcc_res_swimming$PRCC)
colnames(prcc_res_swimming) <- c("prcc", "bias", "std.error", "min.ci", "max.ci")
prcc_res_swimming$parameter <- rownames(prcc_res_swimming)
rownames(prcc_res_swimming) <- NULL
# custom display names
label_map <- c(
  farm_density             = "Farm stocking density",
  target_weight            = "Target bird weight",
  prevalence               = "Initial flock prevalence",
  K                        = "Gut carrying capacity",
  beta.mean                = "Transmission coefficient",
  ed_rate                  = "On-farm decay rate",
  ingestion_rate           = "Bird ingestion rate",
  e_rate                   = "Shedding rate",
  litter_mass              = "Litter mass",
  decay_rate_est           = "Soil decay rate",
  manure_application_rate  = "Manure application rate",
  soil_bulk_density        = "Soil bulk density",
  bacteria_partition_coeff = "Partition coefficient",
  wash_off_fraction        = "Wash-off fraction",
  flow_rate                   = "Flow rate",
  bathing_site_volume          = "Bathing site volume",
  HC_Adult_max             = "Max ingestion volume"
)

prcc_res_swimming$parameter <- factor(
  prcc_res_swimming$parameter,
  levels = prcc_res_swimming$parameter[order(prcc_res_swimming$prcc)]
)

ggplot(prcc_df, aes(x = parameter, y = prcc)) +
  geom_col() +
  geom_errorbar(aes(ymin = min.ci, ymax = max.ci), width = 0.2) +
  coord_flip() +
  scale_x_discrete(labels = label_map) +
  labs(
    x = NULL,
    y = "Partial Rank Correlation\nCoefficient (PRCC)",
    title = ""
  ) +
  theme_minimal(base_size = 12)
