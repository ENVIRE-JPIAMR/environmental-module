## Script to call and execute the SWAT model 
## Outputs the concentration of ESBL E. coli in watershed (CFU/ml)

# Calculate manure runoff in kg/m²
manure_runoff <- data$m_apply[1] * data$runoff[1]

# Convert manure runoff to grams/m²
manure_runoff_grams <- manure_runoff * data$c_factor[1]

# Calculate E.coli load in runoff (CFU/m²)
E_coli_load <- manure_runoff_grams * data$C_barn

# Calculate E.coli reaching the waterbody after transport (CFU/m²)
E_coli_reaching_waterbody <- E_coli_load * data$transport[1]

# Calculate the volume of the waterbody (m³)
volume_waterbody <- data$w_area[1] * data$depth[1]

# Convert volume from m³ to liters (1 m³ = 1000000 mililiters)
volume_waterbody_mililiters <- volume_waterbody * 1000000

# Calculate the initial concentration of E.coli in the waterbody (CFU/mL)
data$C_init <- E_coli_reaching_waterbody / volume_waterbody_mililiters

