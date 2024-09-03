## Script to call and execute the SWAT model 
## Outputs the concentration of ESBL E. coli in watershed (CFU/ml)

#dummy SWAT model
data$C_init <- data$C_barn # assuming no decay given 1ml water = 1g 

