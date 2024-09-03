## Submodule for tap-water drinking
source(here::here("environmental-module/utilities/decay_model.R"))

#computing decay rate
data$k_drinking <- get_decay_rate(data)

#concentration after decay
data$C_decay.drinking <- data$C_init * exp(-data$k_drinking * data$d_dwtp)

#concentration after DWTP processing (all treatments)
data$C_dwtp <-
  data$C_decay.drinking * (1 - data$CFS) * (1 - data$SSF) * (1 - data$RSF) * (1 - data$CF) * 
  (1 - data$C) * (1 / 10**data$UV) * (1 - data$Oz)

#human exposure
data$dose_drinking <- data$C_dwtp * data$w_consum
