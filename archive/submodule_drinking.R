## Submodule for tap-water drinking

#concentration after decay
data$C_decay.drinking <- data$C_init * exp(-data$k * data$d_dwtp)

#concentration after DWTP processing (all treatments)
data$C_dwtp <-
  data$C_decay.drinking * (1 - data$CFS) * (1 - data$SSF) * (1 - data$RSF) * (1 - data$CF) * 
  (1 - data$C) * (1 / 10**data$UV) * (1 - data$Oz)

#human exposure
data$dose_drinking <- data$C_dwtp * data$w_consum
