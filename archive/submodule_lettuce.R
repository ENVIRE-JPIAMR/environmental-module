## Submodule for lettuce harvesting

#concentration after decay
C_decay.lettuce <- matrix(1:Runs, nrow = Runs, ncol = data$n_harvest[1])

for(i in 1:data$n_harvest[1]){
  C_decay.lettuce[, i] <- data$C_init * exp(-data$k * i)
}

#CFU of ESBL E. coli adhered 
C_adhesion <- C_decay.lettuce * data$w_attach

#survival proportion on lettuce
P_decay <- data$f[1] * exp(-data$k1[1] * data$n_harvest[1]:1) + (1 - data$f[1]) * exp(-data$k2[1] * data$n_harvest[1]:1) #in inverse order

#concentration on per gram of lettuce after irrigation
data$C_lettuce <- rowSums(sweep(C_adhesion, MARGIN = 2, STATS = P_decay, FUN = "*"))

#TODO:concentration on per grarm lettuce after postharvest treatment
data$C_lettuce.post <- data$C_lettuce

#concentration on per grarm lettuce after home washing
data$C_lettuce.wash <- data$C_lettuce.post * (1-data$r_wash)

#human exposure
data$dose_lettuce <- data$C_lettuce.wash * data$l_consum

