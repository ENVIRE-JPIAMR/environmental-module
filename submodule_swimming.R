## Submodule for recreational swimming

#concentration after (no) decay in worst case scenario
data$C_decay.swimming <- data$C_init

#concentration in bathing sites
data$C_bath <- (data$C_decay.swimming * data$F)  / data$V
  
#human exposure
if(data$swimmer_type[1] == 1){
  w_ingest <- data$w_ingest.adult
} else {
  w_ingest <- data$w_ingest.child
}

data$dose_swimming <- data$C_bath * w_ingest
