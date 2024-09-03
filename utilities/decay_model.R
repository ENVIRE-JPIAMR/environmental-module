## Script for Mancicni's decay model

get_decay_rate <- function(params) {
  k <-
    (0.8 + 0.006 * params$P_sal) * 1.07 ^ (params$T - 20) + params$IA * (1 - exp(-params$et * params$H)) / (params$et * params$H)
  
  return(k)
}