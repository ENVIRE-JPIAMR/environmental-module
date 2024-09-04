## Script for generating plots

#plot ECDFS
plot_ecdfs <- function(data) {
  # Create a new dataframe with log10 values and a 'pathway' column
  data_long <- data.frame(value = c(
    log10(data$dose_drinking),
    log10(data$dose_swimming),
    log10(data$dose_lettuce)
  ),
  pathway = rep(c("drinking", "swimming", "lettuce"), each = nrow(data)))
  
  # Plot the ECDFs
  ggplot(data_long, aes(x = value, color = pathway)) +
    stat_ecdf(geom = "step") +
    labs(
      title = "ECDF of different pathway exposure dose",
      subtitle = paste(
        "avg. ESBL E. coli conc. in watershed:",
        format(
          mean(data$C_init),
          scientific = TRUE,
          digits = 2
        ),
        "CFU/ml"
      ),
      x = "log10(CFU)",
      y = "ECDF",
      color = "pathway"
    ) +
    scale_color_manual(values = c(
      "drinking" = "blue",
      "swimming" = "green",
      "lettuce" = "red"
    )) +
    theme_minimal()
}