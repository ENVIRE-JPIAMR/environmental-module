## Script for loading libraries

## Set install_libraries to TRUE, or source manually the following
## block, in order to install all the libraries needed.
install_libraries = FALSE

## Install libraries
if (install_libraries)
{
  install.packages(c(
    "readxl",
    "ggplot2",
    "EnvStats",
    "mc2d"
    )
  )
}

####################
## LOAD LIBRARIES ##
####################

library(readxl)       # for reading excel files
library(ggplot2)      # for plotting
library(EnvStats)     # for truncated lognormal dist.
library(mc2d)         # for triangular dist.

## Cleanup
rm("install_libraries")



