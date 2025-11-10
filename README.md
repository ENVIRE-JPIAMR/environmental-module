# environmental-module

This repository provides the implementation of the QMRA model corresponding to the environemtnal model, 
intended to quantify ESBL producing E. coli in from broiler production to humans through the environment.

## Directory layout

This directory contains
- [`environmental_model_complete.R`](./environmental_model_complete.R): the full model script, including farm module, soil module, river module, lettuce exposure module, estimation of GI, UTI and DALY, and sensitivity analysis.
- [`inputs.csv`](./inputs.csv) : the input parameters for the farm module.
- [`sharma.xlsx`](./sharma.xlsx) : supplementary material from Sharma Manan et al (2019), used in the soil module to estimate the E. coli decay rate


