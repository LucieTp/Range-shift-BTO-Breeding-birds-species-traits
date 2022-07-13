# Range-shift-BTO-breeding-birds

This code implements the analysis of the relationship between species traits in and British breeding birds' range shift.

- R code 1 processes data from the BTO, extracts environmental variables, runs the SDMs for the 160 initial species and computes range shift.
- R code 2 computes/extratcs the species traits (biotic, abiotic and dispersal capabilities) and creates the main data file for modelling
- R code 3 runs the analysis (PGLMMs) and exports the results
- R code 4 generates most plots and other minor statistical analysis (wilcox tests etc).

-- utils.R and whois-functions.R are sourced in the scripts 1, 2, 3 and 4.
-- All datasets needed are located in the "data" folder except for the CRU climate which is downloadable from https://catalogue.ceda.ac.uk/uuid/58a8802721c94c66ae45c3baa4d814d0 (CRU v4.01)
