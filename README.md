
British bird's range shifts
===========================
This repository contains the code for assessing which species' traits best explain British Breeding Birds' recent ranges shifts.
This code is in support of Thompson et al. 2023 (under revision in Global Ecology and Biogeography, 2023).


To reproduce the analysis, download the data (see data availability statement in main manuscript and below) and place into a '~/data' folder within your desired directory.
Change working directory to this chosen location on your machine and run code as follows:

- 1.BTO_breeding.atlas.range.shifts.R processes the bird occurrence data, environmental variables, runs species distribution models and computes range shift.
- 2.CreationMainTable_SpeciesTraitsData.R computes/extracts the species traits (biotic, abiotic and dispersal capabilities) and generates the main table used in the subsequent analysis
- 3.Analysis runs the PGLMMs using the tables generated in 2. and exports the outputs into csv format
- 4.Plots generates the figures and other minor statistical analysis (wilcox tests etc).

-- utils.R and whois-functions.R are both annex scripts and contain functions to calculate trophic traits and translate TETRA EU (metaweb) species' IDs to scientific name.
They are sourced in the scripts 1, 2, 3 and 4.
-- All datasets are freely available (see data availability statement in main manuscript for URLs) 

CRU climate data is available at https://catalogue.ceda.ac.uk/uuid/58a8802721c94c66ae45c3baa4d814d0 (CRU v4.01)
