


[![DOI](https://zenodo.org/badge/507051216.svg)](https://zenodo.org/badge/latestdoi/507051216)



British bird's range shifts
===========================
This repository contains the code for assessing the role of species' traits in British breeding bird recent range shifts.
This code is in support of Thompson et al. 2023 - Joint effects of species traits and environmental preferences on range edge shifts of British birds - (Global Ecology and Biogeography).


To reproduce the analysis, download the data (see data availability statement in main manuscript and below) and place into a '~/data' folder within your desired directory.
Change working directory to chosen location on your machine and run code as follows:

- 1.BTO_breeding.atlas.range.shifts.R processes the bird occurrence data, environmental variables and run species distribution models and computes range shift.
- 1.BTO_breeding.atlas.range.shifts.mainland.R processes the bird occurrence data and computes range shift for the mainland only subset.
- 2.CreationMainTable_SpeciesTraitsData.R and 2.CreationMainTable_SpeciesTraitsData.Mainland.R compute and/or extracts the species traits (biotic, abiotic and dispersal capabilities) and generates the main table used in the subsequent analysis
- 3.Analysis.R and 3.Analysis.Mainland.R run the Phylogenetic Linear Mixed Models using the tables generated in 2. and exports the outputs into csv and html format
- 4.Plots generates the figures and other minor statistical analysis (wilcox tests etc) for both the mainland and whole study.


- utils.R and whois-functions.R are both annex scripts and contain functions to calculate trophic traits and translate the metaweb species' IDs to scientific name.
They are sourced in the scripts 1, 2, 3 and 4.
- All datasets are freely available (see data availability statement in main manuscript and below for URLs)

____________________________________________________________________________________________________________________________________________________________________________

## Data availability
- Breeding and wintering bird distributions in Britain and Ireland from citizen science bird atlases (Gillings et al., 2019) can be accessed at https://www.bto.org/our-science/data/what-data-are-available
- Climatic data (Climatic Research Unit (CRU) time series v.4.01 from 1901 to 2016) (Harris et al., 2020) can be accessed at https://catalogue.ceda.ac.uk/uuid/58a8802721c94c66ae45c3baa4d814d0
- Land cover changes in Europe (HILDA version 2.0 from 1900 to 2010) (Fuchs et al., 2015) can be accessed at https://doi.pangaea.de/10.1594/PANGAEA.921846?format=html#download
- Phylogenetic tree (BigBird phylogeny 6714 taxa tree) (Burleigh et al., 2015) can be accessed at https://datadryad.org/stash/dataset/doi:10.5061/dryad.r6b87
- Dispersal capabilities through the global hand-wing index (Sheard et al., 2020) can be accessed at https://zenodo.org/record/3747657
- Species traits (diet composition, body mass - Elton traits) (Wilman et al., 2014) can be accessed at https://figshare.com/articles/Data_Paper_Data_Paper/3559887
- European food web (TETRA EU 1.0) (Maiorano et al., 2020) https://datadryad.org/stash/dataset/doi:10.5061/dryad.jm63xsj7b

