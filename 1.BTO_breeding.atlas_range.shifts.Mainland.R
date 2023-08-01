
### this code is in support of Thompson et al. 2023 - 
### Joint effects of species traits and environmental preferences on range edge shifts of British birds
### published in Global Ecology and biogeography

### it reads in and formats the BTO and computes species' range shifts
### for the mainland study area

#### Authors: Dr Konstans Wells and Lucie Thompson
#### Computational Ecology Lab, Swansea University. UK.
####


################################################################################
library(raster)  
library(ggplot2) 
library(tidyverse)
library(rnaturalearth) # function ne_countries
library(sf)
library(sdm)

#sdm::installAll()
#install.packages("rnaturalearthdata")


setwd("E:/TheseSwansea/TraitStudy/Github")
dir.analysis = "E:/TheseSwansea/TraitStudy/Github"

################################################################################
## Load map

sf_UK  <- ne_states(country = 'United Kingdom', returnclass = "sf")
sf_UK <- subset(sf_UK, !(region %in% c('Northern Ireland','')))
# sf_UK  <- ne_countries(scale = "medium", country = 'United Kingdom', returnclass = "sf")
ggplot(data = sf_UK) + geom_sf()


##########################
## BTO Distribution data

# Breeding and wintering bird distributions in Britain and Ireland from citizen science bird atlases (Gillings et al., 2019) can be accessed at https://www.bto.org/our-science/data/what-data-are-available

BTO_distrib <- read.csv("data/distributions.csv", header=T) # period, sp code, season and grid for GB and Ireland

## BTO Species names data
BTO_sp <- read.csv("data/species_lookup.csv", header=T)

## BTO grid data
BTO_grid <- read.csv("data/grid_square_coordinates_lookup.csv", header=T)

# Subset distribution data to Britain
BTO_distrib <- as_tibble(BTO_distrib) %>% filter(island=="B" & season=="B" & resolution==10)
BTO_distrib$periodN <- paste0("P.", as.numeric(droplevels(as.factor(BTO_distrib$period))))
BTO_distrib$Spec <- paste0("Sp", BTO_distrib$speccode)
BTO_distrib$Pres <- 1



#####################
# Grid and BTO Locations

# Coordinates
Coord10 <- as_tibble(BTO_grid) %>%
  dplyr::filter(resolution==10 & order<5) %>%
  dplyr::group_by(grid) %>%
  dplyr::summarise(long = mean(long), lat = mean(lat))

ngrid <- nrow(Coord10)

# Location data (selecting unique coordinates for centre of each grid cell ?)
Loc10 <- BTO_distrib %>% dplyr::select(grid) %>% distinct(grid)

# Join coordinates to location data
Loc10 <- Loc10 %>%  left_join(Coord10)
# Filter missing coordinates
Loc10 <- Loc10 %>%  filter(!is.na(long) & !is.na(lat))
nloc <- nrow(Loc10)
Loc10_whole = Loc10
# 3867 locations occupied by birds

Loc10 = subset(Loc10, lat < 58.7)

library(plyr)
t.Out = ddply(BTO_distrib, c('speccode','periodN'), summarise, nb.Out = length(grid))

# Join coordinates to occurrence data
BTO_distrib <- BTO_distrib %>%  left_join(Loc10, by='grid')
BTO_distrib <- BTO_distrib %>% filter(!is.na(long) & !is.na(lat))


t.Mainland = ddply(BTO_distrib, c('speccode','periodN'), summarise, nb.Mainland = length(grid))
t = merge(t.Out, t.Mainland, by = c('speccode','periodN'))
plot(t$nb.Out~t$nb.Mainland)

ggplot(data = sf_UK) + geom_sf() +
  geom_point(data = Loc10, aes(x = long, y = lat), size = 1, shape = 23, fill = "darkred")

ggplot(data = sf_UK) + geom_sf() +
  geom_point(data = unique(BTO_distrib[,c('lat','long')]), aes(x = long, y = lat), size = 1, shape = 23, fill = "darkred")

ggplot(data = sf_UK) + geom_sf() +
  geom_point(data = unique(BTO_distrib[which(BTO_distrib$speccode == 71 & BTO_distrib$periodN == "P.1"),c('lat','long',"periodN")]), aes(x = long, y = lat, colour = periodN), size = 1, shape = 23, fill = "darkred")


#####################
# BTO Species


## 3 different periods
period <- unique(BTO_distrib$period)
period

# Nb of grid cells species have been observed
Ngrid <- BTO_distrib %>%
 dplyr::group_by(speccode, periodN) %>%
 dplyr::summarise(ngrid = n()) # group size
Ngrid <- pivot_wider(Ngrid, id_cols = speccode, names_from = periodN, values_from =ngrid)

Spec <- as_tibble(BTO_sp)
Spec  <- Spec %>%  left_join(Ngrid, by = c("speccode" = "speccode"))
Spec  <- Spec %>% filter(P.1>=50 & P.3>=50) 


# Compute distributional core of species (north or south of the mean position of the all 100 km2 grid cells)
Distr.Core <- BTO_distrib %>% dplyr::filter(periodN=="P.1") %>% group_by(speccode) %>% dplyr::summarise(lat_mean.P.1 = mean(lat)) 
Distr.Core$distrib.core = ifelse(Distr.Core$lat_mean.P.1 > mean(Loc10$lat), "north", "south") 
# north south 
# 91 > 93   132 > 129
# slightly more northern species compared to without the Orkney and Shetlands
# because the average latitude changes

Spec <- left_join(Spec, Distr.Core)


### we loose two species when filtering out the Orkney and Shetland Islands. 


#####################
# Species trait data

# Species traits (diet composition, body mass - Elton traits) (Wilman et al., 2014) can be accessed at https://figshare.com/articles/Data_Paper_Data_Paper/3559887
temp <- tempfile()
download.file("https://ndownloader.figshare.com/files/5631081", temp)

ET <- read.table(temp, header = TRUE, fill  = TRUE, quote = "\"", stringsAsFactors = FALSE,sep = "\t")
unlink(temp)

names(ET)[which(names(ET) == 'Scientific')] = 'scientific_name_ET'

ET$scientific_name = Spec$scientific_name[match(ET$scientific_name_ET,Spec$scientific_name)]
names = setdiff(Spec$scientific_name,ET$scientific_name_ET)

## 1 - homogenizing taxonomy

## Tring different methods :
################################################################################
## library for taxonomic queries/synonyms
library(taxize)

syn = synonyms(names,db = "itis")
syn_df = synonyms_df(syn)
syn_df = subset(syn_df, syn_name%in%ET$scientific_name_ET)[,c(".id","syn_name")]
names(syn_df) = c('scientific_name','scientific_name_ET')
ET$scientific_name[match(syn_df$scientific_name_ET,ET$scientific_name_ET)] = syn_df$scientific_name

# then by hand for 11 species 
names = Spec$scientific_name[which(!Spec$scientific_name%in%ET$scientific_name)]

ET$scientific_name[which(ET$scientific_name_ET == 'Anas strepera')] = 'Mareca strepera' # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=9227AE600B1BC4CA&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Tetrao tetrix")] = 'Lyrurus tetrix' # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=D4441CD6E9C993EF&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Parus palustris")] = "Poecile palustris" # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=523763E55D4A8153&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Parus montanus")] = "Poecile montanus" # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=01D4B731BF3D081F&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Parus ater")] = "Periparus ater" #https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=A4EBA919FCAFED5E&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Parus caeruleus")] = "Cyanistes caeruleus" # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=9BE53D340F9A4305&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Carduelis flammea")] = "Acanthis cabaret" # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=9BE53D340F9A4305&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Saxicola torquatus")] = "Saxicola rubicola" # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=0EA8F8B905405FB3&sec=synonyms
ET$scientific_name[which(ET$scientific_name_ET == "Sterna albifrons")] = "Sternula albifrons"  # https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=17310BF0FE9BAC8A&sec=synonyms



names(ET)[which(names(ET) == 'English')] = 'english_name_ET'
Spec = Spec %>%  left_join(ET[,-match(c("SpecID","BLFamilyLatin","BLFamilyEnglish","BLFamSequID"),colnames(ET))], by = c("scientific_name"))

summary(Spec)
Spec = Spec[-which(Spec$scientific_name == 'Corvus cornix'),] # missing values
write.csv(Spec, file="data/SpecTrait_062023_157sp.csv") 

# SpecTrait_122021.csv file with 160 species 
# SpecTrait_062023_151sp.csv analysis restricted to mainland UK
# SpecTrait_062023_157sp < 58.7

nSpec <- nrow(Spec)
spec_speccode <- Spec$speccode





#####################
# 2 - Compute range shift metrics and distance to geographical barriers

## difference in mean location of the 20 most northern or most southern records
library(geosphere)
df_rangeshift <- data.frame(speccode=spec_speccode)

for(i in 1:nSpec){
  lat.all_P.1 <- sort(data.frame(BTO_distrib %>% filter(speccode == spec_speccode[i], periodN=="P.1") %>% dplyr::select(lat))$lat)
  lat.all_P.3 <- sort(data.frame(BTO_distrib %>% filter(speccode == spec_speccode[i], periodN=="P.3") %>% dplyr::select(lat))$lat)

  lat.max20_P.1 <- tail(sort(lat.all_P.1), 20)
  lat.max20_P.3 <- tail(sort(lat.all_P.3), 20)
  
  lat.min20_P.1 <- head(sort(lat.all_P.1), 20)
  lat.min20_P.3 <- head(sort(lat.all_P.3), 20)
  
  # distance from northern boundary at P.1 (distance the species has to expand towards the north)
  df_rangeshift$dist_N_km_max20_P.1[i] <- distGeo(c(3,max(BTO_distrib$lat)), c(3, mean(lat.max20_P.1)))/1000*ifelse(max(Loc10$lat)<mean(lat.max20_P.1), -1, 1)
  
  # distance from southern boundary at P.1 (distance the species has to expand towards the south)
  df_rangeshift$dist_S_km_min20_P.1[i] <- distGeo(c(3,min(BTO_distrib$lat)), c(3, mean(lat.min20_P.1)))/1000*ifelse(min(Loc10$lat)>mean(lat.min20_P.1), -1, 1)
  
  # difference in latitude
  df_rangeshift$shift_max20_P.1.3[i] <- mean(lat.max20_P.3) - mean(lat.max20_P.1)
  df_rangeshift$shift_min20_P.1.3[i] <- mean(lat.min20_P.3) - mean(lat.min20_P.1)
  df_rangeshift$shift_all_P.1.3[i] <- mean(lat.all_P.3) - mean(lat.all_P.1)
  
  # distance
  df_rangeshift$shift_max20_P.1.3_dist_km[i] <- distGeo(c(3,mean(lat.max20_P.3)), c(3, mean(lat.max20_P.1)))/1000*ifelse(df_rangeshift$shift_max20_P.1.3[i]<0, -1, 1)
  df_rangeshift$shift_min20_P.1.3_dist_km[i] <- distGeo(c(3,mean(lat.min20_P.3)), c(3, mean(lat.min20_P.1)))/1000*ifelse(df_rangeshift$shift_min20_P.1.3[i]<0, -1, 1)
  df_rangeshift$shift_all_P.1.3_dist_km[i] <- distGeo(c(3,mean(lat.all_P.3)), c(3, mean(lat.all_P.1)))/1000*ifelse(df_rangeshift$shift_all_P.1.3[i]<0, -1, 1)
 
  df_rangeshift$nrecord_P.1[i] <- Spec$P.1[which(Spec$speccode == spec_speccode[i])]  # length(lat.all_P.1)  
  df_rangeshift$nrecord_P.3[i] <- Spec$P.3[which(Spec$speccode == spec_speccode[i])]  # length(lat.all_P.3)
  
}

## Save

write.csv(x = df_rangeshift, file = "data/df_rangeshift_062023.Mainland1.csv")
# 01.2023 : added velocity of climate change
# 05/2023 : no velocity but removed Ireland
# 06/2023 : Mainland

