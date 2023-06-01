### Code for Thompson et al 2023 

# formatting the BTO data


################################################################################
library(raster)  
library(ggplot2) 
library(tidyverse)
library(rnaturalearth) # function ne_countries
library(sf)
library(sdm)

#sdm::installAll()
#install.packages("rnaturalearthdata")

################################################################################
## Load map

sf_UK  <- ne_states(country = 'United Kingdom', returnclass = "sf")
sf_UK <- subset(sf_UK, !(name %in% c('Orkney','Shetland Islands','Eilean Siar','Isles of Scilly') | region %in% c('Northern Ireland','')))
# sf_UK  <- ne_countries(scale = "medium", country = 'United Kingdom', returnclass = "sf")
ggplot(data = sf_UK) + geom_sf()

## BTO Distribution data
setwd("E:/TheseSwansea/TraitStudy/Github")
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
# 3867 locations occupied by birds

Loc10_sf = st_as_sf(Loc10, coords = c('long', 'lat'), crs = crs(sf_UK))
Loc10_locationNames = st_intersection(Loc10_sf, sf_UK)
Loc10 = subset(Loc10, grid %in% Loc10_locationNames$grid)

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

# EltonTraits v1.0 (Wilman et al. (2014), https://figshare.com/articles/Data_Paper_Data_Paper/3559887)
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

names(ET)[which(names(ET) == 'English')] = 'english_name_ET'
Spec = Spec %>%  left_join(ET[,-match(c("SpecID","BLFamilyLatin","BLFamilyEnglish","BLFamSequID"),colnames(ET))], by = c("scientific_name"))

summary(Spec)
Spec = Spec[-which(Spec$scientific_name == 'Corvus cornix'),] # missing values
write.csv(Spec, file="data/SpecTrait_062023_151sp.csv") 

# SpecTrait_122021.csv file with 160 species 
# SpecTrait_062023_151sp.csv analysis restricted to mainland UK

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

write.csv(x = df_rangeshift, file = "data/df_rangeshift_062023.Mainland.csv")
# 01.2023 : added velocity of climate change
# 05/2023 : no velocity but removed Ireland
# 06/2023 : Mainland


#############
## 3 - Calculating the percentage of coastal occurrences to distinguish marine 
## from non - marine species

Sptraits_Loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode %in% unique(Spec$speccode)),])

# gives the distance between the sp coordinate to the nearest shoreline (sf_UK) 
# if need to do it again: think about only using grid cells from P.1 ?

Sptraits_Loc_sf <- Sptraits_Loc %>% st_as_sf(coords = c('long','lat')) %>%
  st_set_crs(crs(sf_UK))

UK_coast = ne_countries(scale = 'medium', country = 'United Kingdom', returnclass = 'sf')
plot(UK_coast)

#transform from polygon shape to line
UK_line <- st_cast(UK_coast, "MULTILINESTRING")

#calculation of the distance between the coast and our points
d <- st_distance(UK_line, Sptraits_Loc_sf)
df <- data.frame(dist_km = as.vector(d)/1000,
                 st_coordinates(Sptraits_Loc_sf), 
                 Sptraits_Loc_sf[,1:10])

ggplot() + geom_point(data = df[sample(nrow(df), 0.1*nrow(df)),], aes(x = X, y = Y, colour = ifelse(dist_km<10,1,0)), size = .5) 

write.csv(df, "data/DistanceToCoast_062023.csv")

# proportion of points that are within 20km from the shore
prop_marine = df %>%
  dplyr::group_by(speccode = speccode) %>%
  dplyr::summarise(nb_grid = length(speccode), nb_grid_10km = sum(dist_km<10),
            prop_Marine10km = nb_grid_10km/nb_grid*100, 
            nb_grid_20km = sum(dist_km<20),
            prop_Marine20km = nb_grid_20km/nb_grid*100) %>% ungroup()

write.csv(prop_marine,"ProportionMarine_062023.csv")

## check by plotting species that we consider "Marine"
par(mar = c(2, 2, 2, 2))

for (i in unique(prop_marine$speccode[which(prop_marine$prop_Marine20km>75)])){
  print(BTO_sp[which(BTO_sp$speccode == i),] )
  g = ggplot(data = sf_UK) + geom_sf() +
    geom_point(data = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode ==i),]), aes(x = long, y = lat, colour = periodN), size = 1) +
    ggtitle(BTO_sp[which(BTO_sp$speccode == i),]$english_name)
  print(g)
  print(Spec[which(Spec$english_name == BTO_sp[which(BTO_sp$speccode == i),]$english_name),])
}



################################################################################
## 4 - adding environmental covariate (climate and land cover data)

######
# CRU climate data
setwd("E:/TheseSwansea/TraitStudy/Github/Range-shift-BTO-breeding-birds/data/CRU_Climate/")
require(ncdf4)

# Load the CRU TS dataset into R as rasterBrick
tmp <- raster::brick("cru_ts4.04.1901.2019.tmp.dat.nc", varname="tmp")
pre <- raster::brick("cru_ts4.04.1901.2019.pre.dat.nc", varname="pre")

summary(BTO_distrib$period)

# Year for each period to extract environmental data
years_P.1 <- (1968-5):1972
years_P.2 <- (1988-5):1991
years_P.3 <- (2008-5):2011

nP <- 3

# selected CRU climate variables
cru_seas.1_P.1 <- c(paste0("X", years_P.1, ".01.16"), paste0("X", years_P.1, ".02.15"))
cru_seas.1_P.2 <- c(paste0("X", years_P.2, ".01.16"), paste0("X", years_P.2, ".02.15"))
cru_seas.1_P.3 <- c(paste0("X", years_P.3, ".01.16"), paste0("X", years_P.3, ".02.15"))

cru_seas.2_P.1 <- c(paste0("X", years_P.1, ".03.16"), paste0("X", years_P.1, ".04.16"), paste0("X", years_P.1, ".05.16")) 
cru_seas.2_P.2 <- c(paste0("X", years_P.2, ".03.16"), paste0("X", years_P.2, ".04.16"), paste0("X", years_P.2, ".05.16")) 
cru_seas.2_P.3 <- c(paste0("X", years_P.3, ".03.16"), paste0("X", years_P.3, ".04.16"), paste0("X", years_P.3, ".05.16")) 

# temperature for each period and the extent of UK
tmp_seas.1_P.1 <- crop(calc(tmp[[cru_seas.1_P.1]], fun = mean), sf_UK)
tmp_seas.1_P.2 <-  crop(calc(tmp[[cru_seas.1_P.2]], fun = mean), sf_UK)
tmp_seas.1_P.3 <- crop(calc(tmp[[cru_seas.1_P.3]], fun = mean), sf_UK)

tmp_seas.2_P.1 <- crop(calc(tmp[[cru_seas.2_P.1]], fun = mean), sf_UK)
tmp_seas.2_P.2 <-  crop(calc(tmp[[cru_seas.2_P.2]], fun = mean), sf_UK)
tmp_seas.2_P.3 <- crop(calc(tmp[[cru_seas.2_P.3]], fun = mean), sf_UK)

# precipitations
pre_seas.1_P.1 <- crop(calc(pre[[cru_seas.1_P.1]], fun = mean), sf_UK)
pre_seas.1_P.2 <-  crop(calc(pre[[cru_seas.1_P.2]], fun = mean), sf_UK)
pre_seas.1_P.3 <- crop(calc(pre[[cru_seas.1_P.3]], fun = mean), sf_UK)

pre_seas.2_P.1 <- crop(calc(pre[[cru_seas.2_P.1]], fun = mean), sf_UK)
pre_seas.2_P.2 <-  crop(calc(pre[[cru_seas.2_P.2]], fun = mean), sf_UK)
pre_seas.2_P.3 <- crop(calc(pre[[cru_seas.2_P.3]], fun = mean), sf_UK)

# Check similarity between layers
stack_cru <- stack(tmp_seas.1_P.1, tmp_seas.1_P.2, tmp_seas.1_P.3, tmp_seas.2_P.1, tmp_seas.2_P.2, tmp_seas.2_P.3,  pre_seas.1_P.1, pre_seas.1_P.2, pre_seas.1_P.3, pre_seas.2_P.1, pre_seas.2_P.2, pre_seas.2_P.3)
layerStats(stack_cru,'pearson',na.rm=TRUE)

### plotting temperature and precipitation rasters
par(mfrow=c(1,3))
plot(tmp_seas.1_P.1, main = paste0("Winter temperature  ",as.character(period[1])))
plot(tmp_seas.1_P.3, main = paste0("Winter temperature  ",as.character(period[2])))
tmp_seas.1_diff_P1.3 = tmp_seas.1_P.3 - tmp_seas.1_P.1; plot(tmp_seas.1_diff_P1.3, 
                                                             main = paste0("Difference in winter temperature","\n", "between 1968 and 2011"))

plot(pre_seas.1_P.1, main = paste0("Winter precipitation ",as.character(period[1])))
plot(pre_seas.1_P.3, main = paste0("Winter precipitation ",as.character(period[2])))
pre_seas.1_diff_P1.3 = pre_seas.1_P.3 - pre_seas.1_P.1; plot(pre_seas.1_diff_P1.3, 
                                                             main = paste0("Difference in winter precipitation","\n", "between 1968 and 2011"))

plot(pre_seas.2_P.1, main = paste0("Summer precipitation ",as.character(period[1])))
plot(pre_seas.2_P.3, main = paste0("Summer precipitation ",as.character(period[2])))
pre_seas.2_diff_P1.3 = pre_seas.2_P.3 - pre_seas.2_P.1; plot(pre_seas.2_diff_P1.3, 
                                                             main = paste0("Difference in summer precipitation","\n", "between 1968 and 2011"))

plot(tmp_seas.2_P.1, main = paste0("Summer temperature ",as.character(period[1])))
plot(tmp_seas.2_P.3, main = paste0("Summer temperature ",as.character(period[2])))
tmp_seas.2_diff_P1.3 = tmp_seas.2_P.3 - tmp_seas.2_P.1; plot(tmp_seas.2_diff_P1.3, 
                                                             main = paste0("Difference in summer temperature","\n", "between 1968 and 2011"))


plot(values(tmp_seas.1_P.1), values(tmp_seas.1_P.3))

# Convert to dataframe for ggplot
df_tmp_seas.1_P.1 <- as.data.frame(tmp_seas.1_P.1, xy = TRUE)
ggplot() + geom_raster(data = df_tmp_seas.1_P.1, aes(x = x, y = y, fill=layer))


######
# HILDA historical land use (http://www.geo-informatie.nl/fuchs003/#)

setwd("E:/TheseSwansea/TraitStudy/Github/Range-shift-BTO-breeding-birds/HILDA_v2._LandUseChangesEurope")

#lu_P.1_st <- read_stars("./Gross_Final_1km_EU27CH_TIFF/eu27ch1960.tif")
#lu_P.1_st <- st_transform(lu_P.1_st, crs= 4326)
#lu_P.1_st <- st_crop(lu_P.1_st, sf_UK)


lu_P.1 <- raster::raster("./Gross_Final_1km_EU27CH_TIFF/eu27ch1960.tif")
lu_P.2 <- raster::raster("./Gross_Final_1km_EU27CH_TIFF/eu27ch1980.tif")
lu_P.3 <- raster::raster("./Gross_Final_1km_EU27CH_TIFF/eu27ch2000.tif")

lu_P.1 <- projectRaster(lu_P.1, crs = crs(sf_UK))
lu_P.2 <- projectRaster(lu_P.2, crs = crs(sf_UK))
lu_P.3 <- projectRaster(lu_P.3, crs = crs(sf_UK))

lu_P.1 <- crop(lu_P.1, sf_UK)
lu_P.2 <- crop(lu_P.2, sf_UK)
lu_P.3 <- crop(lu_P.3, sf_UK)

values(lu_P.1) <- trunc(values(lu_P.1)/100)
values(lu_P.2) <- trunc(values(lu_P.2)/100)
values(lu_P.3) <- trunc(values(lu_P.3)/100)

# Check similarity between layers
stack_lu <- stack(lu_P.1, lu_P.2, lu_P.3)
layerStats(stack_lu,'pearson',na.rm=TRUE)

lu_diff_P1.3 = lu_P.3 - lu_P.1

par(mfrow=c(1,3))
plot(lu_P.1,main = paste0("Land cover ",as.character(period[1]))); 
plot(lu_P.3, main = paste0("Land cover ",as.character(period[2])))
plot(lu_diff_P1.3, main = "Land cover change")

##   111 - Settlement
##   222 - Cropland
##   333 - Forest
##   444 - Grassland
##   555 - Other Land
##   666 - Water

######
# Extract environmental data for each grid cell

plot(tmp_seas.1_P.1)
points(data = Loc10[ ,c("long", "lat")], lat~long)

Loc10$tmp_seas.1_P.1 <- unlist(raster::extract(tmp_seas.1_P.1, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$tmp_seas.1_P.2 <- unlist(raster::extract(tmp_seas.1_P.2, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$tmp_seas.1_P.3 <- unlist(raster::extract(tmp_seas.1_P.3, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$tmp_seas.2_P.1 <- unlist(raster::extract(tmp_seas.2_P.1, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$tmp_seas.2_P.2 <- unlist(raster::extract(tmp_seas.2_P.2, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$tmp_seas.2_P.3 <- unlist(raster::extract(tmp_seas.2_P.3, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$pre_seas.1_P.1 <- unlist(raster::extract(pre_seas.1_P.1, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$pre_seas.1_P.2 <- unlist(raster::extract(pre_seas.1_P.2, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$pre_seas.1_P.3 <- unlist(raster::extract(pre_seas.1_P.3, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$pre_seas.2_P.1 <- unlist(raster::extract(pre_seas.2_P.1, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$pre_seas.2_P.2 <- unlist(raster::extract(pre_seas.2_P.2, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))
Loc10$pre_seas.2_P.3 <- unlist(raster::extract(pre_seas.2_P.3, data.frame(Loc10[ ,c("long", "lat")]), buffer = 10000))

for(p in 1:nloc){
	vals_P.1 <- unlist(raster::extract(lu_P.1, data.frame(Loc10[p,c("long", "lat")]), buffer = 10000)) 
	Loc10$pForest_P.1[p] <- length(which(vals_P.1==3)) / length(vals_P.1)
	Loc10$pGrass_P.1[p] <- length(which(vals_P.1==4)) / length(vals_P.1)
	Loc10$pCrop_P.1[p] <- length(which(vals_P.1==2)) / length(vals_P.1)
	Loc10$pSettlem_P.1[p] <- length(which(vals_P.1==1)) / length(vals_P.1)
	vals_P.2 <- unlist(raster::extract(lu_P.2, data.frame(Loc10[p,c("long", "lat")]), buffer = 10000)) 
	Loc10$pForest_P.2[p] <- length(which(vals_P.2==3)) / length(vals_P.2)
	Loc10$pGrass_P.2[p] <- length(which(vals_P.2==4)) / length(vals_P.2)
	Loc10$pCrop_P.2[p] <- length(which(vals_P.2==2)) / length(vals_P.2)
	Loc10$pSettlem_P.2[p] <- length(which(vals_P.2==1)) / length(vals_P.2)
	vals_P.3 <- unlist(raster::extract(lu_P.3, data.frame(Loc10[p,c("long", "lat")]), buffer = 10000)) 
	Loc10$pForest_P.3[p] <- length(which(vals_P.3==3)) / length(vals_P.3)
	Loc10$pGrass_P.3[p] <- length(which(vals_P.3==4)) / length(vals_P.3)
	Loc10$pCrop_P.3[p] <- length(which(vals_P.3==2)) / length(vals_P.3)
	Loc10$pSettlem_P.3[p] <- length(which(vals_P.3==1)) / length(vals_P.3)
}


#save(Loc10, file="Loc10.RData")

#load("Loc10.RData")


Loc10_long1 = pivot_longer(Loc10,4:6, names_to = 'tmp_seas1')
Loc10_long2 = pivot_longer(Loc10,7:9, names_to = 'tmp_seas2')
Loc10_long3 = pivot_longer(Loc10,10:12, names_to = 'pre_seas1')
Loc10_long4 = pivot_longer(Loc10,13:15, names_to = 'pre_seas2')

g1 = ggplot(data = Loc10_long1, aes(y = value, x = tmp_seas1)) + geom_boxplot() + ggtitle('Winter - temperature seasonality')
g2 = ggplot(data = Loc10_long2, aes(y = value, x = tmp_seas2)) + geom_boxplot() + ggtitle('Summer - temperature seasonality')
g3 = ggplot(data = Loc10_long4, aes(y = value, x = pre_seas2)) + geom_boxplot() + ggtitle('Summer - precipitation seasonality')
g4 = ggplot(data = Loc10_long3, aes(y = value, x = pre_seas1)) + geom_boxplot() + ggtitle('Winter - precipitation seasonality')

library(ggpubr)
ggarrange(g1,g2,g3,g4)

ggsave('ClimaticVariationAcrossPeriods.jpeg')

##############################################################
## land use and climate change across the three time periods

#################
## FIGURE S3 and S4

par(mfrow = c(1,4))

summary(Loc10)
n = length(Loc10$tmp_seas.1_P.1[!is.na(Loc10$tmp_seas.1_P.1)])

mu <- c(mean(Loc10$tmp_seas.1_P.1, na.rm = T), mean(Loc10$tmp_seas.1_P.2, na.rm = T), mean(Loc10$tmp_seas.1_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$tmp_seas.1_P.1, na.rm = T), sd(Loc10$tmp_seas.1_P.2, na.rm = T),sd(Loc10$tmp_seas.1_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in winter temperature seasonality", names.arg = c("P1", 'P2', "P3"), ylim = c(0, 9))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

mu <- c(mean(Loc10$tmp_seas.2_P.1, na.rm = T), mean(Loc10$tmp_seas.2_P.2, na.rm = T), mean(Loc10$tmp_seas.2_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$tmp_seas.2_P.1, na.rm = T), sd(Loc10$tmp_seas.2_P.2, na.rm = T),sd(Loc10$tmp_seas.2_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in summer temperature seasonality", names.arg = c("P1", 'P2', "P3"), ylim = c(0, 9))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

mu <- c(mean(Loc10$pre_seas.1_P.1, na.rm = T), mean(Loc10$pre_seas.1_P.2, na.rm = T), mean(Loc10$pre_seas.1_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$pre_seas.1_P.1, na.rm = T), sd(Loc10$pre_seas.1_P.2, na.rm = T),sd(Loc10$pre_seas.1_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in winter precipitation seasonality", names.arg = c("P1", 'P2', "P3"), ylim = c(0, 125))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

mu <- c(mean(Loc10$pre_seas.2_P.1, na.rm = T), mean(Loc10$pre_seas.2_P.2, na.rm = T), mean(Loc10$pre_seas.2_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$pre_seas.2_P.1, na.rm = T), sd(Loc10$pre_seas.2_P.2, na.rm = T),sd(Loc10$pre_seas.2_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in summer precipitation seasonality", names.arg = c("P1", 'P2', "P3"), ylim = c(0, 125))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

ggsave('LandUseVariationAcrossPeriods.newplots.jpeg')


#########################
## LAND USE

n = nrow(Loc10)

mu <- c(mean(Loc10$pForest_P.1, na.rm = T), mean(Loc10$pForest_P.2, na.rm = T), mean(Loc10$pForest_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$pForest_P.1, na.rm = T), sd(Loc10$pForest_P.2, na.rm = T),sd(Loc10$pForest_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in forest cover", names.arg = c("P1", 'P2', "P3"), ylim = c(0, .45))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

mu <- c(mean(Loc10$pGrass_P.1, na.rm = T), mean(Loc10$pGrass_P.2, na.rm = T), mean(Loc10$pGrass_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$pGrass_P.1, na.rm = T), sd(Loc10$pGrass_P.2, na.rm = T),sd(Loc10$pGrass_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in grassland cover", names.arg = c("P1", 'P2', "P3"), ylim = c(0,.45))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

mu <- c(mean(Loc10$pCrop_P.1, na.rm = T), mean(Loc10$pCrop_P.2, na.rm = T), mean(Loc10$pCrop_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$pCrop_P.1, na.rm = T), sd(Loc10$pCrop_P.2, na.rm = T),sd(Loc10$pCrop_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in cropland cover", names.arg = c("P1", 'P2', "P3"), ylim = c(0, .45))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

mu <- c(mean(Loc10$pSettlem_P.1, na.rm = T), mean(Loc10$pSettlem_P.2, na.rm = T), mean(Loc10$pSettlem_P.3, na.rm = T))
CIs <-1.96*c(sd(Loc10$pSettlem_P.1, na.rm = T), sd(Loc10$pSettlem_P.2, na.rm = T),sd(Loc10$pSettlem_P.3, na.rm = T))/sqrt(n)
h <- barplot(mu, xlab="Period", ylab="Change in urban area cover", names.arg = c("P1", 'P2', "P3"), ylim = c(0, .45))
arrows(h, mu-CIs, h, mu+CIs, code=3, length=0.05, angle=90)
box()

ggsave('LandUseVariationAcrossPeriods.newplots.jpeg')


################################################################################
## Putting everything together into on table

# Environmental predictor names
env.predictor_names <- c("tmp_seas.1", "tmp_seas.2", "pre_seas.1", "pre_seas.2", "pForest", "pGrass", "pCrop", "pSettlem")
nEnvPred <- length(env.predictor_names)

# Predictor data for different time periods

Pred10_P.1 <- Loc10 %>% select(grid, long, lat, grep("P.1", names(Loc10)))
Pred10_P.2 <- Loc10 %>% select(grid, long, lat, grep("P.2", names(Loc10)))
Pred10_P.3 <- Loc10 %>% select(grid, long, lat, grep("P.3", names(Loc10)))

# Bird occurrence data
spec_occ_P.1 <- BTO_distrib %>% filter( periodN=="P.1" & !is.na(match(speccode, spec_speccode))) %>% select(Spec, grid, Pres)
spec_occ_P.2 <- BTO_distrib %>% filter( periodN=="P.2" & !is.na(match(speccode, spec_speccode))) %>% select(Spec, grid, Pres)
spec_occ_P.3 <- BTO_distrib %>% filter( periodN=="P.3" & !is.na(match(speccode, spec_speccode))) %>% select(Spec, grid, Pres)

spec_occ_wide_P.1 <- pivot_wider(spec_occ_P.1, names_from = Spec, values_from = Pres)  
spec_occ_wide_P.2  <- pivot_wider(spec_occ_P.2, names_from = Spec, values_from = Pres) 
spec_occ_wide_P.3 <- pivot_wider(spec_occ_P.3, names_from = Spec, values_from = Pres) 

# Joint occurence and predictor data
df_datjoint_P.1 <- data.frame(spec_occ_wide_P.1 %>% left_join(Pred10_P.1))
df_datjoint_P.2 <- data.frame(spec_occ_wide_P.2 %>% left_join(Pred10_P.2, by='grid'))
df_datjoint_P.3 <- data.frame(spec_occ_wide_P.3 %>% left_join(Pred10_P.3, by='grid'))

df_datjoint_P.1[,grep("Sp", names(df_datjoint_P.1))][is.na(df_datjoint_P.1[,grep("Sp", names(df_datjoint_P.1))])] <- 0
df_datjoint_P.2[,grep("Sp", names(df_datjoint_P.2))][is.na(df_datjoint_P.2[,grep("Sp", names(df_datjoint_P.2))])] <- 0
df_datjoint_P.3[,grep("Sp", names(df_datjoint_P.3))][is.na(df_datjoint_P.3[,grep("Sp", names(df_datjoint_P.3))])] <- 0

# rows are species presence/absence on grid cells + values for each env variable at each grid cell


par(mfrow = c(1,1))
plot(tmp_seas.1_P.1)
points(data = df_datjoint_P.1[which(df_datjoint_P.1$Sp53 == 1),], lat~long, add = T)
points(data = df_datjoint_P.3[which(df_datjoint_P.1$Sp53 == 1),], lat~long, qdd = T, col = "blue")



#############
# Run SDMs


# sdmData
d_P.1 <- sdmData(~tmp_seas.1_P.1+tmp_seas.2_P.1+pre_seas.1_P.1+pre_seas.2_P.1+pForest_P.1+pGrass_P.1+pCrop_P.1+pSettlem_P.1+coords(long+lat),train=df_datjoint_P.1) 
d_P.2 <- sdmData(~tmp_seas.1_P.2+tmp_seas.2_P.2+pre_seas.1_P.2+pre_seas.2_P.2+pForest_P.2+pGrass_P.2+pCrop_P.2+pSettlem_P.2+coords(long+lat),train=df_datjoint_P.2)  
d_P.3 <- sdmData(~tmp_seas.1_P.3+tmp_seas.2_P.3+pre_seas.1_P.3+pre_seas.2_P.3+pForest_P.3+pGrass_P.3+pCrop_P.3+pSettlem_P.3+coords(long+lat),train=df_datjoint_P.3)  

# bg= list(n=1000, method='gRandom')
# d <- sdmData(formula, train, test, predictors, bg, crs, ...)

for(sp in paste0("Sp",spec_speccode)){

	setwd(paste0(dir.analysis, "/sdm_output"))
	filename <- paste0('model_P.1_',sp,'.sdm')

	if(!file.exists(filename)){
		fo_P.1 <- as.formula(paste(sp,'~ tmp_seas.1_P.1+tmp_seas.2_P.1+pre_seas.1_P.1+pre_seas.2_P.1+pForest_P.1+pGrass_P.1+pCrop_P.1+pSettlem_P.1')) # formula for each species
		##model_P.1 <- sdm(fo_P.1, data=d_P.1, methods=c('glm','brt'),replication='boot',n=10,parallelSettings=list(ncore=2,method='parallel',fork=F))
		model_P.1 <- sdm(fo_P.1, data=d_P.1, methods=c('glm','gam', 'fda', 'svm', 'gbm'),replication='sub', test.p=30, n=10,parallelSettings=list(ncore=2,method='parallel',fork=F))
		setwd(paste0(dir.analysis, "/sdm_output"))
		##setwd("C:/USers/k.l.wells/Kons/BTO_analysis")
		write.sdm(model_P.1, paste0('model_P.1_',sp,'.sdm'), overwrite=TRUE)
	}else{}
	print(paste0(sp, ", Species " , which(paste0("Sp",spec_speccode)==sp), "out of ", length(spec_speccode)))	
}



###########################################################

#Names for SDM variable importance metrics
varimp_names <- c(paste0("AUCtest_", names(Loc10)[-c(1:3)]), paste0("corTest_", names(Loc10)[-c(1:3)]))
# Dataframe SDM variable importance metrics
df_SDM.varimp <- data.frame(array(NA, dim=c(nSpec, (length(varimp_names)+1))))
colnames(df_SDM.varimp) <- c('speccode', varimp_names)
df_SDM.varimp$speccode <- spec_speccode

for(sp in spec_speccode){
	# Prediction/ ensemble from sdm output

	setwd(paste0(dir.analysis, "/sdm_output"))
	
	if(file.exists(paste0('model_P.1', '_', paste0("Sp",sp),'.sdm'))){
		model_P.1 <- read.sdm(paste0('model_P.1', '_', paste0("Sp",sp),'.sdm'))
	    	varimp_P.1_auc <-t(getVarImp(model_P.1)@varImportanceMean$AUCtest)["AUCtest",]
	    	varimp_P.1_cor <-t(getVarImp(model_P.1)@varImportanceMean$corTest)["corTest",]
		df_SDM.varimp[df_SDM.varimp$speccode==sp, paste0("AUCtest_", names(varimp_P.1_auc)) ] <- as.numeric(varimp_P.1_auc)
		df_SDM.varimp[df_SDM.varimp$speccode==sp, paste0("corTest_", names(varimp_P.1_cor)) ] <- as.numeric(varimp_P.1_cor)
 	}else{}

	if(file.exists(paste0('model_P.2', '_', paste0("Sp",sp),'.sdm'))){
		model_P.2 <- read.sdm(paste0('model_P.2', '_', paste0("Sp",sp),'.sdm'))
	    	varimp_P.2_auc <-t(getVarImp(model_P.2)@varImportanceMean$AUCtest)["AUCtest",]
	    	varimp_P.2_cor <-t(getVarImp(model_P.2)@varImportanceMean$corTest)["corTest",]
		df_SDM.varimp[df_SDM.varimp$speccode==sp, paste0("AUCtest_", names(varimp_P.2_auc)) ] <- as.numeric(varimp_P.2_auc)
		df_SDM.varimp[df_SDM.varimp$speccode==sp, paste0("corTest_", names(varimp_P.2_cor)) ] <- as.numeric(varimp_P.2_cor)
 	}else{}

	if(file.exists(paste0('model_P.3', '_', paste0("Sp",sp),'.sdm'))){
		model_P.3 <- read.sdm(paste0('model_P.3', '_', paste0("Sp",sp),'.sdm'))
	    	varimp_P.3_auc <-t(getVarImp(model_P.3)@varImportanceMean$AUCtest)["AUCtest",]
	    	varimp_P.3_cor <-t(getVarImp(model_P.3)@varImportanceMean$corTest)["corTest",]
		df_SDM.varimp[df_SDM.varimp$speccode==sp, paste0("AUCtest_", names(varimp_P.3_auc)) ] <- as.numeric(varimp_P.3_auc)
		df_SDM.varimp[df_SDM.varimp$speccode==sp, paste0("corTest_", names(varimp_P.3_cor)) ] <-as.numeric( varimp_P.3_cor)
 	}else{}
}


## Save
setwd(dir.analysis)
save(df_SDM.varimp, file="df_SDM.varimp.RData")
write.csv(df_SDM.varimp, file="df_SDM.varimp_210707.csv")

load("df_SDM.varimp.RData")



