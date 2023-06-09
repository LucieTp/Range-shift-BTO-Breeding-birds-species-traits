#### 12/2021 Species Traits analysis based on 160 species
#### filtered by GB + >50 grid cells in each time period + Breeding Atlas



####
#### This script implements the analyses accompanying the paper: 'Species
#### trait determine range shifts responses to warming in UK birds'
####
#### Author: Dr Miguel Lurgi
#### Computational Ecology Lab, Swansea University. UK.
####

#### required libraries

require(igraph)
library(rnaturalearth)
library(sf)
require(raster)
require(rgdal)
library(foreign) # for read dbf

#### In first script we ran SDMs for all the bird species using environmental and land use predictor variables

setwd("E:/TheseSwansea/TraitStudy/Github/Range-shift-BTO-breeding-birds")
species_traits <- read.csv('data/SpecTrait_122021.csv', stringsAsFactors = FALSE, row.names = 1)

#### reading the range shift data 
range_shifts <- read.csv('data/df_rangeshift_052023.NoOrkneyNoShetlands.csv', stringsAsFactors = FALSE, row.names = 1)

species_traits <- merge(species_traits, range_shifts)

###########################################################################################################
#### We complement the Elton traits data with topological species traits from the European food web:

#### Here we load the European trophic network (metaweb) that we used in the European
#### bioregions paper to calculate species traits related to their position in the web: 
#### Galiana, N., Barros, C., Braga, J., Ficetola, G.F., Maiorano, L., Thuiller, W., Montoya, J.M. and Lurgi, M. (2021)
#### The spatial scaling of food web structure across European biogeographical regions. 
#### Ecography, 44: 653-664. https://doi.org/10.1111/ecog.05229


metaweb_europe <- read.graph('data/metaweb-europe.graphml', format = 'graphml')
m <- as_adjacency_matrix(metaweb_europe, attr = 'copresence', sparse=FALSE)
m[is.nan(m)] <- 0

#### The metaweb that we use below (from the European bioregions paper) has the species names
#### in a different format. Here we re-format the names to match those in the metaweb.
species_traits$scientific_name_ET_ <- sub(' ', '_', species_traits$scientific_name_ET)

#### and then check that all species in our database are in the metaweb
setdiff(species_traits$scientific_name_ET_, V(metaweb_europe)$name)

#### To make the network calculations more accurate, we subset the European metaweb to only
#### encompass those species found in the UK. To do this, we use the European map and species
#### distributions from the European bioregions paper and overlay on it a layer containing the map
#### of the UK to extract the species only present there


## functions needed to compute network metrics
source("whois-function.R")
source("utils.R")

europeRaster <- raster(x="mask10k/reference_grid_10km.img")
cells_info <- read.dbf('mask10k/reference_grid_10km.img.vat.dbf')

########################################################################################

shape <- ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")
shape <- st_transform(shape, crs(europeRaster))

#### Plot the shape to make sure it corresponds to the UK
cur_ids <- unlist(raster::extract(europeRaster, shape))
cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)

## get the master dtf
load("data/MASTER.bin10000_allhab_tresh0.RData")

region <- data.frame(PAGENAME = master$PAGENAME, SPP = 0)
region$SPP <- 0
region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1

bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
plot(bioregion_raster)

#### from the master table (containing all species distributions) we filter out
#### only those cells falling on the UK map
#master = read.csv2("E:/TheseSwansea/Galiana2021_network-area-europe-master/master.csv")

cur_comm <- master[which(master$PAGENAME %in% cur_codes),]
cur_comm[is.na(cur_comm)] <- 0

row.names(cur_comm) <- cur_comm[,1]
cur_comm <- cur_comm[-1]

#### This removes the species that are not present in any of the cells of the UK
cur_comm <- cur_comm[,-which(colSums(cur_comm) == 0)]
cur_comm <- colSums(cur_comm)
cur_comm[cur_comm > 1] <- 1

#### These are the species in the United Kingdom
sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
uk_species <- unlist(lapply(sppSTRING, function(x){whois(SPPCODE = x)}))

#### once we obtain all the species only present in the UK, we can subset the metaweb 
#### to only consider the metaweb of the UK
m_uk <- m[which(row.names(m) %in% uk_species), which(colnames(m) %in% uk_species)]

# some names don't correspond to the names in the metaweb
species_traits$scientific_name_web <- str_replace(species_traits$scientific_name_ET,' ', '_')
species_traits[which(is.na(match(species_traits$scientific_name_web, row.names(m_uk)))),]$scientific_name_web <- c("Branta_canadensis","Stercorarius_skua","Saxicola_torquata")
# can't find canada goose in the metaweb Branta canadensis

library(stringr)
row.names(m)[which(str_detect(row.names(m),"Branta"))] # not at all in the eu metaweb,


################################################################################
## TROPHIC METRICS - Outdegree / Indegree / diet diversity / Trophic position 
################################################################################

## 1. Number of predators - outdegree / vulberability
species_traits$outdegree <- rowSums(m_uk[match(species_traits$scientific_name_web, row.names(m_uk)),])

## 2. Number of prey - Indegree
species_traits$indegree <- colSums(m_uk[,match(species_traits$scientific_name_web, row.names(m_uk))])
species_traits$normalised_indegree <- species_traits$indegree/dim(m_uk)[1]
#### for the species with 0 in-degree (due to the lack of resolution of species at the basal level in the metaweb)
#### we use information from the diet diversity provided by the BTO species traits:

diet_matrix <- species_traits[which(species_traits$indegree == 0),][15:24]
species_traits$scientific_name_ET_[duplicated(species_traits$scientific_name_ET_)] = "Corvus_cornix"
row.names(diet_matrix) <- species_traits[which(species_traits$indegree == 0),]$scientific_name_ET_
require(vegan)
diet_div <- vegan::diversity(diet_matrix, index='simpson')
species_traits[match(names(diet_div), species_traits$scientific_name_ET_),]$normalised_indegree <- diet_div

## 3. diet diversity from Elton traits (as basal species above) but for all sp
diet_matrix <- species_traits[15:24]
row.names(diet_matrix) <- species_traits$scientific_name_ET_
diet_div <- vegan::diversity(diet_matrix)
species_traits$diet_diversity <- diet_div


## 4. trophic positions for all species according to the metaweb
# subset the EU metaweb for species in the UK
m_uk_160 = as.matrix(m_uk[which(colnames(m_uk)%in%species_traits$scientific_name_web),
                which(rownames(m_uk)%in%species_traits$scientific_name_web)])
trophic_pos <- TrophicPositionsNew(m_uk) ## this function comes from utils.r code
names(trophic_pos) = c("scientific_name_web","trophic_position")
species_traits = merge(species_traits, trophic_pos, all.x = T)

## the canada goose is missing from the metaweb
species_traits$scientific_name_web[is.na(species_traits$trophic_position)]

################################################################################
## ABIOTIC METRICS - Habitat generality / preference 
################################################################################

#### we look at the complement of habitat diversity from the variables identified as important by the SDMs
#### We focus on the covariates concerning habitat specificity: Forest, Grass, Crop, Settlement, for the original distribution
#### of the species (i.e. P1). We do this using the 'Area Under the Curve' (AUC) metric for the habitat covariates identified
#### above and looking into their relative contribution to the prediction of species ranges for the first time period (i.e.
#### historical / original habitat specificity). We use Shannon diversity to quantify habitat generality at the species level

## sdms vars come from script 1.BTO_breeding.atlas_range.shifts.r
sdms_vars <- read.csv('data/df_SDM.varimp_122021.csv', stringsAsFactors = F, row.names = 1)
head(sdms_vars)

var_to_plot <- c('AUCtest_tmp_seas.1_P.1', 'AUCtest_tmp_seas.2_P.1', 'AUCtest_pre_seas.1_P.1', 'AUCtest_pre_seas.2_P.1', 'AUCtest_pForest_P.1', 'AUCtest_pGrass_P.1', 'AUCtest_pCrop_P.1', 'AUCtest_pSettlem_P.1')
par(mfrow=c(3,3))
for(i in var_to_plot){hist(sdms_vars[,i],xlab=i,main=i, breaks =30)}

#### determining habitats for period P.1
habitats <- sdms_vars[match(species_traits$speccode, sdms_vars$speccode), c('AUCtest_tmp_seas.1_P.1', 'AUCtest_tmp_seas.2_P.1', 'AUCtest_pre_seas.1_P.1', 'AUCtest_pre_seas.2_P.1', 'AUCtest_pForest_P.1', 'AUCtest_pGrass_P.1', 'AUCtest_pCrop_P.1', 'AUCtest_pSettlem_P.1')]
#### determining habitats for period P.3
habitats_3 <- sdms_vars[match(species_traits$speccode, sdms_vars$speccode), c('AUCtest_tmp_seas.1_P.3', 'AUCtest_tmp_seas.2_P.3', 'AUCtest_pre_seas.1_P.3', 'AUCtest_pre_seas.2_P.3', 'AUCtest_pForest_P.3', 'AUCtest_pGrass_P.3', 'AUCtest_pCrop_P.3', 'AUCtest_pSettlem_P.3')]

#### We include this information in the species traits database.
species_traits$habitat_gen <- vegan::diversity(habitats)


### habitat generality can also be assessed by looking at the scores of principal components from PCA for each species

library(corrplot)
library(RColorBrewer)


### Climatic
env <- sdms_vars[match(species_traits$speccode, sdms_vars$speccode), c('AUCtest_tmp_seas.1_P.1', 'AUCtest_tmp_seas.2_P.1', 'AUCtest_pre_seas.1_P.1', 'AUCtest_pre_seas.2_P.1')]
row.names(env) <- species_traits$speccode

### Land Cover
lc <- sdms_vars[match(species_traits$speccode, sdms_vars$speccode), c('AUCtest_pForest_P.1', 'AUCtest_pGrass_P.1', 'AUCtest_pCrop_P.1', 'AUCtest_pSettlem_P.1')]
row.names(lc) <- species_traits$speccode

###############################################

### PCA ###
#### first for the environment
pca_env <- prcomp(env, center = TRUE, scale = TRUE)

#### this indicates that the first two components encapsulate 83% of the total variability. So we use those.
summary(pca_env)
biplot(pca_env)
biplot(pca_env,  cex=c(0.2,0.8))
# The following two lines get added to the plot.
abline(h=0, col="grey", lty=4)
abline(v=0, col="grey", lty=4)

library(factoextra)
var_env = get_pca_var(pca_env)
pca_env$rotation

### The Kaiser’s rule (Kaiser-Guttman criterion) is a widely used method to evaluate the maximum number of linear 
### combinations to extract from the data set. According to that rule only those principal components are retained, 
### whose variances exceed 1. The idea behind the Kaiser-Guttman criterion is that any principal with variance less 
### than 1 contains less information than one of the original variables and so is not worth retaining
### The eigenvalues are:
pca_env$sd^2    
#### which suggests that only PC1 and PC2 are relevant to this dataset

#### To obtain the appropriate values for each component we need to apply some transformations
#### as suggested by Hendrickson, A.E., & White, P.O. (1964).

n.axes <- 2
rawLoadings  <- pca_env$rotation[, 1:n.axes] %*% diag(pca_env$sdev, n.axes, n.axes)
pca.scores.x <- scale(pca_env$x[, 1:n.axes]) %*% promax(rawLoadings, m = 4)$rotmat 		# m = 4 recommended by Hendrickson, A.E., & White, P.O. (1964). 

#### once we get the transformed scores, we extract the data for the PC1/2 only as it is the 
#### only one identified as relevant

species_traits$pc1_env <- NA
species_traits$pc2_env <- NA
species_traits[match(rownames(pca_env$x), species_traits$speccode),]$pc1_env <- as.numeric(pca.scores.x[,1])
species_traits[match(rownames(pca_env$x), species_traits$speccode),]$pc2_env <- as.numeric(pca.scores.x[,2])

factorloadings <- data.frame(cor(env, pca.scores.x[,1:2]))
names(factorloadings) = c('PC1', 'PC2')

library(stargazer)
stargazer(var_env$cor, summary = F, type = "html", title = 'Correlation coefficient between the first two components and environmental variables')
stargazer(get_eigenvalue(pca_env), summary = F, type = "html", title = 'Proportion of variance explained by each component - climatic variables')

### The factor loadings (line above) suggest that PC1 is highly negatively correlated to precipitation whereas 
### PC2 is highly negatively correlated to temperature

### We repeat the same procedure for the land cover variables

pca_lc <- prcomp(lc, center = TRUE, scale = TRUE)
var = get_pca_var(pca_lc)
var




#### this indicates that the first two components encapsulate 80% of the total variability. So we use those.
summary(pca_lc)
biplot(pca_lc)
pca_lc$rotation

### Again, using Kaiser's rule, we look at the eigenvalues
#### The eigenvalues are:
pca_lc$sd^2    #### = 1.6320958 0.9452559 0.8406980 0.5819503
#### in this case, only PC1 is relevant to this dataset
#### however, since settlement correlation is very low with PC1, we also take PC2

#### To obtain the appropriate values for each component we need to apply some transformations
#### as suggested by Hendrickson, A.E., & White, P.O. (1964).

n.axes <- 2
rawLoadings  <- pca_lc$rotation[, 1:n.axes] %*% diag(pca_lc$sdev, n.axes, n.axes)
pca.scores.x <- scale(pca_lc$x[, 1:n.axes]) %*% promax(rawLoadings, m = 4)$rotmat 		# m = 4 recommended by Hendrickson, A.E., & White, P.O. (1964). 

#### once we get the transformed scores, we extract the data for the PC1 only as it is the 
#### only one identified as relevant
species_traits$pc1_lc <- NA
species_traits$pc2_lc <- NA
species_traits[match(rownames(pca_lc$x), species_traits$speccode),]$pc1_lc <- pca.scores.x[,1]
species_traits[match(rownames(pca_lc$x), species_traits$speccode),]$pc2_lc <- pca.scores.x[,2]

### These factor loadings suggest PC1 is related to 'good' habitats and PC2 to 'destroyed' habitats
factorloadings <- data.frame(cor(lc, pca.scores.x[,1:2]))
names(factorloadings) = c('PC1', 'PC2')
stargazer(var$cor, summary = F, type = "html", title = 'Correlation coefficient between the first two components and landcover variables')
stargazer(get_eigenvalue(pca_lc), summary = F, type = "html", title = 'Proportion of variance explained by each component - landcover variables')

g1 = fviz_eig(pca_lc, addlabels = TRUE, title = 'Landcover variables')
g2 = fviz_eig(pca_env, addlabels = TRUE, title = 'Climatic variables')
library(ggpubr)
ggarrange(g1,g2)

ggsave('ExplainedVariancePCA.jpeg', dpi = 600)

################################################################################
## MIGRATORY BEHAVIOUR
################################################################################


#### Another ecological trait that might play a role in range shifts is the migratory behaviour
#### of the different species. We fetch these data here and add it to our analyses

migratory_behaviour <- read.table('data/specieslist3_1_migbehav_v1_0.txt', sep='\t', header = TRUE)

species_traits$migratory <- migratory_behaviour[match(species_traits$scientific_name_ET_, migratory_behaviour$IOC3_1_Binomial),]$Migratory_status
# a few species have names that don't correspond to the migratory behaviour dtf so we extract them from the BTO scientific names
species_traits$scientific_name_ <- sub(' ', '_', species_traits$scientific_name)
setdiff(species_traits$scientific_name_, migratory_behaviour$IOC3_1_Binomial)
species_traits$migratory[which(is.na(species_traits$migratory))] <- 
  migratory_behaviour[match(species_traits$scientific_name_[which(is.na(species_traits$migratory))], migratory_behaviour$IOC3_1_Binomial),]$Migratory_status

summary(droplevels(factor(species_traits$migratory)))
summary(factor(species_traits$migrate.strategy.3))

#### among all the migratory categories in this database we chose only the species tagger as 'directional migratory'
#### i.e. those migrating to and fro breeding grounds every year, as migrants
species_traits$migratory_binomial <- 0
species_traits[which(species_traits$migratory == 'directional migratory'),]$migratory_binomial <- 1


species_traits = species_traits[-which(is.na(species_traits$trophic_position)),]
species_traits$core_dist_bin = ifelse(species_traits$distrib.core == "south",1,0)

summary(species_traits)
## 11/01/2022 : changed the indegree and outdegree using the right metaweb (UK not EU)
## 04/02/2022 : changed the trophic position using new function in utils (reran everything)
write.csv(species_traits, "data/SpecTrait_012023_159sp.csv")


################################################################################
## HAND WING INDEX

library(readxl)
HWI = read_excel("data/Global-HWI-master/Dataset HWI 2020-04-10.xlsx")

summary(HWI$HWI)
setdiff(species_traits$scientific_name, HWI$`Species name`)
setdiff(species_traits$scientific_name_ET, HWI$`Species name`)

species_traits$scientific_name_HWI = species_traits$scientific_name
species_traits$scientific_name_HWI[which(!species_traits$scientific_name_HWI %in% HWI$`Species name`)] = species_traits$scientific_name_ET[which(!species_traits$scientific_name_HWI %in% HWI$`Species name`)]

setdiff(species_traits$scientific_name_HWI, HWI$`Species name`)
species_traits$scientific_name_HWI[which(species_traits$scientific_name_HWI == "Phalacrocorax aristotelis")] = "Leucocarbo aristotelis"

species_traits = merge(species_traits, HWI[,c('Species name', 'HWI')], by.x = 'scientific_name_HWI', by.y = 'Species name')
summary(species_traits$HWI)


## 11/01/2022 : changed the indegree and outdegree using the right metaweb (UK not EU)
## 04/02/2022 : changed the trophic position using new function in utils (reran everything)
## 062023 : reran with mainland only and added HWI
write.csv(species_traits, "data/SpecTrait_Full_062023_159sp.csv")


