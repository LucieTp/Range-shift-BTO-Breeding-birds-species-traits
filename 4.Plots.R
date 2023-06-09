################################################################################
### CODE FOR PLOTTING THE MAIN FIGURE, TABLE AND SUPP FIG 1 AND 2
################################################################################


################################################################################
## IMPORTING THE DATA

library(ggplot2)
library(ggpubr) # for ggarrange
library(cowplot)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(stargazer)
library(stringr)

sf_UK  <- ne_countries(scale = "medium", country = 'United Kingdom', returnclass = "sf")

## Main table with all the species traits, range shift etc
setwd("E:/TheseSwansea/TraitStudy/Github")
species_traits <- read.csv('data/SpecTrait_012023_159sp.csv', stringsAsFactors = FALSE, row.names = 1)

# Marine or non marine species : (Marine = >75% of points within 20km of the coastline)
propmarine <- read.csv('data/Speciestraits_ProportionMarineBTOsp_159sp.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits = merge(species_traits, propmarine[,c('speccode','prop_Marine20km')])
species_traits$Marine = ifelse(species_traits$prop_Marine20km>75, 1, 0)

## adding range shift in km 
df.rangeshift <- read.csv('data/df_rangeshift_012023.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits = merge(species_traits, df.rangeshift[,c("speccode","shift_max20_P.1.3_dist_km", "shift_min20_P.1.3_dist_km","dist_N_km_max20_P.1","dist_S_km_min20_P.1")], all.x = T)

## BTO Data (occurrence data per grid cell for all of UK) - used for distribution maps
BTO_grid <- read.csv("data/grid_square_coordinates_lookup.csv", header=T)
BTO_distrib <- read.csv("data/distributions.csv", header=T) # period, sp code, season and grid for GB and Ireland


# Mean coordinates of each grid cell
Coord10 <- as_tibble(BTO_grid) %>%
  filter(resolution==10 & order<5) %>%
  group_by(grid) %>%
  summarise(long = mean(long), lat = mean(lat))

# Location data 
Loc10 <- BTO_distrib %>% dplyr::select(grid) %>% distinct(grid)
# Join coordinates to location data
Loc10 <- Loc10 %>%  left_join(Coord10)
# Filter missing coordinates
Loc10 <- Loc10 %>%  filter(!is.na(long) & !is.na(lat))


# Subset distribution data to Britain
BTO_distrib <- as_tibble(BTO_distrib) %>% filter(island=="B" & season=="B" & resolution==10)
BTO_distrib$periodN <- paste0("P.", as.numeric(droplevels(as.factor(BTO_distrib$period))))
BTO_distrib$Spec <- paste0("Sp", BTO_distrib$speccode)
BTO_distrib$Pres <- 1

# Join coordinates to occurrence data
BTO_distrib <- BTO_distrib %>%  left_join(Loc10, by='grid')
BTO_distrib <- BTO_distrib %>% filter(!is.na(long) & !is.na(lat))

BTO_distrib = merge(Loc10, BTO_distrib, all.y = T)



### wilcoxon test for statistically significant poleward or southward shifts 
### between period 1 and 3

# for example:
# >> for each species we test for significant difference (smaller or larger) in 
# latitude between the 20 leading-edge cells in P1 and in P3 
# >> if the latitude is significantly higher in P3 than P1, then we conclude that 
# the species shifted significantly northward etc 

BTO_distrib$Edge = NA
stable.leading = stable.rear = poleward.leading = poleward.rear = southward.rear = southward.leading = NULL
for (i in unique(species_traits$speccode)){
  for (p in unique(BTO_distrib$periodN)){
    
    loc = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$periodN == p),]
    edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

    BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$periodN == p),'Edge'] = edge
  }
  
  ## leading edge
  
  test = wilcox.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Leading' & BTO_distrib$periodN %in% c('P.1','P.3') ),], alternative = 'less', paired = F)
  if(test$p.value < 0.05 ){
      poleward.leading = c(poleward.leading, i)
    } else {
        test = wilcox.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Leading' & BTO_distrib$periodN %in% c('P.1','P.3') ),], alternative = 'greater', paired = F)
        if(test$p.value < 0.05 ){
          southward.leading = c(southward.leading, i)
        } 
        else {
          stable.leading = c(stable.leading, i)
        }
    }
    
    ## rear edge
    test1 = wilcox.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Rear' & BTO_distrib$periodN %in% c('P.1','P.3') ),], alternative = 'less', paired = F)
    
    if(test1$p.value < 0.05 ){
      poleward.rear = c(poleward.rear, i)
    } else {
      test1 = wilcox.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Rear' & BTO_distrib$periodN %in% c('P.1','P.3') ),], alternative = 'greater', paired = F)
      if(test1$p.value < 0.05 ){
        southward.rear = c(southward.rear, i)
      } 
      else {
        stable.rear = c(stable.rear, i)
      }
    }
}

length(poleward.leading) + length(southward.leading) + length(stable.leading) == 159

edge = poleward.leading
print(paste(length(edge)/159, 'of species shifted their leading edge significantly NORTH by', mean(species_traits[which(species_traits$speccode %in% edge),'shift_max20_P.1.3_dist_km']), 'sd', 
            sd(species_traits[which(species_traits$speccode %in% edge),'shift_max20_P.1.3_dist_km'])))
edge = southward.leading
print(paste(length(edge)/159, 'of species shifted their leading edge significantly SOUTH by', mean(species_traits[which(species_traits$speccode %in% edge),'shift_max20_P.1.3_dist_km']), 'sd', 
            sd(species_traits[which(species_traits$speccode %in% edge),'shift_max20_P.1.3_dist_km'])))
edge = stable.leading
print(paste(length(edge)/159, 'of species shifted did not shift their leading edge significantly'))

edge = poleward.rear
print(paste(length(edge)/159, 'of species shifted their leading edge significantly NORTH by', mean(species_traits[which(species_traits$speccode %in% edge),'shift_min20_P.1.3_dist_km']), 'sd', 
            sd(species_traits[which(species_traits$speccode %in% edge),'shift_min20_P.1.3_dist_km'])))
edge = southward.rear
print(paste(length(edge)/159, 'of species shifted their leading edge significantly SOUTH by', mean(species_traits[which(species_traits$speccode %in% edge),'shift_min20_P.1.3_dist_km']), 'sd', 
            sd(species_traits[which(species_traits$speccode %in% edge),'shift_min20_P.1.3_dist_km'])))
edge = stable.rear
print(paste(length(edge)/159, 'of species shifted did not shift their rear edge significantly'))

### t test are not appropriate given our data
t.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Leading' & BTO_distrib$periodN %in% c('P.1','P.3') ),])

bartlett.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Leading' & BTO_distrib$periodN %in% c('P.1','P.3') ),])
shapiro.test(x$shift[which(x$edge == 'min')])

## Difference in shift between N and S species
# H0: N species shifted at similar rates or faster than S species
wilcox.test(shift_max20_P.1.3_dist_km ~ distrib.core, data = species_traits, alternative = "less") 
wilcox.test(shift_min20_P.1.3_dist_km ~ distrib.core, data = species_traits, alternative = "less")

## looking at the difference in rear and leading-edge shifts of N or S distributed species

par(mfrow = c(2,1))
par(mar = c(2,2,2,2))
boxplot(species_traits[which(species_traits$distrib.core == "south"), 
                       c("shift_max20_P.1.3_dist_km","shift_min20_P.1.3_dist_km")], 
        main = "Southernly-distributed species", col = "red")

boxplot(species_traits[which(species_traits$distrib.core == "north"), 
                       c("shift_max20_P.1.3_dist_km","shift_min20_P.1.3_dist_km")], 
        main = "Northernly-distributed species", col = "blue")

#H0: S species' rear edge shifted at similar rates or faster than S species' leading edge
wilcox.test(species_traits[which(species_traits$distrib.core == "south"),"shift_min20_P.1.3_dist_km"],
            species_traits[which(species_traits$distrib.core == "south"),"shift_max20_P.1.3_dist_km"], alternative = "less")

#H0: N species' rear edge shift is no different from N species' leading edge shift
wilcox.test(species_traits[which(species_traits$distrib.core == "north"),"shift_min20_P.1.3_dist_km"],
            species_traits[which(species_traits$distrib.core == "north"),"shift_max20_P.1.3_dist_km"])

################################################################################

species_traits$PassNonPass = ifelse(species_traits$IOCOrder == "Passeriformes", "Passeriformes","non-Passeriformes")
species_traits$shift_diff = (species_traits$shift_max20_P.1.3_dist_km - species_traits$shift_min20_P.1.3_dist_km)
species_traits$migratory_binomial = factor(species_traits$migratory_binomial)
levels(species_traits$migratory_binomial) = c('Resident','Migrant')

shift_long = pivot_longer(data = species_traits[which(species_traits$Marine == 0),c("speccode","english_name", "distrib.core", "PassNonPass", "IOCOrder", "BodyMass.Value","migratory_binomial","shift_max20_P.1.3_dist_km","shift_min20_P.1.3_dist_km","shift_diff")], cols = c("shift_max20_P.1.3_dist_km","shift_min20_P.1.3_dist_km","shift_diff" ), names_to = "Period",values_to = "Shift")

## subsetting only non coastal species
species_traits = species_traits[which(species_traits$Marine == 0),]


################################################################################
### TABLE S1
################################################################################

names = c('scientific_name','english_name',"Marine",'prop_Marine20km',"dist_S_km_min20_P.1","dist_N_km_max20_P.1",'P.1','lat_mean.P.1','distrib.core','IOCOrder','BodyMass.Value','shift_max20_P.1.3_dist_km','shift_min20_P.1.3_dist_km','shift_diff','outdegree',"pc1_env","pc2_env","pc1_lc" ,"pc2_lc", "normalised_indegree","diet_diversity","trophic_position","habitat_gen")
towrite = species_traits[,names]
names(towrite) = c('Scientific name','English name',"Marine species", 'Proportion of coastal grid cells',"Southern boundary effect","Northern boundary effect",'Number of grid cells in P1','Mean latitude of grid cells in P1','Distribution core','Taxonomy','Body Mass (g)','Shift in leading edge between P1 and P2 (km)','Shift in rear edge between P1 and P2 (km)','Range expansion (km)','Vulnerability',"Association with precipitation seasonality","Association with temperature seasonality","Association with forest and grassland" ,"Association with urban and cropland", "Normalised indegree","Diet diversity","Trophic position","Habitat generality")
towrite[,c(4:8,11:23)] = round(towrite[,c(4:8,11:23)], 2)
write.csv(towrite, 'SpeciesTraitsTable_paper1.csv')

################################################################################
### Outputs of PGLMMs in a nicer format


tb = read.csv2("results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.HWI.csv", sep =",", dec = ".")
tb_km = read.csv2("pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.km.csv", sep =",", dec = ".")

tb[,c("p.value.phylo", "r2", str_subset(names(tb), 'estimate'))] = round(tb[,c("p.value.phylo", "r2", str_subset(names(tb), 'estimate'))],2)

################################################################################
## TABLE 1             #########################################################
################################################################################


tb1 = pivot_wider(tb[which(tb$value %in% c("significance","estimate_mean")),-1], values_from = names(tb)[8:22], names_from = value)
tb1 = tb1[-which(str_detect(tb1$model, "Non passeriformes")),]
tb1$shift = ifelse(tb1$shift == "max", "leading edge shift",ifelse(tb1$shift == "min","rear edge shift", "Expansion/Retraction"))
tb1 = tb1[,-which(str_detect(colnames(tb1),"Loglik|phylo.mean"))]
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = apply(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))], 2, as.numeric)
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = round(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))],2)


test = tb1
test[,c(7:36)] = apply(test[,c(7:36)], 2, as.character)
test[is.na(test)] = ''
for(i in seq(7,36,2)) {test[,i] = apply(test[,c(i, i+1)], 1, function(x) paste(x[1], x[2]))}
test = test[,-str_which(names(test), 'significance')]

stargazer(test, summary = F, type = "html")

################################################################################
## TABLE S3
## same with coef in km
################################################################################

tb1 = pivot_wider(tb_km[which(tb_km$value %in% c("significance","estimate_mean")),-1], values_from = names(tb_km)[8:22], names_from = value)
tb1 = tb1[-which(str_detect(tb1$model, "Non passeriformes")),]
tb1$shift = ifelse(tb1$shift == "max", "leading edge shift",ifelse(tb1$shift == "min","rear edge shift", "Expansion/Retraction"))
tb1 = tb1[,-which(str_detect(colnames(tb1),"Loglik|phylo.mean"))]
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = apply(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))], 2, as.numeric)
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = round(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))],2)


test = tb1
test[,c(7:36)] = apply(test[,c(7:36)], 2, as.character)
test[is.na(test)] = ''
for(i in seq(7,36,2)) {test[,i] = apply(test[,c(i, i+1)], 1, function(x) paste(x[1], x[2]))}
test = test[,-str_which(names(test), 'significance')]

stargazer(test, summary = F, type = "html")

################################################################################
## TABLE S4
## table of results without Orkney and the Shetlands
################################################################################

tb2 = read.csv('E:/TheseSwansea/TraitStudy/Github/Range-shift-BTO-breeding-birds/results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.NoShetlandsNoOrkney.csv')

tb2 = pivot_wider(tb2[which(tb2$value %in% c("significance","estimate_mean")),-1], values_from = names(tb2)[8:22], names_from = value)
tb2 = tb2[-which(str_detect(tb2$model, "Non passeriformes")),]
tb2$shift = ifelse(tb2$shift == "max", "leading edge shift",ifelse(tb2$shift == "min","rear edge shift", "Expansion/Retraction"))
tb2 = tb2[,-which(str_detect(colnames(tb2),"Loglik|phylo.mean"))]
tb2[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb2), 'estimate'))] = apply(tb2[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb2), 'estimate'))], 2, as.numeric)
tb2[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb2), 'estimate'))] = round(tb2[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb2), 'estimate'))],2)


tb2[,c(7:36)] = apply(tb2[,c(7:36)], 2, as.character)
tb2[is.na(tb2)] = ''
for(i in seq(7,36,2)) {tb2[,i] = apply(tb2[,c(i, i+1)], 1, function(x) paste(x[1], x[2]))}
tb2 = tb2[,-str_which(names(tb2), 'significance')]

stargazer(tb2, summary = F, type = "html")


#################################################################################
## FIG 1 and Fig S2 #############################################################
#################################################################################
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)

loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 95 & BTO_distrib$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Northern.dist = ggplot(data = sf_UK) + geom_sf() + geom_point(data = loc, aes(x = long, y = lat), size = 1, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge), size = 1.2) + 
  geom_hline(yintercept = max(BTO_distrib$lat), linetype="dashed") + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed") + 
  ggtitle("Northern species") + ylab("") + xlab("") + theme_classic(base_size =  18) 

loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 272 & BTO_distrib$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Southern.dist = ggplot(data = sf_UK) + geom_sf() + geom_point(data = loc, aes(x = long, y = lat), size = 1, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge), size = 1.2) + 
  geom_hline(yintercept = max(BTO_distrib$lat), linetype="dashed") + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed") + 
  ggtitle("Southern species") + ylab("") + xlab("") + theme_classic(base_size =  20)

### model results

library(tidyverse)
tb[,8:22] = apply(tb[,8:22], 2, as.numeric)

tb2 = pivot_longer(tb, cols = 8:22, names_to = "covariate", values_to = "estimate")
tb2 = pivot_wider(tb2[,2:13], values_from = estimate, names_from = value)
# tb2 = tb2[-which(str_detect(tb2$model, "Non passeriformes")),]
tb2$shift = ifelse(tb2$shift == "max", "leading edge shift",ifelse(tb2$shift == "min","rear edge shift", "Expansion"))

tb2$significance = ifelse(tb2$p.value_mean<0.001,"***",ifelse(tb2$p.value_mean<0.01,"**",ifelse(tb2$p.value_mean<0.05,"*",ifelse(tb2$p.value_mean<0.1,".",""))))
tb2$size = ifelse(tb2$p.value_mean<0.05,0.6,0.5)
tb2$covariate = factor(tb2$covariate)
levels(tb2$covariate) = c("Diet diversity", "Northern boundary effect", "Southern boundary effect","Habitat generality",'Hand wing index',"Log10 body mass", "Normalised number of prey", "Number of predators", "Range size", "Precipitation seasonality", "Association with forest/grasslands","Temperature seasonality","Association with urban areas/croplands","Trophic position","Intercept")

tb2$CI_lo = tb2$estimate_mean - 1.96*tb2$estimate_se
tb2$CI_up = tb2$estimate_mean + 1.96*tb2$estimate_se

tb2$SE_lo = tb2$estimate_mean - tb2$estimate_se
tb2$SE_up = tb2$estimate_mean + tb2$estimate_se

tb2$shift = as.factor(tb2$shift)
tb2$shift = relevel(tb2$shift, "rear edge shift")
levels(tb2$shift)

biogeo = tb2[which(tb2$covariate %in% c("Northern boundary effect", "Southern boundary effect", "Range size", "Intercept")),]
traits = tb2[which(!tb2$covariate %in% c("Northern boundary effect", "Southern boundary effect", "Range size",'Intercept')),]
traits_nodiff = traits[which(traits$shift!='Expansion'),]

### BIOGEOGRAPHICAL COVARIATES

g1 = ggplot(data = biogeo[which(biogeo$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1.5, height = .3) + xlim(min(biogeo[which(biogeo$model == "northern"),'CI_lo']),max(biogeo[which(biogeo$model == "northern"),'CI_up'])) + scale_y_discrete(labels = c('','Intercept','','','Northern boundary effect','','','Range size','','','Southern boundary effect',''))+ ggtitle('Northern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,9.5,3), colour = 'grey', linetype = 'longdash') 
g2 = ggplot(data = biogeo[which(biogeo$model == "southern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1.5, height = .3) + xlim(min(biogeo[which(biogeo$model == "southern"),'CI_lo']),max(biogeo[which(biogeo$model == "southern"),'CI_up'])) + scale_y_discrete(labels = c('','Intercept','','','Northern boundary effect','','','Range size','','','Southern boundary effect',''))+ ggtitle('Southern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,9.5,3), colour = 'grey', linetype = 'longdash')
g3 = ggplot(data = biogeo[which(biogeo$model == "Passeriformes"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1.5, height = .3) + xlim(min(biogeo[which(biogeo$model == "Passeriformes"),'CI_lo']),max(biogeo[which(biogeo$model == "Passeriformes"),'CI_up'])) + scale_y_discrete(labels = c('','Intercept','','','Northern boundary effect','','','Range size','','','Southern boundary effect',''))+ ggtitle('Passeriformes')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,9.5,3), colour = 'grey', linetype = 'longdash')
g4 = ggplot(data = biogeo[which(biogeo$model == "terrestrial"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1.5, height = .3) + xlim(min(biogeo[which(biogeo$model == "terrestrial"),'CI_lo']),max(biogeo[which(biogeo$model == "terrestrial"),'CI_up'])) + scale_y_discrete(labels = c('','Intercept','','','Northern boundary effect','','','Range size','','','Southern boundary effect',''))+ ggtitle('All species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,9.5,3), colour = 'grey', linetype = 'longdash')

# g1 = ggplot(data = biogeo[which(biogeo$model == "northern"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE) + xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))
# g2 = ggplot(data = biogeo[which(biogeo$model == "southern"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))
# g3 = ggplot(data = biogeo[which(biogeo$model == "Passeriformes"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))
# g4 = ggplot(data = biogeo[which(biogeo$model == "terrestrial"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))

library(ggpubr)
## Figure S2 - Estimates from GLMM 
ggarrange(g1,g2,g3,g4)

library(cowplot)
setwd("E:/TheseSwansea/TraitStudy/code_Miguel/Plots")
ggsave2('ModelEstimates_Biogeo.jpeg', scale = 3.5)

### SPECIES TRAITS
labels = c('','Association with forest/grassland','','','Association with urban area/cropland','','','Diet diversity','','','Habitat generality','','','Log 10 body mass','','','Migratory status - Migrant','','','Normalised number of prey','','','Number of predators','','','Precipitation seasonality','','','Temperature seasonality','','','Trophic position','')

g1 = ggplot(data = traits[which(traits$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Northern species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels = labels) + ggtitle('Northern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 15) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,30.5,3), colour = 'grey', linetype = 'longdash')
g2 = ggplot(data = traits[which(traits$model == "southern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Southern species')+ geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up),size = 1, height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels = labels)+ ggtitle('Southern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 15) + theme(axis.ticks = element_blank())  + geom_hline(yintercept =seq(3.5,30.5,3), colour = 'grey', linetype = 'longdash')
g3 = ggplot(data = traits[which(traits$model == "Passeriformes"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Passeriformes') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels =labels)+ ggtitle('Passeriformes')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 15) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,30.5,3), colour = 'grey', linetype = 'longdash')
g4 = ggplot(data = traits[which(traits$model == "terrestrial"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('All species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels =labels)+ ggtitle('All species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 15) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(3.5,30.5,3), colour = 'grey', linetype = 'longdash')

# g1 = ggplot(data = traits[which(traits$model == "northern"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('Northern species') + theme_classic(base_size =  20) + theme(legend.position = 'bottom')
# g2 = ggplot(data = traits[which(traits$model == "southern"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('Southern species')+ theme_classic(base_size =  20) + theme(legend.position = 'bottom')
# g3 = ggplot(data = traits[which(traits$model == "Passeriformes"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('Passeriformes')+ theme_classic(base_size =  20) + theme(legend.position = 'bottom')
# g4 = ggplot(data = traits[which(traits$model == "terrestrial"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('All species')+ theme_classic(base_size =  20) + theme(legend.position = 'bottom')

ggarrange(g1,Northern.dist,g3, g2,Southern.dist,g4, common.legend = T, labels = c('a','','c','b','','d'), font.label = list(size = 18))

setwd("E:/TheseSwansea/TraitStudy/code_Miguel")
## Figure 1 - Estimates from GLMM 
ggsave2('ModelEstimates_Traits_maps.jpeg', scale = 2.2, dpi = 800)


### SPECIES TRAITS - no expansion
labels = c('','','Forest/grassland','','Urban area/cropland','','Body mass','','Diet diversity','','Habitat generality','','Migratory status - Migrant','','Normalised indegree','','Precipitation seasonality','','Temperature seasonality','','Trophic position','','Vulnerability','')
labels = c('','Association with forest/grassland','','Association with urban area/cropland','','Diet diversity','','Habitat generality','','Hand wing index','','Log 10 body mass','','Normalised number of prey','','Number of predators','','Association with precipitation','','Association with temperature','','Trophic position','')


library(ggrepel)
traits_nodiff$shift = relevel(traits_nodiff$shift, "rear edge shift")
traits_nodiff$linesize = ifelse(traits_nodiff$significance %in% c('','.'), 1, 1.5)

g1 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Northern species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels = labels) + ggtitle('Northern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none')
g2 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "southern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Southern species')+ geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels = labels)+ ggtitle('Southern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none')
g3 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "Passeriformes"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Passeriformes') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels =labels)+ ggtitle('Passeriformes')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none')
g4 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "terrestrial"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('All species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels =labels)+ ggtitle('All species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none')

# g1 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + geom_vline(xintercept = 0) + geom_point(size = 2) + geom_text(hjust=0.5, vjust=0.35,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Northern species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1) + xlim(min(traits_nodiff[,'CI_lo']),max(traits_nodiff[,'CI_up'])) + scale_y_discrete(labels = rev(labels), expand = c(0, 1), limits=rev) + ggtitle('Northern species') + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 18) + theme(axis.ticks = element_blank())
# g2 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "southern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + geom_vline(xintercept = 0) + geom_point(size = 2) + geom_text(hjust=0.5, vjust=0.35,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Southern species')+ geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up),  size = 1, height = .1) + xlim(min(traits_nodiff[,'CI_lo']),max(traits_nodiff[,'CI_up'])) + scale_y_discrete(labels = rev(labels), expand = c(0, 1), limits=rev) + ggtitle('Southern species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 18) + theme(axis.ticks = element_blank())
# g3 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "Passeriformes"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance))+ geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + geom_vline(xintercept = 0) + geom_point(size = 2) + geom_text(hjust=0.5, vjust=0.35,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Passeriformes') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1) + xlim(min(traits_nodiff[,'CI_lo']),max(traits_nodiff[,'CI_up'])) + scale_y_discrete(labels = rev(labels), expand = c(0, 1), limits=rev)+ ggtitle('Passeriformes')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 18) + theme(axis.ticks = element_blank())
# g4 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "terrestrial"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + geom_vline(xintercept = 0) + geom_point(size = 2) + geom_text(hjust=0.5, vjust=0.35,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('All species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1) + xlim(min(traits_nodiff[,'CI_lo']),max(traits_nodiff[,'CI_up'])) + scale_y_discrete(labels = rev(labels), expand = c(0, 1), limits=rev)+ ggtitle('All species')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 18) + theme(axis.ticks = element_blank())

ggarrange(g1,Northern.dist,g3, g2,Southern.dist,g4, common.legend = T, labels = c('a','','c','b','','d'), widths = c(1.5,1,1.5), font.label = list(size = 20))


setwd("E:/TheseSwansea/TraitStudy/code_Miguel/Plots")
## Figure 1 - Estimates from GLMM 
ggsave('ModelEstimates_Traits_maps.nodiff.CI.jpeg', scale = 2.2, dpi = 1200, bg = "white")
ggsave('ModelEstimates_Traits_maps.nodiff.CI.pdf', scale = 2.2, dpi = 1200, bg = "white", device = "pdf")


g1bis = ggplot(data = traits_nodiff[which(traits_nodiff$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, label = significance)) + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + geom_vline(xintercept = 0) + geom_point(size = 2) + geom_text(hjust=0.5, vjust=0.35,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Northern species') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up), size = 1, height = .1) + xlim(min(traits_nodiff[,'CI_lo']),max(traits_nodiff[,'CI_up'])) + scale_y_discrete(expand = c(0, 1), limits=rev) 

ggarrange(g1bis, g1)

################################################################################
# non parametric t.test : Wilcoxon rank sum test with continuity correction
wilcox.test(shift_max20_P.1.3_dist_km ~ distrib.core, data = species_traits, alternative = "less") # H0: north shifted less than south
wilcox.test(shift_min20_P.1.3_dist_km ~ distrib.core, data = species_traits, alternative = "less")

x = data.frame(rbind(cbind(shift = species_traits$shift_max20_P.1.3_dist_km, edge = 'max'),
                 cbind(shift = species_traits$shift_min20_P.1.3_dist_km, edge = 'min')))
x$shift = as.numeric(x$shift)
t = t.test(shift ~ edge, data = x, alternative = 'greater', paired = T)
wilcox.test(shift ~ edge, data = x, alternative = 'greater', paired = T)

bartlett.test(shift ~ edge, data = x)
shapiro.test(x$shift[which(x$edge == 'min')])
ggqqplot(x, x = "shift", facet.by = 'edge')
## nope

ggplot(data = x, aes(y = shift, fill = edge)) + geom_boxplot()
ggplot(data = x, aes(x = shift, fill = edge)) + geom_histogram()

wilcox.test(dist_N_km_max20_P.1 ~ distrib.core, data = species_traits, alternative = "less") # H0: north shifted less than south
t.test(P.1 ~ distrib.core, data = species_traits)
t.test(pc1_env ~ distrib.core, data = species_traits)

t.test(P.1 ~ PassNonPass, data = species_traits)
t.test(pc2_lc ~ PassNonPass, data = species_traits)

wilcox.test(shift_max20_P.1.3_dist_km ~ widespread, data = species_traits, alternative = "less")# H0: widespread shifted less than range restricted
wilcox.test(shift_min20_P.1.3_dist_km ~ widespread, data = species_traits, alternative = "less")

wilcox.test(shift_max20_P.1.3_dist_km ~ PassNonPass, data = species_traits, alternative = "less")
wilcox.test(shift_min20_P.1.3_dist_km ~ PassNonPass, data = species_traits, alternative = "less")

wilcox.test(shift_max20_P.1.3_dist_km ~ migratory_binomial, data = species_traits, alternative = "less")
wilcox.test(shift_min20_P.1.3_dist_km ~ migratory_binomial, data = species_traits, alternative = "less")





################################################################################
### Figure S1 - Effect of distance to either boundary 

# color palette
cbp1 <- c("#D9717D","#4DB6D0") 


g1 = ggplot(subset(species_traits, Marine == 0), aes(y = shift_max20_P.1.3_dist_km, x = dist_N_km_max20_P.1)) + geom_point() +
  geom_smooth(method = 'lm', se=T, colour="#D9717D") +
  xlab("Distance between leading edge and the northern boundary") + ggtitle("2") + theme_minimal() + ylab('Leading edge shift (km)')

g2 = ggplot(subset(species_traits, Marine == 0), aes(y = shift_min20_P.1.3_dist_km, x = dist_S_km_min20_P.1)) + geom_point() +
  geom_smooth(method = 'lm', se=T, colour="#4DB6D0") +
  xlab("Distance between rear edge and the southern boundary") + ggtitle("1") + theme_minimal() + ylab('Rear edge shift (km)')

ggarrange(g1,g2)
ggsave2('DistToBoundaries.jpeg', scale = 1.5)



################################################################################
### regression plots

par(mfrow = c(2,2))
plot(data = subset(species_traits, distrib.core == 'north'), 
     shift_max20_P.1.3_dist_km ~ scale(trophic_position), pch = 19, xlab = 'scaled trophic position',
     ylab = 'leading edge shift (km)')
abline(a = tb$X.Intercept.[which(tb$model == 'northern' & tb$shift =='max' & tb$value == 'estimate_mean')], 
       b = tb$scale.trophic_position.[which(tb$model == 'northern' & tb$shift =='max')])

plot(data = subset(species_traits, distrib.core == 'north'), 
     shift_max20_P.1.3_dist_km ~ scale(pc1_lc), pch = 19, xlab = 'scaled land cover',
     ylab = 'leading edge shift (km)')
abline(a = tb$X.Intercept.[which(tb$model == 'northern' & tb$shift =='max' & tb$value == 'estimate_mean')], 
       b = tb$scale.pc1_lc.[which(tb$model == 'northern' & tb$shift =='max')])


plot(data = subset(species_traits, distrib.core == 'north'), 
     shift_max20_P.1.3_dist_km ~ scale(dist_N_km_max20_P.1), pch = 19, xlab = 'scaled distance to northern boundary',
     ylab = 'leading edge shift (km)')
abline(a = tb$X.Intercept.[which(tb$model == 'northern' & tb$shift =='max' & tb$value == 'estimate_mean')], 
       b = tb$scale.dist_N_km_max20_P.1.[which(tb$model == 'northern' & tb$shift =='max')])

plot(data = subset(species_traits, distrib.core == 'north'), 
     shift_max20_P.1.3_dist_km ~ scale(dist_S_km_min20_P.1), pch = 19, xlab = 'scaled distance to southern boundary',
     ylab = 'leading edge shift (km)')
abline(a = tb$X.Intercept.[which(tb$model == 'northern' & tb$shift =='max' & tb$value == 'estimate_mean')], 
       b = tb$scale.dist_S_km_min20_P.1.[which(tb$model == 'northern' & tb$shift =='max')])
