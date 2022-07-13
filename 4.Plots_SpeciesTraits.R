################################################################################
### CODE FOR PLOTTING THE MAIN FIGURE, TABLE AND SUPP FIG 1 AND 2
################################################################################


################################################################################
## IMPORTING THE DATA

library(ggplot2)
library(ggpubr) # for ggarrange
library(cowplot)
library(rnaturalearth)
sf_UK  <- ne_countries(scale = "medium", country = 'United Kingdom', returnclass = "sf")

## Main table with all the species traits, range shift etc
setwd("E:/TheseSwansea/TraitStudy/code_Miguel")
species_traits <- read.csv('E:/TheseSwansea/TraitStudy/code_Miguel/SpecTrait_04022022_159sp.csv', stringsAsFactors = FALSE, row.names = 1)

# Marine or non marine species : (Marine = >75% of points within 20km of the coastline)
propmarine <- read.csv('E:/TheseSwansea/TraitStudy/code_Konstans/Speciestraits_ProportionMarineBTOsp_159sp.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits = merge(species_traits, propmarine[,c('speccode','prop_Marine20km')])
species_traits$Marine = ifelse(species_traits$prop_Marine20km>75, 1, 0)

## adding range shift in km 
df.rangeshift <- read.csv('E:/TheseSwansea/TraitStudy/code_Konstans/df_rangeshift_012022_km_160sp.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits = merge(species_traits, df.rangeshift[,c("speccode","shift_max20_P.1.3_dist_km", "shift_min20_P.1.3_dist_km","dist_N_km_max20_P.1","dist_S_km_min20_P.1")], all.x = T)

## BTO Data (occurrence data per grid cell for all of UK) - used for distribution maps
BTO_grid <- read.csv("E:/TheseSwansea/TraitStudy/code_Konstans/grid_square_coordinates_lookup.csv", header=T)
BTO_distrib <- read.csv("E:/TheseSwansea/TraitStudy/code_Konstans/distributions.csv", header=T) # period, sp code, season and grid for GB and Ireland

library(dplyr)
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

### labelling the leading and rear edge coordinates for mapping
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


t.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Leading' & BTO_distrib$periodN %in% c('P.1','P.3') ),])

bartlett.test(lat ~ periodN, data = BTO_distrib[which(BTO_distrib$speccode == i & BTO_distrib$Edge == 'Leading' & BTO_distrib$periodN %in% c('P.1','P.3') ),])
shapiro.test(x$shift[which(x$edge == 'min')])




################################################################################

species_traits$PassNonPass = ifelse(species_traits$IOCOrder == "Passeriformes", "Passeriformes","non-Passeriformes")
species_traits$shift_diff = (species_traits$shift_max20_P.1.3_dist_km - species_traits$shift_min20_P.1.3_dist_km)
species_traits$migratory_binomial = factor(species_traits$migratory_binomial)
levels(species_traits$migratory_binomial) = c('Resident','Migrant')


library(tidyr)
shift_long = pivot_longer(data = species_traits[which(species_traits$Marine == 0),c("speccode","english_name", "distrib.core", "PassNonPass", "IOCOrder","widespread", "BodyMass.Value","migratory_binomial","shift_max20_P.1.3_dist_km","shift_min20_P.1.3_dist_km","shift_diff")], cols = c("shift_max20_P.1.3_dist_km","shift_min20_P.1.3_dist_km","shift_diff" ), names_to = "Period",values_to = "Shift")

## subsetting only non coastal species
species_traits = species_traits[which(species_traits$Marine == 0),]


################################################################################
### Table to export for supp material 


names = c('scientific_name','english_name','P.1','lat_mean.P.1','distrib.core','IOCOrder','BodyMass.Value','shift_max20_P.1.3_dist_km','shift_min20_P.1.3_dist_km','shift_diff','outdegree','prop_Marine20km',"pc1_env","pc2_env","pc1_lc" ,"pc2_lc", "normalised_indegree","diet_diversity","trophic_position","habitat_gen")
towrite = species_traits[,names]
names(towrite) = c('Scientific name','English name','Number of grid cells in P1','Mean latitude of grid cells in P1','Distribution core','Taxonomy','Body Mass (g)','Shift in leading edge between P1 and P2 (km)','Shift in rear edge between P1 and P2 (km)','Range expansion (km)','Outdegree','Proportion of coastal grid cells',"Association with precipitation seasonality","Association with temperature seasonality","Association with forest and grassland" ,"Association with urban and cropland", "Normalised indegree","Diet diversity","Trophic position","Habitat generality")
write.csv(towrite, 'SpeciesTraitsTable_paper.csv')

################################################################################
### Table of results

library(stargazer)
library(stringr)
library(tidyr)

tb = read.csv2("E:/TheseSwansea/TraitStudy/code_Miguel/pglmm_scaled_terrestrialNS_PCs.csv", sep =",", dec = ".")

# tb[,7:21] = apply(tb[,7:21], 2, as.numeric)
tb[,c("p.value.phylo", "r2", str_subset(names(tb), 'estimate'))] = round(tb[,c("p.value.phylo", "r2", str_subset(names(tb), 'estimate'))],2)

#################################################################################
## TABLE 1 ######################################################################
#################################################################################

tb1 = pivot_wider(tb[which(tb$value %in% c("significance","estimate_mean")),-1], values_from = names(tb)[7:21], names_from = value)
tb1 = tb1[-which(str_detect(tb1$model, "Non passeriformes")),]
tb1$shift = ifelse(tb1$shift == "max", "leading edge shift",ifelse(tb1$shift == "min","rear edge shift", "Expansion/Retraction"))
tb1 = tb1[,-which(str_detect(colnames(tb1),"Loglik|r2phylo"))]
tb1[,c("p.value.phylo", "r2", str_subset(names(tb1), 'estimate'))] = apply(tb1[,c("p.value.phylo", "r2", str_subset(names(tb1), 'estimate'))], 2, as.numeric)
tb1[,c("p.value.phylo", "r2", str_subset(names(tb1), 'estimate'))] = round(tb1[,c("p.value.phylo", "r2", str_subset(names(tb1), 'estimate'))],2)


test = tb1
test[,6:35] = apply(test[,6:35], 2, as.character)
test[is.na(test)] = ''
for(i in seq(6,35,2)) {test[,i] = apply(test[,c(i, i+1)], 1, function(x) paste(x[1], x[2]))}
test = test[,-str_which(names(test), 'significance')]
names(test) = c("model" , "Nb of species","Shift", "p-value phylo signal", "r2", "Intercept ",'Migratory status', 'Diet diversity',"Distance from N tip",
                "Distance from S tip ", "Habitat generality ","log10 Body Mass ", 'Normalised indegree ', "Vulnerability ",
                "Regional coverage " , "Precipitation seasonality ", "Grassland and forest " , "Temperature seasonality ", 'Urban and cropland ',
                "Trophic position " )


stargazer(test, summary = F, type = "html")


#################################################################################
## FIG 1 and Fig S2 #############################################################
#################################################################################

loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 95 & BTO_distrib$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Northern.dist = ggplot(data = sf_UK) + geom_sf() + geom_point(data = loc, aes(x = long, y = lat), size = 1) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge)) + 
  geom_hline(yintercept = max(BTO_distrib$lat), linetype="dashed") + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed") + 
  ggtitle("Northern species") + ylab("") + xlab("") + theme_classic(base_size =  11)


loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 272 & BTO_distrib$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Southern.dist = ggplot(data = sf_UK) + geom_sf() + geom_point(data = loc, aes(x = long, y = lat), size = 1) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge)) + 
  geom_hline(yintercept = max(BTO_distrib$lat), linetype="dashed") + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed") + 
  ggtitle("Southern species") + ylab("") + xlab("") + theme_classic(base_size =  11)

### model results

tb2 = pivot_longer(tb[which(!tb$value %in% c("estimate_sd")),], cols = 7:21, names_to = "covariate", values_to = "estimate")
tb2 = pivot_wider(tb2[,2:12], values_from = estimate, names_from = value)
tb2 = tb2[-which(str_detect(tb2$model, "Non passeriformes")),]
tb2$shift = ifelse(tb2$shift == "max", "leading edge shift",ifelse(tb2$shift == "min","rear edge shift", "Expansion/Retraction"))

tb2$significance = ifelse(tb2$p.value_mean<0.001,"***",ifelse(tb2$p.value_mean<0.01,"**",ifelse(tb2$p.value_mean<0.05,"*",ifelse(tb2$p.value_mean<0.1,".",""))))
tb2$size = ifelse(tb2$p.value_mean<0.05,0.6,0.5)
tb2$covariate = factor(tb2$covariate)
levels(tb2$covariate) = c("Migrant status - Migrant","Diet diversity - Elton", "Northern boundary effect", "Southern boundary effect","Habitat generality","Body mass", "Diet diversity - TETRA EU", "Vulnerability", "Regional coverage", "Precipitation seasonality", "Association with forest/grasslands","Temperature seasonality","Association with urban areas/croplands","Trophic position","Intercept")

biogeo = tb2[which(tb2$covariate %in% c("Northern boundary effect", "Southern boundary effect", "Regional coverage", "Intercept")),]
traits = tb2[which(!tb2$covariate %in% c("Northern boundary effect", "Southern boundary effect", "Regional coverage",'Diet diversity - TETRA EU')),]

### BIOGEOGRAPHICAL COVARIATES
g1 = ggplot(data = biogeo[which(biogeo$model == "northern"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE) + xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))
g2 = ggplot(data = biogeo[which(biogeo$model == "southern"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))
g3 = ggplot(data = biogeo[which(biogeo$model == "Passeriformes"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))
g4 = ggplot(data = biogeo[which(biogeo$model == "terrestrial"),], aes(x = estimate_mean, y = covariate, colour = shift, shape = model, label = significance, alpha = 1/p.value_mean)) + geom_point() + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(biogeo$estimate_mean), max(biogeo$estimate_mean))

library(ggpubr)
## Figure S2 - Estimates from GLMM 
ggarrange(g1,g2,g3,g4)

library(cowplot)
setwd("E:/TheseSwansea/TraitStudy/code_Miguel")
ggsave2('ModelEstimates_Biogeo.jpeg', scale = 3.5)

### SPECIES TRAITS
g1 = ggplot(data = traits[which(traits$model == "northern"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE) + xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('Northern species') + theme_classic(base_size =  15) + theme(legend.position = 'bottom')
g2 = ggplot(data = traits[which(traits$model == "southern"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('Southern species')+ theme_classic(base_size =  15) + theme(legend.position = 'bottom')
g3 = ggplot(data = traits[which(traits$model == "Passeriformes"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('Passeriformes')+ theme_classic(base_size =  15) + theme(legend.position = 'bottom')
g4 = ggplot(data = traits[which(traits$model == "terrestrial"),], aes(x = estimate_mean, y = covariate, colour = shift, label = significance, alpha = 1/p.value_mean)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0, vjust=0,show.legend=FALSE) + scale_alpha(range = c(0.5, 2), guide=FALSE)+ xlab("Estimate") + ylab("") + xlim(min(traits$estimate_mean), max(traits$estimate_mean)) + ggtitle('All species')+ theme_classic(base_size =  15) + theme(legend.position = 'bottom')

ggarrange(g1,Northern.dist,g3, g2,Southern.dist,g4, common.legend = T, labels = c('A','','C','B','','D'))

setwd("E:/TheseSwansea/TraitStudy/code_Miguel")
## Figure 1 - Estimates from GLMM 
ggsave2('ModelEstimates_Traits_maps.jpeg', scale = 2.2, dpi = 800)




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



