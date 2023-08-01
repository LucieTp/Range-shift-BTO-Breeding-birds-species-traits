################################################################################
### This code generates most* of the figures in the paper both for the 
### mainland and whole study area 
################################################################################

### this code is in support of Thompson et al. 2023 - 
### Joint effects of species traits and environmental preferences on range edge shifts of British birds
### published in Global Ecology and biogeography

### *Figures S2.2 and S2.3 are obtained in 1.BTO_breeding.atlas_range.shifts.Mainland.R

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

setwd("F:/TheseSwansea/TraitStudy/Github")

sf_UK  <- ne_states(country = 'United Kingdom', returnclass = "sf")
sf_UK <- subset(sf_UK, !(region %in% c('Northern Ireland')))

## Main table with all the species traits, range shift etc
# for the whole study area
species_traits <- read.csv('data/SpecTrait_Full_062023_159sp.csv', stringsAsFactors = FALSE, row.names = 1)
# for the mainland
species_traits.mainland <- read.csv('data/SpecTrait_Full_062023_156sp.csv', stringsAsFactors = FALSE, row.names = 1)

# Marine or non marine species : (Marine = >75% of points within 20km of the coastline)
propmarine <- read.csv('data/Speciestraits_ProportionMarineBTOsp_159sp.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits = merge(species_traits, propmarine[,c('speccode','prop_Marine20km')])
species_traits$Marine = ifelse(species_traits$prop_Marine20km>75, 1, 0)

# remove marine species
species_traits = subset(species_traits, Marine == 0)
# remove Pintail and Corvus cornix 
species_traits = subset(species_traits, !speccode %in% c(71,910))

# subset the mainland study too
species_traits.mainland = subset(species_traits.mainland, speccode %in% species_traits$speccode)

# we used the distribution core relative to middle of the mainland not the whole study to classify sp
# as Northern or Southern:
species_traits = merge(species_traits[,-which(names(species_traits) == 'distrib.core')], species_traits.mainland[,c('speccode','distrib.core')], by = 'speccode', all.x = T) # keep distribution core from the mainland study to compare results

species_traits$log10BodyMass.Value = log10(species_traits$BodyMass.Value)
species_traits$shift_diff = (species_traits$shift_max20_P.1.3_dist_km - species_traits$shift_min20_P.1.3_dist_km)

species_traits.mainland$log10BodyMass.Value = log10(species_traits.mainland$BodyMass.Value)
species_traits.mainland$shift_diff = (species_traits.mainland$shift_max20_P.1.3_dist_km - species_traits.mainland$shift_min20_P.1.3_dist_km)


################################################################################
## BTO Data (occurrence data per grid cell for all of UK) - used for distribution maps
BTO_grid <- read.csv("data/grid_square_coordinates_lookup.csv", header=T)
BTO_distrib <- read.csv("data/distributions.csv", header=T) # period, sp code, season and grid for GB and Ireland

# Mean coordinates of each grid cell
Coord10 <- as_tibble(BTO_grid) %>%
  dplyr::filter(resolution==10 & order<5) %>%
  dplyr::group_by(grid) %>%
  dplyr::summarise(long = mean(long), lat = mean(lat))

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
BTO_distrib = subset(BTO_distrib, speccode %in% species_traits$speccode)
BTO_distrib = subset(BTO_distrib, period %in% c("1968-72","2008-11"))

BTO_distrib.mainland = subset(BTO_distrib, lat < 58.7)



################################################################################
### Correlation between variables

library(corrplot)
species_traits$log10BodyMass.Value = log(species_traits$BodyMass.Value)
cov = species_traits[,c("log10BodyMass.Value", "pc1_lc" , "pc2_lc" , "pc1_env" , "pc2_env" , "habitat_gen" , "normalised_indegree" , 
                        "diet_diversity" , "nrecord_P.1" , "outdegree" , "trophic_position" , "HWI", 
                        'dist_N_km_max20_P.1','dist_S_km_min20_P.1')]

names(cov) = c('logged body mass', 'Forest and grassland', 'Urban and agricultural', 'Precipitation', 'Temperature', 'Habitat generality',
               'Normalised number of prey','Diet diversity','Range size','Number of prey','Trophic position','Hand-wing index','Northern boundary',
               'Southern boundary')

M <- cor(cov, use = 'complete.obs')


pdf(file = "plots/corrplot.pdf", width = 9, height = 9) 
corrplot(M, type="upper", order="hclust")
dev.off()

################################################################################
## SHIFT DIRECTIONALITY

### wilcoxon test for statistically significant poleward or southward shifts 
### between period 1 and 3

# for example:
# >> for each species we test for significant difference (smaller or larger) in 
# latitude between the 20 leading-edge cells in P1 and in P3 
# >> if the latitude is significantly higher in P3 than P1, then we conclude that 
# the species shifted significantly northward etc 

BTO_distrib.mainland$Edge = NA
stable.leading = stable.rear = poleward.leading = poleward.rear = southward.rear = southward.leading = NULL
for (i in unique(species_traits$speccode)){
  for (p in unique(BTO_distrib.mainland$periodN)){
    
    loc = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$periodN == p),]
    edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

    BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$periodN == p),'Edge'] = edge
  }
  
  ## leading edge
  
  test = wilcox.test(lat ~ periodN, data = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$Edge == 'Leading' & BTO_distrib.mainland$periodN %in% c('P.1','P.3') ),], alternative = 'less', paired = F)
  if(test$p.value < 0.05 ){
      poleward.leading = c(poleward.leading, i)
    } else {
        test = wilcox.test(lat ~ periodN, data = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$Edge == 'Leading' & BTO_distrib.mainland$periodN %in% c('P.1','P.3') ),], alternative = 'greater', paired = F)
        if(test$p.value < 0.05 ){
          southward.leading = c(southward.leading, i)
        } 
        else {
          stable.leading = c(stable.leading, i)
        }
    }
    
    ## rear edge
    test1 = wilcox.test(lat ~ periodN, data = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$Edge == 'Rear' & BTO_distrib.mainland$periodN %in% c('P.1','P.3') ),], alternative = 'less', paired = F)
    
    if(test1$p.value < 0.05 ){
      poleward.rear = c(poleward.rear, i)
    } else {
      test1 = wilcox.test(lat ~ periodN, data = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$Edge == 'Rear' & BTO_distrib.mainland$periodN %in% c('P.1','P.3') ),], alternative = 'greater', paired = F)
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
t.test(lat ~ periodN, data = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$Edge == 'Leading' & BTO_distrib.mainland$periodN %in% c('P.1','P.3') ),])

bartlett.test(lat ~ periodN, data = BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == i & BTO_distrib.mainland$Edge == 'Leading' & BTO_distrib.mainland$periodN %in% c('P.1','P.3') ),])
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

### are the shifts of species on the mainland significantly smaller than those for
### the whole study
wilcox.test(species_traits$shift_max20_P.1.3_dist_km, 
           species_traits.mainland$shift_max20_P.1.3_dist_km[match(species_traits$speccode,species_traits.mainland$speccode)], 
           paired = T)

hist(species_traits.mainland$shift_max20_P.1.3_dist_km[match(species_traits$speccode,species_traits.mainland$speccode)], col = 'blue')
hist(species_traits$shift_max20_P.1.3_dist_km, add = T) 


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

names = c('scientific_name','english_name',"Marine",'prop_Marine20km','lat_mean.P.1','distrib.core','IOCOrder','shift_max20_P.1.3_dist_km','shift_min20_P.1.3_dist_km','shift_diff',"dist_S_km_min20_P.1","dist_N_km_max20_P.1",'P.1','log10BodyMass.Value','outdegree',"pc1_env","pc2_env","pc1_lc" ,"pc2_lc", "normalised_indegree","diet_diversity","trophic_position","habitat_gen",'HWI')
towrite = species_traits[,names]
towrite = merge(towrite, species_traits.mainland[,c('scientific_name','shift_max20_P.1.3_dist_km','shift_min20_P.1.3_dist_km','shift_diff',"dist_S_km_min20_P.1","dist_N_km_max20_P.1",'P.1')], by = 'scientific_name')

names(towrite) = c('Scientific name','English name',"Marine species", 'Proportion of coastal grid cells','Mean latitude of grid cells in P1 (mainland)','Distribution core','Taxonomy','Shift in leading edge between P1 and P3 (km) - mainland','Shift in rear edge between P1 and P3 (km) - mainland','Range expansion (km) - mainland',"Southern boundary effect - mainland","Northern boundary effect - mainland",'Range size - mainland','Body Mass (g)','Number of predators',"Association with precipitation","Association with temperature","Association with forest and grassland" ,"Association with urban and cropland", "Normalised indegree","Diet diversity","Trophic position","Habitat generality",'Hand-wing index','Shift in leading edge between P1 and P3 (km) - whole study','Shift in rear edge between P1 and P3 (km) - whole study','Range expansion (km) - whole study',"Southern boundary effect - mainland","Northern boundary effect - whole study",'Range size - whole study')
towrite[,c(4:5,8:30)] = round(towrite[,c(4:5,8:30)], 2)
write.csv(towrite, 'SpeciesTraitsTableS1_paper.WholeStudy.Mainland.csv')


#################################################################################
## FIG 1  #######################################################################
#################################################################################
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)

################################################################################
### plots for the whole study area:

# NORTHERN SPECIES
loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 95),])
loc$Edge = NA
## create 'Edge' column that classifies the occurrence data into leading or rear edge
loc$Edge[which(loc$period == "1968-72")] = ifelse(loc$lat[which(loc$period == "1968-72")] %in% tail(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "1968-72")] %in% head(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Rear", "Middle"))

loc$Edge[which(loc$period == "2008-11")] = ifelse(loc$lat[which(loc$period == "2008-11")] %in% tail(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "2008-11")] %in% head(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Rear", "Middle"))

## create main plot without legend (we will create legend below to have seperate legends for each element)
Northern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 1.5, shape = 15) +
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed", size = 2) + geom_hline(yintercept = max(Loc10$lat), linetype="dashed", size = 2) + 
  scale_colour_manual(values = c("black", "gray57","#CC6677","#88CCEE")) +
  ggtitle("Northern species - Whole study") + ylab("") + xlab("") + theme_classic(base_size =  23) + xlim(min(loc$long) - 1, max(Coord10$long)) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  theme(legend.position = "none")

# SOUTHERN SPECIES
loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 272),])

loc$Edge = NA
loc$Edge[which(loc$period == "1968-72")] = ifelse(loc$lat[which(loc$period == "1968-72")] %in% tail(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "1968-72")] %in% head(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Rear", "Middle"))

loc$Edge[which(loc$period == "2008-11")] = ifelse(loc$lat[which(loc$period == "2008-11")] %in% tail(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "2008-11")] %in% head(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Rear", "Middle"))

Southern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 1.5, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  geom_hline(yintercept = max(Loc10$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed", size = 2) + 
  ggtitle("Southern species - Whole study ") + ylab("") + xlab("") + theme_classic(base_size =  23) + xlim(min(loc$long) - 1, max(Coord10$long)) + scale_colour_manual(values = c("black", "gray57","#CC6677","#88CCEE")) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  theme(legend.position = "none")

l1 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 2, shape = 15) +
  scale_color_manual(name = "Time period",values = c("black", "gray57","#CC6677","#88CCEE")) + 
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=22), legend.title=element_text(size=22), legend.key.size = unit(1.5, 'cm'))
leg1 <- get_legend(l1)

l2 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") +   
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  scale_color_manual(name = "Range edge",values = c("#CC6677","#88CCEE")) + 
  labs(linetype = 'Time period') +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=22), legend.title=element_text(size=22), legend.key.size = unit(1.5, 'cm'))
leg2 <- get_legend(l2)

library(patchwork)
blank_p <- plot_spacer() + theme_void()
leg12 <- plot_grid(leg1,leg2,blank_p,ncol = 1)


ggarrange(Southern.dist, leg12, Northern.dist, widths = c(1,0.05,1), heights = c(1,0.3,1), nrow = 1)
ggsave('plots/NS.DistChange.WholeStudy.2.0.jpeg', width = 12, height = 7, scale = 2, dpi = 800, bg= 'transparent')


################################################################################
### supplementary figures S2.1
loc66 = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 66),]) # 46, 66 for species in supplementary
loc66$Edge = NA
loc66$Edge[which(loc66$period == "1968-72")] = ifelse(loc66$lat[which(loc66$period == "1968-72")] %in% tail(sort(loc66$lat[which(loc66$period == "1968-72")]), 20), "Leading", 
                                                  ifelse(loc66$lat[which(loc66$period == "1968-72")] %in% head(sort(loc66$lat[which(loc66$period == "1968-72")]), 20), "Rear", "Middle"))
loc66$Edge[which(loc66$period == "2008-11")] = ifelse(loc66$lat[which(loc66$period == "2008-11")] %in% tail(sort(loc66$lat[which(loc66$period == "2008-11")]), 20), "Leading", 
                                                  ifelse(loc66$lat[which(loc66$period == "2008-11")] %in% head(sort(loc66$lat[which(loc66$period == "2008-11")]), 20), "Rear", "Middle"))

loc46 = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 46),]) # 46, 66 for species in supplementary
loc46$Edge = NA
loc46$Edge[which(loc46$period == "1968-72")] = ifelse(loc46$lat[which(loc46$period == "1968-72")] %in% tail(sort(loc46$lat[which(loc46$period == "1968-72")]), 20), "Leading", 
                                                      ifelse(loc46$lat[which(loc46$period == "1968-72")] %in% head(sort(loc46$lat[which(loc46$period == "1968-72")]), 20), "Rear", "Middle"))
loc46$Edge[which(loc46$period == "2008-11")] = ifelse(loc46$lat[which(loc46$period == "2008-11")] %in% tail(sort(loc46$lat[which(loc46$period == "2008-11")]), 20), "Leading", 
                                                      ifelse(loc46$lat[which(loc46$period == "2008-11")] %in% head(sort(loc46$lat[which(loc46$period == "2008-11")]), 20), "Rear", "Middle"))


dist_66 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc66, aes(x = long, y = lat, colour = period), size = 1.5, shape = 15) + 
  stat_ellipse(data = loc66[which(loc66$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  geom_hline(yintercept = max(Loc10$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed", size = 2) + 
  ggtitle("Gadwall - Mareca strepera") + ylab("") + xlab("") + theme_classic(base_size =  23) + xlim(min(loc66$long) - 1, max(Coord10$long)) + scale_colour_manual(values = c("black", "gray57","#CC6677","#88CCEE")) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  theme(legend.position = "none")

dist_46 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc46, aes(x = long, y = lat, colour = period), size = 1.5, shape = 15) + 
  stat_ellipse(data = loc46[which(loc46$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  geom_hline(yintercept = max(Loc10$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed", size = 2) + 
  ggtitle("Mute swan - Cygnus olor") + ylab("") + xlab("") + theme_classic(base_size =  23) + xlim(min(loc46$long) - 1, max(Coord10$long)) + scale_colour_manual(values = c("black", "gray57","#CC6677","#88CCEE")) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  theme(legend.position = "none")


l1 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 2, shape = 15) +
  scale_color_manual(name = "Time period",values = c("black", "gray57","#CC6677","#88CCEE")) + 
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=22), legend.title=element_text(size=22), legend.key.size = unit(1.5, 'cm'))
leg1 <- get_legend(l1)

l2 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") +   
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  scale_color_manual(name = "Range edge",values = c("#CC6677","#88CCEE")) + 
  labs(linetype = 'Time period') +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=22), legend.title=element_text(size=22), legend.key.size = unit(1.5, 'cm'))
leg2 <- get_legend(l2)

library(patchwork)
blank_p <- plot_spacer() + theme_void()
leg12 <- plot_grid(leg1,leg2,blank_p,ncol = 1)

ggarrange(dist_66, leg12, dist_46, widths = c(1,0.05,1), heights = c(1,0.3,1), nrow = 1)

ggsave('plots/NS.DistChange.Gadwall46_Swan66.pdf', scale = 4.5)


################################################################################
### for the mainland study

loc = merge(Loc10, BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == 95),])
loc$Edge[which(loc$period == "1968-72")] = ifelse(loc$lat[which(loc$period == "1968-72")] %in% tail(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "1968-72")] %in% head(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Rear", "Middle"))

loc$Edge[which(loc$period == "2008-11")] = ifelse(loc$lat[which(loc$period == "2008-11")] %in% tail(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "2008-11")] %in% head(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Rear", "Middle"))

Northern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 2, shape = 15) +
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  geom_hline(yintercept = min(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + geom_hline(yintercept = 58.7, linetype="dashed", size = 2) + scale_colour_manual(values = c("black", "gray57","#CC6677","#88CCEE")) +
  ggtitle("Northern species - Mainland") + ylab("") + xlab("") + theme_classic(base_size =  25) + xlim(min(loc$long) - 1, max(Coord10$long)) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  theme(legend.position = "none")

l1 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 2, shape = 15) +
  scale_color_manual(name = "Time period",values = c("black", "gray57","#CC6677","#88CCEE")) + 
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=20))
leg1 <- get_legend(l1)

l2 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") +   
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  scale_color_manual(name = "Range edge",values = c("#CC6677","#88CCEE")) + 
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=20)) 
leg2 <- get_legend(l2)

library(patchwork)
blank_p <- plot_spacer() + theme_void()
leg12 <- plot_grid(leg1,leg2,blank_p,ncol = 1)

final_p <- plot_grid(Northern.dist,leg12, rel_widths = c(1,0.1))
final_p 

# ggsave('plots/NorthernDist.Mainland.MainPlots.Dchange.jpeg', Northern.dist, scale = 2, dpi = 500, bg= 'transparent')

loc = merge(Loc10, BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == 272),])
loc$Edge[which(loc$period == "1968-72")] = ifelse(loc$lat[which(loc$period == "1968-72")] %in% tail(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "1968-72")] %in% head(sort(loc$lat[which(loc$period == "1968-72")]), 20), "Rear", "Middle"))

loc$Edge[which(loc$period == "2008-11")] = ifelse(loc$lat[which(loc$period == "2008-11")] %in% tail(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Leading", 
                                                  ifelse(loc$lat[which(loc$period == "2008-11")] %in% head(sort(loc$lat[which(loc$period == "2008-11")]), 20), "Rear", "Middle"))

Southern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 2, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  geom_hline(yintercept = max(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + 
  ggtitle("Southern species - Mainland") + ylab("") + xlab("") + theme_classic(base_size =  25) + xlim(min(loc$long) - 1, max(Coord10$long)) + scale_colour_manual(values = c("black", "gray57","#CC6677","#88CCEE")) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  theme(legend.position = "none")

l1 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat, colour = period), size = 2, shape = 15) +
  scale_color_manual(name = "Time period",values = c("black", "gray57","#CC6677","#88CCEE")) + 
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=20))
leg1 <- get_legend(l1)

l2 = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") +   
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge, linetype = period), size = 2) +
  scale_color_manual(name = "Range edge",values = c("#CC6677","#88CCEE")) + 
  theme(legend.direction = "vertical", legend.box = "vertical", legend.text=element_text(size=20)) 
leg2 <- get_legend(l2)

library(patchwork)
blank_p <- plot_spacer() + theme_void()
leg12 <- plot_grid(leg1,leg2,blank_p,ncol = 1)

final_p.south <- plot_grid(Southern.dist,leg12, rel_widths = c(1,0.1))
final_p.south

ggarrange(final_p.south, final_p, common.legend = T, widths = c(0.8,0.8))

ggsave('plots/NS.DistChange.Mainland.jpeg', width = 9, scale = 2.5, dpi = 500, bg= 'transparent')


################################################################################
### for generating small png plots with no background in Figure 3 and supp
################################################################################

##################
## Whole study:
loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 95 & BTO_distrib$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Northern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat), size = 2, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge), size = 2) + 
  geom_hline(yintercept = max(BTO_distrib$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed", size = 2) + 
  ggtitle("Northern species") + ylab("") + xlab("") + theme_classic(base_size =  33) + xlim(min(loc$long) - 1, max(Coord10$long)) +
  guides(colour = 'none') +
  theme(
    axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

ggsave('plots/NorthernDist.WholeStudy.png', Northern.dist, scale = 1, dpi = 300, bg= 'transparent')

loc = merge(Loc10, BTO_distrib[which(BTO_distrib$speccode == 272 & BTO_distrib$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Southern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat), size = 2, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge), size = 2) + 
  geom_hline(yintercept = max(BTO_distrib$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib$lat), linetype="dashed", size = 2) + 
  ggtitle("Southern species") + ylab("") + xlab("") + theme_classic(base_size =  33) + xlim(min(loc$long) - 1, max(loc$long)) +
  guides(colour = 'none') +
  theme(
    axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )


ggsave('plots/SouthernDist.WholeStudy.png', Southern.dist, scale = 1, dpi = 300, bg= 'transparent')

##################
## Mainland study

loc = merge(Loc10, BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == 95 & BTO_distrib.mainland$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Northern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat), size = 2, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge), size = 2) + 
  geom_hline(yintercept = max(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + 
  ggtitle("Northern species") + ylab("") + xlab("") + theme_classic(base_size =  33) + xlim(min(loc$long) - 1, max(Coord10$long)) +
  guides(colour = 'none') +
  theme(
    axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

ggsave('plots/NorthernDist.Mainland.png', Northern.dist, scale = 1, dpi = 300, bg= 'transparent')

loc = merge(Loc10, BTO_distrib.mainland[which(BTO_distrib.mainland$speccode == 272 & BTO_distrib.mainland$periodN == "P.1"),])
loc$Edge = ifelse(loc$lat %in% tail(sort(loc$lat), 20), "Leading", ifelse(loc$lat %in% head(sort(loc$lat), 20), "Rear", "Middle"))

Southern.dist = ggplot(data = sf_UK) + geom_sf(colour = "lightgrey") + geom_point(data = loc, aes(x = long, y = lat), size = 2, shape = 15) + 
  stat_ellipse(data = loc[which(loc$Edge %in% c("Leading","Rear")),], aes(x = long, y = lat, colour = Edge), size = 2) + 
  geom_hline(yintercept = max(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + geom_hline(yintercept = min(BTO_distrib.mainland$lat), linetype="dashed", size = 2) + 
  ggtitle("Southern species") + ylab("") + xlab("") + theme_classic(base_size =  33) + xlim(min(loc$long) - 1, max(loc$long)) +
  guides(colour = 'none') +
  theme(
    axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

ggsave('plots/SouthernDist.Mainland.png', Southern.dist, scale = 1, dpi = 300, bg= 'transparent')



################################################################################
### model results 
################################################################################

################################################################################
### Outputs of PGLMMs in a nicer format

tb = read.csv2("results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.WholeStudy.HWI.csv", sep =",", dec = ".")
tb_km = read.csv2("results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.WholeStudy.HWI.km.csv", sep =",", dec = ".")

tb.mainland = read.csv2("results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.Mainland_136sp.HWI.csv", sep =",", dec = ".")
tb_km.mainland = read.csv2("results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.Mainland.HWI.km.csv", sep =",", dec = ".")

## comparing effect sizes of sp traits on leading edge shift for the mainland and whole study

mean(abs(apply(subset(tb, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)))
sd(apply(subset(tb, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric))

mean(abs(apply(subset(tb.mainland, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)))
sd(abs(apply(subset(tb.mainland, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)))

hist(abs(apply(subset(tb, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)))
hist(abs(apply(subset(tb.mainland, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)))

wilcox.test(abs(apply(subset(tb, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)),
            abs(apply(subset(tb.mainland, value == 'estimate_mean' & shift == 'max')[c(9:11,13:22)], 2, as.numeric)), paired = T)

################################################################################
## TABLE S3         
################################################################################

tb1 = pivot_wider(tb.mainland[which(tb.mainland$value %in% c("significance","estimate_mean")),-1], values_from = names(tb.mainland)[8:22], names_from = value)
tb1$shift = ifelse(tb1$shift == "max", "leading edge shift",ifelse(tb1$shift == "min","rear edge shift", "Expansion/Retraction"))
tb1 = tb1[,-which(str_detect(colnames(tb1),"Loglik|phylo.mean"))]
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = apply(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))], 2, as.numeric)
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = round(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))],2)

tb1[,c(7:36)] = apply(tb1[,c(7:36)], 2, as.character)
tb1[is.na(tb1)] = ''
for(i in seq(7,36,2)) {tb1[,i] = apply(tb1[,c(i, i+1)], 1, function(x) paste(x[1], x[2]))}
tb1 = tb1[,-str_which(names(tb1), 'significance')]

stargazer(tb1, summary = F, type = "html")

################################################################################
## TABLE S3
## same with coef in km
################################################################################

tb1 = pivot_wider(tb_km[which(tb_km$value %in% c("significance","estimate_mean")),-1], values_from = names(tb_km)[8:22], names_from = value)
tb1$shift = ifelse(tb1$shift == "max", "leading edge shift",ifelse(tb1$shift == "min","rear edge shift", "Expansion/Retraction"))
tb1 = tb1[,-which(str_detect(colnames(tb1),"Loglik|phylo.mean"))]
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = apply(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))], 2, as.numeric)
tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))] = round(tb1[,c("r2.phylo","p.value.phylo", "R2.rr2", str_subset(names(tb1), 'estimate'))],2)

tb1[,c(7:36)] = apply(tb1[,c(7:36)], 2, as.character)
tb1[is.na(tb1)] = ''
for(i in seq(7,36,2)) {tb1[,i] = apply(tb1[,c(i, i+1)], 1, function(x) paste(x[1], x[2]))}
tb1 = tb1[,-str_which(names(tb1), 'significance')]

stargazer(tb1, summary = F, type = "html")

################################################################################
## TABLE S3
## table of results without Orkney and the Shetlands
################################################################################

tb2 = pivot_wider(tb.mainland[which(tb.mainland$value %in% c("significance","estimate_mean")),-1], values_from = names(tb.mainland)[8:22], names_from = value)
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

################################################################################
## TABLE S3
## table of results without Orkney and the Shetlands in km
################################################################################

tb2 = pivot_wider(tb_km.mainland[which(tb_km.mainland$value %in% c("significance","estimate_mean")),-1], values_from = names(tb_km.mainland)[8:22], names_from = value)
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

################################################################################
## FIGURE 2 and 3

#########################
### whole study
library(tidyverse)
tb[,8:22] = apply(tb[,8:22], 2, as.numeric)

tb2 = pivot_longer(tb, cols = 8:22, names_to = "covariate", values_to = "estimate")
tb2 = pivot_wider(tb2[,2:13], values_from = estimate, names_from = value)
tb2$shift = ifelse(tb2$shift == "max", "leading edge shift",ifelse(tb2$shift == "min","rear edge shift", "Expansion"))

tb2$significance = ifelse(tb2$p.value_mean<0.001,"***",ifelse(tb2$p.value_mean<0.01,"**",ifelse(tb2$p.value_mean<0.05,"*",ifelse(tb2$p.value_mean<0.1,".",""))))
tb2$size = ifelse(tb2$p.value_mean<0.05,0.6,0.5)
tb2$covariate = factor(tb2$covariate)

tb2$CI_lo = tb2$estimate_mean - 1.96*tb2$estimate_se
tb2$CI_up = tb2$estimate_mean + 1.96*tb2$estimate_se

tb2$SE_lo = tb2$estimate_mean - tb2$estimate_se
tb2$SE_up = tb2$estimate_mean + tb2$estimate_se

tb2$shift = as.factor(tb2$shift)
tb2$shift = relevel(tb2$shift, "rear edge shift")
levels(tb2$shift)

biogeo = tb2[which(tb2$covariate %in% c("scale.dist_N_km_max20_P.1.", "scale.dist_S_km_min20_P.1.", "scale.P.1.", "X.Intercept.")),]
biogeo = biogeo[which(biogeo$shift!='Expansion'),]
biogeo = biogeo[which(biogeo$covariate!='X.Intercept.'),]

traits = tb2[which(!tb2$covariate %in% c("scale.dist_N_km_max20_P.1.", "scale.dist_S_km_min20_P.1.", "scale.P.1.", "X.Intercept.")),]
traits_nodiff = traits[which(traits$shift!='Expansion'),]

############################
### mainland
library(tidyverse)
tb.mainland[,8:22] = apply(tb.mainland[,8:22], 2, as.numeric)

tb.mainland2 = pivot_longer(tb.mainland, cols = 8:22, names_to = "covariate", values_to = "estimate")
tb.mainland2 = pivot_wider(tb.mainland2[,2:13], values_from = estimate, names_from = value)
# tb.mainland2 = tb.mainland2[-which(str_detect(tb.mainland2$model, "Non passeriformes")),]
tb.mainland2$shift = ifelse(tb.mainland2$shift == "max", "leading edge shift",ifelse(tb.mainland2$shift == "min","rear edge shift", "Expansion"))

tb.mainland2$significance = ifelse(tb.mainland2$p.value_mean<0.001,"***",ifelse(tb.mainland2$p.value_mean<0.01,"**",ifelse(tb.mainland2$p.value_mean<0.05,"*",ifelse(tb.mainland2$p.value_mean<0.1,".",""))))
tb.mainland2$covariate = factor(tb.mainland2$covariate)

tb.mainland2$CI_lo = tb.mainland2$estimate_mean - 1.96*tb.mainland2$estimate_se
tb.mainland2$CI_up = tb.mainland2$estimate_mean + 1.96*tb.mainland2$estimate_se

tb.mainland2$SE_lo = tb.mainland2$estimate_mean - tb.mainland2$estimate_se
tb.mainland2$SE_up = tb.mainland2$estimate_mean + tb.mainland2$estimate_se

tb.mainland2$shift = as.factor(tb.mainland2$shift)
tb.mainland2$shift = relevel(tb.mainland2$shift, "rear edge shift")
levels(tb.mainland2$shift)

biogeo.mainland = tb.mainland2[which(tb.mainland2$covariate %in% c("scale.dist_N_km_max20_P.1.", "scale.dist_S_km_min20_P.1.", "scale.P.1.", "X.Intercept.")),]
biogeo.mainland = biogeo.mainland[which(biogeo.mainland$shift!='Expansion'),]
biogeo.mainland = biogeo.mainland[which(biogeo.mainland$covariate!='X.Intercept.'),]

traits.mainland = tb.mainland2[which(!tb.mainland2$covariate %in% c("scale.dist_N_km_max20_P.1.", "scale.dist_S_km_min20_P.1.", "scale.P.1.", "X.Intercept.")),]
traits_nodiff.mainland = traits.mainland[which(traits.mainland$shift!='Expansion'),]


##################################
### BIOGEOGRAPHICAL COVARIATES
### FIGURE 2

biogeo.mainland$shift = factor(biogeo.mainland$shift, levels = c('rear edge shift', 'leading edge shift'))
# levels(biogeo$shift) = c('rear edge shift', 'leading edge shift', 'expansion')
biogeo.mainland$linesize = ifelse(biogeo.mainland$significance %in% c(''), 1, 1.5)

labels.biogeo = c('','Northern boundary effect','','Southern boundary effect','','Range size','')

g1 = ggplot(data = biogeo.mainland[which(biogeo.mainland$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  xlab("Estimate")+  ylab('') + ggtitle('a. Northern species - mainland') + xlim(min(biogeo.mainland$CI_lo), max(biogeo.mainland$CI_up)) +
  scale_y_discrete(labels = labels.biogeo) + 
  scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  theme(axis.ticks = element_blank()) + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash') + guides(size = 'none', shape = 'none')

img = png::readPNG("plots/NorthernDist.Mainland (2).png")
img = grid::rasterGrob(img, interpolate=TRUE)

g1 = g1 + annotation_custom(img, xmin=0, xmax=100, ymin=-2, ymax=8) +
  geom_point()

g2 = ggplot(data = biogeo.mainland[which(biogeo.mainland$model == "southern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  xlab("Estimate") + ylab('') +
  geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  scale_y_discrete(labels = labels.biogeo) + 
  ggtitle('b. Southern species - mainland')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  xlim(min(biogeo.mainland$CI_lo), max(biogeo.mainland$CI_up)) +
  theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none', shape = 'none')

img = png::readPNG("plots/SouthernDist.Mainland (2).png")
img = grid::rasterGrob(img, interpolate=TRUE)

g2 = g2 + annotation_custom(img, xmin=-140, xmax=-20, ymin=-3.5, ymax=8) +
  geom_point()

g3 = ggplot(data = biogeo.mainland[which(biogeo.mainland$model == "Passeriformes"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  xlab("Estimate") +  ylab('') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  scale_y_discrete(labels = labels.biogeo) + 
  ggtitle('c. Passeriformes - mainland')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  xlim(min(biogeo.mainland$CI_lo), max(biogeo.mainland$CI_up)) +
  theme(axis.ticks = element_blank(),axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none', shape = 'none')

g4 = ggplot(data = biogeo.mainland[which(biogeo.mainland$model == "terrestrial"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  xlab("Estimate") +  ylab('') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  scale_y_discrete(labels = labels.biogeo) + 
  ggtitle('d. All species - mainland')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  xlim(min(biogeo.mainland$CI_lo), max(biogeo.mainland$CI_up)) +
  theme(axis.ticks = element_blank(),axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none', shape = 'none')

library(ggpubr)
## Figure S2 - Estimates from GLMM 
ggarrange(g1,g3,g2,g4, common.legend = T, widths = c(1.4,1))

library(cowplot)
ggsave2('plots/ModelEstimates_Biogeo.Mainland.pdf', scale = 4.5)


##################################
### BIOGEOGRAPHICAL COVARIATES - whole study region
### FIGURE S2.9

biogeo$shift = factor(biogeo$shift, levels = c('rear edge shift', 'leading edge shift'))
# levels(biogeo$shift) = c('rear edge shift', 'leading edge shift', 'expansion')
biogeo$linesize = ifelse(biogeo$significance %in% c(''), 1, 1.5)

labels.biogeo = c('','Northern boundary effect','','Southern boundary effect','','Range size','')

g1 = ggplot(data = biogeo[which(biogeo$model == "northern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  xlab("Estimate")+  ylab('') + ggtitle('a. Northern species- whole study') + xlim(min(biogeo$CI_lo), max(biogeo$CI_up)) +
  scale_y_discrete(labels = labels.biogeo) + 
  scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  theme(axis.ticks = element_blank()) + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash') + guides(size = 'none', shape = 'none')

g2 = ggplot(data = biogeo[which(biogeo$model == "southern"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  xlab("Estimate") + ylab('') +
  geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  scale_y_discrete(labels = labels.biogeo) + 
  ggtitle('b. Southern species- whole study')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  xlim(min(biogeo$CI_lo), max(biogeo$CI_up)) +
  theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none', shape = 'none')

g3 = ggplot(data = biogeo[which(biogeo$model == "Passeriformes"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  xlab("Estimate") +  ylab('') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  scale_y_discrete(labels = labels.biogeo) + 
  ggtitle('c. Passeriformes- whole study')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  xlim(min(biogeo$CI_lo), max(biogeo$CI_up)) +
  theme(axis.ticks = element_blank(),axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none', shape = 'none')

g4 = ggplot(data = biogeo[which(biogeo$model == "terrestrial"),], aes(x = estimate_mean, y = paste(covariate,shift), colour = shift, shape = model, label = significance)) + 
  geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + 
  xlab("Estimate") +  ylab('') + geom_errorbarh(aes(y = paste(covariate,shift), xmin = CI_lo, xmax = CI_up, size = linesize), height = .3) + 
  scale_y_discrete(labels = labels.biogeo) + 
  ggtitle('d. All species- whole study')+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 25) + 
  xlim(min(biogeo$CI_lo), max(biogeo$CI_up)) +
  theme(axis.ticks = element_blank(),axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,5.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none', shape = 'none')

library(ggpubr)
## Figure S2 - Estimates from GLMM 
ggarrange(g1,g3,g2,g4, common.legend = T, widths = c(1.4,1))

library(cowplot)
ggsave2('plots/ModelEstimates_biogeo-wholestudy.pdf', scale = 4.5)


################################################################################
### SPECIES TRAITS - no expansion
### FIGURE 3

labels = c('','Diet diversity','','Habitat generality','','Hand wing index','','Log 10 body mass','','Normalised number of prey','','Number of predators','','Association with precipitation','','Association with forest/grassland','','Association with temperature','','Association with urban area/cropland','','Trophic position','')


### Whole study area
library(ggrepel)
traits_nodiff$shift = relevel(traits_nodiff$shift, "rear edge shift")
traits_nodiff$linesize = ifelse(traits_nodiff$significance %in% c(''), 1, 1.5)

traits_nodiff$CovShift = paste(traits_nodiff$covariate, traits_nodiff$shift)

# pc1 lc: forest and grassland
# pc2 lc: cropland and urban area
# pc2 env: temperature + winter precipitation
# pc1 env: summer precipiation
g1 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "northern"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Whole study') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels = labels) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none') 
g2 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "southern"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Whole study') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels = labels) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none') 
g3 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "Passeriformes"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Whole study') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels =labels) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none')
g4 = ggplot(data = traits_nodiff[which(traits_nodiff$model == "terrestrial"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Whole study') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits[,'CI_lo']),max(traits[,'CI_up'])) + scale_y_discrete(labels =labels) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash')+ guides(size = 'none')

ggarrange(g1,Northern.dist,g3, g2,Southern.dist,g4, common.legend = T, labels = c('a','','c','b','','d'), widths = c(1.5,1,1.5), font.label = list(size = 20))


### Mainland
traits_nodiff.mainland$shift = relevel(traits_nodiff.mainland$shift, "rear edge shift")
traits_nodiff.mainland$linesize = ifelse(traits_nodiff.mainland$significance %in% c(''), 1, 1.5)

## get the same order as in the whole study plots to compare
traits_nodiff.mainland$CovShift = paste(traits_nodiff$covariate, traits_nodiff$shift)
traits_nodiff.mainland$CovShift = factor(traits_nodiff$CovShift, levels = levels(factor(traits_nodiff$CovShift)))

g1.m = ggplot(data = traits_nodiff.mainland[which(traits_nodiff.mainland$model == "northern"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Mainland') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits_nodiff.mainland[,'CI_lo']),max(traits_nodiff.mainland[,'CI_up']))+ scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank(), axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + guides(size = 'none')
g2.m = ggplot(data = traits_nodiff.mainland[which(traits_nodiff.mainland$model == "southern"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Mainland')+ geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits_nodiff.mainland[,'CI_lo']),max(traits_nodiff.mainland[,'CI_up'])) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank(), axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + guides(size = 'none')
g3.m = ggplot(data = traits_nodiff.mainland[which(traits_nodiff.mainland$model == "Passeriformes"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Mainland') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up, size = linesize), height = .1, alpha = .8) + xlim(min(traits_nodiff.mainland[,'CI_lo']),max(traits_nodiff.mainland[,'CI_up'])) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank(), axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + guides(size = 'none')
g4.m = ggplot(data = traits_nodiff.mainland[which(traits_nodiff.mainland$model == "terrestrial"),], aes(x = estimate_mean, y = CovShift, colour = shift, label = significance)) + geom_point(size = 3) + geom_vline(xintercept = 0) + geom_text(hjust=0.5, vjust=0.2,show.legend=FALSE, size = 10) + xlab("Estimate") + ylab("") + ggtitle('Mainland') + geom_errorbarh(aes(y = CovShift, xmin = CI_lo, xmax = CI_up), size = 1, height = .1, alpha = .8) + xlim(min(traits_nodiff.mainland[,'CI_lo']),max(traits_nodiff.mainland[,'CI_up'])) + scale_colour_manual(values= safe_colorblind_palette) + theme_bw(base_size = 20) + theme(axis.ticks = element_blank(), axis.text.y = element_blank())  + geom_hline(yintercept = seq(2.5,20.5,2), colour = 'grey', linetype = 'longdash') + guides(size = 'none')


### put both mainland and whole study area together into one plot
fig = ggarrange(g1,g1.m,g3, g3.m, g2, g2.m, g4, g4.m, common.legend = T,
          ncol = 4, nrow = 2, widths = c(1.3,0.6,1.3,0.6), 
          labels = c("a. Northern species", "", "c. Passeriformes", "", "b. Southern species", "", "d. All species"),
          font.label = list(size = 20))


## Figure 3:
pdf("plots/ModelEstimates.nodiff.CI.WholeStudy.Mainland.pdf", width = 20, height = 12.5)
fig

## add in insets for northern and southern species
img = png::readPNG("plots/NorthernDist.Mainland (2).png")
img = grid::rasterGrob(img, interpolate=TRUE)
grid::grid.draw(
  grid::grobTree(img,
                 vp = grid::viewport(x = unit(0.05, "npc"), 
                                     y = unit(0.55, "npc"), 
                                     width = unit(0.20, "npc"), 
                                     height = unit(0.30, "npc"))))

img = png::readPNG("plots/SouthernDist.Mainland (2).png")
img = grid::rasterGrob(img, interpolate=TRUE)

grid::grid.draw(
  grid::grobTree(img,
                 vp = grid::viewport(x = unit(0.05, "npc"), 
                                     y = unit(0.1, "npc"), 
                                     width = unit(0.20, "npc"), 
                                     height = unit(0.30, "npc"))))
dev.off()

## Figure 3:
ggsave('plots/ModelEstimates.nodiff.CI.WholeStudy.Mainland.jpeg', scale = 4)



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



