################################################################################
### This code runs the PGLMMs for each subset of species traits 
### (North/South/Passeriformes)
###############################################################################


rm(list = ls())


setwd("E:/TheseSwansea/TraitStudy/code_Miguel")
species_traits <- read.csv('SpecTrait_04022022_159sp.csv', stringsAsFactors = FALSE, row.names = 1)

# Marine or non marine species : (Marine = >75% of points within 20km of the coastline)
propmarine <- read.csv('E:/TheseSwansea/TraitStudy/code_Konstans/Speciestraits_ProportionMarineBTOsp_159sp.csv', stringsAsFactors = FALSE, row.names = 1)

species_traits = merge(species_traits, propmarine[,c('speccode','prop_Marine20km')])
species_traits$Marine = ifelse(species_traits$prop_Marine20km>75, 1, 0)

df.rangeshift <- read.csv('E:/TheseSwansea/TraitStudy/code_Konstans/df_rangeshift_012022_km_160sp.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits = merge(species_traits, df.rangeshift[,c("speccode","shift_max20_P.1.3_dist_km", "shift_min20_P.1.3_dist_km","dist_N_km_max20_P.1","dist_S_km_min20_P.1")], all.x = T)

species_traits$log10BodyMass.Value = log10(species_traits$BodyMass.Value)
species_traits$shift_diff = (species_traits$shift_max20_P.1.3_dist_km - species_traits$shift_min20_P.1.3_dist_km)
summary(species_traits$shift_diff)
# species with very negative values are species that retracted their range
# species with very positive values are species that expanded their ranges 
species_traits$PassNonPass = ifelse(species_traits$IOCOrder == "Passeriformes", "Passeriformes","non-Passeriformes")
species_traits$migratory_binomial = factor(species_traits$migratory_binomial)
levels(species_traits$migratory_binomial) = c('Resident','Migrant')

## Keep only non coastal species
species_traits = species_traits[which(species_traits$Marine == 0),]

library(doBy)
summaryBy(data = species_traits, shift_diff~distrib.core)


###############################################################################

library(phyr)
require(phytools)
require(ape)

tree = ape::read.tree("./BigBirdPhylogeny/BigBird.All.NewNames.6714Taxa.tre/BigBird.All.NewNames.6714Taxa.tre")
# tree <- consensus.edges(tree[1:10])
tree <- tree[[1]]

parseRowNames <- function(x){
  splitted <- strsplit(x, '_')[[1]]
  l <- length(splitted)
  return(paste0(splitted[l-1], '_', splitted[l]))
}

tree$tip.label <- unlist(lapply(tree$tip.label, parseRowNames))
#species_traits$species[!species_traits$species%in%tree$tip.label] = "Acanthis_flammea"
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, species_traits$species))

################################################################################
### LM checking residuals

mod = pglmm(sqrt(shift_min20_P.1.3_dist_km) ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
              (1|species__), data = species_traits, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
mod1 = pglmm(shift_max20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
               (1|species__), data = species_traits, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)

library(DHARMa)
plot(residuals(mod)~fitted(mod))
ggplot(data.frame(residuals = residuals(mod), fitted = fitted(mod)), aes(x = residuals, y = fitted)) + geom_point() + geom_smooth(method = 'lm')
ggplot(data.frame(residuals = residuals(mod1), fitted = fitted(mod1)), aes(x = residuals, y = fitted)) + geom_point() + geom_smooth(method = 'lm')

par(mfrow = c(1,1))
qqnorm(residuals(mod), pch = 1)
qqline(residuals(mod), col = "steelblue", lwd = 2)

qqnorm(residuals(mod1), pch = 1, frame = FALSE)
qqline(residuals(mod1), col = "steelblue", lwd = 2)

hist(species_traits$shift_max20_P.1.3_dist_km)
hist(species_traits$shift_min20_P.1.3_dist_km)

shapiro.test(residuals(mod))
shapiro.test(residuals(mod1))


################################################################################
## TRYING GAM MODELS >> too few points, overfitting, don't know how to include phylo random effect
library(mgcv)
## ALL SPECIES
mod = gam(shift_min20_P.1.3_dist_km ~ s(scale(log10(BodyMass.Value))) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + s(scale(pc1_env)) + s(scale(pc2_env)) + s(scale(pc1_lc)) + s(scale(pc2_lc)) + s(scale(P.1)), data = species_traits)
summary(mod)

mod = gam(shift_max20_P.1.3_dist_km ~ s(scale(log10(BodyMass.Value))) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + s(scale(pc1_env)) + s(scale(pc2_env)) + s(scale(pc1_lc)) + s(scale(pc2_lc)) + s(scale(P.1)), data = species_traits)
summary(mod)

#### SOUTH
mod = gam(shift_min20_P.1.3_dist_km ~ s(scale(log10(BodyMass.Value))) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + s(scale(P.1)), data = subset(species_traits, distrib.core == 'south'))
summary(mod)

mod = gam(shift_max20_P.1.3_dist_km ~ s(scale(log10(BodyMass.Value))) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + s(scale(P.1)), data = subset(species_traits, distrib.core == 'south'))
summary(mod)

#### NORTH
mod = gam(shift_min20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + scale(trophic_position) + s(scale(outdegree)) + scale(pc1_env), data = subset(species_traits, distrib.core == 'north'))
summary(mod)

mod = gam(shift_max20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + scale(trophic_position) + s(scale(outdegree)) + scale(pc1_env), data = subset(species_traits, distrib.core == 'north'))
summary(mod)




################################################################################################################


### FUNCTION FOR RUNNING ALL MODELS IN ONE LOOP
mod_loop = function(dta_list, mod, loop = 1, migrant = FALSE){
  
  library(car)
  model = 0
  all_dta = data.frame() # dataframe to store the data

  for (dta in dta_list){ # list of subsets to run the models for
    
    model = model + 1
    
    print(model)
    
    # variables to save
    estimate.max = NULL; p.value.max = NULL
    estimate.min = NULL; p.value.min = NULL
    estimate.diff = NULL; p.value.diff = NULL
    r2phylo.min = NULL; r2phylo.max = NULL; r2phylo.diff = NULL
    r2sp.min = NULL; r2sp.max = NULL; r2sp.diff = NULL
    Loglik.max = NULL;Loglik.min = NULL;Loglik.diff = NULL
    Loglik.max.lm = NULL; Loglik.min.lm = NULL; Loglik.diff.lm = NULL
    
    
    for (i in 1:loop){
      
      # shuffle the data
      dta1 = dta[sample(nrow(dta),nrow(dta)),]
      
      # model with migratory status as a covariate
      if (migrant == F){ 
        # with pglmm
        mod.min_shift = pglmm(shift_min20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                                (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
        # with lm
        mod.min_shift_lm <- lm(shift_min20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                                 , data = dta1)
        # with pglmm
        mod.max_shift <- pglmm(shift_max20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                                 (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
        
        # with lm
        mod.max_shift_lm <- lm(shift_max20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                                 , data = dta1)
        # with pglmm
        mod.diff_shift = pglmm(shift_diff ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                                (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
        # with lm
        mod.diff_shift_lm <- lm(shift_diff ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                               , data = dta1)
        
        # model without migratory status as a covariate
      }else{
        # with pglmm
        mod.min_shift = pglmm(shift_min20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                                (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
        # with lm
        mod.min_shift_lm <- lm(shift_min20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                                 , data = dta1)
        # with pglmm
        mod.max_shift <- pglmm(shift_max20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                                 (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
        # with lm
        mod.max_shift_lm <- lm(shift_max20_P.1.3_dist_km ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                                 , data = dta1)
        # with pglmm
        mod.diff_shift = pglmm(shift_diff ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                                (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
        # with lm
        mod.diff_shift_lm <- lm(shift_diff ~ scale(log10(BodyMass.Value)) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                               , data = dta1)
        
      }
      # save the outputs
      estimate.min = rbind(estimate.min, mod.min_shift$B)
      p.value.min = rbind(p.value.min, mod.min_shift$B.pvalue)
      r2sp.min = rbind(r2sp.min, mod.min_shift$s2r[[1]])
      r2phylo.min = rbind(r2phylo.min, mod.min_shift$s2r[[2]])
      Loglik.min = rbind(Loglik.min, mod.min_shift$logLik)
      Loglik.min.lm = rbind(Loglik.min.lm, logLik(mod.min_shift_lm))
      
      estimate.max = rbind(estimate.max, mod.max_shift$B)
      p.value.max = rbind(p.value.max, mod.max_shift$B.pvalue)
      r2sp.max = rbind(r2sp.max, mod.max_shift$s2r[[1]])
      r2phylo.max = rbind(r2phylo.max, mod.max_shift$s2r[[2]])
      Loglik.max = rbind(Loglik.max, mod.max_shift$logLik)
      Loglik.max.lm = rbind(Loglik.max.lm, logLik(mod.max_shift_lm))
      
      # expansion
      estimate.diff = rbind(estimate.diff, mod.diff_shift$B)
      p.value.diff = rbind(p.value.diff, mod.diff_shift$B.pvalue)
      r2sp.diff = rbind(r2sp.diff, mod.diff_shift$s2r[[1]])
      r2phylo.diff = rbind(r2phylo.diff, mod.diff_shift$s2r[[2]])
      Loglik.diff = rbind(Loglik.diff, mod.diff_shift$logLik)
      Loglik.diff.lm = rbind(Loglik.diff.lm, logLik(mod.diff_shift_lm))
      
    }
    print(mod[model])
    
    # summary
    test = data.frame(estimate.max = estimate.max, estimate.min = estimate.min, p.value.min = p.value.min, 
                      p.value.max = p.value.max, var.min = rownames(estimate.min), var.max = rownames(estimate.max),
                      model = mod[model], row.names = NULL, r2sp.max = r2sp.max, r2sp.min = r2sp.min, r2phylo.max = r2phylo.max, r2phylo.min = r2phylo.min,
                      nrow = nrow(dta), Loglik.max.lm = Loglik.max.lm, Loglik.min.lm = Loglik.min.lm, Loglik.min = Loglik.min, Loglik.max = Loglik.max,
                      estimate.diff = estimate.diff, p.value.diff = p.value.diff, var.diff = rownames(estimate.diff),Loglik.diff.lm = Loglik.diff.lm, Loglik.diff = Loglik.diff,
                      r2sp.diff = r2sp.diff, r2phylo.diff = r2phylo.diff)
    
    all_dta = rbind(all_dta, test)
  }
  return(all_dta)
}



################################################################################################################

mod = c("terrestrial",'northern','southern',"Passeriformes", "Non passeriformes")
dta_list = list(
  subset(species_traits, Marine == 0), # terrestrial
  subset(species_traits, Marine == 0 & distrib.core == 'north'), 
  subset(species_traits, Marine == 0 & distrib.core == 'south'), 
  subset(species_traits, IOCOrder == "Passeriformes" & Marine == 0), # passeriformes
  subset(species_traits, IOCOrder != "Passeriformes" & Marine == 0) # non passeriformes
  )

all_dta = mod_loop(dta_list, mod = mod)
unique(all_dta$model)


library(dplyr)
library(tidyr)

### putting the results into a nice table

format.res = function(data){
  
  ## 1. rear edge
  x.min <- as_tibble(data) %>% group_by(model, var.min) %>% summarise(p.value_mean = mean(p.value.min), estimate_mean = mean(estimate.min),  estimate_sd = sd(estimate.max), nrow = mean(nrow))
  x.min = merge(x.min, as_tibble(data) %>% group_by(model) %>% summarise(r2phylo.mean = mean(r2phylo.min), Loglik = mean(Loglik.min), Loglik.lm = mean(Loglik.min.lm)))
  print("Rear edge")
  x.min$significance = ifelse(x.min$p.value_mean<0.001,"***",ifelse(x.min$p.value_mean<0.01,"**",ifelse(x.min$p.value_mean<0.05,"*",ifelse(x.min$p.value_mean<0.1,".",""))))
  x.min[which(x.min$p.value_mean<0.1),]
  # round the results
  x.min[,c("p.value_mean","estimate_mean","estimate_sd","r2phylo.mean")] = 
    round(x.min[,c("p.value_mean","estimate_mean","estimate_sd","r2phylo.mean")],3)
  
  x.min_wide = rbind(cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:9)],"estimate_mean")], names_from = var.min, values_from = estimate_mean),value = "estimate_mean"),
                     cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:9)],"estimate_sd")], names_from = var.min, values_from = estimate_sd),value = "estimate_sd"),
                     cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:9)],"p.value_mean")], names_from = var.min, values_from = p.value_mean),value = "p.value_mean"),
                     cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:9)],"significance")], names_from = var.min, values_from = significance),value = "significance"))
  x.min_wide$shift = "min"
  
  ## 1. leading edge
  x.max  <- as_tibble(data) %>% group_by(model, var.max) %>% summarise(p.value_mean = mean(p.value.max), estimate_mean = mean(estimate.max), estimate_sd = sd(estimate.max), nrow = mean(nrow))
  x.max  = merge(x.max , as_tibble(data) %>% group_by(model) %>% summarise(r2phylo.mean = mean(r2phylo.max), Loglik = mean(Loglik.max), Loglik.lm = mean(Loglik.max.lm)))
  x.max$significance = ifelse(x.max$p.value_mean<0.001,"***",ifelse(x.max$p.value_mean<0.01,"**",ifelse(x.max$p.value_mean<0.05,"*",ifelse(x.max$p.value_mean<0.1,".",""))))
  ## print the significant variables
  print("Leading edge")
  x.max[which(x.max$p.value_mean<0.1),]
  # round the results
  x.max[,c("p.value_mean","estimate_mean","estimate_sd","r2phylo.mean")] = 
    round(x.max[,c("p.value_mean","estimate_mean","estimate_sd","r2phylo.mean")],3)
  
  x.max_wide = rbind(cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:9)],"estimate_mean")], names_from = var.max, values_from = estimate_mean),value = "estimate_mean"),
                     cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:9)],"estimate_sd")], names_from = var.max, values_from = estimate_sd),value = "estimate_sd"),
                     cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:9)],"p.value_mean")], names_from = var.max, values_from = p.value_mean),value = "p.value_mean"),
                     cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:9)],"significance")], names_from = var.max, values_from = significance),value = "significance"))
  x.max_wide$shift = "max"
  
  ## 3. expansion
  x.diff  <- as_tibble(data) %>% group_by(model, var.diff) %>%
    summarise(p.value_mean = mean(p.value.diff), estimate_mean = mean(estimate.diff), estimate_sd = sd(estimate.diff), nrow = mean(nrow))
  x.diff  = merge(x.diff , as_tibble(data) %>% group_by(model) %>%
                    summarise(r2phylo.mean = mean(r2phylo.diff), Loglik = mean(Loglik.diff), Loglik.lm = mean(Loglik.diff.lm)))
  x.diff$significance = ifelse(x.diff$p.value_mean<0.001,"***",ifelse(x.diff$p.value_mean<0.01,"**",ifelse(x.diff$p.value_mean<0.05,"*",ifelse(x.diff$p.value_mean<0.1,".",""))))
  print("Leading edge")
  x.diff[which(x.diff$p.value_mean<0.1),]
  x.diff[,c("p.value_mean","estimate_mean","estimate_sd","r2phylo.mean")] = 
    round(x.diff[,c("p.value_mean","estimate_mean","estimate_sd","r2phylo.mean")],3)
  
  x.diff_wide = rbind(cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:9)],"estimate_mean")], names_from = var.diff, values_from = estimate_mean),value = "estimate_mean"),
                      cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:9)],"estimate_sd")], names_from = var.diff, values_from = estimate_sd),value = "estimate_sd"),
                      cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:9)],"p.value_mean")], names_from = var.diff, values_from = p.value_mean),value = "p.value_mean"),
                      cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:9)],"significance")], names_from = var.diff, values_from = significance),value = "significance"))
  x.diff_wide$shift = "diff"
  
  ### combine everything
  x.wide = rbind(x.max_wide,x.min_wide,x.diff_wide)
  
  x.wide$p.value.phylo = pchisq(2*(x.wide$Loglik - x.wide$Loglik.lm), df = 1, lower.tail = FALSE)
  x.wide$r2 = round(1 - exp(-2*(x.wide$Loglik - x.wide$Loglik.lm)/round(x.wide$nrow*0.95)),4)
  
  return(x.wide)
  
}

x.wide = format.res(all_dta)

setwd("E:/TheseSwansea/TraitStudy/code_Miguel")
write.csv(x.wide, "pglmm_scaled_terrestrialNS_PCs_Migrants.csv")


################################################################################
### test of radar chart

library(fmsb)

x.wide.chart_max = x.wide[which(x.wide$value == "estimate_mean"&x.wide$shift == "max"),]
for (i in colnames(x.wide.chart_max[,7:20])){
  x.wide.chart_max[,i] = as.numeric(x.wide.chart_max[,i])
}
max = as.numeric(lapply(x.wide.chart_max[,7:20], FUN = max))
min = as.numeric(lapply(x.wide.chart_max[,7:20], FUN = min))
x.wide.chart_max1 = data.frame(rbind(max,min,x.wide.chart_max[,7:20]))
rownames(x.wide.chart_max1) = c("max","min",paste(x.wide.chart_max$model, x.wide.chart_max$shift))


radarchart(x.wide.chart_max1, pcol = terrain.colors(7))
legend("bottomright",legend = rownames(x.wide.chart_max1)[3:9], col = terrain.colors(7), bty = "n", pch = 20)


x.wide.chart_min = x.wide[which(x.wide$value == "estimate_mean"&x.wide$shift == "min"),]
for (i in colnames(x.wide.chart_min[,7:20])){
  x.wide.chart_min[,i] = as.numeric(x.wide.chart_min[,i])
}
max = as.numeric(lapply(x.wide.chart_min[,7:20], FUN = max))
min = as.numeric(lapply(x.wide.chart_min[,7:20], FUN = min))
x.wide.chart_min1 = data.frame(rbind(max,min,x.wide.chart_min[,7:20]))
rownames(x.wide.chart_min1) = c("max","min",paste(x.wide.chart_min$model, x.wide.chart_min$shift))




radarchart(x.wide.chart_min1, pcol = terrain.colors(7))
legend("bottomright",legend = rownames(x.wide.chart_min1)[3:9], col = terrain.colors(7), bty = "n", pch = 20)

par(mfrow = c(1,2))


radarchart(x.wide.chart_max1, pcol = terrain.colors(7), cglty = 1, cglcol = "gray", plty = 1)
legend("bottomright",legend = rownames(x.wide.chart_max1)[3:9], col = terrain.colors(7), bty = "n", pch = 20)

radarchart(x.wide.chart_min1, pcol = terrain.colors(7), cglty = 1, cglcol = "gray", plty = 1)
legend("bottomright",legend = rownames(x.wide.chart_min1)[3:9], col = terrain.colors(7), bty = "n", pch = 20)


