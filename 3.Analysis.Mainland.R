### Script 3 - Thompson et al 2023



################################################################################
### This code runs the PGLMMs relating species range shift to species traits 
### for three groups of species (North/South/Passeriformes)
###############################################################################


rm(list = ls())


setwd("E:/TheseSwansea/TraitStudy/Github/data")
species_traits <- read.csv('SpecTrait_Full_062023_156sp.csv', stringsAsFactors = FALSE, row.names = 1)
species_traits.full <- read.csv('SpecTrait_012023_159sp.csv', stringsAsFactors = FALSE, row.names = 1)

library(stringr)

species_traits = species_traits[,-which(names(species_traits) %in% c('nrecord_P.1','nrecord_P.3'))]
species_traits.both = merge(species_traits[,c("speccode",'distrib.core',str_subset(names(species_traits), "dist_km"))], 
                            species_traits.full[,c("speccode",'distrib.core',str_subset(names(species_traits), "dist_km"))], 
                            by = 'speccode')

## .x = no Orkney and Shetlands
## .y = full dataset
## leading edge shift are smaller (for the most part, a couple species with
## strange changes)

par(mfrow = c(2,2))
plot(data = species_traits.both, 
     shift_max20_P.1.3_dist_km.x ~ shift_max20_P.1.3_dist_km.y)

## rear edge shifts are identical
plot(data = species_traits.both, 
     shift_min20_P.1.3_dist_km.x ~ shift_min20_P.1.3_dist_km.y)


# Marine or non marine species : (Marine = >75% of points within 20km of the coastline)
propmarine <- read.csv('Speciestraits_ProportionMarineBTOsp_159sp.csv', stringsAsFactors = FALSE, row.names = 1)

species_traits = merge(species_traits, propmarine[,c('speccode','prop_Marine20km')])
species_traits$Marine = ifelse(species_traits$prop_Marine20km>=75, 1, 0)

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


###############################################################################

library(phyr)
require(phytools)
require(ape)

tree = ape::read.tree("BigBird.All.NewNames.6714Taxa.tre/BigBird.All.NewNames.6714Taxa.tre")
# tree <- consensus.edges(tree[1:10])
tree <- tree[[1]]

parseRowNames <- function(x){
  splitted <- strsplit(x, '_')[[1]]
  l <- length(splitted)
  return(paste0(splitted[l-1], '_', splitted[l]))
}

tree$tip.label <- unlist(lapply(tree$tip.label, parseRowNames))

# get matching species names with phylo tree
species_traits$species = species_traits$scientific_name_
sptree_match = read.csv("BigBird.All.NewNames.6714Taxa.tre/PhyloTree_speciesMatch.csv", sep = ";")
names(sptree_match)[1] = "species"
species_traits$species[match(sptree_match$scientific_name_,species_traits$scientific_name_)] = sptree_match$species

pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, species_traits$species))

################################################################################
### LM checking assumptions

mod = pglmm(sqrt(shift_min20_P.1.3_dist_km) ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
              (1|species__), data = species_traits, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
mod1 = pglmm(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
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


library(MASS)

dta = species_traits
dta$log.BodyMass.Value = log10(dta$BodyMass.Value)
mod.rlm = rlm(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1),
          data = dta[which(dta$distrib.core == 'north'),])
summary(mod.rlm)



mod.max_shift <- pglmm(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                         (1|species__), data = dta[which(dta$distrib.core == 'north'),], family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)



################################################################################
## TRYING GAM MODELS >> too few points, overfitting, don't know how to include phylo random effect
library(mgcv)
## ALL SPECIES
mod = gam(shift_min20_P.1.3_dist_km ~ s(scale(log.BodyMass.Value)) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + s(scale(pc1_env)) + s(scale(pc2_env)) + s(scale(pc1_lc)) + s(scale(pc2_lc)) + s(scale(P.1)), data = species_traits)
summary(mod)

mod = gam(shift_max20_P.1.3_dist_km ~ s(scale(log.BodyMass.Value)) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + s(scale(pc1_env)) + s(scale(pc2_env)) + s(scale(pc1_lc)) + s(scale(pc2_lc)) + s(scale(P.1)), data = species_traits)
summary(mod)

#### SOUTH
mod = gam(shift_min20_P.1.3_dist_km ~ s(scale(log.BodyMass.Value)) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + s(scale(P.1)), data = subset(species_traits, distrib.core == 'south'))
summary(mod)

mod = gam(shift_max20_P.1.3_dist_km ~ s(scale(log.BodyMass.Value)) + s(scale(habitat_gen)) + s(scale(normalised_indegree)) + s(scale(diet_diversity)) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + s(scale(trophic_position)) + s(scale(outdegree)) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + s(scale(P.1)), data = subset(species_traits, distrib.core == 'south'))
summary(mod)

#### NORTH
mod = gam(shift_min20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + scale(trophic_position) + s(scale(outdegree)) + scale(pc1_env), data = subset(species_traits, distrib.core == 'north'))
summary(mod)

mod = gam(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + s(scale(dist_N_km_max20_P.1)) + s(scale(dist_S_km_min20_P.1))  +  migratory_binomial  + scale(trophic_position) + s(scale(outdegree)) + scale(pc1_env), data = subset(species_traits, distrib.core == 'north'))
summary(mod)




################################################################################################################


hist(species_traits$shift_max20_P.1.3_dist_km)

### FUNCTION FOR RUNNING ALL MODELS IN ONE LOOP
pglmm_shift = function(dta_list, mod){
  
  library(car)
  model = 0
  all_dta = data.frame() # dataframe to store the data
  
  for (dta in dta_list){ # list of subsets to run the models for
    
    model = model + 1
    
    print(model)
    
    # variables to save
    estimate.max = NULL; p.value.max = NULL;std.Error.max = NULL
    estimate.min = NULL; p.value.min = NULL;std.Error.min = NULL
    estimate.diff = NULL; p.value.diff = NULL;std.Error.diff = NULL
    phylo.min = NULL; phylo.max = NULL; phylo.diff = NULL # species__ effect
    sp.min = NULL; sp.max = NULL; sp.diff = NULL # species effect 
    Loglik.max = NULL;Loglik.min = NULL;Loglik.diff = NULL # log likelihood with phylo signal
    Loglik.max.lm = NULL; Loglik.min.lm = NULL; Loglik.diff.lm = NULL # log likelihood without phylo signal
    R2.min = NULL; R2.max = NULL; R2.diff = NULL # overall R2 of pglmm
    
    # dta$log.trophic_position = log(dta$trophic_position)
    dta$log.BodyMass.Value = log10(dta$BodyMass.Value)
    # dta$dist_S_km_min20_P.1 = log(dta$dist_S_km_min20_P.1)
    # dta$dist_N_km_max20_P.1 = log(dta$dist_N_km_max20_P.1)
    
    # model with migratory status as a covariate
    
    
    # with pglmm
    mod.min_shift = pglmm(shift_min20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                            (1|species__), data = dta, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
    # with lm
    mod.min_shift_lm <- lm(shift_min20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                           , data = dta)
    # with pglmm
    mod.max_shift <- pglmm(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                             (1|species__), data = dta, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
    # with lm
    mod.max_shift_lm <- lm(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                           , data = dta)
    # with pglmm
    mod.diff_shift = pglmm(shift_diff ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                             (1|species__), data = dta, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
    # with lm
    mod.diff_shift_lm <- lm(shift_diff ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1)  +  migratory_binomial  + scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1)
                            , data = dta)
    
    
    # save the outputs
    estimate.min = rbind(estimate.min, mod.min_shift$B)
    p.value.min = rbind(p.value.min, mod.min_shift$B.pvalue)
    std.Error.min = rbind(std.Error.min, mod.min_shift$B.se)
    sp.min = rbind(sp.min, mod.min_shift$s2r[[1]])
    phylo.min = rbind(phylo.min, mod.min_shift$s2r[[2]])
    Loglik.min = rbind(Loglik.min, mod.min_shift$logLik)
    Loglik.min.lm = rbind(Loglik.min.lm, logLik(mod.min_shift_lm))
    R2.min = rbind(R2.min, rr2::R2(mod.min_shift)[[1]])
    
    estimate.max = rbind(estimate.max, mod.max_shift$B)
    p.value.max = rbind(p.value.max, mod.max_shift$B.pvalue)
    std.Error.max = rbind(std.Error.max, mod.max_shift$B.se)
    sp.max = rbind(sp.max, mod.max_shift$s2r[[1]])
    phylo.max = rbind(phylo.max, mod.max_shift$s2r[[2]])
    Loglik.max = rbind(Loglik.max, mod.max_shift$logLik)
    Loglik.max.lm = rbind(Loglik.max.lm, logLik(mod.max_shift_lm))
    R2.max = rbind(R2.max, rr2::R2(mod.max_shift)[[1]])
    
    # expansion
    estimate.diff = rbind(estimate.diff, mod.diff_shift$B)
    p.value.diff = rbind(p.value.diff, mod.diff_shift$B.pvalue)
    std.Error.diff = rbind(std.Error.max, mod.diff_shift$B.se)
    sp.diff = rbind(sp.diff, mod.diff_shift$s2r[[1]])
    phylo.diff = rbind(phylo.diff, mod.diff_shift$s2r[[2]])
    Loglik.diff = rbind(Loglik.diff, mod.diff_shift$logLik)
    Loglik.diff.lm = rbind(Loglik.diff.lm, logLik(mod.diff_shift_lm))
    R2.diff = rbind(R2.diff, rr2::R2(mod.diff_shift)[[1]])
    
    
    print(mod[model])
    
    # summary
    test = data.frame(estimate.max = estimate.max, estimate.min = estimate.min, p.value.min = p.value.min, 
                      p.value.max = p.value.max, var.min = rownames(estimate.min), var.max = rownames(estimate.max),
                      model = mod[model], row.names = NULL, sp.max = sp.max, sp.min = sp.min, phylo.max = phylo.max, phylo.min = phylo.min,
                      nrow = nrow(dta), Loglik.max.lm = Loglik.max.lm, Loglik.min.lm = Loglik.min.lm, Loglik.min = Loglik.min, Loglik.max = Loglik.max,
                      estimate.diff = estimate.diff, p.value.diff = p.value.diff, var.diff = rownames(estimate.diff),Loglik.diff.lm = Loglik.diff.lm, Loglik.diff = Loglik.diff,
                      sp.diff = sp.diff, phylo.diff = phylo.diff, se.diff = std.Error.diff, se.min = std.Error.min, se.max = std.Error.max,
                      R2.min = R2.min, R2.max = R2.max, R2.diff = R2.diff)
    
    all_dta = rbind(all_dta, test)
    
    ###############################
    ## plotting
    
    list.var = c('log.BodyMass.Value', 'habitat_gen', 'normalised_indegree', 'diet_diversity', 'dist_N_km_max20_P.1', 'dist_S_km_min20_P.1', 
                 'trophic_position', 'outdegree', 'pc1_env', 'pc2_env', 'pc1_lc', 'pc2_lc', 'P.1')
    
    for (shift in c('shift_max20_P.1.3_dist_km','shift_min20_P.1.3_dist_km','shift_diff')){
      
      dta$Y = dta[,shift]
      if (shift == 'shift_max20_P.1.3_dist_km'){mod.p = mod.max_shift_lm} else if (shift =='shift_min20_P.1.3_dist_km'){mod.p = mod.min_shift_lm} else {mod.p = mod.diff_shift_lm}
      
      coef.all = mod.p$coefficients
      pvalue = summary(mod.p)$coefficients[,4]
      
      names.all = names(coef.all)
      
      names = names.all[which(pvalue<=0.1)]
      
      if("(Intercept)" %in% names){names = names[-1]}
      if("migratory_binomialMigrant" %in% names){names = names[-which(names == "migratory_binomialMigrant")]}
      
      par(mfrow = c(3,2))
      
      if(length(names)>0){
        for(co in 1:length(names)){          
          print(names[co])
          jpeg(paste0('E:/TheseSwansea/TraitStudy/Github/plots/mcPlots.LinearModels.',mod[model],'.',shift,'.',names[co],".jpg"), width = 800, height = 600, quality = 1000)
          car::mcPlot(mod.p, variable = names[co], overlaid = F, col.marginal='black', col.conditional= 'green4', title = F, new = F)
          # The second conditional plot is the added-variable plot of e(Y|Z) versus e(X|Z) where e(a|b) means the Pearson residuals from the regression of a on b.
          dev.off()
        }
      }
    }
  }
  return(all_dta)
}




### putting the results into a nice table

format.res = function(data){
  
  library(dplyr)
  library(tidyr)
  ## 1. rear edge
  x.min <- as_tibble(data) %>% group_by(model, var.min) %>% summarise(p.value_mean = mean(p.value.min), estimate_mean = mean(estimate.min),  estimate_se = mean(se.min), nrow = mean(nrow), R2.rr2 = mean(R2.min))
  x.min = merge(x.min, as_tibble(data) %>% group_by(model) %>% summarise(phylo.mean = mean(phylo.min), Loglik = mean(Loglik.min), Loglik.lm = mean(Loglik.min.lm)))
  print("Rear edge")
  x.min$significance = ifelse(x.min$p.value_mean<0.001,"***",ifelse(x.min$p.value_mean<0.01,"**",ifelse(x.min$p.value_mean<0.05,"*",ifelse(x.min$p.value_mean<0.1,".",""))))
  x.min[which(x.min$p.value_mean<0.1),]
  # round the results
  x.min[,c("p.value_mean","estimate_mean","estimate_se","phylo.mean")] = 
    round(x.min[,c("p.value_mean","estimate_mean","estimate_se","phylo.mean")],3)
  
  x.min_wide = rbind(cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:10)],"estimate_mean")], names_from = var.min, values_from = estimate_mean),value = "estimate_mean"),
                     cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:10)],"estimate_se")], names_from = var.min, values_from = estimate_se),value = "estimate_se"),
                     cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:10)],"p.value_mean")], names_from = var.min, values_from = p.value_mean),value = "p.value_mean"),
                     cbind(pivot_wider(x.min[,c(colnames(x.min)[c(1,2,6:10)],"significance")], names_from = var.min, values_from = significance),value = "significance"))
  x.min_wide$shift = "min"
  
  ## 1. leading edge
  x.max  <- as_tibble(data) %>% group_by(model, var.max) %>% summarise(p.value_mean = mean(p.value.max), estimate_mean = mean(estimate.max),  estimate_se = mean(se.max), nrow = mean(nrow), R2.rr2 = mean(R2.max))
  x.max  = merge(x.max , as_tibble(data) %>% group_by(model) %>% summarise(phylo.mean = mean(phylo.max), Loglik = mean(Loglik.max), Loglik.lm = mean(Loglik.max.lm)))
  x.max$significance = ifelse(x.max$p.value_mean<0.001,"***",ifelse(x.max$p.value_mean<0.01,"**",ifelse(x.max$p.value_mean<0.05,"*",ifelse(x.max$p.value_mean<0.1,".",""))))
  ## print the significant variables
  print("Leading edge")
  x.max[which(x.max$p.value_mean<0.1),]
  # round the results
  x.max[,c("p.value_mean","estimate_mean","estimate_se","phylo.mean")] = 
    round(x.max[,c("p.value_mean","estimate_mean","estimate_se","phylo.mean")],3)
  
  x.max_wide = rbind(cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:10)],"estimate_mean")], names_from = var.max, values_from = estimate_mean),value = "estimate_mean"),
                     cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:10)],"estimate_se")], names_from = var.max, values_from = estimate_se),value = "estimate_se"),
                     cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:10)],"p.value_mean")], names_from = var.max, values_from = p.value_mean),value = "p.value_mean"),
                     cbind(pivot_wider(x.max[,c(colnames(x.max)[c(1,2,6:10)],"significance")], names_from = var.max, values_from = significance),value = "significance"))
  x.max_wide$shift = "max"
  
  ## 3. expansion
  x.diff  <- as_tibble(data) %>% group_by(model, var.diff) %>%
    summarise(p.value_mean = mean(p.value.diff), estimate_mean = mean(estimate.diff), estimate_se = mean(se.diff), nrow = mean(nrow), R2.rr2 = mean(R2.diff))
  x.diff  = merge(x.diff , as_tibble(data) %>% group_by(model) %>%
                    summarise(phylo.mean = mean(phylo.diff), Loglik = mean(Loglik.diff), Loglik.lm = mean(Loglik.diff.lm)))
  x.diff$significance = ifelse(x.diff$p.value_mean<0.001,"***",ifelse(x.diff$p.value_mean<0.01,"**",ifelse(x.diff$p.value_mean<0.05,"*",ifelse(x.diff$p.value_mean<0.1,".",""))))
  print("Leading edge")
  x.diff[which(x.diff$p.value_mean<0.1),]
  x.diff[,c("p.value_mean","estimate_mean","estimate_se","phylo.mean")] = 
    round(x.diff[,c("p.value_mean","estimate_mean","estimate_se","phylo.mean")],3)
  
  x.diff_wide = rbind(cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:10)],"estimate_mean")], names_from = var.diff, values_from = estimate_mean),value = "estimate_mean"),
                      cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:10)],"estimate_se")], names_from = var.diff, values_from = estimate_se),value = "estimate_se"),
                      cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:10)],"p.value_mean")], names_from = var.diff, values_from = p.value_mean),value = "p.value_mean"),
                      cbind(pivot_wider(x.diff[,c(colnames(x.diff)[c(1,2,6:10)],"significance")], names_from = var.diff, values_from = significance),value = "significance"))
  x.diff_wide$shift = "diff"
  
  ### combine everything
  x.wide = rbind(x.max_wide,x.min_wide,x.diff_wide)
  
  x.wide$p.value.phylo = pchisq(2*(x.wide$Loglik - x.wide$Loglik.lm), df = 1, lower.tail = FALSE)
  x.wide$r2.phylo = round(1 - exp(-2*(x.wide$Loglik - x.wide$Loglik.lm)/round(x.wide$nrow*0.95)),4)
  
  return(x.wide)
  
}

################################################################################################################

mod = c("terrestrial",'northern','southern',"Passeriformes")
dta_list = list(
  subset(species_traits, Marine == 0 & speccode != 71), # terrestrial
  subset(species_traits, Marine == 0 & distrib.core == 'north' & speccode != 71),
  subset(species_traits, Marine == 0 & distrib.core == 'south' & speccode != 71),
  subset(species_traits, IOCOrder == "Passeriformes" & Marine == 0 & speccode != 71) # passeriformes
)

all_dta = pglmm_shift(dta_list, mod = mod)
unique(all_dta$model)


x.wide = format.res(all_dta)

write.csv(x.wide, "E:/TheseSwansea/TraitStudy/Github/results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.NoOrkneyNoShetlands_136sp.csv")





### convert coefficients into km
dta1 = subset(species_traits, Marine == 0 & distrib.core == 'north')
dta1 = subset(species_traits, IOCOrder == "Passeriformes" & Marine == 0)

## converting in units without scaling 
mod.max_shift_km <- pglmm(shift_max20_P.1.3_dist_km ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1 +
                         (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
mod.max_shift <- pglmm(shift_max20_P.1.3_dist_km ~ scale(log.BodyMass.Value) + scale(habitat_gen) + scale(normalised_indegree) + scale(diet_diversity) + scale(dist_N_km_max20_P.1) + scale(dist_S_km_min20_P.1) + migratory_binomial +  scale(trophic_position) + scale(outdegree) + scale(pc1_env) + scale(pc2_env) + scale(pc1_lc) + scale(pc2_lc) + scale(P.1) +
                         (1|species__), data = dta1, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)

summary(mod.max_shift_km)
summary(mod.max_shift)

21.7/sd(dta1$trophic_position)
71.1/sd(dta1$outdegree)

ggplot(data = dta1, aes(x = habitat_gen, y = shift_diff)) + geom_point() + geom_smooth()


hist(species_traits$shift_max20_P.1.3_dist_km[which(species_traits$distrib.core == "north")])


###################################################################################
### RUNNING MODELS WITHOUT SCALING TO GET COEFFICIENTS IN KM


pglmm_shift_km = function(dta_list, mod){
  
  library(car)
  model = 0
  all_dta = data.frame() # dataframe to store the data
  
  for (dta in dta_list){ # list of subsets to run the models for
    
    model = model + 1
    
    print(model)
    
    # variables to save
    estimate.max = NULL; p.value.max = NULL;std.Error.max = NULL
    estimate.min = NULL; p.value.min = NULL;std.Error.min = NULL
    estimate.diff = NULL; p.value.diff = NULL;std.Error.diff = NULL
    phylo.min = NULL; phylo.max = NULL; phylo.diff = NULL # species__ effect
    sp.min = NULL; sp.max = NULL; sp.diff = NULL # species effect 
    Loglik.max = NULL;Loglik.min = NULL;Loglik.diff = NULL # log likelihood with phylo signal
    Loglik.max.lm = NULL; Loglik.min.lm = NULL; Loglik.diff.lm = NULL # log likelihood without phylo signal
    R2.min = NULL; R2.max = NULL; R2.diff = NULL # overall R2 of pglmm
    
    # dta$log.trophic_position = log(dta$trophic_position)
    dta$log.BodyMass.Value = log10(dta$BodyMass.Value)
    # dta$dist_S_km_min20_P.1 = log(dta$dist_S_km_min20_P.1)
    # dta$dist_N_km_max20_P.1 = log(dta$dist_N_km_max20_P.1)
    
    # model with migratory status as a covariate
    
    
    # with pglmm
    mod.min_shift = pglmm(shift_min20_P.1.3_dist_km ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1 +
                            (1|species__), data = dta, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
    # with lm
    mod.min_shift_lm <- lm(shift_min20_P.1.3_dist_km ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1
                           , data = dta)
    # with pglmm
    mod.max_shift <- pglmm(shift_max20_P.1.3_dist_km ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1 +
                             (1|species__), data = dta, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
    # with lm
    mod.max_shift_lm <- lm(shift_max20_P.1.3_dist_km ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1
                           , data = dta)
    # with pglmm
    mod.diff_shift = pglmm(shift_diff ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1 +
                             (1|species__), data = dta, family = "gaussian", cov_ranef = list(species = pruned_tree), REML = F)
    # with lm
    mod.diff_shift_lm <- lm(shift_diff ~ log.BodyMass.Value + habitat_gen + normalised_indegree + diet_diversity + dist_N_km_max20_P.1 + dist_S_km_min20_P.1  +  migratory_binomial  + trophic_position + outdegree + pc1_env + pc2_env + pc1_lc + pc2_lc + P.1
                            , data = dta)
    
    
    # save the outputs
    estimate.min = rbind(estimate.min, mod.min_shift$B)
    p.value.min = rbind(p.value.min, mod.min_shift$B.pvalue)
    std.Error.min = rbind(std.Error.min, mod.min_shift$B.se)
    sp.min = rbind(sp.min, mod.min_shift$s2r[[1]])
    phylo.min = rbind(phylo.min, mod.min_shift$s2r[[2]])
    Loglik.min = rbind(Loglik.min, mod.min_shift$logLik)
    Loglik.min.lm = rbind(Loglik.min.lm, logLik(mod.min_shift_lm))
    R2.min = rbind(R2.min, rr2::R2(mod.min_shift)[[1]])
    
    estimate.max = rbind(estimate.max, mod.max_shift$B)
    p.value.max = rbind(p.value.max, mod.max_shift$B.pvalue)
    std.Error.max = rbind(std.Error.max, mod.max_shift$B.se)
    sp.max = rbind(sp.max, mod.max_shift$s2r[[1]])
    phylo.max = rbind(phylo.max, mod.max_shift$s2r[[2]])
    Loglik.max = rbind(Loglik.max, mod.max_shift$logLik)
    Loglik.max.lm = rbind(Loglik.max.lm, logLik(mod.max_shift_lm))
    R2.max = rbind(R2.max, rr2::R2(mod.max_shift)[[1]])
    
    # expansion
    estimate.diff = rbind(estimate.diff, mod.diff_shift$B)
    p.value.diff = rbind(p.value.diff, mod.diff_shift$B.pvalue)
    std.Error.diff = rbind(std.Error.max, mod.diff_shift$B.se)
    sp.diff = rbind(sp.diff, mod.diff_shift$s2r[[1]])
    phylo.diff = rbind(phylo.diff, mod.diff_shift$s2r[[2]])
    Loglik.diff = rbind(Loglik.diff, mod.diff_shift$logLik)
    Loglik.diff.lm = rbind(Loglik.diff.lm, logLik(mod.diff_shift_lm))
    R2.diff = rbind(R2.diff, rr2::R2(mod.diff_shift)[[1]])
    
    
    print(mod[model])
    
    # summary
    test = data.frame(estimate.max = estimate.max, estimate.min = estimate.min, p.value.min = p.value.min, 
                      p.value.max = p.value.max, var.min = rownames(estimate.min), var.max = rownames(estimate.max),
                      model = mod[model], row.names = NULL, sp.max = sp.max, sp.min = sp.min, phylo.max = phylo.max, phylo.min = phylo.min,
                      nrow = nrow(dta), Loglik.max.lm = Loglik.max.lm, Loglik.min.lm = Loglik.min.lm, Loglik.min = Loglik.min, Loglik.max = Loglik.max,
                      estimate.diff = estimate.diff, p.value.diff = p.value.diff, var.diff = rownames(estimate.diff),Loglik.diff.lm = Loglik.diff.lm, Loglik.diff = Loglik.diff,
                      sp.diff = sp.diff, phylo.diff = phylo.diff, se.diff = std.Error.diff, se.min = std.Error.min, se.max = std.Error.max,
                      R2.min = R2.min, R2.max = R2.max, R2.diff = R2.diff)
    
    all_dta = rbind(all_dta, test)
    
    ###############################
    ## plotting
    
    for (shift in c('shift_max20_P.1.3_dist_km','shift_min20_P.1.3_dist_km','shift_diff')){
      
      dta$Y = dta[,shift]
      if (shift == 'shift_max20_P.1.3_dist_km'){mod.p = mod.max_shift} else if (shift =='shift_min20_P.1.3_dist_km'){mod.p = mod.min_shift} else {mod.p = mod.diff_shift}
      
      coef = mod.p$B
      pvalue = mod.p$B.pvalue
      
      int = coef[1]
      
      coef = coef[which(pvalue<=0.05)]
      names = rownames(pvalue)[which(pvalue<=0.05)]
      
      
      library(stringr)
      names = str_remove_all(names, c('scale\\('))
      names = str_remove_all(names, c('\\)'))
      names = str_remove_all(names, c('log\\('))
      if("(Intercept" %in% names){names = names[-1]}
      if("migratory_binomialMigrant" %in% names){names = names[-which(names == "migratory_binomialMigrant")]}
      
      par(mfrow = c(2,2))
      
      if(length(names)>0){
        for(co in 1:length(names)){
          print(names[co])
          plot(y = dta[,'Y'], x = dta[,names[co]], pch = 20, main = paste(mod[model]), cex = 1,
               ylab = shift, xlab = names[co])
          # legend("bottomright", legend = d$ID.new.Bioregion, col = c(1:44), pch = 20, bty = "n")
          abline(int, coef[co])
          abline(0,0, col = 'pink')
        }
      }
    }
  }
  return(all_dta)
}



################################################################################################################


all_dta_km = pglmm_shift_km(dta_list, mod = mod)

x.wide_km = format.res(all_dta_km)

write.csv(x.wide_km, "results/pglmm_scaled_terrestrialNS_PCs_Std.Error.R2.NoShetlandsNoOrkney.km.csv")



