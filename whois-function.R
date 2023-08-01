#################################################
# Spatial networks WHOIS FUNCTION
# 10 KM
# 28/06/2016 Joao Braga
#################################################
# rm(list = ls())

#################################################
# Species codes and species names
# Function to Identify a spp by the code
whois <- function(SPPCODE = NULL, SPPNAME = NULL) {
  # Function to Identify a spp by the code
  
  if(is.null(SPPCODE) & is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  if(!is.null(SPPCODE) & !is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  
  SppID <- read.table(file = "F:/TheseSwansea/TraitStudy/Github/data/SppID.txt", header = TRUE, stringsAsFactors = F)
  
  if(length(SPPCODE) > 1){
    SPPCODE <- paste0(SPPCODE, "$", collapse = "|")
  }
  if(length(SPPNAME) > 1){
    SPPNAME <- paste0(SPPNAME, "$", collapse = "|")
  }
  
  if(!is.null(SPPCODE))    who <- SppID[which(SppID$ID==SPPCODE),]$SPPname
  if(!is.null(SPPNAME))    who <- SppID[grepl(pattern = SPPNAME,x = SppID$SPPname),]
  
  return(who)
}


whoisGBIF <- function(SPPCODE = NULL, SPPNAME = NULL) {
  
  library(stringr)
  # Function to Identify a spp by the code
  
  if(is.null(SPPCODE) & is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  if(!is.null(SPPCODE) & !is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  
  SppID <- read.csv(file = "E:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/GBIFSpecies_BackboneTaxonomy_minfo.csv", header = TRUE, stringsAsFactors = F)
  SPPCODE = as.numeric(str_replace(SPPCODE,'X',''))
  
  if(length(SPPCODE) > 1){
    SPPCODE <- paste0(SPPCODE, "$", collapse = "|")
  }
  if(length(SPPNAME) > 1){
    SPPNAME <- paste0(SPPNAME, "$", collapse = "|")
  }
  
  if(!is.null(SPPCODE))    who <- SppID[which(SppID$usageKey==SPPCODE),]$masterSpecies
  if(!is.null(SPPNAME))    who <- SppID[grepl(pattern = SPPNAME,x = SppID$masterSpecies),]
  
  return(who)
}

# Function to transform spp distribution (it can be network properties per pixel) into raster
fun.dbf2raster <- function(SPPPA, mask.dir = NULL){
  # SPPPA must be in this format - first colmun with CELL ID (PAGENALE) and second column with the value (to plot)
  if(is.null(mask.dir)) stop("Must specify the mask file directory")
  require(raster)
  require(foreign)
  
  maskID <- read.dbf(list.files(path = mask.dir, full.names = TRUE, pattern = ".img.vat.dbf$"))
  maskk <- raster(x = list.files(path = mask.dir, full.names = TRUE, pattern = ".img$"))
  
  spp <- maskID
  spp$val <- NA
  spp$PageName <- as.character(spp$PageName)
  row.names(spp) <- spp$PageName
  
  SPPPA$PAGENAME <- as.character(SPPPA$PAGENAME)
  SPPPA[,2] <- as.numeric(as.character(SPPPA[,2]))
  row.names(SPPPA) <- SPPPA$PAGENAME
  
  cellID <- as.character(SPPPA$PAGENAME)
  if( nrow(spp[cellID,]) != nrow(SPPPA[cellID,])) stop("Cell IDs do not match")
  spp <- spp[cellID,] 
  spp$val <- SPPPA[,2]
  
  xx <- raster::values(maskk)
  
  if( length(xx[!is.na(xx)]) != nrow(spp) ) stop("Mask size inadequate")
  xx[!is.na(xx)] <- spp$val    #[xx[!is.na(xx)]]
  
  raster::values(maskk) <- xx
  return(maskk)
}
