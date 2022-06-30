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
  
  SppID <- read.table(file = "E:/TheseSwansea/Galiana2021_network-area-europe-master/SppID.txt", header = TRUE, stringsAsFactors = F)
  
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

