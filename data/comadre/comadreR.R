library(tidyverse)
library(dplyr)
library(magrittr)


subsetDB <- function(db,sub){
  e <- substitute(sub)
  r <- eval(e, db$metadata, parent.frame())
  subsetID <- (1:length(r))[r & !is.na(r)]
  
  # First make a copy of the database.
  ssdb <- db
  
  # Subset the sub-parts of the database
  ssdb$metadata <- ssdb$metadata[subsetID,]
  ssdb$mat <- ssdb$mat[subsetID]
  ssdb$matrixClass <- ssdb$matrixClass[subsetID]
  
  # Version information is retained, but modified as follows.
  if("version" %in% names(ssdb)){
    ssdb$version$Version <- paste(ssdb$version$Version," - subset created on ",format(Sys.time(), "%b_%d_%Y"),sep="")
    ssdb$version$DateCreated <- paste(ssdb$version$DateCreated," - subset created on ",format(Sys.time(), "%b_%d_%Y"),sep="")
    ssdb$version$NumberAcceptedSpecies <- length(unique(ssdb$metadata$SpeciesAccepted))
    ssdb$version$NumberStudies <- length(unique(ssdb$metadata$SpeciesAuthor))
    ssdb$version$NumberMatrices <- length(ssdb$mat)
  }
  
  return(ssdb)
}

getwd()
setwd("/Users/dcsuh/Desktop/datasets/COMPADRE")

# How to make a data frame with only Species names?
compadreSpecies <- subsetDB(comadre, SpeciesAccepted == SpeciesAccepted)

x <- subsetDB(comadre, MatrixCriteriaOntogeny == "Yes")
x$metadata$SpeciesAccepted
x$metadata$CommonName


n=21
x$metadata$SpeciesAccepted[n]
x$metadata$CommonName[n]
x$mat[[n]]$matA
x$matrixClass[[n]]$MatrixClassAuthor


eid <- read_csv("/Users/dcsuh/Desktop/datasets/EID2/SpeciesInteractions_EID2.csv")
mammal <- eid %>% select(`Carrier classification`, Carrier, Cargo, `Cargo classification`)
mammal %<>% filter(`Carrier classification` == "Mammal")


modeOfLife <- read_csv("/Users/dcsuh/Desktop/fossorial/Supplementary_data.txt")
