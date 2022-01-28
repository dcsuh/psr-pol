#Comparison of GMPD and Healy et al dataset
#Annakate Schatz

library(plyr)
library(dplyr)
library(magrittr)
library(tidyverse) #ggplot2, purrr, readr, tidyr

#load GMPD
g<-read_csv("data/source/gmpd_full.csv")[,-29]

#load Healy et al data
h<-read_csv("data/healy_axis_analysis.csv")
#replace underscore with space in species names
h$species<-gsub("_"," ",h$species)

#compare species names
sum(unique(g$HostCorrectedName) %in% unique(h$species))
  length(unique(g$HostCorrectedName))
  length(unique(h$species))
#27 of 462 GMPD hosts are included in Healy's 121 hosts

#limit GMPD to those hosts in Healy's dataset (and records with parasites identified to species)
gh<-g %>% filter(HostCorrectedName %in% h$species & HasBinomialName=="yes")

#summarize hosts
table(gh$Group) #3 host groups
table(gh$HostOrder) #3 host orders
table(gh$HostFamily) #14 host families
table(gh$HostEnvironment) #3 host environments
  sum(gh$HostEnvironment=="terrestrial")/nrow(gh)

#summarize parasites
length(unique(gh$ParasiteCorrectedName)) #661 parasites

table(gh$ParType) #5/7 are good, fungi and prions not so much

length(unique(gh$ParPhylum)) #26 parasite phyla
sum(is.na(gh$ParPhylum)) #5 NA records

length(unique(gh$ParClass)) #44 parasite classes
sum(is.na(gh$ParClass)) #8 NA records

sum(is.na(gh$Prevalence))
sum(!is.na(gh$Prevalence))/nrow(gh)

sum(is.na(gh$Intensity))/nrow(gh)
unique(gh[!is.na(gh$Intensity),10])

write_csv(gh, "data/gmpd_subset.csv", col_names=T)


