##Script originator: Daniel suh
##Date: 10/12/20

#Parasite Pace of Life project with Park Lab
#look at EID2 and compare with database from Kevin Healy


library(tidyverse)
library(magrittr)


eid <- read_csv("data/source/eid2_full.csv")
healy <- read_csv("data/healy_axis_analysis.csv")

eid_hosts <- eid %>% select(Carrier) %>% distinct()
healy %<>% mutate(Carrier = gsub("_", " ", species)) %>% mutate_all(., tolower)
healy_hosts <- healy %>% select(Carrier) %>% distinct() 

match <- which(healy_hosts$Carrier %in% eid_hosts$Carrier)
matched_hosts <- healy_hosts$Carrier[match]
healy_matched_hosts <- healy[which(healy$Carrier %in% matched_hosts),]

#this is a dataframe with all host-parasite associations from EID2 that match hosts from the healy dataset
eid_matched_hosts <- eid[which(eid$Carrier %in% matched_hosts),]

#this is a datafame with all unique hosts and their PSR calculated from data from EID2
eid_PSR <- eid_matched_hosts %>% group_by(Carrier) %>% summarize(PSR=n())

which(matched_hosts %in% eid_PSR$Carrier)

healy_matched_hosts %<>% mutate(Carrier = gsub("_", " ", species)) %>% mutate_all(., tolower)

#this is a dataframe with all host data from Healy data (including multiple sources for hosts) and PSR for each host
PSR_traits <- full_join(healy_matched_hosts, eid_PSR, by = "Carrier")
PSR_traits %>% filter(., species!="homo_sapiens") %>% filter(., species!="canis_lupus") %>% group_by(species) %>% summarize(PSR=PSR) %>% distinct() %>% ggplot(.,aes(x=PSR)) +
  geom_histogram(binwidth = 1)

#write_csv(eid_matched_hosts, "data/eid2_subset.csv")
