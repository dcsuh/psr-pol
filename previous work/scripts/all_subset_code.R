#Healy host-parasite lists
#Annakate Schatz

#Compile the full list of host-parasite associations from all sources for which data is available from Healy et al

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(tidyverse)
library(taxize)

##### load data #####

#host-parasite data to filter (re-doing it altogether for consistency)
ban<-read_csv("data/source/banesh_full.csv")
bro<-read_csv("data/source/brook_full.csv")
eid<-read_csv("data/source/eid2_full.csv")
gmpd<-read_csv("data/source/gmpd_full.csv")[,-29]
nhml<-read_csv("data/nhml_subset.csv") #note that this is already filtered through the download process
oli<-read_csv("data/source/olival_full.csv")

#Healy et al life history data against which to compare host species
healy<-read_csv("data/healy_axis_analysis.csv")

##### prep host data and subset #####

#set common column name for host species
    #edit contents to use spaces rather than underscores and sentence case
ban1<-ban %>% mutate(host=Host.species)
bro1<-bro %>% mutate(host=gsub("_"," ",hHostNameFinal))
eid1<-eid %>% mutate(host=str_to_sentence(Carrier))
gmpd1<-gmpd %>% mutate(host=HostCorrectedName)
nhml1<-nhml %>% mutate(host=Host)
oli1<-oli %>% mutate(host=gsub("_"," ",hHostNameFinal))
healy1<-healy %>% mutate(host=gsub("_"," ",species))

#filter by host species
ban.sub<-ban1 %>% filter(host %in% healy1$host)
bro.sub<-bro1 %>% filter(host %in% healy1$host)
eid.sub<-eid1 %>% filter(host %in% healy1$host)
gmpd.sub<-gmpd1 %>% filter(host %in% healy1$host)
nhml.sub<-nhml1 %>% filter(host %in% healy1$host)
oli.sub<-oli1 %>% filter(host %in% healy1$host)

##### edit host data #####

#specific edit for Banesh data
ban.sub %<>% filter(Def.int=="def") #only definitive hosts

## host orders ##

#for datasets with host order data: set common column name
    #edit contents to sentence case
bro.sub1<-bro.sub %>% mutate(hostOrder=str_to_sentence(hOrder))
gmpd.sub1<-gmpd.sub %>% mutate(hostOrder=HostOrder)

#for datasets without host order data: obtain from NCBI, add column
need.orders<-unique(c(ban.sub$host, eid.sub$host, nhml.sub$host, oli.sub$host))
ncbi.orders<-tax_name(sci=need.orders, get="order", db="ncbi", messages=F)
ncbi.orders1<-tibble(ncbi.orders[,-1])
colnames(ncbi.orders1)<-c("host","hostOrder")

ban.sub1<-left_join(ban.sub, ncbi.orders1, by="host")
eid.sub1<-left_join(eid.sub, ncbi.orders1, by="host")
nhml.sub1<-left_join(nhml.sub, ncbi.orders1, by="host")
oli.sub1<-left_join(oli.sub, ncbi.orders1, by="host")

##### edit parasite data #####

#set common column name for parasite species
  #edit contents to use spaces rather than underscores and sentence case
ban.sub2<-ban.sub1 %>% mutate(parasite=Parasite.species)
bro.sub2<-bro.sub1 %>% mutate(parasite=gsub("_"," ",vVirusNameCorrected))
eid.sub2<-eid.sub1 %>% mutate(parasite=str_to_sentence(Cargo))
gmpd.sub2<-gmpd.sub1 %>% mutate(parasite=ParasiteCorrectedName)
nhml.sub2<-nhml.sub1 %>% mutate(parasite=Parasite)
oli.sub2<-oli.sub1 %>% mutate(parasite=gsub("_"," ",vVirusNameCorrected))

#add column and/or set common column name for parasite type
    #edit contents to lowercase
ban.sub3<-ban.sub2 %>% mutate(parType="helminth")
bro.sub3<-bro.sub2 %>% mutate(parType="virus")
eid.sub3<-eid.sub2 %>% mutate(parType=tolower(`Cargo classification`)) %>%
  filter(parType %in% c("arthropod", "bacteria", "helminth", "protozoa", "virus")) #drop extra parasite types
gmpd.sub3<-gmpd.sub2 %>% mutate(parType=tolower(ParType)) %>% 
  filter(parType!="fungus" & parType!="prion") #drop parasite types with too little data (12 fungi, 6 prions)
nhml.sub3<-nhml.sub2 %>% mutate(parType=paraType)
oli.sub3<-oli.sub2 %>% mutate(parType="virus")

#drop parasite records not identified to species
ban.sub4<-ban.sub3 %>%
  filter(grepl(" ",parasite)==T) %>% #remove parasite species without spaces, since it means they're only named to genus
  filter(grepl("\\.",parasite)==F) %>%  #remove parasite species with ".", since it means "spp." i.e. genus only
  filter(grepl(" .$",parasite)==F) #remove parasite species with "[space][single char][end line] (except viruses!)

bro.sub4<-bro.sub3 %>% 
  filter(grepl(" ",parasite)==T) %>% 
  filter(grepl("\\.",parasite)==F)

eid.not.vir<-eid.sub3 %>% filter(parType!="virus") %>% 
  filter(grepl(" ",parasite)==T) %>%
  filter(grepl("\\.",parasite)==F) %>%
  filter(grepl(" .$",parasite)==F)
eid.vir<-eid.sub3 %>% filter(parType=="virus") %>% 
  filter(grepl(" ",parasite)==T) %>%
  filter(grepl("\\.",parasite)==F)
eid.sub4<-rbind(eid.not.vir, eid.vir)

gmpd.sub4<-gmpd.sub3 %>%
  filter(HasBinomialName=="yes")

nhml.sub4<-nhml.sub3 %>% #other filtering steps have already been done
  filter(grepl(" larva",parasite)==F) %>% #drop larval records not identified to species
  filter(grepl(" larvae",parasite)==F)

oli.sub4<-oli.sub3 %>%
  filter(grepl(" ",parasite)==T) %>% 
  filter(grepl("\\.",parasite)==F)

##### merge to get unique host-parasite associations #####

#reduce dataframes to relevant columns and add column for data source
ban.sub5<-ban.sub4 %>% select(host, hostOrder, parasite, parType)
bro.sub5<-bro.sub4 %>% select(host, hostOrder, parasite, parType)
eid.sub5<-eid.sub4 %>% select(host, hostOrder, parasite, parType)
gmpd.sub5<-gmpd.sub4 %>% select(host, hostOrder, parasite, parType)
nhml.sub5<-nhml.sub4 %>% select(host, hostOrder, parasite, parType)
oli.sub5<-oli.sub4 %>% select(host, hostOrder, parasite, parType)

#combine source dataframes to get all unique host-parasite associations
hp.assoc<-rbind(ban.sub5, bro.sub5, eid.sub5, gmpd.sub5, nhml.sub5, oli.sub5) %>% #14142 rows total
  group_by(host) %>%
  distinct(parasite, .keep_all=T)

#add source data
hp.assoc %<>% ungroup() %>%
  mutate(banesh=ifelse(hp.assoc$host %in% ban.sub5$host & hp.assoc$parasite %in% ban.sub5$parasite, 1, 0),
         brook=ifelse(hp.assoc$host %in% bro.sub5$host & hp.assoc$parasite %in% bro.sub5$parasite, 1, 0),
         eid=ifelse(hp.assoc$host %in% eid.sub5$host & hp.assoc$parasite %in% eid.sub5$parasite, 1, 0),
         gmpd=ifelse(hp.assoc$host %in% gmpd.sub5$host & hp.assoc$parasite %in% gmpd.sub5$parasite, 1, 0),
         nhml=ifelse(hp.assoc$host %in% nhml.sub5$host & hp.assoc$parasite %in% nhml.sub5$parasite, 1, 0),
         oli=ifelse(hp.assoc$host %in% oli.sub5$host & hp.assoc$parasite %in% oli.sub5$parasite, 1, 0))

#check how many hosts we have overall that are in Healy's data
sum(unique(hp.assoc$host) %in% healy1$host)

write_csv(hp.assoc, "data/healy_host_par.csv")


