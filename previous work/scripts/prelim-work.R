#Preliminary work: data organization and basic visualizations
#Annakate Schatz, David Vasquez, John Vinson
#11/11/20

##### load #####

#packages
library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(tidyverse)

#Healy et al life history data
healy<-read_csv("data/healy_axis_analysis.csv")

#collected host-parasite associations
hp<-read_csv("data/healy_host_par.csv")

#function used later to replace NaN with NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

##### summarize parasite species richness for each host species #####

#create dataframe with PSR per host
hp.psr<-data.frame(table(hp$host))
colnames(hp.psr)<-c("host", "PSR")

#add PSR proportion for each parasite type
hosts<-hp.psr$host
pars<-c("arthropod", "bacteria", "helminth", "protozoa", "virus")
hp.psr1<-hp.psr %>% mutate(prop.arth=NA, prop.bact=NA, prop.helm=NA, prop.prot=NA, prop.vir=NA)
for (i in 1:length(hosts)) {
  hp.sub<-hp[hp$host==hosts[i],]
  temp<-c()
  for (j in 1:length(pars)) {
    temp[j]<-nrow(hp.sub[hp.sub$parType==pars[j],])/nrow(hp.sub)
  }
  hp.psr1[hp.psr1$host==hosts[i],3:7]<-as.numeric(temp)
}

##### merge host PSR with Healy data #####

### Healy prep ###

#add host column for merging and drop extra columns
healy1<-healy %>% mutate(host=gsub("_"," ",species), matrix_size=NULL, surv_95=NULL, surv_99=NULL, surv_sd=NULL)

#to collapse multiple populations into one host record, we will take the mean of numeric lh traits
#non-numeric columns have the same entries across populations so we can just carry over that data

#set column types to ensure means work for numeric lh traits
healy1[,c(3:12,14:15,17)]<-as.numeric(healy1[,c(3:12,14:15,17)])
healy1[,c(1:2,13,16,18:20)]<-as.character(healy1[,c(1:2,13,16,18:20)])

#take means of numeric lh traits and carry over non-numeric lh traits
hosts<-unique(healy1$host)
healy.means<-NULL
for (i in 1:length(hosts)) {
  healy.sub<-data.frame(healy1[healy1$host==hosts[i],])
  healy.means.temp<-healy.sub[1,1:2]
  for (j in 3:ncol(healy.sub)) {
    if (class(healy.sub[1,j])=="numeric") {
      healy.means.temp<-data.frame(healy.means.temp, mean(healy.sub[,j],na.rm=T))
    }
    if (class(healy.sub[1,j])=="character") {
      healy.means.temp<-data.frame(healy.means.temp, healy.sub[1,j])
    }
  }
  healy.means<-rbind(healy.means, healy.means.temp)
}

#add column names and replace NaN values with NA
colnames(healy.means)<-colnames(healy1)
healy.means[is.nan(healy.means)]<-NA

### merge ###
hp.healy<-left_join(hp.psr1, healy.means, by="host")

##### preliminary visualizations #####

#histogram of PSR
hist(hp.healy$PSR, breaks=30)

#parasite type proportions vs PSR
hp.healy1<-hp.healy %>% 
  pivot_longer(cols=c(3:7), names_to="type", values_to="prop")
ggplot(data=hp.healy1) +
  geom_point(aes(x=PSR, y=prop, color=type)) +
  xlim(0,400)

ggplot(data=hp.healy1, aes(x=type, y=prop)) +
  geom_boxplot() +
  facet_wrap(~met_type) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

#PSR vs numeric life history traits
hp.healy2<-hp.healy %>% 
  pivot_longer(cols=c(10:19, 21:22, 24), names_to="lh_trait_num") %>% 
  filter(PSR<400)
ggplot(data=hp.healy2) +
  geom_point(aes(x=value, y=PSR)) +
  # geom_smooth(aes(x=value, y=PSR)) +
  facet_wrap(~lh_trait_num, scales="free_x") #+
  # ylim(0,400)

#possibly try this against parasite types (ecto vs endo, macro vs micro)

#PSR vs non-numeric life history traits
hp.healy3<-hp.healy %>% 
  pivot_longer(cols=c(9, 20, 23, 25), names_to="lh_trait_char")
ggplot(data=hp.healy3) +
  geom_boxplot(aes(x=value, y=log(PSR))) +
  facet_wrap(~lh_trait_char, scales="free_x") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) #+
  # ylim(0,400)

##### save data #####

save(hp.psr1, healy.means, hp.healy, file="data/prelim-work.Rda")
