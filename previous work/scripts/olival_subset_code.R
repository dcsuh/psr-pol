#John Vinson

require(dplyr)

healy.data = read.csv("data/healy_axis_analysis.csv", header = TRUE)

olival.data = read.csv("data/source/olival_full.csv", header = TRUE)

hosts.inter = intersect(healy.data$species, olival.data$hHostNameFinal)

olival.data.sub = filter(olival.data, hHostNameFinal %in% hosts.inter)

write.csv(olival.data.sub, "data/olival_subset.csv")
