## Daniel Suh
## August 24, 2021
## Playing around with data from Becker & Han 2020 GBE
## https://doi.org/10.1111/geb.13256

library(tidyverse)
library(magrittr)
library(here)

data <- read_csv(here("data/Becker Han 2020 GEB_gpf & BRT data.csv"))
#trait data for 4691 avian species
#competence data for 183 species from meta-analysis 
#(only used studies with xenodiagnostic testing or larvae from wild birds that tested positive)
#183 sampled for bbsl competence
##92 not competent
##91 competent
#4508 unsampled for competence
#study effort approximated using Web of Science citation counts per species

comp_sampled <- data %>% filter(!is.na(.$bcomp))
