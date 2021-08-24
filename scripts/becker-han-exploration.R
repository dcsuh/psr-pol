## Daniel Suh
## August 24, 2021
## Playing around with data from Becker & Han 2020 GBE

library(tidyverse)
library(magrittr)
library(here)

data <- read_csv(here("data/Becker Han 2020 GEB_gpf & BRT data.csv"))

dat1 <- data %>% select(X1,genus,common,inames,bcomp) %>% distinct()

