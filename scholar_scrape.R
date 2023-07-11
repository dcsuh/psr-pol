## web scrape google scholar for citation counts based off search query

library(here)
library(tidyverse)
library(magrittr)
library(rvest)


names <- read_csv(here("data/host_names.csv"))

names$results <- NA

#names$species <- tolower(names$species) 
names$species <- str_replace(string = names$species, pattern = " ", replacement = "+")


for (i in 1:nrow(names)){
  tmp <- pull(names[i,1])
  
  link <- str_replace(string="https://scholar.google.com/scholar?hl=en&as_sdt=0%2C11&q=REPLACEME&btnG=", pattern="REPLACEME", replacement=tmp)
  
  page <- read_html(link)
  
  summaries_css <- page %>% html_elements(css = ".gs_ab_mdw")
  
  text <- capture.output(print(summaries_css[2]))
  
  line <- stringr::str_extract(text, "About.+results")
  
  line_a <- str_replace(line, ",", "")
  
  line_b <- line_a[2]
  
  number <- as.numeric(str_extract(line_b, "[:digit:]+"))

  names$results[i] <- number
}

write_csv(names, here("data/google_refs.csv"))
