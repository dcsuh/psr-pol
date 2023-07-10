## web scrape google scholar for citation counts based off search query

library(here)
library(tidyverse)
library(magrittr)
library(rvest)


names <- tibble(species = c("Mustela nigripes", "Mus muscula"), results = NA)

names$species <- tolower(names$species) 
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

link <- "https://scholar.google.com/scholar?hl=en&as_sdt=0%2C11&q=mustela+nigripes&btnG="

page <- read_html(link)

summaries_css <- page %>% html_elements(css = ".gs_ab_mdw")

text <- capture.output(print(summaries_css[2]))

line <- stringr::str_extract(text, "About.+results")

line <- str_replace(line, ",", "")

line <- line[2]

number <- as.numeric(str_extract(line, "[:digit:]+"))