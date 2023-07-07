## RSelenium method to scrape citation counts from WOS

## Currently set up to work with the RSelenium docker image
## https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-docker.html 




library(RSelenium)
library(plyr); library(dplyr)

checkForServer()
startServer()

remDr <- remoteDriver()


remDr <- remoteDriver(remoteServerAddr = "localhost", port = 4445L, browserName = "chrome'")

remDr$open()

## replace the link below with a fresh advanced search link
remDr$navigate("https://apps.webofknowledge.com/WOS_AdvancedSearch_input.do?SID=1EwJfvSCrzymRmVr9i5&product=WOS&search_mode=AdvancedSearch")


# Loop to get citation data for each page of results, each iteration will save a txt file,
# I used selectorgadget to check the css ids, they might be different for you.
#Get the host names, loop through the web of science page using them
# Note: query below is a vector of host names to search through Web of Science.
# The default is to search NAMES + "parasite"

query <- c('mus musculus', 'gorilla gorilla')


ret <- vector()


for(i in 1:length(query)){
  
  # find where to enter text, and then enter text
  webElem <- remDr$findElement(using = 'css', value = "div textarea")
  webElem$sendKeysToElement(list(paste("TS=(", query[i], "AND 'parasite')", sep=' ')))
  
  # wait a few seconds
  Sys.sleep(10)
  
  # click on 'Search'
  webElem2 <- remDr$findElement(using = "css", "span input")
  webElem2$clickElement()
  
  # wait a few seconds
  Sys.sleep(10)
  
  # attempt to snag the citation count, if it doesn't exist, clear the element
  err <- tryCatch(unlist(remDr$findElement(using='css', '.newErrorHead')), error = function(e){})
  if(!is.null(err)){remDr$findElement(using = 'css', value = "div textarea")$clearElement()}
  
  # extract the citation counts from the html page and put count into vector element (ret)
  tab <- html(unlist(remDr$getCurrentUrl()))
  ret[i] <- trimws(tab %>% html_node(paste("#set_",i,"_div",sep='')) %>% html_text())
  
  print(i)
  print(ret[i])
  
}
