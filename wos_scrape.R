## Scrape web of science citations
## source code from Tad Dallas: https://gist.github.com/taddallas/cbf0589df1fa9051d774b2d282b6e635


# RSelenium method to scrape citation counts from WOS

**Updates**: RSelenium package has switched some of the code around, causing breaking changes to this old gist. Web of Science has also changed their layout, causing breaking changes. The current code is a hacky patch to allow for some WoS interactivity and data acquistion. 


### Currently set up to work with the RSelenium docker image

library(RSelenium)
library(plyr); library(dplyr)

driver <- rsDriver(port = 4545L, 
                   browser = "firefox", 
                   version = "latest"
)

 # https://resulumit.com/teaching/scrp_workshop.html#190

browser <- driver$client
server <- driver$server

browser$navigate("https://www.webofscience.com/wos/woscc/advanced-search")

query <- c('mus musculus', 'gorilla gorilla')

ret <- vector()

for(i in 1:length(query)){
  
  # find where to enter text, and then enter text
  webElem <- browser$findElement(using = 'css', value = "textarea.search-criteria-input")
  webElem$sendKeysToElement(list(paste("TS=('", query[i], "' AND 'parasite')", sep=' '), "\uE007")) 
  
  Sys.sleep(5)
  browser$goBack()
  Sys.sleep(5)
  
  webElem <- browser$findElement(using = 'css', value = "textarea.search-criteria-input")
  webElem$clearElement()
  
  print(i)
  
}

webElem <- browser$findElement(using = 'css', value = "app-history-entries-list.ng-star-inserted")
ret <- webElem$getElementText() 

driver[['server']]$stop()

