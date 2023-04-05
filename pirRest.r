library(httr)
library(jsonlite)

api_key <- "XXX"
url <- "https://research.pirche.com"

token_response <- POST(paste0(url, "/portal/rest/oauth/token"), body = list('grant_type'= 'authorization_code', 'client_id'= 'api-client', 'code'= api_key), encode = "form")
access_token <- paste("bearer", content(token_response)$access_token)
print(access_token)

patient = list(id = "123P", population = "NMDP EUR haplotypes (2007)", glString = "A*01+A*29^B*08+B*44^C*07+C*16^DRB1*03+DRB1*07^DQB1*02+DQB1*02")
donors = list(list(id = "987D", population = "NMDP EUR haplotypes (2007)", glString = "A*02+A*29^B*08+B*44^C*07+C*16^DRB1*03+DRB1*07^DQB1*02+DQB1*02"))

pirche_response <- POST(paste0(url, "/pirche/rest/sot/api/match"), body = list(patient = patient, donors = donors), encode = "json", add_headers(Authorization = access_token))
pirche_data <- content(pirche_response)

print(pirche_data$pircheII$'987D'$sum)
