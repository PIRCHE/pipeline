# -------------------------------------------------------------------------
# Script: pirche-sot-score-requests.R
# Author: Hilary Mehler
# Event: 2026 ASHI Clinical Informatics Workshop
# Date: February 24, 2026
#
#    1. Obtain PIRCHE-T2 and PIRCHE-B scores using PIRCHE REST API.
#    2. Open TxPredictor in browser with HLA typings pre-filled.
#
#    NOTE: This script accepts molecular HLA typings only
#
#    See R Notebook for a complete walkthrough of this script (pirche-sot-wkshp.Rmd)
#
# License: MIT (See main repository LICENSE file for details)
# -------------------------------------------------------------------------#


library(dplyr)
library(readr)
library(httr2)
library(immunogenetr)


##### CONFIGURATION #####
input_file <- "test_data.csv" # HLA typings grouped by case
demo_case_id <- "2"
client_id <- "api-client"
pipeline <- "Frost-1.1_IMGT-3.54"
api_key <- "___API-KEY___" # Please reach out to support@pirche.com - We're happy to help with your integration!

# GET Access Token
cat("Authenticating with PIRCHE Server...\n")

req_auth <- request("https://research.pirche.com/portal/rest/oauth/token") %>% 
  req_headers("Accept" = "application/json") %>% 
  req_body_form(
    code       = api_key,
    grant_type = "authorization_code",
    client_id  = client_id
  )

resp_auth <- req_perform(req_auth)
auth_token <- resp_body_json(resp_auth)$access_token

cat("...Access Token acquired successfully.\n")

##### Process HLA Cases #####
raw_data <- read_csv(input_file)  #import csv
case_data <- raw_data %>% filter(case_id == demo_case_id) #filter by case_id

if(nrow(case_data) == 0) stop("Aborted. Case ID not found.") #stop if no data

print(paste("Data Loaded for: ", demo_case_id))
print(case_data) # Print raw typings

## Create GL Strings using immunogenetr ##
# Function collects input data for individual (ID, HLA GL String)
build_hla_block <- function(row_data) {
  
  # Define columns
  hla_cols <- c("A1", "A2", "B1", "B2", "C1", "C2", 
                "DRB1_1", "DRB1_2", "DRB345_1", "DRB345_2", 
                "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", 
                "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2")
  
  # Filter to columns containing alleles
  active_cols <- hla_cols[
    !is.na(row_data[hla_cols]) & 
      row_data[hla_cols] != ""
  ]
  
  # Ensure that immunogenetr identifies low res molecular 
  ## (seems to look for a colon; if not found string is converted to serology)
  row_mod <- row_data
  for(col in active_cols) {
    val <- row_mod[[col]]
    if(!grepl(":", val)) {
      row_mod[[col]] <- paste0(val, ":99")
    }
  }
  
  # Create GL String
  gl_list <- immunogenetr::HLA_columns_to_GLstring(
    data = row_mod,
    HLA_typing_columns = all_of(active_cols)
  )
  
  clean_gl <- gsub(":99", "", gl_list[[1]])
  
  list(
    id       = row_data$person_id,
    glString = clean_gl
  )
}

# Build data objects for patient & donors using build_hla_block function
patient_data <- build_hla_block(case_data %>% filter(person_type == "patient"))

donor_rows  <- case_data %>% filter(person_type == "donor")
donor_data_list <- list()

for(i in 1:nrow(donor_rows)) {
  donor_data_list[[i]] <- build_hla_block(donor_rows[i, ])
}

# Final Payload (complete data / JSON object)
payload <- list(
  databaseVersion = pipeline,
  patient = patient_data,
  donors  = donor_data_list,
  hlaValueFormatVersion = "v2"
)

print(payload)

##### Call REST API - SOT and SNOW endpoints (PIRCHE-T2, PIRCHE-B) #####
base_url  <- "https://research.pirche.com/pirche/rest/sot"
endpoint_sot  <- "api/match"       # PIRCHE-T2
endpoint_snow <- "snow/api/match"  # PIRCHE-B

# Function takes endpoint, case data & calls REST API using httr2
call_pirche_api <- function(endpoint_path, payload) {
  full_url <- paste0(base_url, "/", endpoint_path)
  request(full_url) %>% 
    req_headers(
      "Authorization" = paste("Bearer", auth_token),
      "Content-Type"  = "application/json;charset=utf-8",
      "Accept"        = "application/json"
    ) %>% 
    req_body_json(payload) %>% 
    req_perform() %>% 
    resp_body_json()
}

cat("\n--- Calling PIRCHE API ---\n")
cat("\nCalling SOT Endpoint (PIRCHE-T2 scores)...\n")
json_sot  <- call_pirche_api(endpoint_sot, payload)

cat("Calling SNOW Endpoint (PIRCHE-B scores)...\n")
json_snow <- call_pirche_api(endpoint_snow, payload)

##### Format and Print Results #####
summary_rows  <- list() 
pirche_details   <- list()
snow_details  <- list()

# Function formats numbers to 2 decimals, keeping NA/NULL distinct
fmt <- function(x) {
  if(is.null(x)) return("NULL")
  if(is.na(x))   return("NA")
  return(sprintf("%.2f", as.numeric(x)))
}


donor_ids <- names(json_sot$pircheII) 

for(id in donor_ids) {
  
  pirche_dat  <- json_sot$pircheII[[id]]
  snow_dat <- json_snow$snow[[id]]
  
  # Summary Table - Total Scores
  summary_rows[[length(summary_rows)+1]] <- tibble(
    Case      = demo_case_id,
    Donor_ID  = id,
    PIRCHE_Sum   = fmt(pirche_dat$sum),
    SNOW_Sum  = fmt(snow_dat$sum)
  )
  
  # Detailed Table - PIRCHE-T2
  row_sot <- as_tibble(pirche_dat) %>%
    mutate(across(everything(), fmt)) %>%
    mutate(Donor_ID = id, .before = 1) %>%
    rename(Total = sum)
  
  pirche_details[[length(pirche_details)+1]] <- row_sot
  
  # Detailed Table - PIRCHE-B
  row_snow <- as_tibble(snow_dat) %>%
    mutate(across(everything(), fmt)) %>%
    mutate(Donor_ID = id, .before = 1) %>%
    rename(Total = sum)
  
  snow_details[[length(snow_details)+1]] <- row_snow
}


# Print Result Tables  
cat("\n--- Summary - Total Scores ---\n")
print(bind_rows(summary_rows))
  
df_sot <- bind_rows(pirche_details)
locus_cols <- sort(setdiff(names(df_sot), c("Donor_ID", "Total"))) 
df_sot <- df_sot %>% select(Donor_ID, Total, all_of(locus_cols)) 
  
cat("\n--- Detailed Output: PIRCHE-T2 (sot endpoint) ---\n")
print(df_sot)

df_snow <- bind_rows(snow_details)
locus_cols_snow <- sort(setdiff(names(df_snow), c("Donor_ID", "Total")))
df_snow <- df_snow %>% select(Donor_ID, Total, all_of(locus_cols_snow))  
  
cat("\n--- Detailed Output: PIRCHE-B (snow endpoint) ---\n")
print(df_snow)


##### Function sends case information to web application #####

# (script narrows to first two donors to limit the url length...
# ... user may cache donor information to include more donors) 
open_tx_predictor <- function(payload, auth_token) {
  
  weblink_base_url <- "https://research.pirche.com/pirche/#/pirche/landing"
  
  # Function prepares string for URL
  encode_url <- function(gl_string) {
    clean <- gsub("HLA-", "", gl_string) # Remove 'HLA-'
    return(curl::curl_escape(clean)) # Percent-encode special chars using curl_escape from httr2
  }
  
  # Get Patient
  pat_gl_encoded <- encode_url(payload$patient$glString)
  pat_id_encoded <- curl::curl_escape(payload$patient$id)
  
  # Get Donors (limit to 2)
  target_donors <- head(payload$donors, 2)
  
  donor_gls_encoded <- c()
  donor_ids_encoded <- c()
  
  for(i in 1:length(target_donors)) {
    current_donor <- target_donors[[i]]
    
    clean_gl <- encode_url(current_donor$glString)
    donor_gls_encoded[i] <- clean_gl
    
    clean_id <- curl::curl_escape(current_donor$id)
    donor_ids_encoded[i] <- clean_id
  }
  
  # Join donors with semicolon
  final_donor_gl_string <- paste(donor_gls_encoded, collapse = ";")
  final_donor_id_string <- paste(donor_ids_encoded, collapse = ";")
  
  # Construct URL
  full_url <- paste0(
    weblink_base_url,
    "?patientGLString=", pat_gl_encoded,
    "&patientId=",       pat_id_encoded,
    "&donorGLString=",   final_donor_gl_string,
    "&donorId=",         final_donor_id_string,
    "&authorization=",   auth_token 
  )
  
  cat("\n--- Opening TxPredictor in Default Browser ---\n")
  print(full_url)
  browseURL(full_url)
}

# Call function to open TxPredictor in browser
open_tx_predictor(payload = payload, auth_token = auth_token)

