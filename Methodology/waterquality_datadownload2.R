#Water quality data download for Needs Assessment and Aquifer Risk Map

#Methodology write-up: https://gispublic.waterboards.ca.gov/portal/home/item.html?id=70feb9f4b00f4b3384a9a0bf89f9f18a

#Updated 5/25/2021 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

#Components:
#1. create depth filter to identify state small and domestic well depth water quality 
#2. download water quality data from GAMA GIS from
  #Division of Drinking Water
  #Department of Water Resources
  #USGS NWIS
  #GAMA projects
  #ILRP, LOCALGW, and WB_DAIRY datasets from Geotracker
#3. standardize and clean water quality data (standardize chemicals, handle outliers and non-detects)

#load required libraries
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('rgdal')) install.packages('rgdal'); library('rgdal')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('janitor')) install.packages('janitor'); library('janitor')
if (!require('httr')) install.packages('httr'); library('httr')
if (!require('data.table')) install.packages('data.table'); library('data.table')

#define special functions
RLfun <- function(rl, mcl) {  #always make adjusted RL-results smaller than Reporting Limit
  ifelse(rl/sqrt(2) >= mcl/2, mcl/2, rl/sqrt(2))
}
my_max <- function(x) {  #grab maxes even when there are NAs in the list                                          
  ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
}
my_min <- function(x) {  #grab mins even when there are NAs in the list                                          
  ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
}
st_drop_geometry <- function(x) {
  if(inherits(x,"sf")) {
    x <- st_set_geometry(x, NULL)
    class(x) <- 'data.frame'
  }
  return(x)
}

#define chemicals of interest
chem_list <- read.csv("Datasets/NA_allchems.csv", stringsAsFactors = F) %>% 
  filter(CCT == "MCL" | chemVVL %in% c("CU", "PB", "NNSM", "CR6")) %>%
  filter(!chemVVL %in% c("ASBESTOS", "COLIFORM", "FCOLIFORM", "RN-222"))

#1. create depth filter ####
#import groundwater unit boundaries
gu_link <- "https://pubs.usgs.gov/ds/796/downloads/ds796_GIS.zip"
temp <- tempfile()
temp2 <- tempfile()
download.file(gu_link, temp)
unzip(zipfile = temp, exdir = temp2)
GU <- st_read(paste0(temp2, "/ds796_GIS/CA_Groundwater_Units.shp")) %>% st_transform(crs = 4326)
GU_table <- st_drop_geometry(GU)
unlink(temp)
unlink(temp2)

#import MTRS boundaries and domestic/public well depths per section
url <- list(hostname = "gis.water.ca.gov/arcgis/rest/services",
            scheme = "https",
            path = "Environment/i07_WellCompletionReports/FeatureServer/1/query",
            query = list(where = "1=1",
                         outFields = "*",
                         returnGeometry = "true",
                         f = "geojson")) %>%
  setattr("class","url")
request <- build_url(url)
MTRS <- st_read(request)

#join mtrs to groundwater units and remove duplicate sections (that are duplicated across county line splits of PLSS sections)
#define columns that must be numeric
numeric_cols <- c("DomWellCount", "DomWellDepthAvg", "DomWellDepthMin", "DomWellDepthMax",
                  "PubWellCount", "PubWellDepthAvg", "PubWellDepthMin", "PubWellDepthMax")

#identify columns to be deleted
delete <- c("COUNTY_CD", "WCRFolderLink")

#join mtrs to groundwater units, then drop geometry
MTRS_GU <- st_join(MTRS, GU) %>% st_drop_geometry()

#remove duplicated mtrs items
MTRS_GU <- MTRS_GU[, !(colnames(MTRS_GU) %in% delete), drop=FALSE] %>% distinct() %>%
  mutate_at(numeric_cols, as.numeric)

# remove depth zeros (sections without any wells with depth information) (some sections have domestic wells but no depth data associated
#with them, so I am filtering on depth, not count)
clean_d <- MTRS_GU %>% filter(!DomWellDepthMax == 0) #for domestic wells
clean_p <- MTRS_GU %>% filter(!PubWellDepthMax == 0) #for public wells

#summarize by basin and display depth data:
#calculate average maximum and minimum well depths, and standard deviation. If there 2 or less sections
#in the basin, the basin max/min domestic is the average max domestic+150 ft/average min domestic-150 ft.
depth_dom <- clean_d %>% dplyr::group_by(GU_ID) %>% 
  dplyr::summarize(dom_section_n = n(),
                   dom_well_n = sum(DomWellCount), 
                   avgmax_dom = mean(DomWellDepthMax), 
                   avgmaxsd_dom = sd(DomWellDepthMax),
                   avgmin_dom = mean(DomWellDepthMin), 
                   avgminsd_dom = sd(DomWellDepthMin)) %>% 
  dplyr::mutate(maxdomestic = ifelse(dom_section_n <= 2, avgmax_dom + 3*150, avgmax_dom + 3*avgmaxsd_dom),
                mindomestic = ifelse(dom_section_n <= 2, avgmin_dom - 3*150, avgmin_dom - 3*avgminsd_dom))

#calculate average maximum public well depths. Same as above - if 2 or less sections add 150 to average max.
depth_pub <- clean_p %>% dplyr::group_by(GU_ID) %>% 
  dplyr::summarize(pub_section_n = n(),
                   pub_well_n = sum(PubWellCount),
                   avgmax_pub = mean(PubWellDepthMax), 
                   avgmaxsd_pub = sd(PubWellDepthMax)) %>% 
  dplyr::mutate(maxpublic = ifelse(pub_section_n <= 2, avgmax_pub + 3*150, avgmax_pub + 3*avgmaxsd_pub))

#merge public and domestic columns together by groundwater unit, calculate % difference of domestic/public bottom
#check if basins are same (within 10%) or different (not within 10%)
#check if public wells are shallower than domestic wells
alldata <- full_join(depth_dom, depth_pub, by = c("GU_ID")) %>% 
  mutate(percent_max = (abs(maxpublic - maxdomestic)/((maxpublic + maxdomestic)/2)*100),
         comparison = ifelse(is.na(avgmax_dom), NA,
                             ifelse(is.na(avgmax_pub), NA,
                                    ifelse(percent_max <= 10, "same", "different"))),
         public_shallow = ifelse(maxpublic < maxdomestic, "shallower", "deeper"),
         use_public_as_dom = ifelse(is.na(public_shallow) & !is.na(maxdomestic), "yes", 
                                    ifelse(public_shallow == "shallower", "yes", 
                                           ifelse(comparison == "same", "yes", "no"))))
rm(clean_d, clean_p, depth_dom, depth_pub)

#add full list of units (including ones that have no groundwater depth filter)
GU_depths <- left_join(GU_table, alldata)
rm(GU_table, alldata)

#write file to working directory
write.table(GU_depths, paste0("well_depth_gwunits", Sys.Date(), ".txt"), sep = "\t", row.names = F)

#2. download water quality data ####
#identify wells that are within domestic depth filter
#download location data
temp <- tempfile()
download.file('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_location_construction_v2.zip', temp)
loc <- readr::read_tsv(unz(temp, 'gama_location_construction_v2.txt'), guess_max = 100000) %>% clean_names()
unlink(temp)

#convert well locations to sf object and join with GU and MTRS locations
wells_int <- st_as_sf(loc, coords = c('gm_longitude', 'gm_latitude'), crs = 4326) %>%
  st_join(., left = TRUE, GU["GU_ID"]) %>% 
  st_join(., left = TRUE, MTRS["MTRS"]) %>% 
  st_drop_geometry() %>%
  select("gm_well_id", "gm_well_category", "gm_well_depth_ft",
         "gm_top_depth_of_screen_ft", "gm_bottom_depth_of_screen_ft", "GU_ID", "gm_dataset_name", "MTRS") %>%
  left_join(., select(GU_depths, GU_ID, avgmax_dom, maxdomestic, mindomestic, use_public_as_dom), by = c("GU_ID")) %>%
  dplyr::filter(!is.na(GU_ID))

#define numeric columns and replace 0 with NULL
numeric_cols <- c("gm_well_depth_ft", "gm_top_depth_of_screen_ft", "gm_bottom_depth_of_screen_ft")
wells_int <- wells_int %>% mutate_at(numeric_cols, as.numeric)
wells_int$gm_well_depth_ft[wells_int$gm_well_depth_ft == 0] <- NA                           #replace 0s with NA
wells_int$gm_top_depth_of_screen_ft[wells_int$gm_top_depth_of_screen_ft == 0] <- NA
wells_int$gm_bottom_depth_of_screen_ft[wells_int$gm_bottom_depth_of_screen_ft == 0] <- NA
wells_int <- mutate(wells_int, depth = ifelse(is.na(wells_int$gm_well_depth_ft),          #determine numeric depth, if possible
                                       ifelse(is.na(wells_int$gm_top_depth_of_screen_ft) | is.na(wells_int$gm_bottom_depth_of_screen_ft), NA, 
                                                     wells_int$gm_bottom_depth_of_screen_ft), 
                                                     wells_int$gm_well_depth_ft))
#filter on 1) if domestic, 2) if numeric, apply numeric filter OR if depth unknown apply unit use
allwells <- wells_int %>% mutate(domdepth = ifelse(gm_well_category == "Domestic" | gm_dataset_name == "GAMA_DOM", "yes",
                                                   ifelse(gm_dataset_name == "WB_CLEANUP", "no",
                                                   ifelse(is.na(depth),
                                                          ifelse(is.na(use_public_as_dom) | use_public_as_dom == "yes", "yes", "no"),
                                                          ifelse(!is.na(maxdomestic) & depth <= maxdomestic & depth > mindomestic, "yes", "no")))))
rm(wells_int)
lat_long <- loc %>% select(gm_well_id, gm_longitude, gm_latitude)
allwells <- left_join(allwells, lat_long)
write.table(allwells, paste0("allwells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
domwells <- allwells %>% filter(domdepth == "yes")

#retrieve water quality data
#define column types
col_sty <- cols(
  GM_DATASET_NAME = col_character(),
  GM_WELL_CATEGORY = col_character(),
  GM_DATA_SOURCE = col_character(),
  GM_WELL_ID = col_character(),
  GM_CHEMICAL_VVL = col_character(),
  GM_CHEMICAL_NAME = col_character(),
  GM_RESULT_MODIFIER = col_character(),
  GM_RESULT = col_double(),
  GM_RESULT_UNITS = col_character(),
  GM_SAMP_COLLECTION_DATE = col_character(),
  GM_REPORTING_LIMIT = col_double(),
  GM_LATITUDE = col_double(),
  GM_LONGITUDE = col_double(),
  GM_WELL_DEPTH_FT = col_double(),
  GM_TOP_DEPTH_OF_SCREEN_FT = col_double(),
  GM_BOTTOM_DEPTH_OF_SCREEN_FT = col_double(),
  GM_CAS_NUMBER = col_character(),
  GM_ALTWELL_ID1 = col_character(),
  GM_ALTWELL_ID2 = col_character(),
  GM_ALTWELL_ID3 = col_character(),
  SRC_CHEMICAL = col_character(),
  SRC_RESULT_MODIFIER = col_character(),
  SRC_RESULT = col_double(),
  SRC_RESULT_UNITS = col_character(),
  SRC_SAMP_COLLECTION_DATE = col_character(),
  SRC_SAMP_COLLECTION_TIME = col_character(),
  SRC_REPORTING_LIMIT = col_double(),
  SRC_ANALYTICAL_METHOD = col_character(),
  SRC_LAB_NOTE = col_character(),
  SRC_LATITUDE = col_double(),
  SRC_LONGITUDE = col_double(),
  SRC_DATUM = col_character(),
  SRC_WELL_DEPTH_FT = col_double(),
  SRC_TOP_DEPTH_OF_SCREEN_FT = col_double(),
  SRC_BOTTOM_DEPTH_OF_SCREEN_FT = col_double()
)
#define function to download data
dwd_wq_tbl <- function(link, filename, col_str, dom_wells, mychems) {
  temp <- tempfile()
  download.file(link, temp)
  mytable <- readr::read_tsv(unz(temp, filename), quote = "", col_types = col_str) %>% clean_names()
  unlink(temp)
  mytable$gm_samp_collection_date <- lubridate::mdy(mytable$gm_samp_collection_date)
  mytable <- mytable %>% filter(gm_samp_collection_date >= ymd("2001-05-25"),
                                gm_well_id %in% dom_wells$gm_well_id,
                                gm_chemical_vvl %in% mychems$chemVVL)
  return(mytable)
}

#download various datasets
ddw <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_ddw_statewide_v2.zip',
                  'gama_ddw_statewide_v2.txt', col_sty, domwells, chem_list)
dwr <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_dwr_statewide_v2.zip',
                  'gama_dwr_statewide_v2.txt', col_sty, domwells, chem_list)
gamadom <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_gama_dom_statewide_v2.zip',
                      'gama_gama_dom_statewide_v2.txt', col_sty, domwells, chem_list)
gamausgs <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_gama_usgs_statewide_v2.zip',
                       'gama_gama_usgs_statewide_v2.txt', col_sty, domwells, chem_list)
localgw <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_localgw_statewide_v2.zip',
                      'gama_localgw_statewide_v2.txt', col_sty, domwells, chem_list)
wbdairy <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_wb_dairy_statewide_v2.zip',
                      'gama_wb_dairy_statewide_v2.txt', col_sty, domwells, chem_list)
usgsnwis <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_usgs_nwis_statewide_v2.zip',
                       'gama_usgs_nwis_statewide_v2.txt', col_sty, domwells, chem_list)
wbilrp <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_wb_ilrp_statewide_v2.zip',
                     'gama_wb_ilrp_statewide_v2.txt', col_sty, domwells, chem_list)

#compile datasets
domwqdata <- rbind(ddw, dwr, gamadom, gamausgs, localgw, wbdairy, usgsnwis, wbilrp)
rm(ddw, dwr, gamadom, gamausgs, localgw, wbdairy, usgsnwis, wbilrp)

#remove duplicated data (duplicated between GAMA USGS and USGS NWIS datasets)
duplicated_wells <- read.csv("Datasets/NWIS_dupes2.csv", stringsAsFactors = F)
domwqdata <- domwqdata %>% filter(!gm_well_id == "Vanadium Study",
                                  !gm_well_id %in% duplicated_wells$Duplicate_NWIS,
                                  !gm_well_id == "USGS-415953121292901")
domwells_wq <- domwqdata %>% select(gm_well_id) %>% distinct() %>% left_join(., domwells)

#save water quality data and well information
write.table(domwqdata, paste0("domwqdata_", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
write.table(domwells_wq, paste0("domwells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

#3. standardize and clean data ####
domwqdata <- read.table("domwqdata_2021-05-25.txt", stringsAsFactors = F, header = T)
#downsize data for cleaning
domwqdata <- domwqdata %>% select(gm_well_id, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date,
                                  gm_reporting_limit)
#read dates as dates
domwqdata$gm_samp_collection_date <- lubridate::ymd(domwqdata$gm_samp_collection_date)

#convert Nitrate/Nitrite combined to Nitrate
#all nitrogen measurements
nitrogen <- domwqdata %>% dplyr::filter(gm_chemical_vvl %in% c("NO3N", "NO2", "NO3NO2N")) %>% dplyr::mutate(uID = paste0(gm_well_id, gm_samp_collection_date))
#NO3NO2N measurements
combined <- nitrogen %>% dplyr::filter(gm_chemical_vvl == "NO3NO2N") %>% dplyr::select(uID, gm_samp_collection_date) %>% distinct() 
#nitrate (NO3N) measurements
nitrate <- nitrogen %>% dplyr::filter(gm_chemical_vvl == "NO3N") %>% dplyr::select(uID) %>% distinct()
#nitrite (NO2) measurements
nitrite <- nitrogen %>% dplyr::filter(gm_chemical_vvl == "NO2") %>% dplyr::select(uID, gm_samp_collection_date) %>% distinct()
#measurements with NO3NO2N but without NO3N (to be substituted or converted)
nonitrate <- setdiff(combined$uID, nitrate$uID)
#substitution measurements (meas with NO3NO2N but without NO2; just use NO3NO2N as NO3N)
sub_data <- nitrogen %>% dplyr::filter(uID %in% setdiff(nonitrate, nitrite$uID)) %>% 
  select(gm_well_id, gm_result, gm_chemical_vvl, gm_samp_collection_date) %>% mutate(newCHEMICAL = "NO3N")
#conversion measurements (meas with NO3NO2N but with NO2; subtract NO2 from No3NO2N to get NO3N)
con_data <- nitrogen %>% dplyr::filter(uID %in% setdiff(nonitrate, setdiff(nonitrate, nitrite$uID))) %>%  
  select(gm_well_id, gm_result, gm_chemical_vvl, gm_samp_collection_date) %>% 
  pivot_wider(names_from = gm_chemical_vvl, values_from = gm_result) %>% 
  mutate(gm_result = NO3NO2N - NO2, newCHEMICAL = "NO3N", gm_chemical_vvl = "NO3NO2N")
con_data$NO2 <- NULL
con_data$NO3NO2N <- NULL
#recombine adjusted nitrogen meas. with full data
n_add <- rbind(con_data, sub_data)         
domwqdata <- left_join(domwqdata, n_add, by = c("gm_well_id", "gm_samp_collection_date", "gm_chemical_vvl")) %>%  
  mutate(adjCHEMICAL = ifelse(is.na(newCHEMICAL), gm_chemical_vvl, newCHEMICAL),
         gm_result = ifelse(is.na(newCHEMICAL), gm_result.x, gm_result.y))
domwqdata$newCHEMICAL <- NULL
domwqdata$gm_chemical_vvl <- NULL
colnames(domwqdata)[colnames(domwqdata)=="adjCHEMICAL"] <- "gm_chemical_vvl"
domwqdata$gm_result.x <- NULL
domwqdata$gm_result.y <- NULL
domwqdata <- domwqdata2 %>% dplyr::filter(!gm_chemical_vvl == "NO3NO2N") #delete duplicated rows
rm(nitrate, nitrogen, nitrite, n_add, nonitrate, sub_data, combined, con_wide, con_data)

#standardize Qualifiers to either "<" or "="
domwqdata$gm_result_modifier[is.na(domwqdata$gm_result_modifier)] <- "="                        
domwqdata <- domwqdata %>% dplyr::filter(!gm_result_modifier %in% c("S", "V", "M")) %>%
  dplyr::mutate(gm_result_modifier = ifelse(gm_result_modifier == "<" | gm_result_modifier == "ND", "<", "="))

#standardize Reporting Limits/Results (fix Null results, results == 0, and missing reporting limits)
domwqdata$gm_reporting_limit <- as.numeric(domwqdata$gm_reporting_limit)
finishedchems <- c()         #hold chemicals as they are cleaned
additions <- data.frame()    #holds rows with "newRL"s
chemicals <- unique(domwqdata$gm_chemical_vvl)     #list of chemicals that need to be fixed
domwqdata <- dplyr::arrange(domwqdata, gm_samp_collection_date)
chemicals <- sort(chemicals)
for (ch in chemicals) {
  uRL <- domwqdata %>% dplyr::filter(is.na(gm_reporting_limit) | gm_reporting_limit == 0, gm_result <= 0, gm_chemical_vvl == ch) #list of measurements with unknown RL
  kRL <- domwqdata %>% dplyr::filter(gm_chemical_vvl == ch, !is.na(gm_reporting_limit), !gm_reporting_limit == 0)             #list of measurements with known RL
  if(dim(kRL)[1] == 0) {   #if chemical has no RLs listed, take implied non-detects
    kRL <- domwqdata %>% dplyr::filter(gm_chemical_vvl == ch, gm_result_modifier == "<", !gm_result == 0) %>% 
      mutate(gm_reporting_limit = gm_result) 
  }
  if(dim(uRL)[1] > 0 & dim(kRL)[1] > 0) {
    indx_b <- neardate(uRL$gm_chemical_vvl, kRL$gm_chemical_vvl, uRL$gm_samp_collection_date, kRL$gm_samp_collection_date, best = "prior") #get index of closest RL prior DATE
    indx_a <- neardate(uRL$gm_chemical_vvl, kRL$gm_chemical_vvl, uRL$gm_samp_collection_date, kRL$gm_samp_collection_date, best = "after") #get index of closest RL after DATE
    estRLs <- kRL[ifelse(is.na(indx_b), indx_a, indx_b), c("gm_reporting_limit", "gm_samp_collection_date")]              #assign newRLs and date of newRL, first looking at priors and then afters, if no priors exist
    uRL$estRL <- estRLs$gm_reporting_limit                     #add new RLs to measurement
    uRL$estRLdate <- estRLs$gm_samp_collection_date            #add date new RL came from, to identify potential mistakes/keep track of method
    print(ifelse(anyNA(uRL$estRL), nrow(uRL %>% filter(is.na(estRL))), paste0("Clear for ", ch))) #check if RLs have been standardized, if not it prints number of problem rows
    additions <- rbind(additions, uRL)                         #accumulate table of unknown RLs (with new RLs) for joining later
    finishedchems <- append(finishedchems, ch)                 #keep track of successful chemicals
    rm(indx_a, indx_b, estRLs)
  } else {
    print(paste0("Clear for ", ch))
  }
  rm(uRL, kRL)
} #assigns a RL to non-detects without one by looking at closest earlier RL of the same chemical
domwqdata <- left_join(domwqdata, additions)
rm(additions)

#compress RLs into a single column, if applicable
domwqdata <- domwqdata %>% mutate(RL = ifelse(!is.na(gm_reporting_limit) & gm_reporting_limit > 0, gm_reporting_limit, 
                                        ifelse(!is.na(estRL), estRL, gm_result)),
               RLmethod = ifelse(!is.na(gm_reporting_limit) & gm_reporting_limit > 0, "reported", 
                          ifelse(!is.na(estRL), "estimated", "result"))) %>%
  left_join(., chem_list, by = c("gm_chemical_vvl" = "chemVVL"))      #join with chemical MCLs
domwqdata$chemNAME <- NULL

#calculate MCL index and track detections/non-detections, outliers
domwqdata <- domwqdata %>% mutate(detection = ifelse(gm_result_modifier == "<", "non-detect", 
                                    ifelse(gm_result <= 0, "non-detect",        #negative radioactive measurements are non-detects
                                    ifelse(RLmethod == "reported" & gm_result <= gm_reporting_limit, "non-detect", "detect"))),
                                  RES = ifelse(detection == "non-detect", RLfun(RL, CCV), gm_result),
                                  MCLindex = RES/CCV)                                            #calculate index
md_max <- summarize(group_by(domwqdata %>% dplyr::filter(detection == "detect"), gm_chemical_vvl), avg = mean(RES), sd = sd(RES)) %>% 
  mutate(max = (avg + (10*sd))) %>% select(gm_chemical_vvl, max)               #calculates outlier maximum
domwqdata <- left_join(domwqdata, md_max, by = "gm_chemical_vvl") %>% 
  mutate(outlier_status = ifelse(is.na(max), "no outliers", 
                                 ifelse(RES > max & detection == "detect", "outlier", "standard")))

#clean table and output results
domwqdata <- domwqdata %>% select(gm_well_id, gm_samp_collection_date, gm_chemical_vvl,
                                  RES, RL, RLmethod, detection, outlier_status, MCLindex)
write.table(domwqdata, paste0("domwqdata_clean", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
