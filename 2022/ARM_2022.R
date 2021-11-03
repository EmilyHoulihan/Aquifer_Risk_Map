#Water quality data download for Needs Assessment and Aquifer Risk Map

#Methodology write-up: https://gispublic.waterboards.ca.gov/portal/home/item.html?id=70feb9f4b00f4b3384a9a0bf89f9f18a

#Updated 6/24/2021 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

#Components####
#1. Load libraries, data tables, and define functions
#2. Create depth filter to identify state small and domestic well depth water quality 
#3. Download water quality data from GAMA GIS from
  #Division of Drinking Water
  #Department of Water Resources
  #USGS NWIS
  #GAMA projects
  #ILRP
  #LocalGW
  #GeoTracker cleanup sites
#4. Standardize and clean water quality data (standardize chemicals, handle outliers and non-detects)
#5. Calculate point risk scores (wells)
# 5.1 Calculate well averages
# 5.2 Calculate count of recent samples above/close to MCL
# 5.3 Determine point risk
# 5.4 "HR2W" domestic wells
#6. Calculate square mile section risk (MTRS)
# 6.1 Calculate source section average/recent results
# 6.2 Calculate neighbor section average/recent results 
# 6.3 Determine section risk
# 6.4 Write individual chemical tables and shapefiles
#7. Summarize section data by block group

#1. Load required libraries, tables, and define special functions####
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('rgdal')) install.packages('rgdal'); library('rgdal')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('janitor')) install.packages('janitor'); library('janitor')
if (!require('httr')) install.packages('httr'); library('httr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('sf')) install.packages('sf'); library('sf')

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
chem_list <- read.csv("Datasets/NA_allchems_updated.csv", stringsAsFactors = F) %>% 
  filter(CCT == "MCL" | chemVVL %in% c("CU", "PB", "NNSM", "CR6")) %>%
  filter(!chemVVL %in% c("ASBESTOS", "COLIFORM", "FCOLIFORM", "RN-222"))
#get domestic well data information
oswcr_domdata <- read.csv("Datasets/oswcr_dom_counts2021-09-24.csv", stringsAsFactors = F) %>% distinct()
#oswcr_domdata2 <- read.csv("Datasets/CALGW_MTRS_with_Domestic_Well_Counts.csv", stringsAsFactors = F)

#2. Create depth filter ####
#import groundwater unit boundaries
gu_link <- "https://pubs.usgs.gov/ds/796/downloads/ds796_GIS.zip"
temp <- tempfile()
temp2 <- tempfile()
download.file(gu_link, temp)
unzip(zipfile = temp, exdir = temp2)
GU <- st_read(paste0(temp2, "/ds796_GIS/CA_Groundwater_Units.shp")) %>% st_transform(crs = 4326) %>% st_make_valid()
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
 MTRS <- st_read(request) %>% st_transform(crs = 4326) #%>% select(MTRS) %>% distinct()
#MTRS <- st_read("Datasets/i07_WellReportStatsBySection_1009.shp") %>% st_transform(crs = 4326) #%>% select(MTRS) %>% distinct()

#join mtrs to groundwater units, then drop geometry
MTRS_GU <- st_join(MTRS, GU) %>% st_drop_geometry()
 
#join mtrs to groundwater units and remove duplicate sections (that are duplicated across county line splits of PLSS sections)
#define columns that must be numeric
numeric_cols <- c("DomWellCount", "DomWellDepthAvg", "DomWellDepthMin", "DomWellDepthMax",
                  "PubWellCount", "PubWellDepthAvg", "PubWellDepthMin", "PubWellDepthMax")

#identify columns to be deleted
delete <- c("COUNTY_CD", "WCRFolderLink")

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
write.table(GU_depths, paste0("Tables/well_depth_gwunits", Sys.Date(), ".txt"), sep = "\t", row.names = F)

#3. Download water quality data ####
#GU_depths <- read.table("Tables/well_depth_gwunits2021-09-08.txt", stringsAsFactors = F, header = T)

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
  left_join(., select(GU_depths, GU_ID, avgmax_dom, maxdomestic, mindomestic, use_public_as_dom), by = c("GU_ID")) #%>%
  #dplyr::filter(!is.na(GU_ID))

#define numeric columns and replace 0 with NULL
numeric_cols <- c("gm_well_depth_ft", "gm_top_depth_of_screen_ft", "gm_bottom_depth_of_screen_ft")
wells_int <- wells_int %>% mutate_at(numeric_cols, as.numeric)
wells_int$gm_well_depth_ft[wells_int$gm_well_depth_ft == 0] <- NA                           #replace 0s with NA
wells_int$gm_top_depth_of_screen_ft[wells_int$gm_top_depth_of_screen_ft == 0] <- NA
wells_int$gm_bottom_depth_of_screen_ft[wells_int$gm_bottom_depth_of_screen_ft == 0] <- NA
wells_int <- mutate(wells_int, depth = ifelse(is.na(wells_int$gm_well_depth_ft),          #determine numeric depth, if possible
                                       ifelse(is.na(wells_int$gm_bottom_depth_of_screen_ft),
                                       ifelse(is.na(wells_int$gm_top_depth_of_screen_ft), NA, 
                                                     wells_int$gm_top_depth_of_screen_ft), 
                                                     wells_int$gm_bottom_depth_of_screen_ft),
                                                     wells_int$gm_well_depth_ft))
#filter on 1) if domestic, 2) if numeric, apply numeric filter OR if depth unknown apply unit use
allwells <- wells_int %>% 
  mutate(domdepth = ifelse(gm_well_category == "DOMESTIC" | gm_dataset_name == "GAMA_DOM", "yes", #keep "domestic" wells
                    ifelse(gm_dataset_name == "WB_CLEANUP",                                       #look at Geotracker wells
                           ifelse(is.na(depth), "no",                                             #remove all GT wells without numeric depth
                           ifelse(!is.na(maxdomestic) & depth <= maxdomestic & depth > mindomestic & depth > 100, "yes", "no")), #min 100 ft depth
                    ifelse(is.na(use_public_as_dom), "yes",                                       #if basin has no domestic well depth data, use all wq wells
                      ifelse(is.na(depth),                                                        #for non-GT wells without numeric depth
                        ifelse(use_public_as_dom == "yes", "yes", "no"),                           #use public/domestic comparison by basin
                        ifelse(depth <= maxdomestic & depth > mindomestic, "yes",                 #for non-GT wells with numeric depth, use numeric depth filter
                    "no"))))))
rm(wells_int)
lat_long <- loc %>% select(gm_well_id, gm_longitude, gm_latitude)
allwells <- left_join(allwells, lat_long)
write.table(allwells, paste0("Tables/allwells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

#allwells <- read.table("Tables/allwells2021-09-08.txt", stringsAsFactors = F, header = T)
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
  mytable <- mytable %>% filter(gm_samp_collection_date >= ymd("2001-01-01"), gm_samp_collection_date <= ymd("2021-01-01"),
                                gm_well_id %in% dom_wells$gm_well_id,
                                gm_chemical_vvl %in% mychems$chemVVL)
  return(mytable)
}

#download various datasets
ddw <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_ddw_statewide_v2.zip',
                  'gama_ddw_statewide_v2.txt', col_sty, domwells, chem_list)
ddw <- ddw %>% filter(!str_detect(gm_altwell_id1, "SPRING") | !str_detect(gm_altwell_id2, "SPRING")) #remove ~240 spring sites

dwr <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_dwr_statewide_v2.zip',
                  'gama_dwr_statewide_v2.txt', col_sty, domwells, chem_list)

gamadom <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_gama_dom_statewide_v2.zip',
                      'gama_gama_dom_statewide_v2.txt', col_sty, domwells, chem_list)

gamausgs <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_gama_usgs_statewide_v2.zip',
                       'gama_gama_usgs_statewide_v2.txt', col_sty, domwells, chem_list)

localgw <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_localgw_statewide_v2.zip',
                      'gama_localgw_statewide_v2.txt', col_sty, domwells, chem_list)

usgsnwis <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_usgs_nwis_statewide_v2.zip',
                       'gama_usgs_nwis_statewide_v2.txt', col_sty, domwells, chem_list)

wbilrp <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_wb_ilrp_statewide_v2.zip',
                     'gama_wb_ilrp_statewide_v2.txt', col_sty, domwells, chem_list)

#compile datasets
domwqdata <- rbind(ddw, dwr, gamadom, gamausgs, localgw, usgsnwis, wbilrp)
rm(ddw, dwr, gamadom, gamausgs, localgw, usgsnwis, wbilrp)

#remove duplicated data (duplicated between GAMA USGS and USGS NWIS datasets)
duplicated_wells <- read.csv("Datasets/NWIS_dupes2.csv", stringsAsFactors = F)
domwqdata <- domwqdata %>% filter(!gm_well_id == "Vanadium Study",
                                  !gm_well_id %in% duplicated_wells$Duplicate_NWIS,
                                  !gm_well_id == "USGS-415953121292901")
domwells_wq <- domwqdata %>% select(gm_well_id) %>% distinct() %>% left_join(., domwells)

#save water quality data and well information
write.table(domwqdata, paste0("Tables/wqdatanoGT_", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
write.table(domwells_wq, paste0("Tables/wqdatawellsnoGT", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

#Geotracker data (too large to download at once, must get by county)
cl <- loc %>% select(gm_gis_county) %>% distinct()
cl$gm_gis_county <- tolower(cl$gm_gis_county)                       #make lower case and remove NAs
cl <- cl %>% filter(!is.na(gm_gis_county), !gm_gis_county == "no county found", !gm_gis_county == "orange county")
cl <- cl[order(cl$gm_gis_county),]
edf <- data.frame()
for (x in 1:nrow(cl)) {
  baselink <- "https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_wb_cleanup_"
  county <- gsub(" ", "", cl[x,1], fixed = TRUE)
  countyedflink <- paste0(baselink, county, "_v2.zip")
  newdat <- dwd_wq_tbl(countyedflink, paste0("gama_wb_cleanup_", county, "_v2.txt"), col_sty, domwells, chem_list)
  edf <- rbind(edf, newdat)
  print(cl[x,1])     #to keep track of progress!
  rm(newdat)
}
rm(cl, col_sty)
edfwqwells <- edf %>% select(gm_well_id) %>% distinct() %>% left_join(., domwells)
#save water quality data and well information
write.table(edf, paste0("Tables/wqdataGT_", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
write.table(edfwqwells, paste0("Tables/wqdatawellsGT", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

#4. Standardize and clean data ####
domwqdata <- read.table("Tables/wqdatanoGT_2021-09-08.txt", stringsAsFactors = F, header = T) %>% 
  select(gm_well_id, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit)
edf <- read.table("Tables/wqdataGT_2021-09-08.txt", stringsAsFactors = F, header = T) %>% 
  select(gm_well_id, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit)
domwqdata <- rbind(domwqdata, edf)

edfwqwells <- read.table("Tables/wqdatawellsGT2021-09-08.txt", stringsAsFactors = F, header = T)
domwells_wq <- read.table("Tables/wqdatawellsnoGT2021-09-08.txt", stringsAsFactors = F, header = T)
domwells_wq <- rbind(edfwqwells, domwells_wq)

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
domwqdata <- domwqdata %>% dplyr::filter(!gm_chemical_vvl == "NO3NO2N") #delete duplicated rows
rm(nitrate, nitrogen, nitrite, n_add, nonitrate, sub_data, combined, con_data)

#standardize Qualifiers to either "<" or "="
domwqdata$gm_result_modifier[is.na(domwqdata$gm_result_modifier)] <- "="                        
domwqdata <- domwqdata %>% dplyr::filter(!gm_result_modifier %in% c("S", "V", "M")) %>%
  dplyr::mutate(gm_result_modifier = ifelse(is.na(gm_result_modifier), "=",
                                     ifelse(gm_result_modifier == "<" | gm_result_modifier == "ND", "<", "=")))

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
domwqdata <- domwqdata %>% 
  mutate(detection = ifelse(gm_result_modifier == "<", "non-detect",
                     ifelse(is.na(gm_result), "non-detect", 
                     ifelse(gm_result <= 0, "non-detect",        #negative radioactive measurements are non-detects
                     ifelse(RLmethod == "reported" & gm_result <= gm_reporting_limit, "non-detect", "detect")))),
                                  RES = ifelse(detection == "non-detect", RLfun(RL, CCV), gm_result),
                                  MCLindex = RES/CCV)                                            #calculate index
md_max <- summarize(group_by(domwqdata %>% dplyr::filter(detection == "detect"), gm_chemical_vvl), avg = mean(RES), sd = sd(RES)) %>% 
  mutate(max = (avg + (10*sd))) %>% select(gm_chemical_vvl, max)               #calculates outlier maximum
domwqdata <- left_join(domwqdata, md_max, by = "gm_chemical_vvl") %>% 
  mutate(outlier_status = ifelse(is.na(max), "no outliers", 
                                 ifelse(RES > max & detection == "detect", "outlier", "standard")))

#clean table and output results
write.table(domwqdata, paste0("Tables/wqdata_clean", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
domwqdata_sm <- domwqdata %>% select(gm_well_id, gm_samp_collection_date, gm_chemical_vvl,
                                  RES, RL, RLmethod, detection, outlier_status, MCLindex) %>%
  left_join(., select(domwells_wq, gm_well_id, gm_dataset_name, MTRS)) #attach MTRS for easy filtering later
write.table(domwqdata_sm, paste0("Tables/wqdata_small_", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)
write.table(domwells_wq, paste0("Tables/wqdatawells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

#5. Calculate point risk scores (wells)####
domwqdata_sm <- read.table("Tables/wqdata_small_2021-09-08.txt", stringsAsFactors = F, header = T)
wellinfo <- read.table("Tables/wqdatawells2021-09-08.txt", stringsAsFactors = F, header = T) %>% select(gm_well_id, gm_latitude, gm_longitude, gm_well_category, gm_dataset_name, MTRS)

domwqdata_sm$gm_samp_collection_date <- lubridate::ymd(domwqdata_sm$gm_samp_collection_date) #read dates as dates
domwqdata_sm <- domwqdata_sm %>% 
  filter(!outlier_status == "outlier", !is.na(MCLindex)) %>% #remove outliers and samples without values (?)
  select(gm_well_id, gm_samp_collection_date, gm_chemical_vvl, MCLindex) %>% #select only necessary columns to calculate averages
  mutate(year = year(gm_samp_collection_date), #prepare for grouping samples by year
         MCLindex = ifelse(gm_chemical_vvl == "CR6", MCLindex*2, MCLindex)) #CR6 adjustment - change to 10 ug/L (original has it at 20 ccv)
#5.1 calculate well averages ####
domwellavg <- domwqdata_sm %>% 
  group_by(gm_chemical_vvl, gm_well_id, year) %>% #group samples by chemical, well, and year
  summarize(yr_avg = mean(MCLindex))%>%
  group_by(gm_chemical_vvl, gm_well_id) %>% #group samples by chemical, well
  summarize(well_avg = mean(yr_avg))

#5.2 calculate recent results####
domwellre_avg <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2016-01-01"), MCLindex > 1) %>% 
  group_by(gm_well_id, gm_chemical_vvl) %>%
  summarize(re_avg_over = mean(MCLindex))
domwellrecent <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2016-01-01")) %>%
  mutate(re_status = ifelse(MCLindex > 1, 'over',
                            ifelse(MCLindex > 0.8, 'close', 'under'))) %>% 
  group_by(gm_chemical_vvl, gm_well_id, re_status) %>% count() %>%
  pivot_wider(names_from = re_status, names_prefix = "re_", values_from = n) %>%
  left_join(., domwellre_avg)
domwellrecent[is.na(domwellrecent)] <- 0

#5.3 Calculate point risk####
pointdata <- full_join(domwellavg, domwellrecent)
pointdata[is.na(pointdata)] <- 0
pointdata <- pointdata %>% mutate(avg_overMCL = ifelse(well_avg > 1 & re_avg_over > 0, (well_avg + re_avg_over)/2,
                                                       ifelse(well_avg > 1 & re_avg_over == 0, well_avg,
                                                              ifelse(!well_avg > 1 & re_avg_over > 0, re_avg_over, 0)))) %>%
  left_join(., wellinfo)
rm(domwqdata, domwqdata_sm, domwellavg, domwellrecent)

pdrisk <- pointdata %>% group_by(gm_well_id) %>% arrange(-avg_overMCL) %>%
  summarize(PRF1 = n_distinct(gm_chemical_vvl[well_avg > 1 | re_over > 0]),
            PRF2 = n_distinct(gm_chemical_vvl[well_avg <= 1 & well_avg > 0.8 | re_close > 0 & re_over == 0]),
            PRF3 = mean(avg_overMCL[avg_overMCL > 0], na.rm = T),
            PL1 = paste(gm_chemical_vvl[well_avg > 1 | re_over > 0], collapse = "; "),
            PL2 = paste(gm_chemical_vvl[well_avg <= 1 & well_avg > 0.8 | re_close > 0 & re_over == 0], collapse = "; ")) %>% 
  mutate(risk = ifelse(PRF1 > 0, "high", ifelse(PRF1 == 0 & PRF2 > 0, "med", "low"))) %>%
  left_join(., wellinfo)

write.csv(pointdata, paste0("Tables/pdatadetailed", Sys.Date(), ".csv"), row.names = F)
write.csv(pdrisk, paste0("Tables/pdata_risk", Sys.Date(), ".csv"), row.names = F)

#5.4 "HR2W" domestic wells####
pdrisk <- read.csv("Tables/pdata_risk2021-09-08.csv", stringsAsFactors = F)
pdrisk_dom <- pdrisk %>% filter(gm_well_category == "DOMESTIC")
write.csv(pdrisk_dom, "Shapefiles/hr2w_arm_domestic.csv", row.names = F)

#6. Calculates square mile section risk####
#load point data tables
pointdata <- read.csv("Tables/pdatadetailed2021-09-08.csv", stringsAsFactors = F)

#6.1 Calculate source section average/recent results####
#load file that contains neighbor reference data
temp <- tempfile()
download.file("https://gispublic.waterboards.ca.gov/portal/sharing/rest/content/items/05f49f6f8bd24ac38143bdda402ba098/data", temp)
neighborsec1 <- read.csv(temp)
unlink(temp)
temp <- tempfile()
download.file("https://gispublic.waterboards.ca.gov/portal/sharing/rest/content/items/a72861be691643359e62e78d2c1a8347/data", temp)
neighborsec2 <- read.csv(temp)
unlink(temp)

neighborsec <- rbind(neighborsec1, neighborsec2) %>% distinct() %>% filter(!MTRS_1 == MTRS_2)
colnames(neighborsec) <- c("MTRS", "MTRS_2")

#group point data by MTRS (source data first)
source_sections <- pointdata %>% left_join(., select(wellinfo, gm_well_id, MTRS)) %>%
  group_by(gm_chemical_vvl, MTRS) %>%
  summarize(SecDet = mean(well_avg),
            re_under = sum(re_under),
            re_close = sum(re_close),
            re_over = sum(re_over),
            re_avg_over = mean(re_avg_over[re_avg_over > 0])) %>%
  mutate(method = "source")

#6.2 Calculate neighbor section average/recent results####
neighbor_sections <- left_join(source_sections, neighborsec) %>%
  group_by(gm_chemical_vvl, MTRS_2) %>%
  summarize(SecDet = mean(SecDet),
            re_under = mean(re_under),
            re_close = mean(re_close),
            re_over = mean(re_over),
            re_avg_over = mean(re_avg_over, na.rm = T)) %>%
  mutate(method = "neighbor")
colnames(neighbor_sections)[colnames(neighbor_sections)== "MTRS_2"] <- "MTRS"
#for section/chemical pairs that do not have source data, use neighbor data
neighbor_sections_add <- anti_join(neighbor_sections, source_sections, by = c("gm_chemical_vvl", "MTRS"))
section_data <- rbind(source_sections, neighbor_sections_add) %>% filter(!is.na(MTRS))
section_data[is.na(section_data)] <- 0
rm(neighbor_sections, neighbor_sections_add, neighborsec, neighborsec1, neighborsec2, source_sections,
   pdrisk, pdrisk_dom, pointdata)

#calculate average magnitude over MCL
section_data <- section_data %>% mutate(avg_overMCL = ifelse(SecDet > 1 & re_avg_over > 0, (SecDet + re_avg_over)/2,
                                                             ifelse(SecDet > 1 & re_avg_over == 0, SecDet,
                                                                    ifelse(!SecDet > 1 & re_avg_over > 0, re_avg_over, 0))),
                                        WQriskbin = ifelse(SecDet > 1 | re_over > 0, "high",
                                                           ifelse(SecDet > 0.8 | re_close > 0, "med", "low")))
section_data[is.na(section_data)] <- 0

#6.3 Determine section risk####
sdrisk <- section_data %>% 
  group_by(MTRS) %>% arrange(-SecDet) %>%
  summarize(SRF1 = n_distinct(gm_chemical_vvl[WQriskbin == "high"]),
            SRF2 = n_distinct(gm_chemical_vvl[WQriskbin == "med"]),
            SRF3 = mean(avg_overMCL[avg_overMCL > 0]),
            SL1 = paste(gm_chemical_vvl[WQriskbin == "high"], collapse = "; "),
            SL2 = paste(gm_chemical_vvl[WQriskbin == "med"], collapse = "; "))
sdrisk[is.na(sdrisk)] <- 0
sdrisk <- sdrisk %>% mutate(WQriskbins = ifelse(SRF1 > 0, "high",
                                                ifelse(SRF2 > 0 & SRF1 == 0, "med", "low")))

#load state small water system count per section and domestic well count per section - join with wq risk
ss_mtrs <- read.csv("Datasets/SWS_rcac_MTRS_table.csv", stringsAsFactors = F) %>% 
  select(sysname, mtrs) %>%  distinct() %>% group_by(mtrs) %>%
  summarize(ssws_count = n())
sdrisk <- sdrisk %>% 
  full_join(., oswcr_domdata) %>%
  full_join(., select(oswcr_domdata2, MTRS, Count.of.Domestic.Wells)) %>%
  full_join(ss_mtrs, by = c("MTRS" = "mtrs"))
sdrisk$ssws_count[is.na(sdrisk$ssws_count)] <- 0
sdrisk$dom_wells1970na[is.na(sdrisk$dom_wells1970na)] <- 0

write.csv(section_data, paste0("Tables/sectiondata_detailed", Sys.Date(), ".csv"), row.names = F)
write.csv(sdrisk, paste0("Tables/sectiondata_risk", Sys.Date(), ".csv"), row.names = F)

#MTRS <- st_read("Datasets/i07_WellReportStatsBySection_1009.shp") %>% st_transform(crs = 4326) %>% select(MTRS) %>% distinct()
sdrisk_mtrs <- left_join(MTRS, sdrisk)
sdrisk_mtrs[is.na(sdrisk_mtrs)] <- 0
st_write(sdrisk_mtrs, paste0("Shapefiles/arm_sectionrisk", Sys.Date(), ".shp"))

#6.4 write individual chemical shapefiles####
section_data_chems <- section_data
section_data_chems[is.na(section_data_chems)] <- 0
indiv_chems <- c("NO3N", "AS", "CR6", "TCPR123", "U")
for (x in 1:length(indiv_chems)) {
  indiv_table <- section_data_chems %>% filter(gm_chemical_vvl == indiv_chems[x])
  write.csv(indiv_table, paste0("Tables/section_risk_", as.character(indiv_chems[x]), ".csv"), row.names = F, na = "")
  indiv_shp <- left_join(MTRS, indiv_table)
  st_write(indiv_shp, paste0("Shapefiles/", as.character(indiv_chems[x]), Sys.Date(), ".shp"))
}
rm(indiv_shp, indiv_table)

#7. summarize section data by block group####
sdrisk <- read.csv("Tables/sectiondata_risk2021-09-08.csv", stringsAsFactors = F) %>% select(MTRS, WQriskbins)

mtrs_cbg <- st_read("Datasets/MTRS_2018bg_intersect.shp") %>% st_drop_geometry() #%>% 
  filter(!MTRS %in% c("SALTONSEA", "LAKETAHOE", "BAY/DELTA")) %>%
  group_by(MTRS) %>%
  mutate(tot_area = sum(ShapeArea),
         per_area = ShapeArea/tot_area) %>%
  filter(!TRACTCE %in% c(990000, 990100, 990200, 990300)) #remove ocean census blocks
cbg_wq <- left_join(mtrs_cbg, sdrisk) %>% 
  left_join(., oswcr_domdata) %>%
  mutate(per_domwell = dom_wells1970na*per_area,
         #per_domwell = Count.of.Domestic.Wells*per_area,
         WQriskbins = ifelse(is.na(WQriskbins), "nodata", WQriskbins)) %>%
  group_by(GEOID) %>%
  summarize(pa_hr = sum(ShapeArea[WQriskbins == "high"])/sum(ShapeArea),
            pa_mr = sum(ShapeArea[WQriskbins == "med"])/sum(ShapeArea),
            pa_lr = sum(ShapeArea[WQriskbins == "low"])/sum(ShapeArea),
            pa_nr = sum(ShapeArea[WQriskbins == "nodata"])/sum(ShapeArea),
            dw_hr = sum(per_domwell[WQriskbins == "high"]),
            dw_mr = sum(per_domwell[WQriskbins == "med"]),
            dw_lr = sum(per_domwell[WQriskbins == "low"]),
            dw_nr = sum(per_domwell[WQriskbins == "nodata"]),
            dw_tot = sum(per_domwell))

#small detour for state small water systems
ssws <- read.csv("Datasets/ssws_cbg.csv", stringsAsFactors = F) %>%
  select(primary_id, sysname, city, COUNTYFP, GEOID)
ss_mtrs <- read.csv("Datasets/SWS_rcac_MTRS_table.csv", stringsAsFactors = F) %>%
  select(primary_id, sysname, city, mtrs)
ssws <- full_join(ssws, ss_mtrs) %>%
  left_join(., select(sdrisk, MTRS, WQriskbins), by = c("mtrs" = "MTRS")) %>%
  mutate(WQriskbins = ifelse(is.na(WQriskbins), "nr", 
                             ifelse(WQriskbins == "high", "hr",
                                    ifelse(WQriskbins == "med", "mr", "lr"))))

cbg_wq_ssws <- ssws %>% group_by(GEOID, WQriskbins) %>%
  count() %>%
  pivot_wider(names_from = WQriskbins, values_from = n, values_fill = 0, names_prefix = "ssws_") %>%
  mutate(ssws_tot = ssws_lr + ssws_nr + ssws_hr + ssws_mr,
         GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                        paste0("0", GEOID)))
#join ssws with domestic well counts by block groups
cbg_wq <- left_join(cbg_wq, cbg_wq_ssws)
cbg_wq[c("ssws_lr", "ssws_nr", "ssws_hr", "ssws_mr", "ssws_tot")][is.na(cbg_wq[c("ssws_lr", "ssws_nr", "ssws_hr", "ssws_mr", "ssws_tot")])] <- 0

#join with mhi data
mhi <- read.dbf("Datasets/2018bgmhi_dwr.dbf", as.is = T) %>% select(GEOID10, MHI18, COUNTYFP10) %>%
  mutate(DAC_status = ifelse(MHI18 == 0, "no data",
                             ifelse(MHI18 < 42737, "SDAC",
                                    ifelse(MHI18 < 56982, "DAC", "none"))))
cbg_wq_mhi <- left_join(cbg_wq, mhi, by = c("GEOID" = "GEOID10"))

#join with cdag data
cdag_cbg_scores <- read.csv("WaterShortageRisk/cdag_cbg_adj_09022021.csv", stringsAsFactors = F) %>% 
  select(GEOID, RCsum_rescale_NEW) %>%
  mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                        paste0("0", GEOID)))
cbg_wq_mhi_cdag <- left_join(cbg_wq_mhi, cdag_cbg_scores) %>%
  mutate(CDAG_bin = ifelse(is.na(RCsum_rescale_NEW), "no_data",
                           ifelse(RCsum_rescale_NEW > 55, "upper",
                                  ifelse(RCsum_rescale_NEW > 28, "middle", "lower"))))

ssws <- ssws %>% left_join(., cdag_cbg_scores)
write.csv(ssws, "Tables/ssws_wqrisk.csv", row.names = F)

#county names
fips_link <- 'https://www2.census.gov/programs-surveys/popest/geographies/2019/all-geocodes-v2019.xlsx'
temp <- tempfile()
download.file(url = fips_link, destfile = temp, mode = 'wb')
fips_codes <- read_xlsx(temp, skip = 4)
fips_codes <- clean_names(fips_codes) 
fips_codes <- fips_codes %>% dplyr::filter(state_code_fips == "06", !county_code_fips == "000") %>% 
  dplyr::select(county_code_fips, area_name_including_legal_statistical_area_description)
colnames(fips_codes)[colnames(fips_codes) == "area_name_including_legal_statistical_area_description"] <- "county"
cbg_wq_mhi_cdag <- left_join(cbg_wq_mhi_cdag, fips_codes, by = c("COUNTYFP10" = "county_code_fips")) %>%
  filter(!is.na(county))

#prep for ArcGIS
round_cols <- c("dw_hr", "dw_mr", "dw_lr", "dw_nr", "dw_tot",
                "RCsum_rescale_NEW")
cbg_wq_mhi_cdag <- cbg_wq_mhi_cdag %>% 
  mutate_at(round_cols, round, 2)
cbg_wq_mhi_cdag$RCsum_rescale_NEW[is.na(cbg_wq_mhi_cdag$RCsum_rescale_NEW)] <- -99
write.csv(cbg_wq_mhi_cdag, paste0("Tables/cbg_wq_mhi_cdag", Sys.Date(), ".csv"), row.names = F)
