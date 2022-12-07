#Water quality data download for 2023 Aquifer Risk Map

#Updated 12/72022 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

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
if (!require('readxl')) install.packages('readxl'); library('readxl')

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

#2. Create well depth filter ####
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
write.table(GU_depths, paste0("Datasets/well_depth_gwunits", Sys.Date(), ".txt"), sep = "\t", row.names = F)

#3. Download water quality data ####
GU_depths <- read.table("Datasets/well_depth_gwunits2022-08-19.txt", stringsAsFactors = F, header = T)

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
wells_int[numeric_cols][wells_int[numeric_cols] == 0] <- NA
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
                                         ifelse(!is.na(maxdomestic) & depth <= maxdomestic & depth > mindomestic & depth > 100, "yes", "no")), #min 100 ft depth for WB cleanup wells
                                  ifelse(is.na(use_public_as_dom), "yes",                                       #if basin has no domestic well depth data, use all wq wells
                                         ifelse(is.na(depth),                                                        #for non-GT wells without numeric depth
                                                ifelse(use_public_as_dom == "yes", "yes", "no"),                           #use public/domestic comparison by basin
                                                ifelse(depth <= maxdomestic & depth > mindomestic, "yes",                 #for non-GT wells with numeric depth, use numeric depth filter
                                                       "no"))))))
rm(wells_int)
lat_long <- loc %>% select(gm_well_id, gm_longitude, gm_latitude)
allwells <- left_join(allwells, lat_long)
write.table(allwells, paste0("Datasets/allwells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

allwells <- read.table("Datasets/allwells2022-08-23.txt", stringsAsFactors = F, header = T)
domwells <- allwells %>% filter(domdepth == "yes")
wells_unk_area <- allwells %>% filter(domdepth == "no", !gm_dataset_name == "WB_CLEANUP")
possiblewells <- rbind(domwells, wells_unk_area)

#retrieve water quality data
#define column types
col_sty <- cols(
  GM_DATASET_NAME = col_character(),GM_WELL_CATEGORY = col_character(),GM_DATA_SOURCE = col_character(),
  GM_WELL_ID = col_character(),GM_CHEMICAL_VVL = col_character(),GM_CHEMICAL_NAME = col_character(),
  GM_RESULT_MODIFIER = col_character(),GM_RESULT = col_double(),GM_RESULT_UNITS = col_character(),
  GM_SAMP_COLLECTION_DATE = col_character(),GM_REPORTING_LIMIT = col_double(),GM_LATITUDE = col_double(),
  GM_LONGITUDE = col_double(),GM_WELL_DEPTH_FT = col_double(),GM_TOP_DEPTH_OF_SCREEN_FT = col_double(),
  GM_BOTTOM_DEPTH_OF_SCREEN_FT = col_double(),GM_CAS_NUMBER = col_character(),GM_ALTWELL_ID1 = col_character(),
  GM_ALTWELL_ID2 = col_character(),GM_ALTWELL_ID3 = col_character(),SRC_CHEMICAL = col_character(),
  SRC_RESULT_MODIFIER = col_character(),SRC_RESULT = col_double(),SRC_RESULT_UNITS = col_character(),
  SRC_SAMP_COLLECTION_DATE = col_character(),SRC_SAMP_COLLECTION_TIME = col_character(),SRC_REPORTING_LIMIT = col_double(),
  SRC_ANALYTICAL_METHOD = col_character(),SRC_LAB_NOTE = col_character(),SRC_LATITUDE = col_double(),
  SRC_LONGITUDE = col_double(),SRC_DATUM = col_character(),SRC_WELL_DEPTH_FT = col_double(),
  SRC_TOP_DEPTH_OF_SCREEN_FT = col_double(),SRC_BOTTOM_DEPTH_OF_SCREEN_FT = col_double()
)
#define function to download data
dwd_wq_tbl <- function(link, filename, col_str, my_wells, mychems) {
  temp <- tempfile()
  download.file(link, temp)
  mytable <- readr::read_tsv(unz(temp, filename), quote = "", col_types = col_str) %>% clean_names()
  unlink(temp)
  mytable$gm_samp_collection_date <- lubridate::mdy(mytable$gm_samp_collection_date)
  mytable <- mytable %>% filter(gm_samp_collection_date >= ymd("2002-01-01"),
                                gm_well_id %in% my_wells$gm_well_id,
                                gm_chemical_vvl %in% mychems$chemVVL)
  return(mytable)
}

#download various datasets for domestic depth wells
ddw <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_ddw_statewide_v2.zip',
                  'gama_ddw_statewide_v2.txt', col_sty, possiblewells, chem_list)
ddw <- ddw %>% filter(!str_detect(gm_altwell_id1, "SPRING") | !str_detect(gm_altwell_id2, "SPRING")) #remove ~240 spring sites
ddw_wells <- ddw %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
fwrite(ddw, "Datasets/rawgw/ddw2022-09-29.csv")
fwrite(ddw_wells, "Datasets/rawgw/ddw_wells2022-09-29.csv")

dwr <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_dwr_statewide_v2.zip',
                  'gama_dwr_statewide_v2.txt', col_sty, possiblewells, chem_list)
dwr_wells <- dwr %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
fwrite(dwr, "Datasets/rawgw/dwr2022-09-29.csv")
fwrite(dwr_wells, "Datasets/rawgw/dwr_wells2022-09-29.csv")

gamadom <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_gama_dom_statewide_v2.zip',
                      'gama_gama_dom_statewide_v2.txt', col_sty, possiblewells, chem_list)
gamadom_wells <- gamadom %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
fwrite(gamadom, "Datasets/rawgw/gamadom2022-09-29.csv")
fwrite(gamadom_wells, "Datasets/rawgw/gamadom_wells2022-09-29.csv")

#gamausgs <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_gama_usgs_statewide_v2.zip',
#                       'gama_gama_usgs_statewide_v2.txt', col_sty, domwells, chem_list)

localgw <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_localgw_statewide_v2.zip',
                      'gama_localgw_statewide_v2.txt', col_sty, possiblewells, chem_list)
localgw_wells <- localgw %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
fwrite(localgw, "Datasets/rawgw/localgw2022-09-29.csv")
fwrite(localgw_wells, "Datasets/rawgw/localgw_wells2022-09-29.csv")

usgsnwis <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_usgs_nwis_statewide_v2.zip',
                       'gama_usgs_nwis_statewide_v2.txt', col_sty, possiblewells, chem_list)
usgsnwis_wells <- usgsnwis %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
fwrite(usgsnwis, "Datasets/rawgw/usgsnwis2022-09-29.csv")
fwrite(usgsnwis_wells, "Datasets/rawgw/usgsnwis_wells2022-09-29.csv")

wbilrp <- dwd_wq_tbl('https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_wb_ilrp_statewide_v2.zip',
                     'gama_wb_ilrp_statewide_v2.txt', col_sty, possiblewells, chem_list)
wbilrp_wells <- wbilrp %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
fwrite(wbilrp, "Datasets/rawgw/wbilrp2022-09-29.csv")
fwrite(wbilrp_wells, "Datasets/rawgw/wbilrp_wells2022-09-29.csv")

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
  newdat <- dwd_wq_tbl(countyedflink, paste0("gama_wb_cleanup_", county, "_v2.txt"), col_sty, possiblewells, chem_list)
  edf <- rbind(edf, newdat)
  print(cl[x,1])     #to keep track of progress!
  rm(newdat)
}
rm(cl, col_sty)
edfwqwells <- edf %>% select(gm_well_id) %>% distinct() %>% left_join(., possiblewells)
#save water quality data and well information
fwrite(edfwqwells, "Datasets/rawgw/edf_wells2022-09-29.csv")
fwrite(edf, "Datasets/rawgw/edf2022-09-29.csv")

#compile datasets
ddw <- fread("Datasets/rawgw/ddw2022-09-29.csv")
ddw_wells <- fread("Datasets/rawgw/ddw_wells2022-09-29.csv")
dwr <- fread("Datasets/rawgw/dwr2022-09-29.csv")
dwr_wells <- fread("Datasets/rawgw/dwr_wells2022-09-29.csv")
gamadom <- fread("Datasets/rawgw/gamadom2022-09-29.csv")
gamadom_wells <- fread("Datasets/rawgw/gamadom_wells2022-09-29.csv")
localgw <- fread("Datasets/rawgw/localgw2022-09-29.csv")
localgw_wells <- fread("Datasets/rawgw/localgw_wells2022-09-29.csv")
usgsnwis <- fread("Datasets/rawgw/usgsnwis2022-09-29.csv")
usgsnwis$gm_well_id <- as.character(usgsnwis$gm_well_id)
usgsnwis_wells <- fread("Datasets/rawgw/usgsnwis_wells2022-09-29.csv")
usgsnwis_wells$gm_well_id <- as.character(usgsnwis_wells$gm_well_id)
wbilrp <- fread("Datasets/rawgw/wbilrp2022-09-29.csv")
wbilrp_wells <- fread("Datasets/rawgw/wbilrp_wells2022-09-29.csv")
edf <- fread("Datasets/rawgw/edf2022-09-29.csv")
edf_wells <- fread("Datasets/rawgw/edf_wells2022-09-29.csv")

domwqdata <- rbind(ddw, dwr, gamadom, localgw, usgsnwis, wbilrp, edf)
domwqwells <- rbind(ddw_wells, dwr_wells, gamadom_wells, localgw_wells, usgsnwis_wells, wbilrp_wells, edf_wells)

#remove duplicated data (not needed anymore because not downloading GAMA_USGS
domwqdata <- domwqdata %>% filter(!gm_well_id == "415953121292901") #not in california
domwqwells <- domwqwells %>% filter(!gm_well_id == "415953121292901")

#save water quality data and well information
fwrite(domwqdata, paste0("Datasets/wqdata_", Sys.Date(), ".csv"))
fwrite(domwqwells, paste0("Datasets/wqdata_wells_", Sys.Date(), ".csv"))

#4. Standardize and clean data ####
domwqdata <- fread("Datasets/wqdata_2022-09-29.csv") %>% 
  select(gm_well_id, gm_dataset_name, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit)
domwells_wq <- fread("Datasets/wqdata_wells_2022-09-29.csv")

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
#remove duplicate entries in GAMA (multiple locations for one well ID)
duplicates <- c("CA5800868_001_001", "AGL020003332-MAIN WELL", "AGW080010455-WELL_#24")
domwells_wq_clean <- domwells_wq[!duplicated(domwells_wq$gm_well_id),]
domwqdata_clean <- domwqdata %>% distinct()

#clean table and output results
fwrite(domwqdata_clean, paste0("Datasets/wqdata_clean", Sys.Date(), ".csv"))
domwqdata_sm <- domwqdata_clean %>% select(gm_well_id, gm_dataset_name, gm_samp_collection_date, gm_chemical_vvl,
                                     RES, RL, RLmethod, detection, outlier_status, MCLindex) %>%
  left_join(., select(domwells_wq_clean, gm_well_id, domdepth, MTRS)) #attach MTRS for easy filtering later
fwrite(domwqdata_sm, paste0("Datasets/wqdata_small_", Sys.Date(), ".csv"))
fwrite(domwells_wq_clean, paste0("Datasets/wqdata_wellsclean_", Sys.Date(), ".csv"))

#5. Load data tables####
domwqdata_sm <- fread("Datasets/wqdata_small_2022-09-29.csv")
wellinfo <- fread("Datasets/wqdata_wellsclean_2022-09-29.csv") %>% 
  select(gm_well_id, domdepth, gm_latitude, gm_longitude, gm_well_category, gm_dataset_name, MTRS)

#6. Calculate point risk scores (wells)####
domwqdata_sm$gm_samp_collection_date <- lubridate::ymd(domwqdata_sm$gm_samp_collection_date) #read dates as dates
domwqdata_sm <- domwqdata_sm %>% 
  filter(!outlier_status == "outlier", !is.na(MCLindex)) %>% #remove outliers and samples without values (?)
  select(gm_well_id, gm_samp_collection_date, gm_chemical_vvl, MCLindex) %>% #select only necessary columns to calculate averages
  mutate(year = year(gm_samp_collection_date), #prepare for grouping samples by year
         MCLindex = ifelse(gm_chemical_vvl == "CR6", MCLindex*2, MCLindex)) #CR6 adjustment - change to 10 ug/L (original has it at 20 ccv)
#6.1 calculate well averages ####
domwellavg <- domwqdata_sm %>% 
  group_by(gm_chemical_vvl, gm_well_id, year) %>% #group samples by chemical, well, and year
  summarize(yr_avg = mean(MCLindex))%>%
  group_by(gm_chemical_vvl, gm_well_id) %>% #group samples by chemical, well
  summarize(well_avg = mean(yr_avg))

#6.2 calculate recent results####
domwellre_avg <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2017-01-01"), MCLindex > 1) %>% 
  group_by(gm_well_id, gm_chemical_vvl) %>%
  summarize(re_avg_over = mean(MCLindex))
domwellrecent <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2017-01-01")) %>%
  mutate(re_status = ifelse(MCLindex > 1, 'over',
                            ifelse(MCLindex > 0.8, 'close', 'under'))) %>% 
  group_by(gm_chemical_vvl, gm_well_id, re_status) %>% count() %>%
  pivot_wider(names_from = re_status, names_prefix = "re_", values_from = n) %>%
  left_join(., domwellre_avg)
domwellrecent[is.na(domwellrecent)] <- 0

#6.3 Calculate point risk####
pointdata <- full_join(domwellavg, domwellrecent)
pointdata[is.na(pointdata)] <- 0
pointdata <- pointdata %>% mutate(avg_overMCL = ifelse(well_avg > 1 & re_avg_over > 0, (well_avg + re_avg_over)/2,
                                                       ifelse(well_avg > 1 & re_avg_over == 0, well_avg,
                                                              ifelse(!well_avg > 1 & re_avg_over > 0, re_avg_over, 0)))) %>%
  left_join(., wellinfo)
rm(domwqdata_sm, domwellavg, domwellrecent)

pdrisk <- pointdata %>% group_by(gm_well_id) %>% arrange(-avg_overMCL) %>%
  summarize(PRF1 = n_distinct(gm_chemical_vvl[well_avg > 1 | re_over > 0]),
            PRF2 = n_distinct(gm_chemical_vvl[well_avg <= 1 & well_avg > 0.8 | re_close > 0 & re_over == 0]),
            PRF3 = mean(avg_overMCL[avg_overMCL > 0], na.rm = T),
            PL1 = paste(gm_chemical_vvl[well_avg > 1 | re_over > 0], collapse = "; "),
            PL2 = paste(gm_chemical_vvl[well_avg <= 1 & well_avg > 0.8 | re_close > 0 & re_over == 0], collapse = "; ")) %>% 
  mutate(risk = ifelse(PRF1 > 0, "high", ifelse(PRF1 == 0 & PRF2 > 0, "med", "low"))) %>%
  left_join(., wellinfo)

fwrite(pointdata, paste0("Tables/Point Data (wells)/pdatadetailed", Sys.Date(), ".csv"))
fwrite(pdrisk, paste0("Tables/Point Data (wells)/pdata_risk", Sys.Date(), ".csv"))

#7. Calculates square mile section risk####
#load point data tables
pointdata <- fread("Tables/Point Data (wells)/pdatadetailed2022-09-29.csv")

#7.1 Calculate source section average/recent results####
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
source_sections <- pointdata %>% 
  filter(domdepth == "yes") %>%
  left_join(., select(wellinfo, gm_well_id, MTRS)) %>%
  group_by(gm_chemical_vvl, MTRS) %>%
  summarize(SecDet = mean(well_avg),
            re_under = sum(re_under),
            re_close = sum(re_close),
            re_over = sum(re_over),
            re_avg_over = mean(re_avg_over[re_avg_over > 0])) %>%
  mutate(method = "source")

#7.2 Calculate neighbor section average/recent results####
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
rm(neighbor_sections, neighbor_sections_add, source_sections)

#calculate average magnitude over MCL
section_data <- section_data %>% mutate(avg_overMCL = ifelse(SecDet > 1 & re_avg_over > 0, (SecDet + re_avg_over)/2,
                                                             ifelse(SecDet > 1 & re_avg_over == 0, SecDet,
                                                                    ifelse(!SecDet > 1 & re_avg_over > 0, re_avg_over, 0))),
                                        WQriskbin = ifelse(SecDet > 1 | re_over > 0, "high",
                                                           ifelse(SecDet > 0.8 | re_close > 0, "med", "low")))
section_data[is.na(section_data)] <- 0

#7.3 Determine section risk####
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

fwrite(section_data, paste0("Tables/Square Mile Sections/sectiondata_detailed", Sys.Date(), ".csv"))
fwrite(sdrisk, paste0("Tables/Square Mile Sections/sectiondata_risk", Sys.Date(), ".csv"))

#7.4 get location data for state smalls and domestic wells####
#dwr domestic well records data download process:
  #connect to https://services.arcgis.com/aa38u6OgfNoCkTJ6/ArcGIS/rest/services/i07_WellCompletionReports_Exported_gdb/FeatureServer
  #save subset of well completion records using following parameters
    #B118WellUse == "Domestic"
    #Date Work Ended > 12/31/1969
    #Record Type == “WellCompletion/New/Production or Monitoring/NA”
  #export as .csv and save in "Datasets" folder
dwr_domesticwellrecords <- st_read("Datasets/i07_domestic_1970_completion_9_28_2022.shp")
MTRS <- st_read("DWrisk_newcount/DWrisk_newcount.shp")
dwr_domesticwellrecords_mtrs <- st_join(dwr_domesticwellrecords, MTRS["MTRS"], left = F) %>%
  st_drop_geometry()
colnames(dwr_domesticwellrecords_mtrs)[colnames(dwr_domesticwellrecords_mtrs)=="CountyName"] <- "NAME"

#dwr_domesticwellrecords <- read.csv("Datasets/i07_export_domestic_1970_new_validloc.csv")
#colnames(dwr_domesticwellrecords)[colnames(dwr_domesticwellrecords)=="CountyName"] <- "County"
mtrs_county <- read.csv("Tables/MTRS_Counties_new.csv") %>%
  #remove split objects
  select(MTRS, NAME, NAMELSAD) %>% distinct()
domwc_known <- dwr_domesticwellrecords_mtrs %>% filter(!NAME %in% c("", "Unknown"), !is.na(NAME)) %>% 
  group_by(MTRS.y, NAME) %>% count(name = "DWR_dom") %>% ungroup()
domwc_lost <- dwr_domesticwellrecords_mtrs %>% filter(NAME %in% c("", "Unknown") | is.na(NAME)) %>%
  group_by(MTRS.y) %>% count(name = "DWR_dom_nocounty") %>%
  left_join(., select(mtrs_county, MTRS, NAME), by = c("MTRS.y" = "MTRS")) %>%
  filter(!is.na(NAME)) %>% ungroup()#removes six domestic wells without a valid county
domwc_lost[duplicated(domwc_lost$MTRS.y),]
domwc <- domwc_known %>% full_join(., domwc_lost)
domwc$DWR_dom[is.na(domwc$DWR_dom)] <- 0
domwc$DWR_dom_nocounty[is.na(domwc$DWR_dom_nocounty)] <- 0
domwc <- domwc %>% group_by(MTRS.y, NAME) %>%
  mutate(DWR_dom_count = sum(DWR_dom, DWR_dom_nocounty)) %>%
  select(-c("DWR_dom", "DWR_dom_nocounty"))
colnames(domwc)[colnames(domwc)=="MTRS.y"] <- "MTRS"
fwrite(domwc, "Tables/domesticwells_mtrs_county_10_3_2022.csv")

#get state small water system count per section
ss_mtrs <- read.csv("Datasets/SSWS_mtrs_county_ddw2022.csv", stringsAsFactors = F) %>%
  distinct() %>% group_by(MTRS, County) %>%
  summarize(ssws_count = n())
fwrite(ss_mtrs, "Tables/ssws_mtrs_county.csv")

#put it all together (wq risk, domwellcounts, sswscounts):
section_data <- fread("Tables/Square Mile Sections/sectiondata_detailed2022-09-29.csv")
#section_data_expanded <- fread("Tables/Square Mile Sections/sectiondata_detailed_expanded2022-09-29.csv")
sdrisk <- fread("Tables/Square Mile Sections/sectiondata_risk2022-09-29.csv") %>%
  mutate(method = "depth_filter")
#sdrisk_expanded <- fread("Tables/Square Mile Sections/sectiondata_risk_expanded2022-09-29.csv") %>%
#  mutate(method = "no_depth_filter")
sdrisk_combined <- sdrisk
#sdrisk_combined <- rbind(sdrisk, sdrisk_expanded)
#note - some of these mtrs from dwr file are not in California!! so they will not appear when joined with mtrs_county
domwc <- fread("Tables/domesticwells_mtrs_county_10_3_2022.csv")
ss_mtrs <- fread("Tables/ssws_mtrs_county.csv")
mtrs_county <- read.csv("Tables/MTRS_Counties_new.csv") %>%
  #remove split objects
  select(COUNTYFP, NAME, NAMELSAD, MTRS) %>% distinct()

sdrisk_wells <- mtrs_county %>% select(MTRS, NAME, NAMELSAD) %>%
  left_join(., sdrisk_combined, by = "MTRS") %>%
  left_join(., domwc, by = c("MTRS", "NAME")) %>%
  left_join(., ss_mtrs, by = c("MTRS", "NAME")) %>%
  filter(!MTRS %in% c("", "BAY/DELTA", "SALTONSEA", "LAKETAHOE"))
sdrisk_wells$DWR_dom_count[is.na(sdrisk_wells$DWR_dom_count)] <- 0
sdrisk_wells$ssws_count[is.na(sdrisk_wells$ssws_count)] <- 0
sdrisk_wells$WQriskbins[is.na(sdrisk_wells$WQriskbins)] <- "unknown"
#sdrisk_wells <- sdrisk_wells %>% 
#  mutate(WQrisk_1 = ifelse(method == "no_depth_filter", "UNK", WQriskbins))
#sdrisk_wells$WQrisk_1[is.na(sdrisk_wells$WQrisk_1)] <- "UNK"

colnames(sdrisk_wells)[colnames(sdrisk_wells)=="WQriskbins"] <- "WQr_2023"
#colnames(sdrisk_wells)[colnames(sdrisk_wells)=="WQrisk_1"] <- "WQr_2023a"

arm_2021 <- fread("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2020/sectiondata_risk_20202021-06-14.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(arm_2021) <- c("MTRS", "WQr_2021")
arm_2022 <- fread("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_project/Webtool/tables/Square Mile Sections/sectiondata_risk2021-11-10.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(arm_2022) <- c("MTRS", "WQr_2022")

sdrisk_wells <- sdrisk_wells %>%
  left_join(., arm_2022) %>%
  left_join(., arm_2021)
sdrisk_wells$WQr_2022[is.na(sdrisk_wells$WQr_2022)] <- "unknown"
sdrisk_wells$WQr_2021[is.na(sdrisk_wells$WQr_2021)] <- "unknown"
sdrisk_wells$WQr_2023[sdrisk_wells$WQr_2023 == "med"] <- "medium"
sdrisk_wells$WQr_2022[sdrisk_wells$WQr_2022 == "med"] <- "medium"
sdrisk_wells$WQr_2021[sdrisk_wells$WQr_2021 == "med"] <- "medium"
fwrite(sdrisk_wells, "Tables/Square Mile Sections/ARM_2022-10-20.csv", na = "")

sdrisk_wells <- fread("Tables/Square Mile Sections/ARM_2022-10-20.csv") %>% select(-method)
MTRS <- st_read("Datasets/MTRS_Counties_new.shp") 
sdrisk_mtrs <- left_join(MTRS, sdrisk_wells)
sdrisk_mtrs[is.na(sdrisk_mtrs)] <- 0
st_write(sdrisk_mtrs, paste0("Tables/arm_sectionrisk", Sys.Date(), ".shp"))

#7.5 write individual contaminant shapefiles
MTRS <- st_read("Datasets/MTRS_Counties_new.shp")
indivchems <- c("NO3N", "AS", "CR6", "U", "TCPR123")
for (x in 1:length(indivchems)) {
  chem_sd <- section_data %>% filter(gm_chemical_vvl == indivchems[x])
  chem_sd_mtrs <- left_join(MTRS, chem_sd) %>% filter(!is.na(MTRS), !MTRS %in% c("", "BAY/DELTA", "SALTONSEA", "LAKETAHOE"))
  colnames(chem_sd_mtrs)[colnames(chem_sd_mtrs)=="WQriskbin"] <- "WQr_2023"
  chem_sd_mtrs$WQr_2023[is.na(chem_sd_mtrs$WQr_2023)] <- "unknown"
  chem_sd_mtrs$WQr_2023[chem_sd_mtrs$WQr_2023 == "med"] <- "medium"
  chem_sd_mtrs[is.na(chem_sd_mtrs)] <- 0
  st_write(chem_sd_mtrs, paste0("Tables/", indivchems[x], "_sectionrisk", Sys.Date(), ".shp"))
}
