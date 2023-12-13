#Water quality data download for 2024 Aquifer Risk Map

#Methodology write-up:

#Updated 12/8/2023 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

#0 Table of Contents####
# 1 Load required libraries, tables, and define special functions
# 2 Create well depth filter
# 3 Download water quality data
# 4 Standardize and clean data
  # 4.1 Special chemicals
  # 4.2 Standardize reporting limits/results
# 5 Load data tables
# 6 Calculate point risk scores
  # 6.1 Calculate well averages
  # 6.2 Calculate recent results
  # 6.3 Calculate point risk
# 7 Calculate square mile section risk
  # 7.1 Calculate source section average/recent results
  # 7.2 Calculate neighbor section average/recent results
  # 7.3 Calculate section risk
  # 7.4 Get location data for SSWS and DW
  # 7.5 Write shapefile for all section data
  # 7.6 Write shapefiles for individual contaminants
# 8 Calculate demographic statistics (block group/census tract)

#1. Load required libraries, tables, and define special functions####
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
if (!require('survival')) install.packages('survival'); library('survival')
#if (!require('rgdal')) install.packages('rgdal'); library('rgdal')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('janitor')) install.packages('janitor'); library('janitor')
#if (!require('httr')) install.packages('httr'); library('httr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
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
#https://gamagroundwater.waterboards.ca.gov/gama/translation_table.asp
chem_list <- read.csv("Reference_files/GAMAchemicalmapping.csv", stringsAsFactors = F) %>% clean_names() %>%
  filter(comparison_concentration_type == "MCL" | chemical_vvl %in% c("CU", "PB", "NNSM", "CR6", "PFOA", "PFOS", "PFBSA", "PFHXSA")) %>%
  filter(!chemical_vvl %in% c("ASBESTOS", "COLIFORM", "FCOLIFORM", "RN-222")) %>%
  select(chemical_vvl, chemical_name, units, comparison_concentration_value, comparison_concentration_type)

chem_list$comparison_concentration_value[chem_list$chemical_vvl == "CR6"] <- 10

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
#https://utility.arcgis.com/usrsvcs/servers/08dce2ff86434cfe8582ccdd2333b7fa/rest/services/Environment/i07_WellCompletionReports/MapServer
#https://gis.water.ca.gov/arcgis/rest/services/Environment/i07_WellCompletionReports/FeatureServer/1
url <- list(hostname = "gis.water.ca.gov/arcgis/rest/services",
            scheme = "https",
            path = "Environment/i07_WellCompletionReports/FeatureServer/1/query",
            query = list(where = "1=1",
                         outFields = "*",
                         returnGeometry = "true",
                         f = "geojson")) %>%
  setattr("class","url")
request <- build_url(url)
i07 <- st_read(request) %>% st_transform(crs = 4326) #%>% select(MTRS) %>% distinct()
#MTRS <- st_read("Datasets/i07_WellReportStatsBySection_1009.shp") %>% st_transform(crs = 4326) #%>% select(MTRS) %>% distinct()

#join mtrs to groundwater units, then drop geometry
MTRS_GU <- st_join(i07, GU) %>% st_drop_geometry()

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
write.table(GU_depths, paste0("Reference_files/well_depth_gwunits", Sys.Date(), ".txt"), sep = "\t", row.names = F)

#3. Download water quality data ####
GU_depths <- read.table("Reference_files/well_depth_gwunits2023-08-01.txt", stringsAsFactors = F, header = T)

col_sty <- cols(
  GM_DATASET_NAME = col_character(),
  GM_WELL_CATEGORY = col_character(),
  GM_DATA_SOURCE = col_character(),
  GM_WELL_ID = col_character(),
  GM_CHEMICAL_VVL = col_character(),
  GM_CHEMICAL_NAME = col_character(),
  GM_RESULT_MODIFIER = col_character(),
  GM_RESULT = col_double(),
  GM_CHEMICAL_UNITS = col_character(),
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
  SRC_CHEMICAL_UNITS = col_character(),
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

col_sty_loc <- cols(
  GM_DATASET_NAME = col_character(),
  GM_WELL_CATEGORY = col_character(),
  GM_DATA_SOURCE = col_character(),
  GM_WELL_ID = col_character(),
  GM_LATITUDE = col_double(),
  GM_LONGITUDE = col_double(),
  GM_WELL_DEPTH_FT = col_double(),
  GM_TOP_DEPTH_OF_SCREEN_FT = col_double(),
  GM_BOTTOM_DEPTH_OF_SCREEN_FT = col_double(),
  GM_GIS_COUNTY = col_character(),
  GM_GIS_DWR_BASIN = col_character(),
  GM_GIS_REGIONAL_BOARD = col_character(),
  GM_GIS_SENATE_DISTRICT = col_double(),
  GM_GIS_HVA = col_character(),
  GM_GIS_GAMA_STUDY_AREA = col_character(),
  GM_GIS_ASSEMBLY_DISTRICT = col_double(),
  GM_GIS_GSA = col_character(),
  GM_GIS_DWR_REGION = col_character()
)

#identify wells that are within domestic depth filter
#download location data
temp <- tempfile()
download.file('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_location_construction_v2.zip', temp)
loc <- readr::read_tsv(unz(temp, 'gama_location_construction_v2.txt'), col_types = col_sty_loc) %>% clean_names() %>% distinct() #%>%
unlink(temp)

#well "AGW080016134-14805\xa0S.\xa0M" and "L10007696901-VLT\xdbFLDNON" have unusual characters. Address this here in locations and later in wq download
#loc[loc$gm_well_id == "AGW080016134-14805\xa0S.\xa0M", 4] <- "AGW080016134-14805"
#loc[loc$gm_well_id == "L10007696901-VLT\xdbFLDNON", 4] <- "L10007696901-VLT"

#create PWSID for DDW wells
ddw_pwsid <- loc %>% filter(gm_dataset_name == "DDW") %>%
  select(gm_dataset_name, gm_well_id) %>% 
  mutate(PWSID = ifelse(gm_dataset_name == "DDW", substr(gm_well_id, 1, 9), NA))
loc <- left_join(loc, ddw_pwsid)

#identify state small water system wells (not available as a well category on GAMA, but I use the list of SSWS I get from DDW)
ssws_list <- readxl::read_excel("Data/SSWS_10_17_23.xlsx")

#get MTRS shapefile
MTRS <- st_read("Reference_files/MTRS_joined_clean.shp") %>% st_transform(crs = 4326) %>% st_make_valid() %>% distinct()

#convert well locations to sf object and join with GU and MTRS locations
wells_int <- st_as_sf(loc, coords = c('gm_longitude', 'gm_latitude'), crs = 4326) %>%
  st_join(., left = TRUE, GU["GU_ID"]) %>% 
  st_join(., left = TRUE, MTRS["MTRS"]) %>% 
  mutate(gm_latitude = sf::st_coordinates(.)[,2],
         gm_longitude = sf::st_coordinates(.)[,1]) %>%
  st_drop_geometry() %>%
  select("gm_well_id", "gm_well_category", "gm_latitude", "gm_longitude", "gm_well_depth_ft",
         "gm_top_depth_of_screen_ft", "gm_bottom_depth_of_screen_ft", "GU_ID", "gm_dataset_name", "gm_data_source", "PWSID", "MTRS") %>%
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
  mutate(domdepth = ifelse(gm_well_category == "DOMESTIC" | gm_dataset_name == "GAMA_DOM" | PWSID %in% ssws_list$PWSID, "yes", #keep "domestic" and "state small water system" wells
                           ifelse(gm_well_category == "MONITORING",                                       #special filter for monitoring wells
                                  ifelse(is.na(depth), "no",                                             #remove all MON wells without numeric depth
                                         ifelse(!is.na(maxdomestic) & depth <= maxdomestic & depth > mindomestic & depth > 100, "yes", "no")), #min 100 ft depth for WB cleanup wells
                                  ifelse(is.na(use_public_as_dom), "yes",                                       #if basin has no domestic well depth data, use all wq wells
                                         ifelse(is.na(depth),                                                        #for non-MON wells without numeric depth
                                                ifelse(use_public_as_dom == "yes", "yes", "no"),                           #use public/domestic comparison by basin
                                                ifelse(depth <= maxdomestic & depth > mindomestic, "yes",                 #for non-GT wells with numeric depth, use numeric depth filter
                                                       "no"))))))
write.table(allwells, paste0("Data/allwells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

allwells <- read.table("Data/allwells2023-10-25.txt", stringsAsFactors = F, header = T)
domwells <- allwells %>% filter(domdepth == "yes")

#retrieve water quality data
#define function to download data
dwd_wq_tbl <- function(link, filename, col_str, my_wells, mychems) {
  temp <- tempfile()
  download.file(link, temp)
  mytable <- readr::read_tsv(unz(temp, filename), quote = "", col_types = col_str) %>% clean_names()
  unlink(temp)
  mytable$gm_samp_collection_date <- lubridate::mdy(mytable$gm_samp_collection_date)
  mytable <- mytable %>% filter(gm_samp_collection_date >= ymd("2003-01-01"),
                                gm_well_id %in% my_wells$gm_well_id,
                                gm_chemical_vvl %in% mychems$chemical_vvl)
  return(mytable)
}

#download various datasets for domestic depth wells
#too large, download manually and read from file 'https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_ddw_statewide_v2.zip'
ddw <- fread("Reference_files/gama_ddw_statewide_v2.txt") %>% clean_names()
ddw$gm_samp_collection_date <- lubridate::mdy(ddw$gm_samp_collection_date)
ddw <- ddw %>% filter(gm_samp_collection_date >= ymd("2003-01-01"),
                              gm_well_id %in% domwells$gm_well_id,
                              gm_chemical_vvl %in% chem_list$chemical_vvl)
ddw <- ddw %>% filter(!str_detect(gm_altwell_id1, "SPRING") | !str_detect(gm_altwell_id2, "SPRING")) #remove ~240 spring sites
ddw_wells <- ddw %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(ddw, paste0("Data/ddw", Sys.Date(), ".csv"))
fwrite(ddw_wells, paste0("Data/ddw_wells", Sys.Date(), ".csv"))
rm(ddw, ddw_wells)

dwr <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_dwr_statewide_v2.zip',
                  'gama_dwr_statewide_v2.txt', col_sty, domwells, chem_list)
dwr_wells <- dwr %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(dwr, paste0("Data/dwr", Sys.Date(), ".csv"))
fwrite(dwr_wells, paste0("Data/dwr_wells", Sys.Date(), ".csv"))
rm(dwr, dwr_wells)

dpr <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_dpr_statewide_v2.zip',
                  'gama_dpr_statewide_v2.txt', col_sty, domwells, chem_list)
dpr_wells <- dpr %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(dpr, paste0("Data/dpr", Sys.Date(), ".csv"))
fwrite(dpr_wells, paste0("Data/dpr_wells", Sys.Date(), ".csv"))
rm(dpr, dpr_wells)

gamadom <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_gama_dom_statewide_v2.zip',
                      'gama_gama_dom_statewide_v2.txt', col_sty, domwells, chem_list)
gamadom_wells <- gamadom %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(gamadom, paste0("Data/gamadom", Sys.Date(), ".csv"))
fwrite(gamadom_wells, paste0("Data/gamadom_wells", Sys.Date(), ".csv"))
rm(gamadom, gamadom_wells)

gama_localgw <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_gama_localgw_statewide_v2.zip',
                           'gama_gama_localgw_statewide_v2.txt', col_sty, domwells, chem_list)
gama_localgw_wells <- gama_localgw %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(gama_localgw, paste0("Data/gama_localgw", Sys.Date(), ".csv"))
fwrite(gama_localgw_wells, paste0("Data/gama_localgw_wells", Sys.Date(), ".csv"))
rm(gama_localgw, gama_localgw_wells)

gamaspstudy <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_gama_sp-study_statewide_v2.zip',
                          'gama_gama_sp-study_statewide_v2.txt', col_sty, domwells, chem_list)
gamaspstudy_wells <- gamaspstudy %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(gamaspstudy, paste0("Data/gamaspstudy", Sys.Date(), ".csv"))
fwrite(gamaspstudy_wells, paste0("Data/gamaspstudy_wells", Sys.Date(), ".csv"))
rm(gamaspstudy, gamaspstudy_wells)

gamausgs <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_gama_usgs_statewide_v2.zip',
                       'gama_gama_usgs_statewide_v2.txt', col_sty, domwells, chem_list)
gamausgs_wells <- gamausgs %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(gamausgs, paste0("Data/gamausgs", Sys.Date(), ".csv"))
fwrite(gamausgs_wells, paste0("Data/gamausgs_wells", Sys.Date(), ".csv"))


localgw <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_localgw_statewide_v2.zip',
                      'gama_localgw_statewide_v2.txt', col_sty, domwells, chem_list)
localgw_wells <- localgw %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(localgw, paste0("Data/localgw", Sys.Date(), ".csv"))
fwrite(localgw_wells, paste0("Data/localgw_wells", Sys.Date(), ".csv"))
rm(localgw, localgw_wells)

usgsnwis <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_usgs_nwis_statewide_v2.zip',
                       'gama_usgs_nwis_statewide_v2.txt', col_sty, domwells, chem_list)
#remove results that are a duplication in gama_usgs dataset
usgsnwis_small <- usgsnwis %>% mutate(temp_id = gm_altwell_id1)
gamausgs_small <- gamausgs %>%  mutate(temp_id = gm_well_id)
usgsnwis_unique <- anti_join(usgsnwis_small, gamausgs_small, by = c("temp_id", "gm_chemical_vvl", "gm_samp_collection_date")) %>% select(-temp_id)
usgsnwis_wells <- usgsnwis_unique %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(usgsnwis_unique, paste0("Data/usgsnwis", Sys.Date(), ".csv"))
fwrite(usgsnwis_wells, paste0("Data/usgsnwis_wells", Sys.Date(), ".csv"))
rm(usgsnwis, usgsnwis_small, gamausgs_small, usgsnwis_unique, usgsnwis_wells)
rm(gamausgs, gamausgs_wells)

wbilrp <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_wb_ilrp_statewide_v2.zip',
                     'gama_wb_ilrp_statewide_v2.txt', col_sty, domwells, chem_list)
wbilrp_wells <- wbilrp %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(wbilrp, paste0("Data/wbilrp", Sys.Date(), ".csv"))
fwrite(wbilrp_wells, paste0("Data/wbilrp_wells", Sys.Date(), ".csv"))
rm(wbilrp, wbilrp_wells)

wrd <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_wrd_statewide_v2.zip',
                  'gama_wrd_statewide_v2.txt', col_sty, domwells, chem_list)
wrd_wells <- wrd %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(wrd, paste0("Data/wrd", Sys.Date(), ".csv"))
fwrite(wrd_wells, paste0("Data/wrd_wells", Sys.Date(), ".csv"))
rm(wrd, wrd_wells)

ucd_no3 <- dwd_wq_tbl('https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_ucd_no3_statewide_v2.zip',
                      'gama_ucd_no3_statewide_v2.txt', col_sty, domwells, chem_list)
#remove values that may be a duplication of other datasets
ucd_no3_unique <- ucd_no3 %>% filter(gm_data_source %in% c("UCD_NO3 - KECO", "UCD_NO3 - MOCO", "UCD_NO3 - MCEM",
                                                            "UCD_NO3 - FRCO", "UCD_NO3 - TCEHS", "UCD_NO3 - KICO",
                                                            "UCD_NO3 - TUCO", "UCD_NO3 - KRB98"))
ucd_no3_wells <- ucd_no3_unique %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
fwrite(ucd_no3_unique, paste0("Data/ucdno3", Sys.Date(), ".csv"))
fwrite(ucd_no3_wells, paste0("Data/ucdno3_wells", Sys.Date(), ".csv"))
rm(ucd_no3, ucd_no3_unique, ucd_no3_wells)

#download 'https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_wb_cleanup_statewide_v2.zip'
edf <- fread("Reference_files/gama_wb_cleanup_statewide_v2.txt") %>% clean_names() %>%
  filter(gm_well_id %in% domwells$gm_well_id,
         gm_chemical_vvl %in% chem_list$chemical_vvl)
edf$gm_samp_collection_date <- lubridate::mdy(edf$gm_samp_collection_date)
edf <- edf %>% filter(gm_samp_collection_date >= ymd("2003-01-01"))

edfwqwells <- edf %>% select(gm_well_id, gm_dataset_name) %>% distinct() %>% left_join(., domwells)
#save water quality data and well information
fwrite(edfwqwells, paste0("Data/edf_wells", Sys.Date(), ".csv"))
fwrite(edf, paste0("Data/edf", Sys.Date(), ".csv"))
rm(edf, edfwqwells)

#compile datasets
date_dwld <- c("2023-10-25")

ddw <- fread(paste0("Data/ddw", date_dwld, ".csv"))
ddw_wells <- fread(paste0("Data/ddw_wells", date_dwld, ".csv"))
dpr <- fread(paste0("Data/dpr", date_dwld, ".csv"))
dpr_wells <- fread(paste0("Data/dpr_wells", date_dwld, ".csv"))
dwr <- fread(paste0("Data/dwr", date_dwld, ".csv"))
dwr_wells <- fread(paste0("Data/dwr_wells", date_dwld, ".csv"))
gamadom <- fread(paste0("Data/gamadom", date_dwld, ".csv"))
gamadom_wells <- fread(paste0("Data/gamadom_wells", date_dwld, ".csv"))
gama_localgw <- fread(paste0("Data/gama_localgw", date_dwld, ".csv"))
gama_localgw_wells <- fread(paste0("Data/gama_localgw_wells", date_dwld, ".csv"))
gamaspstudy <- fread(paste0("Data/gamaspstudy", date_dwld, ".csv"))
gamaspstudy_wells <- fread(paste0("Data/gamaspstudy_wells", date_dwld, ".csv"))
gamausgs <- fread(paste0("Data/gamausgs", date_dwld, ".csv"))
gamausgs$gm_altwell_id1 <- as.character(gamausgs$gm_altwell_id1)
gamausgs_wells <- fread(paste0("Data/gamausgs_wells", date_dwld, ".csv"))
localgw <- fread(paste0("Data/localgw", date_dwld, ".csv"))
localgw_wells <- fread(paste0("Data/localgw_wells", date_dwld, ".csv"))
usgsnwis <- fread(paste0("Data/usgsnwis", date_dwld, ".csv"))
usgsnwis_wells <- fread(paste0("Data/usgsnwis_wells", date_dwld, ".csv"))
wbilrp <- fread(paste0("Data/wbilrp", date_dwld, ".csv"))
wbilrp_wells <- fread(paste0("Data/wbilrp_wells", date_dwld, ".csv"))
wrd <- fread(paste0("Data/wrd", date_dwld, ".csv"))
wrd_wells <- fread(paste0("Data/wrd_wells", date_dwld, ".csv"))
ucd_no3 <- fread(paste0("Data/ucdno3", date_dwld, ".csv"))
ucd_no3_wells <- fread(paste0("Data/ucdno3_wells", date_dwld, ".csv"))
edf <- fread(paste0("Data/edf", date_dwld, ".csv"))
edf_wells <- fread(paste0("Data/edf_wells", date_dwld, ".csv"))

domwqdata <- rbind(ddw, dpr, dwr, gamadom, gama_localgw, gamaspstudy, gamausgs, localgw, usgsnwis, wbilrp, wrd, ucd_no3, edf)
domwells_wq <- rbind(ddw_wells, dpr_wells, dwr_wells, gamadom_wells, gama_localgw_wells, gamaspstudy_wells, gamausgs_wells,
                    localgw_wells, usgsnwis_wells, wbilrp_wells, wrd_wells, ucd_no3_wells, edf_wells)
rm(ddw, dpr, dwr, gamadom, gama_localgw, gamaspstudy, gamausgs, localgw, usgsnwis, wbilrp, wrd, ucd_no3, edf)
rm(ddw_wells, dpr_wells, dwr_wells, gamadom_wells, gama_localgw_wells, gamaspstudy_wells, gamausgs_wells, localgw_wells, usgsnwis_wells, wbilrp_wells, wrd_wells, ucd_no3_wells, edf_wells)

#remove well not in CA
domwqdata <- domwqdata %>% filter(!gm_well_id == "USGS-415953121292901") #not in california
domwells_wq <- domwells_wq %>% filter(!gm_well_id == "USGS-415953121292901")

#save water quality data and well information
fwrite(domwqdata, paste0("Data/wqdata_", Sys.Date(), ".csv"))
fwrite(domwells_wq, paste0("Data/wqdata_wells_", Sys.Date(), ".csv"))

#4. Standardize and clean data ####
domwqdata <- fread(paste0("Data/wqdata_", date_dwld, ".csv")) %>% 
  select(gm_well_id, gm_dataset_name, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit) %>%
  #deal with wells with multiple lat/longs
  distinct()
domwells_wq <- fread(paste0("Data/wqdata_wells_", date_dwld, ".csv"))
dupe_locs <- domwells_wq %>% group_by(gm_well_id, gm_dataset_name) %>% count() %>% filter(n > 1)
#remove duplicate entries in GAMA (multiple locations for one well ID)
domwells_wq_clean <- domwells_wq %>%
  distinct(gm_well_id, gm_dataset_name, .keep_all= TRUE)

#read dates as dates
domwqdata$gm_samp_collection_date <- lubridate::ymd(domwqdata$gm_samp_collection_date)

#identify non-quantified results
domwqdata <- domwqdata %>%
  #convert result values to simplified form (null, zero, or nonzero)
  mutate(simple_res = case_when(
    is.na(gm_result) ~ "NL",
    gm_result == "" ~ "NL",
    gm_result == 0 ~ "zero",
    !gm_result == 0 ~ "nonzero",
    TRUE ~ "HELP")) %>%
  mutate(new_modifier = case_when(
    #modifier 1.1 "NULL"
    is.na(gm_result_modifier) & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    is.na(gm_result_modifier) & simple_res %in% c("nonzero") ~ "Q",
    #modifier 1 "-, =, P, X, '', ' '"
    gm_result_modifier %in% c("-", "=", "P", "X", "", " ") & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_result_modifier %in% c("-", "=", "P", "X", "", " ") & simple_res %in% c("nonzero") ~ "Q",
    #modifier 2 "<, F, Y, ND, M, <=, OS, DN"
    gm_result_modifier %in% c("<", "F", "Y", "ND", "M", "<=", "OS", "DN") ~ "NQ_LT",
    #modifier 3 ">, >="
    gm_result_modifier %in% c(">", ">=") ~ "GNQ",
    #modifier 4 "E, S, A, TI"
    gm_result_modifier %in% c("E", "S", "A", "TI") & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_result_modifier %in% c("E", "S", "A", "TI") & simple_res %in% c("nonzero") & gm_result <= gm_reporting_limit ~ "NQ_LT",
    gm_result_modifier %in% c("E", "S", "A", "TI")  & simple_res %in% c("nonzero") ~ "GNQ",
    #modifier 5 "N"
    gm_dataset_name == "DDW" & gm_result_modifier == "N" & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_dataset_name == "DDW" & gm_result_modifier == "N" & simple_res %in% c("nonzero") ~ "Q",
    gm_dataset_name == "DPR" & gm_result_modifier == "N" ~ "NQ",
    #modifier 6 "V"
    gm_dataset_name == "DPR" & gm_result_modifier == "V" & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_dataset_name == "DPR" & gm_result_modifier == "V" & simple_res %in% c("nonzero") ~ "Q",
    #modifier 7 "R, U" GAMA/NWIS only
    gm_dataset_name %in% c("USGS_NWIS", "GAMA_USGS") & gm_result_modifier %in% c("R", "U") ~ "NQ_LT",
    TRUE ~ "REJECT")
  )
domwqdata <- domwqdata %>% filter(!gm_result_modifier %in% c("I", "Q", "NA", "RF", "SU", "DU", "IN", "NR", "PA")) %>%
  filter(!new_modifier == "REJECT")

#4.1 special chemicals####
#NITROGEN: convert Nitrate/Nitrite combined to Nitrate
no3n <- domwqdata %>% filter(gm_chemical_vvl %in% c("NO3N")) %>% select(-gm_chemical_vvl)
colnames(no3n) <- c("gm_well_id", "gm_dataset_name", "NO3N_result_modifier", "NO3N_result", "gm_samp_collection_date", "NO3N_reporting_limit",
                     "NO3N_simple_res", "NO3N_new_modifier")
no2n <- domwqdata %>% filter(gm_chemical_vvl %in% c("NO2")) %>% select(-gm_chemical_vvl)
colnames(no2n) <- c("gm_well_id", "gm_dataset_name", "NO2_result_modifier", "NO2_result", "gm_samp_collection_date", "NO2_reporting_limit",
                     "NO2_simple_res", "NO2_new_modifier")
no3no2n <- domwqdata %>% filter(gm_chemical_vvl %in% c("NO3NO2N")) %>% select(-gm_chemical_vvl)
colnames(no3no2n) <- c("gm_well_id", "gm_dataset_name", "NO3NO2N_result_modifier", "NO3NO2N_result", "gm_samp_collection_date", "NO3NO2N_reporting_limit",
                     "NO3NO2N_simple_res", "NO3NO2N_new_modifier")
#join no3n, no2n and no3no2n files
nitrogen <- full_join(no3n, no2n) %>% 
  full_join(., no3no2n) %>% distinct() %>%
  filter(!is.na(NO3NO2N_new_modifier), is.na(NO3N_new_modifier)) %>%
  #for each row
  rowwise() %>%
  #create new result, rl, and modifier columns
  mutate(gm_result = case_when(
    #if there is no no2n, use no3no2n value
    is.na(NO2_new_modifier) ~ NO3NO2N_result,
    #if there is no2n value, subtract no2n from no3no2n
    !is.na(NO2_new_modifier) ~ NO3NO2N_result - NO2_result,
    TRUE ~ -99),
  gm_reporting_limit = case_when(
    is.na(NO2_new_modifier) ~ NO3NO2N_reporting_limit,
    !is.na(NO2_new_modifier) ~ max(NO3NO2N_reporting_limit, NO2_reporting_limit, na.rm = T),
    TRUE ~ -99),
  #never the case where NO2 is quantified but NO3NO2N isn't (but NO3NO2N can be quantified and NO2 not) so always take NO3NO2N modifier
  gm_result_modifier = NO3NO2N_result_modifier,
  new_modifier = NO3NO2N_new_modifier) %>%
  mutate(gm_chemical_vvl = "NO3N",
         simple_res = NA) %>%
  select(gm_well_id, gm_dataset_name, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit, simple_res, new_modifier) %>% ungroup()
domwqdata <- rbind(domwqdata, nitrogen) %>% filter(!gm_chemical_vvl == "NO3NO2N")
rm(nitrogen, no2n, no3n, no3no2n)

#RADIUM: combine radium-226 and radium-228 (if both sampled on the same well on same date) (mcl is based on combined total)
ra226 <- domwqdata %>% filter(gm_chemical_vvl %in% c("RA-226")) %>% select(-gm_chemical_vvl)
colnames(ra226) <- c("gm_well_id", "gm_dataset_name", "RA226_result_modifier", "RA226_result", "gm_samp_collection_date", "RA226_reporting_limit",
                     "RA226_simple_res", "RA226_new_modifier")
ra228 <- domwqdata %>% filter(gm_chemical_vvl %in% c("RA-228")) %>% select(-gm_chemical_vvl)
colnames(ra228) <- c("gm_well_id", "gm_dataset_name", "RA228_result_modifier", "RA228_result", "gm_samp_collection_date", "RA228_reporting_limit",
                     "RA228_simple_res", "RA228_new_modifier")
#join ra226 and ra228 files
ra226_228 <- full_join(ra226, ra228) %>% distinct() %>%
  #for each row
  rowwise() %>%
  #create new result, rl, and modifier columns
  mutate(gm_result = case_when(
    #if there is no ra228 value, use ra226 value
    is.na(RA228_new_modifier) ~ RA226_result,
    #if there is no ra226 value, use ra228 value
    is.na(RA226_new_modifier) ~ RA228_result,
    #if ra226 and ra228 values are not quantified, take the maximum non-quantified result
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" ~ max(RA226_result, RA228_result, na.rm = T),
    #if ra226 is not quantified but ra228 is quantified, use ra228 result
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_result,
    #same but flipped
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_result,
    #if ra226 and ra228 values are quantified, take the sum as the result
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" ~ sum(RA226_result, RA228_result, na.rm = T),
    TRUE ~ -99
  ),
  gm_reporting_limit = case_when(
    is.na(RA228_new_modifier) ~ RA226_reporting_limit,
    is.na(RA226_new_modifier) ~ RA228_reporting_limit,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" ~ max(RA226_reporting_limit, RA228_reporting_limit, na.rm = T),
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" ~ max(RA226_reporting_limit, RA228_reporting_limit, na.rm = T),
    TRUE ~ -99
  ),
  gm_result_modifier = case_when(
    is.na(RA228_new_modifier) ~ RA226_result_modifier,
    is.na(RA226_new_modifier) ~ RA228_result_modifier,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" ~ "<",
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_result_modifier,
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_result_modifier,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" ~ "=",
    TRUE ~ "HELP"
  ),
  new_modifier = case_when(
    is.na(RA228_new_modifier) ~ RA226_new_modifier,
    is.na(RA226_new_modifier) ~ RA228_new_modifier,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" ~ "NQ_LT",
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_new_modifier,
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_new_modifier,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" ~ "Q",
    TRUE ~ "HELP"
  ))

ra226_228_2 <- ra226_228 %>%
  mutate(gm_chemical_vvl = "RA-226-228",
         simple_res = NA) %>%
  select(gm_well_id, gm_dataset_name, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit, simple_res, new_modifier)
domwqdata <- rbind(domwqdata, ra226_228_2) %>% filter(!gm_chemical_vvl %in% c("RA-226", "RA-228"))
rm(ra226, ra228, ra226_228, ra226_228_2)

#edit chem_list to include combined radium vvl
chem_list[102, 1] <- "RA-226-228"
chem_list[102, 2] <- "Combined Radium 226 and Radium 228"
chem_list[102, 3] <- "pCi/L"
chem_list[102, 4] <- 5
chem_list[102, 5] <- "MCL"

chem_list2 <- chem_list %>% filter(!chemical_vvl %in% c("NO3NO2N", "RA-226", "RA-228"))

#4.2 standardize Reporting Limits/Results (fix Null results, results == 0, and missing reporting limits)####
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
  #join with chemical MCLs
  left_join(., select(chem_list2, chemical_vvl, comparison_concentration_value), by = c("gm_chemical_vvl" = "chemical_vvl"))      


#calculate MCL index and track detections/non-detections, outliers
domwqdata <- domwqdata %>% 
  #calculate numeric value for all results (with caveat for non-detects that are greater than CCV)
  mutate(RES = ifelse(new_modifier %in% c("NQ_LT", "GNQ"), RLfun(RL, comparison_concentration_value), gm_result),
         #calculate MCL index
         MCLindex = RES/comparison_concentration_value)                                            
md_max <- summarize(group_by(domwqdata %>% dplyr::filter(new_modifier == "Q"), gm_chemical_vvl), avg = mean(RES), sd = sd(RES)) %>% 
  #calculate outlier maximum
  mutate(max = (avg + (10*sd))) %>% select(gm_chemical_vvl, max)               
write.csv(md_max, paste0("Data/outliers", Sys.Date(), ".csv"), row.names = F)
domwqdata <- left_join(domwqdata, md_max, by = "gm_chemical_vvl") %>% 
  mutate(outlier_status = ifelse(is.na(max), "no outliers", 
                                 ifelse(RES > max & new_modifier == "Q", "outlier", "standard")))

#clean table and output results
fwrite(distinct(domwqdata), paste0("Data/wqdata_clean", Sys.Date(), ".csv"))
domwqdata_sm <- domwqdata %>% select(gm_well_id, gm_dataset_name, gm_samp_collection_date, gm_chemical_vvl,
                                           RES, RL, RLmethod, new_modifier, outlier_status, MCLindex) %>%
  left_join(., select(domwells_wq_clean, gm_well_id, gm_dataset_name, domdepth, MTRS)) #attach MTRS for easy filtering later
fwrite(domwqdata_sm, paste0("Data/wqdata_small_", Sys.Date(), ".csv"))
fwrite(domwells_wq_clean, paste0("Data/wqdata_wellsclean_", Sys.Date(), ".csv"))

#5. Load data tables####
domwqdata_sm <- fread(paste0("Data/wqdata_small_", date_dwld, ".csv"))
wellinfo <- fread(paste0("Data/wqdata_wellsclean_", date_dwld, ".csv")) %>% 
  select(gm_well_id, domdepth, gm_latitude, gm_longitude, gm_well_category, gm_dataset_name, MTRS)

#6. Calculate point risk scores (wells)####
domwqdata_sm$gm_samp_collection_date <- lubridate::ymd(domwqdata_sm$gm_samp_collection_date) #read dates as dates
domwqdata_sm <- domwqdata_sm %>% 
  filter(!outlier_status == "outlier", !is.na(MCLindex)) %>% #remove outliers and samples without values (?)
  select(gm_well_id, gm_dataset_name, gm_samp_collection_date, gm_chemical_vvl, MCLindex) %>% #select only necessary columns to calculate averages
  mutate(year = year(gm_samp_collection_date)) #prepare for grouping samples by year
#6.1 calculate well averages ####
domwellavg <- domwqdata_sm %>% 
  group_by(gm_chemical_vvl, gm_well_id, gm_dataset_name, year) %>% #group samples by chemical, well, and year
  summarize(yr_avg = mean(MCLindex))%>%
  group_by(gm_chemical_vvl, gm_well_id, gm_dataset_name) %>% #group samples by chemical, well
  summarize(well_avg = mean(yr_avg))

#6.2 calculate recent results####
domwellre_avg <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2018-01-01"), MCLindex > 1) %>% 
  group_by(gm_well_id, gm_dataset_name, gm_chemical_vvl) %>%
  summarize(re_avg_over = mean(MCLindex))
domwellrecent <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2018-01-01")) %>%
  mutate(re_status = ifelse(MCLindex > 1, 'over',
                            ifelse(MCLindex > 0.8, 'close', 'under'))) %>% 
  group_by(gm_chemical_vvl, gm_well_id, gm_dataset_name, re_status) %>% count() %>%
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

pdrisk <- pointdata %>% group_by(gm_well_id, gm_dataset_name) %>% arrange(-avg_overMCL) %>%
  summarize(PRF1 = n_distinct(gm_chemical_vvl[well_avg > 1 | re_over > 0]),
            PRF2 = n_distinct(gm_chemical_vvl[well_avg <= 1 & well_avg > 0.8 | re_close > 0 & re_over == 0]),
            PRF3 = mean(avg_overMCL[avg_overMCL > 0], na.rm = T),
            PL1 = paste(gm_chemical_vvl[well_avg > 1 | re_over > 0], collapse = "; "),
            PL2 = paste(gm_chemical_vvl[well_avg <= 1 & well_avg > 0.8 | re_close > 0 & re_over == 0], collapse = "; ")) %>% 
  mutate(risk = ifelse(PRF1 > 0, "high", ifelse(PRF1 == 0 & PRF2 > 0, "medium", "low"))) %>%
  left_join(., wellinfo)

fwrite(pointdata, paste0("Tables/pdatadetailed.csv"))
fwrite(pdrisk, paste0("Tables/pdata_risk.csv"))

#7. Calculates square mile section risk####
#load point data tables
pointdata <- fread("Tables/pdatadetailed.csv")
wellinfo <- fread(paste0("Data/wqdata_wellsclean_", date_dwld, ".csv")) %>% 
  select(gm_well_id, domdepth, gm_latitude, gm_longitude, gm_well_category, gm_dataset_name, MTRS) %>% distinct()

#7.1 Calculate source section average/recent results####
#load file that contains neighbor reference data - not working for some reason - can contact Emily Houlihan for reference file
#temp <- tempfile()
#download.file("https://gispublic.waterboards.ca.gov/portal/sharing/rest/content/items/05f49f6f8bd24ac38143bdda402ba098/data", temp)
#neighborsec1 <- read.csv(temp)
#unlink(temp)
#temp <- tempfile()
#download.file("https://gispublic.waterboards.ca.gov/portal/sharing/rest/content/items/a72861be691643359e62e78d2c1a8347/data", temp)
#neighborsec2 <- read.csv(temp)
#unlink(temp)

#neighborsec <- rbind(neighborsec1, neighborsec2) %>% distinct() %>% filter(!MTRS_1 == MTRS_2)
#colnames(neighborsec) <- c("MTRS", "MTRS_2")


neighborsec <- fread("Reference_files/MTRS_source_neighbors.csv")

#group point data by MTRS (source data first)
source_sections <- pointdata %>% 
  filter(domdepth == "yes") %>% #annoyingly some datasets use the same well id, so when my original wq data pull grabbed "well id's in the domestic depth well list" it got some extras. they are missing the location join though (lat/long and dom depth)
  left_join(., select(wellinfo, gm_well_id, gm_dataset_name, MTRS)) %>%
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
section_data <- section_data %>% mutate(avg_overMCL = ifelse(SecDet > 1 & re_over >= 1, (SecDet + re_avg_over)/2,
                                                             ifelse(SecDet > 1 & re_over < 1, SecDet,
                                                                    ifelse(!SecDet > 1 & re_over >= 1, re_avg_over, 0))),
                                        WQriskbin = ifelse(SecDet > 1 | re_over >= 1, "high",
                                                           ifelse(SecDet > 0.8 | re_close >= 1 | re_over > 0, "medium", "low")))
section_data[is.na(section_data)] <- 0

#7.3 Determine section risk####
sdrisk <- section_data %>% 
  group_by(MTRS) %>% arrange(-SecDet) %>%
  summarize(SRF1 = n_distinct(gm_chemical_vvl[WQriskbin == "high"]),
            SRF2 = n_distinct(gm_chemical_vvl[WQriskbin == "medium"]),
            SRF3 = mean(avg_overMCL[avg_overMCL > 0]),
            SL1 = paste(gm_chemical_vvl[WQriskbin == "high"], collapse = "; "),
            SL2 = paste(gm_chemical_vvl[WQriskbin == "medium"], collapse = "; "),
            source_high = +("source" %in% method[WQriskbin == "high"]),
            neighbor_high = +("neighbor" %in% method[WQriskbin == "high"]),
            source_med = +("source" %in% method[WQriskbin == "medium"]),
            neighbor_med = +("neighbor" %in% method[WQriskbin == "medium"]),
            source_low = +("source" %in% method[WQriskbin == "low"]),
            neighbor_low = +("neighbor" %in% method[WQriskbin == "low"])) %>%
  mutate(method = ifelse(source_high == 1, "source",
                         ifelse(neighbor_high == 1, "neighbor",
                                ifelse(source_med == 1, "source",
                                       ifelse(neighbor_med == 1, "neighbor",
                                              ifelse(source_low == 1, "source", 
                                                     ifelse(neighbor_low == 1, "neighbor", "help"))))))) %>%
  select(MTRS, SRF1, SRF2, SRF3, SL1, SL2, method)
sdrisk[is.na(sdrisk)] <- 0
sdrisk <- sdrisk %>% mutate(WQriskbin = ifelse(SRF1 > 0, "high",
                                                ifelse(SRF2 > 0 & SRF1 == 0, "medium", "low")))

fwrite(section_data, paste0("Tables/sectiondata_detailed.csv"))
fwrite(sdrisk, paste0("Tables/sectiondata_risk.csv"))

#7.4 get location data for state smalls and domestic wells####
#dwr domestic well records data download process:
#connect to https://services.arcgis.com/aa38u6OgfNoCkTJ6/arcgis/rest/services/i07_WellCompletionReports_Exported_v2_gdb/FeatureServer
#save subset of well completion records using following parameters
#B118WellUse == "Domestic"
#Date Work Ended > 12/31/1969
#Record Type == “WellCompletion/New/Production or Monitoring/NA”
#export as .csv and save in "Data" folder
dwr_domesticwellrecords <- st_read("Data/i07_domestic_1970_completin_10_16_23.shp") %>% st_transform(crs = 4326)
MTRS <- st_read("Data/MTRS.shp") %>% st_transform(crs = 4326)
dwr_domesticwellrecords_mtrs <- st_join(dwr_domesticwellrecords, MTRS["MTRS"], left = F) %>%
  st_drop_geometry()
colnames(dwr_domesticwellrecords_mtrs)[colnames(dwr_domesticwellrecords_mtrs)=="CountyName"] <- "NAME"
mtrs_county <- MTRS %>% st_drop_geometry() %>%
  #remove split objects
  select(MTRS, NAME, NAMELSA) %>% distinct()
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
fwrite(domwc, "Tables/domesticwells_mtrs_county_10_17_2023.csv")

#get state small water system count per section
ss_mtrs <- read.csv("Data/SSWS_10_17_23_INTmtrs.csv", stringsAsFactors = F) 
mismatch_county_ssws <- ss_mtrs %>% filter(!tolower(NAME) == tolower(COUNTY), !is.na(Lat)) %>% select(PWSID, SYSTEM_NAME, COUNTY, REGULATING_AGENCY, Lat, Long, NAME, MTRS)
write.csv(mismatch_county_ssws, "ssws_county_mismatch.csv", row.names = F)

#put it all together (wq risk, domwellcounts, sswscounts):
section_data <- fread("Tables/sectiondata_detailed.csv")
sdrisk <- fread("Tables/sectiondata_risk.csv")

domestics <- fread("Reference_files/domesticwells_mtrs_county_10_17_2023.csv")
ssws_all <- fread("Data/SSWS_10_17_23_INTmtrs.csv")
#note - some SSWS do not have a lat/long, and some have a lat/long that is outside the listed regulating county
#for this analysis I am leaving those SSWS off of the location count but will include them in statewide summary
ssws <- ssws_all %>% filter(acc_loc == "yes") %>% 
  group_by(MTRS, NAME) %>%
  summarize(ssws_sum = n())

colnames(sdrisk)[colnames(sdrisk)=="WQriskbin"] <- "WQr_2024"
statesmall_points <- ssws_all %>% left_join(., sdrisk)
statesmall_points$WQr_2024[is.na(statesmall_points$WQr_2024)] <- "unknown"
statesmall_points <- statesmall_points %>% select(-c("CURRENT_STATUS", "OID_", "Join_Count", "TARGET_FID"))
#remove inaccurate SSWS locations after export
write.csv(statesmall_points, "Tables/SSWS_2024ARM.csv", row.names = F, na = "")


v2021 <- fread("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2020/sectiondata_risk_20202021-06-14.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(v2021) <- c("MTRS", "WQr_2021")
v2022 <- fread("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_project/Webtool/tables/Square Mile Sections/sectiondata_risk2021-11-10.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(v2022) <- c("MTRS", "WQr_2022")
v2023_whole <- fread("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2023/Tables/Square Mile Sections/sectiondata_risk2022-09-29.csv") %>% distinct()
colnames(v2023_whole) <- paste(colnames(v2023_whole),"2023",sep="_")
v2023 <- fread("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2023/Tables/Square Mile Sections/sectiondata_risk2022-09-29.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(v2023) <- c("MTRS", "WQr_2023")

#7.5 write shapefiles####
colnames(sdrisk)[colnames(sdrisk)=="WQriskbin"] <- "WQr_2024"
MTRS <- st_read("Reference_files/MTRS_joined_clean.shp") %>% st_transform(crs = 4326) %>% st_make_valid() %>% distinct()
sdrisk_mtrs <- left_join(MTRS, sdrisk) %>% 
  left_join(., domestics) %>% 
  left_join(., ssws) %>% 
  left_join(., v2021) %>%
  left_join(., v2022) %>%
  left_join(., v2023) %>%
  filter(!is.na(MTRS), !MTRS %in% c("", "BAY/DELTA", "SALTONSEA", "LAKETAHOE"))
sdrisk_mtrs$WQr_2024[is.na(sdrisk_mtrs$WQr_2024)] <- "unknown"
sdrisk_mtrs$WQr_2023[is.na(sdrisk_mtrs$WQr_2023)] <- "unknown"
sdrisk_mtrs$WQr_2022[is.na(sdrisk_mtrs$WQr_2022)] <- "unknown"
sdrisk_mtrs$WQr_2021[is.na(sdrisk_mtrs$WQr_2021)] <- "unknown"
sdrisk_mtrs$WQr_2023[sdrisk_mtrs$WQr_2023 == "med"] <- "medium"
sdrisk_mtrs$WQr_2022[sdrisk_mtrs$WQr_2022 == "med"] <- "medium"
sdrisk_mtrs$WQr_2021[sdrisk_mtrs$WQr_2021 == "med"] <- "medium"
sdrisk_mtrs[is.na(sdrisk_mtrs)] <- 0
st_write(sdrisk_mtrs, paste0("Shapefiles/arm_sectionrisk.shp"))
write.csv(st_drop_geometry(sdrisk_mtrs), "Tables/arm_sectionrisk.csv", row.names = F)

ssws_nolocations <- ssws_all %>% filter(acc_loc == "no") %>%
  group_by(COUNTY) %>%
  summarize(ssws_sum = n()) %>%
  mutate(NAME = str_to_title(COUNTY),
         NAMELSAD = paste0(str_to_title(COUNTY), " County"),
         COUNTYFP = NA,
         MTRS = NA,
         GSA_ID = NA,
         GSA_Name = NA,
         GSP_ID = NA,
         Basin_Numb = NA,
         Basin_Subb = NA,
         Basin_Name = NA,
         Basin_Su_1 = NA,
         SRF1 = NA,
         SRF2 = NA,
         SRF3 = NA,
         SL1 = NA,
         SL2 = NA,
         method = NA,
         WQr_2024 = "unknown",
         DWR_dom_count = NA,
         WQr_2021 = "unknown",
         WQr_2022 = "unknown",
         WQr_2023 = "unknown") %>%
  select(-COUNTY)
sdrisk_mtrs_table <- sdrisk_mtrs %>% st_drop_geometry() %>%
  rbind(., ssws_nolocations)
sdrisk_mtrs_table %>% group_by(WQr_2024) %>%
  summarize(sections = n_distinct(MTRS),
            domesticwells = sum(DWR_dom_count, na.rm = T),
            ssws = sum(ssws_sum, na.rm = T))
county_summary <- sdrisk_mtrs_table %>%group_by(WQr_2024, NAME) %>%
  summarize(sections= n_distinct(MTRS),
            domesticwells = sum(DWR_dom_count, na.rm = T),
            ssws = sum(ssws_sum, na.rm = T))
write.csv(county_summary, "county_summary.csv", row.names = F)

domestics_mtrs <- domestics %>% group_by(MTRS) %>%
  summarize(DWR_dom_count = sum(DWR_dom_count, na.rm = T))
topcontam_dw <- section_data %>% filter(WQriskbin == "high") %>%
  left_join(.,domestics_mtrs) %>%
  group_by(gm_chemical_vvl) %>%
  summarize(hr_dwc = sum(DWR_dom_count, na.rm = T)) %>% arrange(-hr_dwc)
fwrite(topcontam_dw, "topcontam2024_dw.csv")

ssws_mtrs <- ssws %>% group_by(MTRS) %>%
  summarize(ssws_sum = sum(ssws_sum, na.rm = T))
topcontam_ssws <- section_data %>% filter(WQriskbin == "high") %>%
  left_join(., ssws_mtrs) %>%
  group_by(gm_chemical_vvl) %>%
  summarize(hr_ssws = sum(ssws_sum, na.rm = T)) %>% arrange(-hr_ssws)
fwrite(topcontam_ssws, "topcontam2024_ssws.csv")


#7.6 write individual contaminant shapefiles####
MTRS <- st_read("Reference_files/MTRS_joined_clean.shp") %>% st_transform(crs = 4326) %>% st_make_valid() %>% distinct()
indivchems <- c("NO3N", "AS", "CR6", "U", "TCPR123")
for (x in 1:length(indivchems)) {
  chem_sd <- section_data %>% filter(gm_chemical_vvl == indivchems[x])
  chem_sd_mtrs <- left_join(MTRS, chem_sd) %>% filter(!is.na(MTRS), !MTRS %in% c("", "BAY/DELTA", "SALTONSEA", "LAKETAHOE"))
  colnames(chem_sd_mtrs)[colnames(chem_sd_mtrs)=="WQriskbin"] <- "WQr_2024"
  chem_sd_mtrs$WQr_2024[is.na(chem_sd_mtrs$WQr_2024)] <- "unknown"
  chem_sd_mtrs$WQr_2024[chem_sd_mtrs$WQr_2024 == "med"] <- "medium"
  chem_sd_mtrs[is.na(chem_sd_mtrs)] <- 0
  st_write(chem_sd_mtrs, paste0("Shapefiles/", indivchems[x], "_sectionrisk", Sys.Date(), ".shp"))
}

#8.0 calculate demographic statistics####
sdrisk_mtrs <- read.csv("Tables/arm_sectionrisk.csv")

#summarize results by block group
cbg2021_mtrs <- read.csv("Exported Maps/MTRS_cbg2021_areaintersect.csv") %>% 
  select(MTRS, NAME, GEOID, Area_mi, Int_Area_mi) %>%
  mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                        paste0("0", GEOID)),
         per_area = Int_Area_mi/Area_mi)
ssws_geoid <- read.csv("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ssws2024_bg.csv") %>%
  filter(!is.na(Lat)) %>%
  mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                        paste0("0", GEOID)))
sdrisk_census_ssws <- ssws_geoid %>% group_by(GEOID) %>%
  summarize(hr_ssws = sum(WQr_2024 == "high"),
            tot_ssws = n())
sdrisk_census <- full_join(cbg2021_mtrs, sdrisk_mtrs) %>%
  mutate(per_dw = DWR_dom_count*per_area)
sdrisk_census_dw <- sdrisk_census %>% group_by(GEOID) %>%
  summarize(hr_dw = sum(per_dw[WQr_2024 == "high"]),
            tot_dw = sum(per_dw))

cbg_summary <- full_join(sdrisk_census_ssws, sdrisk_census_dw)
cbg_summary[c("hr_ssws", "tot_ssws")][is.na(cbg_summary[c("hr_ssws", "tot_ssws")])] <- 0
ct_summary <- cbg_summary %>% mutate(GEOID_tract = substr(GEOID, 1, 11)) %>%
  group_by(GEOID_tract) %>% summarize(HR_ssws = sum(hr_ssws),
                                      TOT_ssws = sum(tot_ssws),
                                      HR_dw = sum(hr_dw),
                                      TOT_dw = sum(tot_dw))

#get race/mhi data by block group and tract
bgrace <- read.csv("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ACSDT5Y2022.B03002-Data.csv")
colstokeep <- c("GEO_ID", "B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E")
bgrace2 <- bgrace[,(names(bgrace) %in% colstokeep)]
bgrace2 <- bgrace2 %>% mutate_at(c("B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E"), as.numeric) %>% 
                       mutate(B03002_008E_009E = B03002_008E + B03002_009E,
                              GEOID = substr(GEO_ID, 10, 21)) %>% 
                       filter(!GEO_ID == "Geography") %>% 
                       select(GEOID, B03002_001E, B03002_003E, B03002_012E, B03002_004E, B03002_005E, B03002_006E, B03002_007E, B03002_008E_009E)

bgmhi <- read.csv("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ACSDT5Y2022.B19013-Data.csv")
CAmhi_2022 <- 91905
bgmhi2 <- bgmhi %>% mutate_at(c("B19013_001E"), as.numeric) %>%
                    mutate(GEOID = substr(GEO_ID, 10, 21),
                           DAC_status = ifelse(is.na(B19013_001E), "no data",
                                        ifelse(B19013_001E <= 0.6*CAmhi_2022, "SDAC",
                                        ifelse(B19013_001E <= 0.8*CAmhi_2022, "DAC", "none")))) %>%
                    select(GEOID, B19013_001E, DAC_status) %>%
                    filter(!GEOID == "")

bg <- full_join(bgmhi2, bgrace2)
write.csv(bg, "C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ACS2022_blockgroup.csv", row.names = F, na = "")

ctrace <- read.csv("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ACSDT5Y2022.B03002-Data (2).csv")
colstokeep <- c("GEO_ID", "B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E")
ctrace2 <- ctrace[,(names(ctrace) %in% colstokeep)]
ctrace2 <- ctrace2 %>% mutate_at(c("B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E"), as.numeric) %>% 
  mutate(B03002_008E_009E = B03002_008E + B03002_009E,
         GEOID_tract = substr(GEO_ID, 10, 20)) %>% 
  filter(!GEO_ID == "Geography") %>% 
  select(GEOID_tract, B03002_001E, B03002_003E, B03002_012E, B03002_004E, B03002_005E, B03002_006E, B03002_007E, B03002_008E_009E)

ctmhi <- read.csv("C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ACSDT5Y2022.B19013-Data (2).csv")
ctmhi2 <- ctmhi %>% mutate_at(c("B19013_001E"), as.numeric) %>%
  mutate(GEOID_tract = substr(GEO_ID, 10, 20),
         DAC_status = ifelse(is.na(B19013_001E), "no data",
                             ifelse(B19013_001E <= 0.6*CAmhi_2022, "SDAC",
                                    ifelse(B19013_001E <= 0.8*CAmhi_2022, "DAC", "none")))) %>%
  select(GEOID_tract, B19013_001E, DAC_status) %>%
  filter(!GEOID_tract == "")

ct <- full_join(ctmhi2, ctrace2)
write.csv(ct, "C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_2024/Census Data/ACS2022_censustract.csv", row.names = F, na = "")

#join with dw/ssws data
tract_export <- left_join(ct, ct_summary)
tract_export[c("HR_ssws", "TOT_ssws", "HR_dw", "TOT_dw")][is.na(tract_export[c("HR_ssws", "TOT_ssws", "HR_dw", "TOT_dw")])] <- 0
write.csv(tract_export, "Tables/tract.csv")

bg_export <- left_join(bg, cbg_summary)
bg_export[c("hr_ssws", "tot_ssws", "hr_dw", "tot_dw")][is.na(bg_export[c("hr_ssws", "tot_ssws", "hr_dw", "tot_dw")])] <- 0
write.csv(bg_export, "Tables/blockgroup.csv")

