#Water quality data download for 2025 Aquifer Risk Map

#Methodology write-up:

#Updated 11/12/2024 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

setwd("C:/Users/EHoulihan/OneDrive - Water Boards/Aquifer Risk Map/ARM 2025")

#1. Load required libraries, tables, and define special functions####
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('janitor')) install.packages('janitor'); library('janitor')
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
#https://gis.water.ca.gov/arcgis/rest/services/Environment/i07_WellCompletionReports/FeatureServer
MTRS <- st_read("Reference_files/i07_WellReportStatsBySection.shp") %>% st_transform(crs = 4326) #%>% select(MTRS) %>% distinct()

#join mtrs to groundwater units, then drop geometry
MTRS_GU <- st_join(MTRS, GU) %>% st_drop_geometry()

#join mtrs to groundwater units and remove duplicate sections (that are duplicated across county line splits of PLSS sections)
#define columns that must be numeric
numeric_cols <- c("DomWellCou", "DomWellDep", "DomWellD_1", "DomWellD_2",
                  "PubWellCou", "PubWellDep", "PubWellD_1", "PubWellD_2")

#identify columns to be deleted
delete <- c("COUNTY_CD", "WCRFolderL")

#remove duplicated mtrs items
MTRS_GU <- MTRS_GU[, !(colnames(MTRS_GU) %in% delete), drop=FALSE] %>% distinct() %>%
  mutate_at(numeric_cols, as.numeric)

# remove depth zeros (sections without any wells with depth information) (some sections have domestic wells but no depth data associated
#with them, so I am filtering on depth, not count)
clean_d <- MTRS_GU %>% filter(!DomWellD_2 == 0) #for domestic wells
clean_p <- MTRS_GU %>% filter(!PubWellD_2 == 0) #for public wells

#summarize by basin and display depth data:
#calculate average maximum and minimum well depths, and standard deviation. If there 2 or less sections
#in the basin, the basin max/min domestic is the average max domestic+150 ft/average min domestic-150 ft.
depth_dom <- clean_d %>% dplyr::group_by(GU_ID) %>% 
  dplyr::summarize(dom_section_n = n(),
                   dom_well_n = sum(DomWellCou), 
                   avgmax_dom = mean(DomWellD_2), 
                   avgmaxsd_dom = sd(DomWellD_2),
                   avgmin_dom = mean(DomWellD_1), 
                   avgminsd_dom = sd(DomWellD_1)) %>% 
  dplyr::mutate(maxdomestic = ifelse(dom_section_n <= 2, avgmax_dom + 3*150, avgmax_dom + 3*avgmaxsd_dom),
                mindomestic = ifelse(dom_section_n <= 2, avgmin_dom - 3*150, avgmin_dom - 3*avgminsd_dom))

#calculate average maximum public well depths. Same as above - if 2 or less sections add 150 to average max.
depth_pub <- clean_p %>% dplyr::group_by(GU_ID) %>% 
  dplyr::summarize(pub_section_n = n(),
                   pub_well_n = sum(PubWellCou),
                   avgmax_pub = mean(PubWellD_2), 
                   avgmaxsd_pub = sd(PubWellD_2)) %>% 
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
GU_depths <- read.table("Reference_files/well_depth_gwunits2024-08-22.txt", stringsAsFactors = F, header = T)
#column style####
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
  GM_REPORTING_LIMIT = col_double(),
  GM_SAMP_COLLECTION_DATE = col_character(),
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
  SRC_REPORTING_LIMIT = col_double(),
  SRC_SAMP_COLLECTION_DATE = col_character(),
  SRC_SAMP_COLLECTION_TIME = col_character(),
  SRC_ANALYTICAL_METHOD = col_character(),
  SRC_LAB_NOTE = col_character(),
  SRC_LATITUDE = col_double(),
  SRC_LONGITUDE = col_double(),
  SRC_DATUM = col_character(),
  SRC_WELL_DEPTH_FT = col_double(),
  SRC_TOP_DEPTH_OF_SCREEN_FT = col_double(),
  SRC_BOTTOM_DEPTH_OF_SCREEN_FT = col_double(),
  SRC_WELL_CATEGORY = col_character()
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

#identify wells that are within domestic depth filter####
#download location data
loc <- fread(archive::archive_extract("Reference_files/gama_location_construction_v2.7z"), colClasses = c("GM_WELL_ID" = "character"), quote = F) %>% 
  clean_names() %>% mutate(gm_well_id = str_remove_all(gm_well_id, '"'))
#loc <- fread("https://documents.geotracker.waterboards.ca.gov/gama/data_download/gama_location_construction_v2.zip", colClasses = c("GM_WELL_ID" = "character"), quote =  F) %>% 
#  clean_names() %>% distinct() #%>%

#well "AGW080016134-14805\xa0S.\xa0M" and "L10007696901-VLT\xdbFLDNON" have unusual characters. Address this here in locations and later in wq download
#loc[loc$gm_well_id == "AGW080016134-14805\xa0S.\xa0M", 4] <- "AGW080016134-14805"
#loc[loc$gm_well_id == "L10007696901-VLT\xdbFLDNON", 4] <- "L10007696901-VLT"

#create PWSID for DDW wells
ddw_pwsid <- loc %>% filter(gm_dataset_name == "DDW") %>%
  select(gm_dataset_name, gm_well_id) %>% 
  mutate(PWSID = ifelse(gm_dataset_name == "DDW", substr(gm_well_id, 1, 9), NA))
loc <- left_join(loc, ddw_pwsid)

#identify state small water system wells (not available as a well category on GAMA, but I use the list of SSWS I get from DDW)
ssws_list <- readxl::read_excel("C:/Users/EHoulihan/OneDrive - Water Boards/Aquifer Risk Map/ARM 2024/Data/SSWS_10_17_23.xlsx")

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

ind <- duplicated(wells_int[,c(1,9)])
view(wells_int[ind,])
#for wells with duplicate locations, keep first instance
wells_int <- wells_int[!ind,]

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
allwells <- allwells %>%
  mutate(uniquewellid = paste0(gm_dataset_name, gm_well_id))
write.table(allwells, paste0("Reference_files/allwells", Sys.Date(), ".txt"), sep = "\t", row.names = F, col.names = T)

#retrieve water quality data####
allwells <- read.table("Reference_files/allwells2024-10-08.txt", stringsAsFactors = F, header = T)
domwells <- allwells %>% filter(domdepth == "yes")
#define function to download data
#download various datasets for domestic depth wells
read7zip <- function(filename) {
  mytable <- fread(archive::archive_extract(paste0("Reference_files/", filename, ".7z")), colClasses = c("GM_WELL_ID" = "character"), quote = F) %>% 
    clean_names() %>% mutate(gm_well_id = str_remove_all(gm_well_id, '"'))
  mytable$gm_samp_collection_date <- lubridate::mdy(mytable$gm_samp_collection_date)
  mytable <- mytable %>% mutate(uniquewellid = paste0(gm_dataset_name, gm_well_id)) %>%
    filter(gm_samp_collection_date >= ymd("2000-01-01"),
           uniquewellid %in% domwells$uniquewellid,
           gm_chemical_vvl %in% chem_list$chemical_vvl)
  return(mytable)
}

#updated yearly (or less)
dpr <- read7zip("gama_dpr_statewide_v2")
dwr <- read7zip("gama_dwr_statewide_v2")
gama_dom <- read7zip("gama_gama_dom_statewide_v2")
gama_spstudy <- read7zip("gama_gama_sp-study_statewide_v2")
wrd <- read7zip("gama_wrd_statewide_v2")
ucdno3 <- read7zip("gama_ucd_no3_statewide_v2") %>% 
  #only keep ucd source data that is unique (not replicated in other datasets)
  filter(gm_data_source %in% c("UCD_NO3 - KECO", "UCD_NO3 - MOCO", "UCD_NO3 - MCEM", "UCD_NO3 - FRCO", "UCD_NO3 - TCEHS", "UCD_NO3 - KICO", "UCD_NO3 - TUCO", "UCD_NO3 - KRB98"))
wq1 <- rbind(dpr, dwr, gama_dom, gama_spstudy, wrd, ucdno3)
saveRDS(wq1, "wq1.RDS")

#updated monthly
ddw <- read7zip("gama_ddw_statewide_v2") %>% 
  #remove spring sites
  filter(!str_detect(gm_altwell_id1, "SPRING") | !str_detect(gm_altwell_id2, "SPRING"))
gama_localgw <- read7zip("gama_gama_localgw_statewide_v2")
gama_usgs <- read7zip("gama_gama_usgs_statewide_v2")
gama_usgs$gm_altwell_id1 <- as.character(gama_usgs$gm_altwell_id1)
localgw <- read7zip("gama_localgw_statewide_v2")
usgs_nwis <- read7zip("gama_usgs_nwis_statewide_v2")
#remove results that are a duplication in gama_usgs dataset
usgsnwis_small <- usgs_nwis %>% mutate(temp_id = gm_altwell_id1)
gamausgs_small <- gama_usgs %>%  mutate(temp_id = gm_well_id)
usgsnwis_unique <- anti_join(usgsnwis_small, gamausgs_small, by = c("temp_id", "gm_chemical_vvl", "gm_samp_collection_date")) %>% select(-temp_id)
wb_ilrp <- read7zip("gama_wb_ilrp_statewide_v2")
wq2 <- rbind(ddw, gama_localgw, gama_usgs, localgw, usgsnwis_unique, wb_ilrp)
saveRDS(wq2, "wq2.RDS")

wb_cleanup <- read7zip("gama_wb_cleanup_statewide_v2")
wq3 <- wb_cleanup
saveRDS(wq3, "wq3.RDS")

wq1 <- readRDS("wq1.RDS")
wq2 <- readRDS("wq2.RDS")
wq3 <- readRDS("wq3.RDS")

domwqdata <- rbind(wq1, wq2, wq3) 
#remove duplicates created by multiple well locations
domwqdata <- domwqdata %>% select(gm_dataset_name, gm_well_category, gm_well_id, gm_chemical_vvl, gm_result_modifier,
                                  gm_result, gm_reporting_limit, gm_samp_collection_date, src_result_modifier, uniquewellid) %>%
  distinct()
#remove well not in CA
domwqdata <- domwqdata %>% filter(!gm_well_id == "USGS-415953121292901") #not in california

#4. Standardize and clean data ####
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
    #group 1 (NULL modifier, Q/NQ_LT across all datasets)
    is.na(src_result_modifier) & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    is.na(src_result_modifier) & simple_res %in% c("nonzero") ~ "Q",
    #group 1.1 (modifiers that are always Q/NQ_LT across all datasets)
    src_result_modifier %in% c("", " ", "=", "E", "P", "X", "S", "A", "TI", "Estimated", "J") & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    src_result_modifier %in% c("", " ", "=", "E", "P", "X", "S", "A", "TI", "Estimated", "J") & simple_res %in% c("nonzero") ~ "Q",
    #group 2 (modifiers that are always NQ_LT across all datasets)
    src_result_modifier %in% c("<", "F", "Y", "ND", "M", "<=", "OS", "DN", "Not Detected", "Detected Not Quantified", "-") ~ "NQ_LT",
    #group 3 (modifiers that are always GNQ across all datasets)
    src_result_modifier %in% c(">", ">=", "Present Above Quantification Limit") ~ "GNQ",
    #group 4 (modifiers that have different meaning in different datasets)
    #group 4.1 "N" (Q/NQ_LT in newDDW, NQ_LT in DPR)
    gm_dataset_name == "DDW" & src_result_modifier == "N" & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_dataset_name == "DDW" & src_result_modifier == "N" & simple_res %in% c("nonzero") ~ "Q",
    gm_dataset_name == "DPR" & src_result_modifier == "N" ~ "NQ_LT",
    #group 4.2 "R" (Reject in DPR, NQ_LT in GAMA_USGS)
    gm_dataset_name %in% c("GAMA_USGS") & src_result_modifier %in% c("R") ~ "NQ_LT",
    #group 4.3 "U" (Q/NQ_LT in DPR, NQ_LT in GAMA_USGS)
    gm_dataset_name == "DPR" & src_result_modifier %in% c("U") & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_dataset_name == "DPR" & src_result_modifier %in% c("U") & simple_res %in% c("nonzero") ~ "Q",
    gm_dataset_name == "GAMA_USGS" & src_result_modifier %in% c("U") ~ "NQ_LT",
    #group 4.4 "V" (Q/NQ_LT in DPR, Reject in GAMA_USGS)
    gm_dataset_name == "DPR" & src_result_modifier %in% c("V") & simple_res %in% c("NL", "zero") ~ "NQ_LT",
    gm_dataset_name == "DPR" & src_result_modifier %in% c("V") & simple_res %in% c("nonzero") ~ "Q",
    TRUE ~ "REJECT")
  )
check <- domwqdata %>% group_by(src_result_modifier, gm_dataset_name, new_modifier) %>% count()
view(check)
unique(domwqdata$new_modifier)
domwqdata <- domwqdata %>%
  filter(!new_modifier == "REJECT")

#special chemicals####
#NITROGEN: convert Nitrate/Nitrite combined to Nitrate
#for the conversion, only keep highest values if multiple samples were taken on the same day (but leave duplicates in domwqdata)
no3n <- domwqdata %>% filter(gm_chemical_vvl %in% c("NO3N")) %>% select(-gm_chemical_vvl)
colnames(no3n) <- c("gm_dataset_name", "gm_well_category", "gm_well_id", "NO3N_result_modifier", "NO3N_result", "NO3N_reporting_limit", "gm_samp_collection_date", 
                    "NO3N_src_result_modifier", "uniquewellid", "NO3N_simple_res", "NO3N_new_modifier")
no2n <- domwqdata %>% filter(gm_chemical_vvl %in% c("NO2")) %>% select(-gm_chemical_vvl)
colnames(no2n) <- c("gm_dataset_name", "gm_well_category", "gm_well_id", "NO2_result_modifier", "NO2_result", "NO2_reporting_limit", "gm_samp_collection_date", 
                    "NO2_src_result_modifier", "uniquewellid", "NO2_simple_res", "NO2_new_modifier")
no3no2n <- domwqdata %>% filter(gm_chemical_vvl %in% c("NO3NO2N")) %>% select(-gm_chemical_vvl)
colnames(no3no2n) <- c("gm_dataset_name", "gm_well_category", "gm_well_id", "NO3NO2N_result_modifier", "NO3NO2N_result", "NO3NO2N_reporting_limit", "gm_samp_collection_date", 
                       "NO3NO2N_src_result_modifier", "uniquewellid", "NO3NO2N_simple_res", "NO3NO2N_new_modifier")
#join no3n, no2n and no3no2n files
nitrogen <- left_join(no3no2n, no3n) %>% 
  left_join(., no2n) %>% distinct() #%>%
  #filter table so only instances where no3no2n results are available without a corresponding no3n result on the same day are kept
nitrogen <- nitrogen %>% 
  filter(is.na(NO3N_new_modifier)) %>%
  #create new result, rl, and modifier columns
  mutate(gm_result = case_when(
    #if there is no no2n, use no3no2n value
    is.na(NO2_new_modifier) ~ NO3NO2N_result,
    #if there is no2n value, subtract no2n from no3no2n
    !is.na(NO2_new_modifier) & NO2_result <= NO3NO2N_result ~ NO3NO2N_result - NO2_result,
    !is.na(NO2_new_modifier) & NO2_result > NO3NO2N_result ~ 0,
    TRUE ~ -99),
  gm_reporting_limit = case_when(
    is.na(NO2_new_modifier) ~ NO3NO2N_reporting_limit,
    !is.na(NO2_new_modifier) & NO3NO2N_reporting_limit >= NO2_reporting_limit ~ NO3NO2N_reporting_limit,
    !is.na(NO2_new_modifier) & NO3NO2N_reporting_limit < NO2_reporting_limit ~ NO2_reporting_limit,
    TRUE ~ NA),
  #never the case where NO2 is quantified but NO3NO2N isn't (but NO3NO2N can be quantified and NO2 not) so always take NO3NO2N modifier
  gm_result_modifier = NO3NO2N_result_modifier,
  src_result_modifier = NO3NO2N_src_result_modifier, 
  new_modifier = NO3NO2N_new_modifier) %>%
  mutate(gm_chemical_vvl = "NO3N",
         simple_res = NA) %>%
  select(gm_dataset_name, gm_well_category, gm_well_id, gm_chemical_vvl, gm_result_modifier, gm_result, gm_reporting_limit, 
         gm_samp_collection_date, src_result_modifier, uniquewellid, simple_res, new_modifier) %>% ungroup()
domwqdata2 <- rbind(domwqdata, nitrogen) %>% filter(!gm_chemical_vvl == "NO3NO2N")
rm(nitrogen, no2n, no3n, no3no2n)

#RADIUM: combine radium-226 and radium-228 (if both sampled on the same well on same date) (mcl is based on combined total)
ra226 <- domwqdata2 %>% filter(gm_chemical_vvl %in% c("RA-226")) %>% select(-gm_chemical_vvl)
colnames(ra226) <- c("gm_dataset_name", "gm_well_category", "gm_well_id", "RA226_result_modifier", "RA226_result", "RA226_reporting_limit", "gm_samp_collection_date", 
                     "RA226_src_result_modifier", "uniquewellid", "RA226_simple_res", "RA226_new_modifier")
ra228 <- domwqdata2 %>% filter(gm_chemical_vvl %in% c("RA-228")) %>% select(-gm_chemical_vvl)
colnames(ra228) <- c("gm_dataset_name", "gm_well_category", "gm_well_id", "RA228_result_modifier", "RA228_result", "RA228_reporting_limit", "gm_samp_collection_date", 
                     "RA228_src_result_modifier", "uniquewellid", "RA228_simple_res", "RA228_new_modifier")
#join ra226 and ra228 files
ra226_228 <- full_join(ra226, ra228) %>% distinct() %>%
  #create new result, rl, and modifier columns
  mutate(gm_result = case_when(
    #if there is no ra228 value, use ra226 value
    is.na(RA228_new_modifier) ~ RA226_result,
    #if there is no ra226 value, use ra228 value
    is.na(RA226_new_modifier) ~ RA228_result,
    #if ra226 and ra228 values are not quantified, take the maximum non-quantified result
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" & RA226_result > RA228_result ~ RA226_result,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" & RA226_result <= RA228_result ~ RA228_result,
    #if ra226 is not quantified but ra228 is quantified, use ra228 result
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_result,
    #same but flipped
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_result,
    #if ra226 and ra228 values are quantified, take the sum as the result
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" ~ RA226_result + RA228_result,
    TRUE ~ -99
  ),
  gm_reporting_limit = case_when(
    is.na(RA228_new_modifier) ~ RA226_reporting_limit,
    is.na(RA226_new_modifier) ~ RA228_reporting_limit,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" & RA226_result > RA228_result ~ RA226_reporting_limit,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" & RA226_result <= RA228_result ~ RA228_reporting_limit,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" & is.na(RA226_reporting_limit) & !is.na(RA228_reporting_limit) ~ RA228_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" & !is.na(RA226_reporting_limit) & is.na(RA228_reporting_limit) ~ RA226_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" & is.na(RA226_reporting_limit) & is.na(RA228_reporting_limit) ~ NA,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" & RA226_reporting_limit > RA228_reporting_limit ~ RA226_reporting_limit,
    RA226_new_modifier == "Q" & RA228_new_modifier == "Q" & RA226_reporting_limit <= RA228_reporting_limit ~ RA228_reporting_limit,
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
  src_result_modifier = case_when(
    is.na(RA228_new_modifier) ~ RA226_src_result_modifier,
    is.na(RA226_new_modifier) ~ RA228_src_result_modifier,
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "NQ_LT" ~ "<",
    RA226_new_modifier == "NQ_LT" & RA228_new_modifier == "Q" ~ RA228_src_result_modifier,
    RA226_new_modifier == "Q" & RA228_new_modifier == "NQ_LT" ~ RA226_src_result_modifier,
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
#some days have more than one ra226 or ra228 values, which multiples the total number of values when they are joined. reduce to one sample per day.
#ra226_228_one <- ra226_228 %>%
#  group_by(gm_well_id, gm_dataset_name, gm_samp_collection_date) %>%
#  #define average quantified results or maximum non-quantified result
#  summarize(q_res = mean(gm_result[new_modifier == "Q"]),
#            q_rl = max(gm_reporting_limit[new_modifier == "Q"]),
#            nq_res = max(gm_result[new_modifier == "NQ_LT"]),
#            nq_rl = max(gm_reporting_limit[new_modifier == "NQ_LT"])) %>%
#  #if quantified average exists, use that. if not, use the non-quantifed maximum
#  mutate(gm_result = ifelse(!q_res == -Inf & !is.na(q_res), q_res, nq_res),
#         gm_reporting_limit = ifelse(!q_res == -Inf & !is.na(q_res), q_rl, nq_rl),
#         new_modifier = ifelse(!q_res == -Inf & !is.na(q_res), "Q", "NQ_LT")) %>%
#  mutate(gm_chemical_vvl = "RA-226-228",
#         simple_res = NA) %>%
#  select(gm_well_id, gm_dataset_name, gm_chemical_vvl, gm_result_modifier, gm_result, gm_samp_collection_date, gm_reporting_limit, simple_res, new_modifier)
ra226_228_2 <- ra226_228 %>%
  mutate(gm_chemical_vvl = "RA-226-228",
         simple_res = NA) %>%
  select(gm_dataset_name, gm_well_category, gm_well_id, gm_chemical_vvl, gm_result_modifier, gm_result, gm_reporting_limit, 
         gm_samp_collection_date, src_result_modifier, uniquewellid, simple_res, new_modifier)
domwqdata2 <- rbind(domwqdata2, ra226_228_2) %>% filter(!gm_chemical_vvl %in% c("RA-226", "RA-228"))
rm(ra226, ra228, ra226_228, ra226_228_2)

#edit chem_list to include combined radium vvl
chem_list[102, 1] <- "RA-226-228"
chem_list[102, 2] <- "Combined Radium 226 and Radium 228"
chem_list[102, 3] <- "pCi/L"
chem_list[102, 4] <- 5
chem_list[102, 5] <- "MCL"

chem_list2 <- chem_list %>% filter(!chemical_vvl %in% c("NO3NO2N", "RA-226", "RA-228"))

#standardize Reporting Limits/Results (fix Null results, results == 0, and missing reporting limits)####
domwqdata2$gm_reporting_limit <- as.numeric(domwqdata2$gm_reporting_limit)
finishedchems <- c()         #hold chemicals as they are cleaned
additions <- data.frame()    #holds rows with "newRL"s
chemicals <- unique(domwqdata2$gm_chemical_vvl)     #list of chemicals that need to be fixed
domwqdata2 <- dplyr::arrange(domwqdata2, gm_samp_collection_date)
chemicals <- sort(chemicals)
for (ch in chemicals) {
  uRL <- domwqdata2 %>% dplyr::filter(is.na(gm_reporting_limit) | gm_reporting_limit == 0, gm_result <= 0, gm_chemical_vvl == ch) #list of measurements with unknown RL
  kRL <- domwqdata2 %>% dplyr::filter(gm_chemical_vvl == ch, !is.na(gm_reporting_limit), !gm_reporting_limit == 0)             #list of measurements with known RL
  if(dim(kRL)[1] == 0) {   #if chemical has no RLs listed, take implied non-detects
    kRL <- domwqdata2 %>% dplyr::filter(gm_chemical_vvl == ch, gm_result_modifier == "<", !gm_result == 0) %>% 
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
domwqdata2 <- left_join(domwqdata2, additions)
rm(additions)
domwqdata2 <- domwqdata2 %>% distinct()

#compress RLs into a single column, if applicable
domwqdata2 <- domwqdata2 %>% mutate(RL = ifelse(!is.na(gm_reporting_limit) & gm_reporting_limit > 0, gm_reporting_limit, 
                                              ifelse(!is.na(estRL), estRL, gm_result)),
                                  RLmethod = ifelse(!is.na(gm_reporting_limit) & gm_reporting_limit > 0, "reported", 
                                                    ifelse(!is.na(estRL), "estimated", "result"))) %>%
  #join with chemical MCLs
  left_join(., select(chem_list2, chemical_vvl, comparison_concentration_value), by = c("gm_chemical_vvl" = "chemical_vvl"))      


#calculate MCL index and track detections/non-detections, outliers
domwqdata2 <- domwqdata2 %>% 
  #calculate numeric value for all results (with caveat for non-detects that are greater than CCV)
  mutate(RES = ifelse(new_modifier %in% c("NQ_LT"), RLfun(RL, comparison_concentration_value), gm_result),
         #calculate MCL index
         MCLindex = RES/comparison_concentration_value)                                            
md_max <- summarize(group_by(domwqdata2 %>% dplyr::filter(new_modifier == "Q"), gm_chemical_vvl), avg = mean(RES), sd = sd(RES)) %>% 
  #calculate outlier maximum
  mutate(max = (avg + (10*sd))) %>% select(gm_chemical_vvl, max) %>% 
  left_join(., select(chem_list2, chemical_vvl, comparison_concentration_value), by = c("gm_chemical_vvl" = "chemical_vvl")) %>%
  mutate(max = ifelse(max < comparison_concentration_value | is.na(max), comparison_concentration_value*2, max))

write.csv(md_max, paste0("Reference_files/outliers", Sys.Date(), ".csv"), row.names = F)
domwqdata2 <- left_join(domwqdata2, select(md_max, -comparison_concentration_value), by = "gm_chemical_vvl") %>% 
  mutate(outlier_status = ifelse(is.na(max), "no outliers", 
                                 ifelse(RES > max & new_modifier == "Q", "outlier", "standard")))

#clean table and output results
domwqdata2 <- domwqdata2 %>% select(gm_well_id, gm_dataset_name, gm_samp_collection_date, gm_chemical_vvl,
                                           RES, RL, RLmethod, new_modifier, outlier_status, MCLindex)
saveRDS(domwqdata2, "wqdata.RDS")

#5. Load data tables####
domwqdata_sm <- readRDS("wqdata.RDS")
wellinfo <- fread("Reference_files/allwells2024-10-08.txt") %>% 
  select(gm_well_id, domdepth, gm_latitude, gm_longitude, gm_well_category, gm_dataset_name, MTRS)

#6. Calculate point risk scores (wells)####
domwqdata_sm <- domwqdata_sm %>% 
  left_join(., select(wellinfo, gm_well_id, gm_dataset_name, MTRS)) %>%
  filter(!outlier_status == "outlier", !is.na(MCLindex)) %>% #remove outliers and samples without values (?)
  select(gm_well_id, gm_dataset_name, gm_samp_collection_date, gm_chemical_vvl, MCLindex) %>% #select only necessary columns to calculate averages
  mutate(year = year(gm_samp_collection_date)) %>% #prepare for grouping samples by year
  filter(year >= 2004) #define date range

#6.1 calculate well averages ####
domwellavg <- domwqdata_sm %>% 
  group_by(gm_chemical_vvl, gm_well_id, gm_dataset_name, year) %>% #group samples by chemical, well, and year
  summarize(yr_avg = mean(MCLindex))%>%
  group_by(gm_chemical_vvl, gm_well_id, gm_dataset_name) %>% #group samples by chemical, well
  summarize(well_avg = mean(yr_avg))

#6.2 calculate recent results####
domwellre_avg <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2019-01-01"), MCLindex > 1) %>% 
  group_by(gm_well_id, gm_dataset_name, gm_chemical_vvl) %>%
  summarize(re_avg_over = mean(MCLindex))
domwellrecent <- domwqdata_sm %>% filter(gm_samp_collection_date >= ymd("2019-01-01")) %>%
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

fwrite(pointdata, paste0("Tables/pdatadetailed2025.csv"))
fwrite(pdrisk, paste0("Tables/pdata_risk2025.csv"))

#7. Calculates square mile section risk####
#load point data tables
pointdata <- fread("Tables/pdatadetailed2025.csv")
#wellinfo <-  fread("Reference_files/allwells2024-10-08.txt") %>% 
#  select(gm_well_id, domdepth, gm_latitude, gm_longitude, gm_well_category, gm_dataset_name, MTRS)

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

fwrite(section_data, paste0("Tables/sectiondata_detailed2025.csv"))
fwrite(sdrisk, paste0("Tables/sectiondata_risk2025.csv"))

#7.4 get location data for state smalls and domestic wells####
#dwr domestic well records data download process:
#connect to https://services.arcgis.com/aa38u6OgfNoCkTJ6/arcgis/rest/services/i07_WellCompletionReports_Exported_v2_gdb/FeatureServer
#save subset of well completion records using following parameters
#B118WellUse == "Domestic"
#Date Work Ended > 12/31/1969
#Record Type == “WellCompletion/New/Production or Monitoring/NA”
#export as .csv and save in "Data" folder
dwr_domesticwellrecords <- st_read("Reference_files/i07_domestic_1970_completin_08_29_24.shp") %>% st_transform(crs = 4326) %>% st_make_valid()
MTRS <- st_read("Reference_files/MTRS_joined_clean.shp") %>% st_transform(crs = 4326) %>% st_make_valid()
dwr_domesticwellrecords_mtrs <- st_join(dwr_domesticwellrecords, MTRS["MTRS"], left = F) %>%
  st_drop_geometry()
colnames(dwr_domesticwellrecords_mtrs)[colnames(dwr_domesticwellrecords_mtrs)=="CountyName"] <- "NAME"
mtrs_county <- MTRS %>% st_drop_geometry() %>%
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
fwrite(domwc, "Reference_files/domesticwells_mtrs_county_9_5_2024.csv")

#get state small water system count per section
#ss_mtrs <- read.csv("Reference_files/SSWS_2025.csv", stringsAsFactors = F) 
#mismatch_county_ssws <- ss_mtrs %>% filter(!tolower(NAME) == tolower(COUNTY), !is.na(Lat)) %>% select(PWSID, SYSTEM_NAME, COUNTY, REGULATING_AGENCY, Lat, Long, NAME, MTRS)
#write.csv(mismatch_county_ssws, "ssws_county_mismatch.csv", row.names = F)

#put it all together (wq risk, domwellcounts, sswscounts):
section_data <- fread("Tables/sectiondata_detailed2025.csv")
sdrisk <- fread("Tables/sectiondata_risk2025.csv")

domestics <- fread("Reference_files/domesticwells_mtrs_county_9_5_2024.csv")
#ssws_all <- fread("Reference_files/SSWS_2025.csv")
#ssws_deactivated <- readxl::read_excel("Reference_files/DataExport_11_01_2024.xlsx") 
#ssws_deactivated <- ssws_deactivated %>% filter(`CURRENT SAFER STATUS` == "Deactivated")
#ssws_all <- ssws_all %>% filter(!PWSID %in% ssws_deactivated$PWSID)

#ssws_newlocs <- readxl::read_excel("20241105_SSWS-Clearinghouse.xlsx")
#ssws_newlocs <- ssws_newlocs %>% filter(!LATITUDE == "NULL") %>% select(PWSID, LATITUDE, LONGITUDE)
#ssws_all <- left_join(ssws_all, ssws_newlocs)
#ssws_all <- ssws_all %>% mutate(LAT = ifelse(!is.na(LATITUDE), LATITUDE, LAT),
#                                LONG = ifelse(!is.na(LONGITUDE), LONGITUDE, LONG))
ssws_final <- read.csv("ssws_mtrsjoin_11-8-24.csv")

#note - some SSWS do not have a lat/long, and some have a lat/long that is outside the listed regulating county
#for this analysis I am leaving those SSWS off of the location count but will include them in statewide summary
ssws <- ssws_final %>% filter(!MTRS == "") %>% 
  group_by(MTRS, NAME) %>%
  summarize(ssws_sum = n())

colnames(sdrisk)[colnames(sdrisk)=="WQriskbin"] <- "WQr_2025"
#statesmall_points <- ssws_all %>% left_join(., sdrisk)
#statesmall_points$WQr_2024[is.na(statesmall_points$WQr_2024)] <- "unknown"
#statesmall_points <- statesmall_points %>% select(-c("CURRENT_STATUS", "OID_", "Join_Count", "TARGET_FID"))
#statesmall_points$WQr_2025[is.na(statesmall_points$WQr_2025)] <- "unknown"
#write.csv(statesmall_points, "Tables/SSWS_2025ARM.csv", row.names = F, na = "")

#accidentally was using an incorrect version of the 2021 map for previous year comparisons, now corrected
v2021 <- fread("C:/Users/EHoulihan/OneDrive - Water Boards/Aquifer Risk Map/ARM 2021/webtool/tables/sectiondata_arm.csv") %>%
  select(MTRS, s_risk) %>% distinct()
colnames(v2021) <- c("MTRS", "WQr_2021")
v2022 <- fread("C:/Users/EHoulihan/OneDrive - Water Boards/Aquifer Risk Map/ARM 2022/Webtool/tables/Square Mile Sections/sectiondata_risk2021-11-10.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(v2022) <- c("MTRS", "WQr_2022")
v2023 <- fread("C:/Users/EHoulihan/OneDrive - Water Boards/Aquifer Risk Map/ARM 2023/ARM_2023/Tables/Square Mile Sections/sectiondata_risk2022-09-29.csv") %>%
  select(MTRS, WQriskbins) %>% distinct()
colnames(v2023) <- c("MTRS", "WQr_2023")
v2024 <- fread("C:/Users/EHoulihan/OneDrive - Water Boards/Aquifer Risk Map/ARM 2024/Tables/sectiondata_risk.csv") %>%
  select(MTRS, WQriskbin) %>% distinct()
colnames(v2024) <- c("MTRS", "WQr_2024")

#write shapefiles
#colnames(sdrisk)[colnames(sdrisk)=="WQriskbin"] <- "WQr_2025"
MTRS <- st_read("Reference_files/MTRS_joined_clean.shp") %>% st_transform(crs = 4326) %>% st_make_valid() %>% distinct()
sdrisk_mtrs <- left_join(MTRS, sdrisk) %>% 
  left_join(., domestics) %>% 
  left_join(., ssws) %>% 
  left_join(., v2021) %>%
  left_join(., v2022) %>%
  left_join(., v2023) %>%
  left_join(., v2024) %>%
  filter(!is.na(MTRS), !MTRS %in% c("", "BAY/DELTA", "SALTONSEA", "LAKETAHOE"))
sdrisk_mtrs$WQr_2025[is.na(sdrisk_mtrs$WQr_2025)] <- "unknown"
sdrisk_mtrs$WQr_2024[is.na(sdrisk_mtrs$WQr_2024)] <- "unknown"
sdrisk_mtrs$WQr_2023[is.na(sdrisk_mtrs$WQr_2023)] <- "unknown"
sdrisk_mtrs$WQr_2022[is.na(sdrisk_mtrs$WQr_2022)] <- "unknown"
sdrisk_mtrs$WQr_2021[sdrisk_mtrs$WQr_2021 == ""] <- "unknown"
sdrisk_mtrs$WQr_2023[sdrisk_mtrs$WQr_2023 == "med"] <- "medium"
sdrisk_mtrs$WQr_2022[sdrisk_mtrs$WQr_2022 == "med"] <- "medium"
sdrisk_mtrs$WQr_2021[sdrisk_mtrs$WQr_2021 == "med"] <- "medium"
sdrisk_mtrs[is.na(sdrisk_mtrs)] <- 0
st_write(sdrisk_mtrs, paste0("Tables/arm_sectionrisk.shp"), append=F)
write.csv(st_drop_geometry(sdrisk_mtrs), "Tables/arm_sectionrisk.csv", row.names = F)

#7.5 write individual contaminant shapefiles####
MTRS <- st_read("Reference_files/MTRS_joined_clean.shp") %>% st_transform(crs = 4326) %>% st_make_valid() %>% distinct()
indivchems <- c("NO3N", "AS", "CR6", "U", "TCPR123")
for (x in 1:length(indivchems)) {
  chem_sd <- section_data %>% filter(gm_chemical_vvl == indivchems[x])
  chem_sd_mtrs <- left_join(MTRS, chem_sd) %>% filter(!is.na(MTRS), !MTRS %in% c("", "BAY/DELTA", "SALTONSEA", "LAKETAHOE"))
  colnames(chem_sd_mtrs)[colnames(chem_sd_mtrs)=="WQriskbin"] <- "WQr_2025"
  chem_sd_mtrs$WQr_2025[is.na(chem_sd_mtrs$WQr_2025)] <- "unknown"
  chem_sd_mtrs$WQr_2025[chem_sd_mtrs$WQr_2025 == "med"] <- "medium"
  chem_sd_mtrs[is.na(chem_sd_mtrs)] <- 0
  st_write(chem_sd_mtrs, paste0("Tables/", indivchems[x], "_sectionrisk", Sys.Date(), ".shp"))
}

#8.0 calculate ARM statistics####
sdrisk_mtrs <- read.csv("Tables/arm_sectionrisk.csv")
sdrisk_ssws <- read.csv("Tables/SSWS_2025ARM.csv")

#summarize results by block group
mtrs_cbg <- read.csv("Reference_files/census/MTRS_cbg2021_areaintersect.csv") %>%
  select(MTRS, NAME, GEOID, Area_mi, Int_Area_mi) %>%
  mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                        paste0("0", GEOID)),
         per_area = Int_Area_mi/Area_mi)
sdrisk_mtrs_cbg <- full_join(sdrisk_mtrs, mtrs_cbg) %>%
  mutate(per_dw = DWR_dom_count*per_area) %>% group_by(GEOID) %>%
  summarize(hr_dw = sum(per_dw[WQr_2025 == "high"]),
            tot_dw = sum(per_dw))
  
ssws_cbg <- read.csv("Reference_files/census/ssws_bg.csv") %>%
  mutate(GEOID = ifelse(grepl("^0...........", geoid), geoid,
                        paste0("0", geoid)))
sdrisk_ssws_cbg <- full_join(sdrisk_ssws, ssws_cbg) %>% group_by(GEOID) %>%
  summarize(hr_ssws = sum(WQ_2025 == "high"),
            tot_ssws = n())
cbg_summary <- full_join(sdrisk_mtrs_cbg, sdrisk_ssws_cbg) %>% filter(!is.na(GEOID))
cbg_summary[c("hr_ssws", "tot_ssws")][is.na(cbg_summary[c("hr_ssws", "tot_ssws")])] <- 0
ct_summary <- cbg_summary %>% mutate(GEOID_tract = substr(GEOID, 1, 11)) %>%
  group_by(GEOID_tract) %>% summarize(HR_ssws = sum(hr_ssws),
                                      TOT_ssws = sum(tot_ssws),
                                      HR_dw = sum(hr_dw),
                                      TOT_dw = sum(tot_dw))

#get race/mhi data by block group and tract
bgrace <- read.csv("Reference_files/census/ACSDT5Y2023.B03002-Data.csv")
colstokeep <- c("GEO_ID", "B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E")
bgrace2 <- bgrace[,(names(bgrace) %in% colstokeep)]
bgrace2 <- bgrace2 %>% mutate_at(c("B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E"), as.numeric) %>% 
                       mutate(B03002_008E_009E = B03002_008E + B03002_009E,
                              GEOID = substr(GEO_ID, 10, 21)) %>% 
                       filter(!GEO_ID == "Geography") %>% 
                       select(GEOID, B03002_001E, B03002_003E, B03002_012E, B03002_004E, B03002_005E, B03002_006E, B03002_007E, B03002_008E_009E)

bgmhi <- read.csv("Reference_files/census/ACSDT5Y2023.B19013-Data.csv")
#CAmhi_2022 <- 91905
#from ACS 5 year 2023 CA estimate
CAmhi_2023 <- 96334
bgmhi2 <- bgmhi %>% mutate_at(c("B19013_001E"), as.numeric) %>%
                    mutate(GEOID = substr(GEO_ID, 10, 21),
                           DAC_status = ifelse(is.na(B19013_001E), "no data",
                                        ifelse(B19013_001E <= 0.6*CAmhi_2023, "SDAC",
                                        ifelse(B19013_001E <= 0.8*CAmhi_2023, "DAC", "none")))) %>%
                    select(GEOID, B19013_001E, DAC_status) %>%
                    filter(!GEOID == "")

bg <- full_join(bgmhi2, bgrace2)
write.csv(bg, "Tables/ACS2023_blockgroup.csv", row.names = F, na = "")

ctrace <- read.csv("Reference_files/census/ACSDT5Y2023.B03002-Data (2).csv")
colstokeep <- c("GEO_ID", "B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E")
ctrace2 <- ctrace[,(names(ctrace) %in% colstokeep)]
ctrace2 <- ctrace2 %>% mutate_at(c("B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_006E", "B03002_007E", "B03002_008E", "B03002_009E", "B03002_012E"), as.numeric) %>% 
  mutate(B03002_008E_009E = B03002_008E + B03002_009E,
         GEOID_tract = substr(GEO_ID, 10, 20)) %>% 
  filter(!GEO_ID == "Geography") %>% 
  select(GEOID_tract, B03002_001E, B03002_003E, B03002_012E, B03002_004E, B03002_005E, B03002_006E, B03002_007E, B03002_008E_009E)

ctmhi <- read.csv("Reference_files/census/ACSDT5Y2023.B19013-Data (2).csv")
ctmhi2 <- ctmhi %>% mutate_at(c("B19013_001E"), as.numeric) %>%
  mutate(GEOID_tract = substr(GEO_ID, 10, 20),
         DAC_status = ifelse(is.na(B19013_001E), "no data",
                             ifelse(B19013_001E <= 0.6*CAmhi_2023, "SDAC",
                                    ifelse(B19013_001E <= 0.8*CAmhi_2023, "DAC", "none")))) %>%
  select(GEOID_tract, B19013_001E, DAC_status) %>%
  filter(!GEOID_tract == "")

ct <- full_join(ctmhi2, ctrace2)
write.csv(ct, "Tables/ACS2023_censustract.csv", row.names = F, na = "")

#join with dw/ssws data
tract_export <- left_join(ct, ct_summary)
tract_export[c("HR_ssws", "TOT_ssws", "HR_dw", "TOT_dw")][is.na(tract_export[c("HR_ssws", "TOT_ssws", "HR_dw", "TOT_dw")])] <- 0
write.csv(tract_export, "Tables/tract.csv")

bg_export <- left_join(bg, cbg_summary)
bg_export[c("hr_ssws", "tot_ssws", "hr_dw", "tot_dw")][is.na(bg_export[c("hr_ssws", "tot_ssws", "hr_dw", "tot_dw")])] <- 0
write.csv(bg_export, "Tables/blockgroup.csv")

