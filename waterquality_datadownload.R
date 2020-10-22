#Water quality data download for Needs Assessment and Aquifer Risk Map

#Methodology write-up: https://gispublic.waterboards.ca.gov/portal/home/item.html?id=70feb9f4b00f4b3384a9a0bf89f9f18a

#Updated 10/9/2020 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

#Components:
#1. download water quality data from GAMA GIS from
  #Division of Drinking Water
  #Department of Water Resources
  #USGS
  #GAMA projects
  #Geotracker - ILRP and Local Projects ONLY
#2. create depth filter to identify state small and domestic well depth water quality
#3. apply depth filter to wells with water quality data
#4. standardize and clean water quality data (standardize chemicals, handle outliers and non-detects)

#set working directory
if (!dir.exists("~/R_projects")) dir.create("~/R_projects", FALSE); setwd("~/R_projects")

#load required libraries
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('rgdal')) install.packages('rgdal'); library('rgdal')
if (!require('sf')) install.packages('sf'); library('sf')

#define special functions
fix_pun <- function(x) { #fix punctuation issues
  x_names <- make.names(names(x))
  x_names <- gsub(x = x_names, pattern = '\\.\\.', replacement = '_') 
  x_names <- gsub(x = x_names, pattern = '\\.', replacement = '_')
  x_names <- gsub(x = x_names, pattern = '\\_$', replacement = '')
  x_names
}
NA_clean <- function(x) { #filter well data for Needs Assessment (date/chemical)
  x$DATE <- lubridate::mdy(x$DATE)
  df <- x %>% filter(CHEMICAL %in% chem_list$chemVVL, 
                 DATE >= ymd("2000-01-01"), DATE <= ymd("2020-01-01"))
  return(df)
}
wells_from_data <- function(x) {
  wells <- x %>% select(WELL_ID) %>% distinct()
  left_join(wells, loc)
}
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
chem_list <- read.csv("Rfiles/NA_allchems.csv", stringsAsFactors = F)        
chem_list <- chem_list %>% filter(CCT == "MCL" | chemVVL %in% c("CU", "PB", "NNSM", "CR6")) %>%
  filter(!chemVVL %in% c("ASBESTOS", "COLIFORM", "FCOLIFORM", "RN-222"))

#1. download data ####
#EI well category dataset - identify domestic wells
well_source <- read.csv("Rfiles/gama_raw_well_export.csv", stringsAsFactors = F)
well_domestics <- well_source %>% filter(Category == "Domestic")

#well location information (the well categories here are incorrect- the EI file above is used to identify well type)
loc_link <- 'https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_location_construction_gis.zip'
temp <- tempfile()
download.file(loc_link, temp)
loc <- readr::read_tsv(unz(temp, 'gama_location_construction_gis.txt'), guess_max = 100000)
unlink(temp)
names(loc) <- fix_pun(loc)
loc <- left_join(loc, well_source) %>% 
  select(WELL_ID, LATITUDE, LONGITUDE, WELL_TYPE, Category, COUNTY,
         WELL_DEPTH_FT, TOP_OF_SCREEN_FT, SCREEN_LENGTH_FT, SOURCE)
#remove EDF wells without a well type
edf_rm <- loc %>% filter(SOURCE == "EDF", is.na(Category))
loc <- anti_join(loc, edf_rm)

#DHS/DDW water quality data/wells
ddw_link <- 'https://geotracker.waterboards.ca.gov/gama/data_download/gama_ddw_statewide.zip'
temp <- tempfile()
download.file(ddw_link, temp)
DDW <- readr::read_tsv(unz(temp, 'gama_ddw_statewide.txt'), guess_max = 1000)
unlink(temp)
names(DDW) <- fix_pun(DDW)                                              #fix punctuation in column headers
DDW <- NA_clean(DDW)                                                    #only take chems of interest, only look at last 20 years of data
DDW <- DDW[!grepl("SPRING", DDW$OTHER_NAMES),]                          #remove springs from DHS dataset
write.csv(DDW, paste0("DDW_data_", Sys.Date(), ".csv"), row.names = F)
DDW_wells <- wells_from_data(DDW)                                       #make list of unique wells and join to location data
DDW_wells <- DDW_wells %>% filter(!is.na(LATITUDE))                     

#DWR water quality data/wells
dwr_link <- 'https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_dwr_statewide.zip'
temp <- tempfile()
download.file(dwr_link, temp)
DWR <- readr::read_tsv(unz(temp, 'gama_dwr_statewide.txt'), guess_max = 1000)
unlink(temp)
names(DWR) <- fix_pun(DWR)                                              #fix punctuation in column headers
DWR <- NA_clean(DWR)                                                    #only take chems of interest, only look at last 20 years of data
write.csv(DWR, paste0("DWR_data_", Sys.Date(), ".csv"), row.names = F)
DWR_wells <- wells_from_data(DWR)                                       #make list of unique wells and join to location data

#EDF (AGLAND/LOCAL GW PROJECTS) water quality data/wells
#note - the EDF is too large to download statewide, so this is a workaround where data is downloaded and filtered by county
GT_wells_AG_Local <- loc %>% filter(SOURCE == "AGLAND" | SOURCE == "LOCALGW" | Category == "Domestic")
#join GeoTracker wells with location data to determine unique counties with GT data
cl <- GT_wells_AG_Local %>% select(COUNTY) %>% distinct()
cl$COUNTY <- tolower(cl$COUNTY)                       #make lower case and remove NAs
cl <- cl %>% filter(!is.na(COUNTY), !COUNTY == "no county found")
cl <- cl[order(cl$COUNTY),]

#create empty dataframe to host downloaded GT/EDF data from gama site
edf <- data.frame()
col_sty <- cols(
  `WELL ID` = col_character(),
  RESULTS = col_double(),
  CHEMICAL = col_character(),
  DATE = col_character(),
  UNITS = col_character(),
  QUALIFER = col_character(),
  RL = col_character(),
  LATITUDE = col_double(),
  LONGITUDE = col_double(),
  `WELL TYPE` = col_character(),
  `WELL DEPTH (FT)` = col_character(),
  `TOP OF SCREEN (FT)` = col_character(),
  `SCREEN LENGTH (FT)` = col_character(),
  SOURCE = col_character(),
  `SOURCE NAME` = col_character(),
  `OTHER NAMES` = col_character()
)
for (x in 1:nrow(cl)) {
  baselink <- "https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_edf_"
  county <- gsub(" ", "", cl[x,1], fixed = TRUE)
  countyedflink <- paste0(baselink, county, ".zip")
  temp <- tempfile()
  download.file(url = countyedflink, destfile = temp, method = 'curl')
  newdat <- readr::read_tsv(file = unzip(zipfile = temp), quote = '', col_types = col_sty, skip_empty_rows = T)
  unlink(temp)
  names(newdat) <- fix_pun(newdat)
  newdat <- newdat %>% filter(WELL_ID %in% loc$WELL_ID)
  newdat <- NA_clean(newdat)
  edf <- rbind(edf, newdat)
  print(cl[x,1])     #to keep track of progress!
  rm(newdat)
}
write.csv(edf, paste0("EDF_AG_LG_data_", Sys.Date(), ".csv"), row.names = F)
edf_wells <- wells_from_data(edf)                                 #make list of unique wells and join to location data

#USGS, NWIS, GAMA DOMESTIC water quality data/wells
gama_link <- 'https://geotracker.waterboards.ca.gov/gama/data_download/gama_gama_statewide.zip'
usgs_link <- 'https://geotracker.waterboards.ca.gov/gama/data_download/gama_usgs_statewide.zip'
usgsnwis_link <- 'https://gamagroundwater.waterboards.ca.gov/gama/data_download/gama_usgsnwis_statewide.zip' #gama_usgsnew_statewide
ca_data <- data.frame()
link_list <- c(gama_link, usgs_link, usgsnwis_link)
filenames <- c("gama_gama_statewide.txt","gama_usgs_statewide.txt", "gama_usgsnew_statewide.txt")
for (n in 1:length(link_list)) {
  temp <- tempfile()
  download.file(link_list[n], temp)
  mytable <- readr::read_tsv(unz(temp, filenames[n]), guess_max = 10000)
  unlink(temp)
  names(mytable) <- fix_pun(mytable)
  mytable$WELL_ID <- as.character(mytable$WELL_ID)                      #force characters
  mytable$SOURCE_NAME <- as.character(mytable$SOURCE_NAME)
  mytable$OTHER_NAMES <- as.character(mytable$OTHER_NAMES)
  ca_data <- rbind(ca_data, mytable)                                    #grow big data table (continue appending new data)
  rm(mytable)
}
unique(ca_data$SOURCE)                                                  #check to make sure all datasets were downloaded
ca_data <- ca_data %>% filter(!WELL_ID == "Vanadium Study")             #remove Vanadium Study wells (mistake - duplicate WELLIDs)
USGS <- NA_clean(ca_data)                                               #only take chems of interest, only look at last 20 years of data
dupes <- read.csv("Rfiles/NWIS_dupes2.csv", stringsAsFactors = F)           #duplicates (from GAMA)
USGS <- USGS %>% filter(!WELL_ID %in% dupes$Duplicate_NWIS,
                         !WELL_ID == "USGS-415953121292901")
write.csv(USGS, paste0("USGS_GAMA_NWIS_data_", Sys.Date(), ".csv"), row.names = F)
rm(ca_data)
USGS_wells <- wells_from_data(USGS)                                     #make list of unique wells and join to location data

#compile
wells <- rbind(DDW_wells, DWR_wells, edf_wells, USGS_wells)
wqdata <- rbind(DDW, DWR, edf, ca_data)
#write.csv(wells, paste0("wells", Sys.Date(), ".csv"), row.names = F)

#2. create depth filter ####
#import groundwater unit boundaries
gu_link <- "https://pubs.usgs.gov/ds/796/downloads/ds796_GIS.zip"
temp <- tempfile()
temp2 <- tempfile()
download.file(gu_link, temp)
unzip(zipfile = temp, exdir = temp2)
GU <- readOGR(dsn = paste0(temp2, "/ds796_GIS"), layer = "CA_Groundwater_Units")
unlink(temp)
unlink(temp2)

#import MTRS boundaries (no interface with ArcGIS Online via R, so must be downloaded separately)
#####DomWellD_1 is min, DomWellD_2 is max, DomWellDep is average
MTRS <- readOGR(dsn = "~/R_projects", 
                layer = "i07_WellReportStatsBySection_1009")

#define coordinate system for points
WGS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#transform groundwater unit spatial object to correct coordinate system, and create sf object
GU_wgs <- spTransform(GU, CRS(WGS))
GU_wgs_sf <- st_as_sf(GU_wgs)
MTRS_wgs <- spTransform(MTRS, CRS(WGS))
MTRS_wgs_sf <- st_as_sf(MTRS_wgs)
rm(MTRS, MTRS_wgs, GU, GU_wgs)

MTRS_GU <- st_join(MTRS_wgs_sf, GU_wgs_sf)

#remove duplicate sections (that are duplicated across county line splits of PLSS sections)
depth <- st_drop_geometry(MTRS_GU)
GU <- st_drop_geometry(GU_wgs_sf)
delete <- c("COUNTY_CD", "WCRFolderL", "AreaKm")
depth <- depth[, !(colnames(depth) %in% delete), drop=FALSE] %>% distinct()
numeric_cols <- c("DomWellCou", "DomWellDep", "DomWellD_1", "DomWellD_2",
                 "PubWellCou", "PubWellDep", "PubWellD_1", "PubWellD_2")
depth <- depth %>% mutate_at(numeric_cols, as.numeric)

# remove depth zeros (sections without any wells with depth information) (some sections have domestic wells but no depth data associated
# with them, so I am filtering on depth, not count)
clean_d <- depth %>% filter(!DomWellD_2 == 0) #for domestic wells
clean_p <- depth %>% filter(!PubWellD_2 == 0) #for public wells

#double check (should be zero)
length(which(clean_d$DomWellD_2 == 0))
length(which(clean_p$PubWellD_2 == 0))

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

#add full list of units (including ones that have no groundwater depth filter)
GU_depths <- left_join(GU, alldata)

#write file to working directory
write.table(GU_depths, "well_depth_units.txt", sep = "\t", row.names = F)

#create table with only essential fields
GU_depths_slim <- GU_depths %>% select(GU_ID, GU_Name, maxdomestic, use_public_as_dom)
GU_depths_slim$maxdomestic <- round(GU_depths_slim$maxdomestic, digits = 2)
colnames(GU_depths_slim)[colnames(GU_depths_slim) == "maxdomestic"] <- "max_domestic_depth_ft"
write.table(GU_depths_slim, "depth_filter_output.txt", sep = "\t", row.names = F)

#3. apply depth filter ####
#convert well list to sf object
pnts_sf <- st_as_sf(wells, coords = c('LONGITUDE', 'LATITUDE'), crs = CRS(WGS))

#determine which points are in which groundwater units and MTRS areas
pnts_int <- st_join(pnts_sf, left = TRUE, GU_wgs_sf["GU_ID"])
pnts_int <- st_join(pnts_int, left = TRUE, MTRS_wgs_sf["MTRS"])
wells_int <- st_drop_geometry(pnts_int)

GU_depths_sm <- GU_depths %>%
  select("GU_ID", "avgmax_dom", "maxdomestic", "mindomestic", "use_public_as_dom")
wells_int <- wells_int %>%
  select("WELL_ID", "WELL_TYPE", "WELL_DEPTH_FT",
         "TOP_OF_SCREEN_FT", "SCREEN_LENGTH_FT", "GU_ID", "SOURCE", "MTRS")

wells_int <- left_join(wells_int, GU_depths_sm, by = c("GU_ID"))
wells_int <- wells_int %>% dplyr::filter(!is.na(GU_ID))                       #remove 14 wells that are not in a GW unit

numeric_cols <- c("WELL_DEPTH_FT", "TOP_OF_SCREEN_FT", "SCREEN_LENGTH_FT")
wells_int <- wells_int %>% mutate_at(numeric_cols, as.numeric)
wells_int$WELL_DEPTH[wells_int$WELL_DEPTH == 0] <- NA                           #replace 0s with NA
wells_int$TOP_OF_SCR[wells_int$TOP_OF_SCR == 0] <- NA
wells_int$SCREEN_LEN[wells_int$SCREEN_LEN == 0] <- NA
wells_int <- mutate(wells_int, depth = ifelse(is.na(wells_int$WELL_DEPTH),          #determine numeric depth, if possible
                               ifelse(is.na(wells_int$TOP_OF_SCR) | is.na(wells_int$SCREEN_LEN), NA, (wells_int$TOP_OF_SCR + wells_int$SCREEN_LEN)), 
                               wells_int$WELL_DEPTH))
#filter on 1) if domestic, 2) if numeric, apply numeric filter OR if depth unknown apply unit use
allwells <- wells_int %>% mutate(domesticdepth1 = ifelse(WELL_TYPE == "DOMESTIC" | SOURCE == "GAMA" | is.na(use_public_as_dom), "yes", #domestic wells get a YES
                                                    ifelse(!is.na(depth) & !depth == 0, 
                                                           ifelse(depth <= maxdomestic & depth >= mindomestic, "yes", "no"), #numeric wells get pushed to this criteria, get a yes or no
                                                    ifelse(use_public_as_dom == "yes", "yes", "no"))))
#dom_wells <- filtered_wells %>% filter(domesticdepth == "yes")
write.table(allwells, paste0("allwells", Sys.Date(), ".txt"), sep = "\t", row.names = F)

#4. standardize and clean data ####
md_sm <- wqdata %>% dplyr::filter(WELL_ID %in% wells$WELL_ID) %>%              #filter out unused data
  dplyr::select(WELL_ID, RESULTS, CHEMICAL, DATE, QUALIFER, RL, SOURCE)
md_sm$DATE <- lubridate::ymd(md_sm$DATE)                                #read dates as dates

#convert Nitrate/Nitrite combined to Nitrate
nitrogen <- md_sm %>% dplyr::filter(CHEMICAL %in% c("NO3N", "NO2", "NO3NO2N")) %>%    #nitrogen measurements
  dplyr::mutate(uID = paste0(WELL_ID, DATE))
combined <- nitrogen %>% dplyr::filter(CHEMICAL == "NO3NO2N") %>%                     #NO3NO2N measurements
  dplyr::select(uID, DATE, SOURCE) %>% distinct() 
nitrate <- nitrogen %>% dplyr::filter(CHEMICAL == "NO3N") %>%                         #nitrate (NO3N) measurements
  dplyr::select(uID, SOURCE) %>% distinct()
nitrite <- nitrogen %>% dplyr::filter(CHEMICAL == "NO2") %>%                          #nitrite (NO2) measurements
  dplyr::select(uID, DATE, SOURCE) %>% distinct()
nonitrate <- setdiff(combined$uID, nitrate$uID)                                       #measurmements without nitrate (to be substituted or converted)
sub_data <- nitrogen %>% dplyr::filter(uID %in% setdiff(nonitrate, nitrite$uID)) %>%  #substitution measurements
  select(WELL_ID, RESULTS, CHEMICAL, DATE) %>% mutate(newCHEMICAL = "NO3N")
con_data <- nitrogen %>% dplyr::filter(uID %in% setdiff(nonitrate, setdiff(nonitrate, nitrite$uID)))  #conversion measurements
con_wide <- con_data %>% select(WELL_ID, RESULTS, CHEMICAL, DATE) %>%                                 #conversion wide format
  spread(CHEMICAL, RESULTS) %>% mutate(RESULTS = NO3NO2N - NO2, newCHEMICAL = "NO3N", CHEMICAL = "NO3NO2N")
con_wide$NO2 <- NULL
con_wide$NO3NO2N <- NULL
n_add <- rbind(con_wide, sub_data)                                      #recombine adjusted nitrogen meas. with full data
md_sm <- left_join(md_sm, n_add, by = c("WELL_ID", "DATE", "CHEMICAL")) %>%  
  mutate(adjCHEMICAL = ifelse(is.na(newCHEMICAL), CHEMICAL, newCHEMICAL),
         RESULTS = ifelse(is.na(newCHEMICAL), RESULTS.x, RESULTS.y))
md_sm$newCHEMICAL <- NULL
md_sm$CHEMICAL <- NULL
colnames(md_sm)[colnames(md_sm)=="adjCHEMICAL"] <- "CHEMICAL"
md_sm$RESULTS.x <- NULL
md_sm$RESULTS.y <- NULL
md_c <- md_sm %>% dplyr::filter(!CHEMICAL == "NO3NO2N")                 #delete dupes
rm(nitrate, nitrogen, nitrite, n_add, nonitrate)
rm(md_sm)

#standardize Qualifiers to either "<" or "="
md_c$QUALIFER[is.na(md_c$QUALIFER)] <- "="                        
md_c <- dplyr::mutate(md_c, newXMOD = ifelse(QUALIFER == "<" | QUALIFER == "ND", "<", "="))

#standardize RL/Results (fix NAs, RESULTS == 0, and missing reporting limits)
md_c$RL <- as.numeric(md_c$RL)
finishedchems <- c()                                              #hold chemicals as they are cleaned
additions <- data.frame()                                         #holds rows with "newRL"s
chemicals <- unique(md_c$CHEMICAL)                                #list of chemicals that need to be fixed
md_c2 <- dplyr::arrange(md_c, DATE)
chemicals <- sort(chemicals)
for (ch in chemicals) {                                                 #or for (ch in chemicals_2try) {
  uRL <- md_c2 %>% dplyr::filter(is.na(RL), RESULTS <= 0, CHEMICAL == ch) #list of measurements with unknown RL
  kRL <- md_c2 %>% dplyr::filter(CHEMICAL == ch, !is.na(RL))             #list of measurements with known RL
  if(dim(kRL)[1] == 0) {                                                #if chemical has no RLs listed, take implied non-detects
    kRL <- md_c2 %>% dplyr::filter(CHEMICAL == ch, newXMOD == "<", !RESULTS == 0) %>% mutate(RL = RESULTS) 
  }
  if(dim(uRL)[1] > 0 & dim(kRL)[1] > 0) {
    indx_b <- neardate(uRL$CHEMICAL, kRL$CHEMICAL, uRL$DATE, kRL$DATE, best = "prior") #get index of closest RL prior DATE
    indx_a <- neardate(uRL$CHEMICAL, kRL$CHEMICAL, uRL$DATE, kRL$DATE, best = "after") #get index of closest RL after DATE
    estRLs <- kRL[ifelse(is.na(indx_b), indx_a, indx_b), c("RL", "DATE")]              #assign newRLs and date of newRL, first looking at priors and then afters, if no priors exist
    uRL$estRL <- estRLs$RL                                                             #add new RLs to measurement
    uRL$estRLdate <- estRLs$DATE                                                       #add date new RL came from, to identify potential mistakes/keep track of method
    print(ifelse(anyNA(uRL$estRL), nrow(uRL %>% filter(is.na(estRL))), paste0("Clear for ", ch))) #check if RLs have been standardized, if not it prints number of problem rows
    additions <- rbind(additions, uRL)                                                 #accumulate table of unknown RLs (with new RLs) for joining later
    finishedchems <- append(finishedchems, ch)                                         #keep track of successful chemicals
    rm(indx_a, indx_b, estRLs)
  } else {
    print(paste0("Clear for ", ch))
  }
  rm(uRL, kRL)
} #assigns a RL to non-detects without one by looking at closest earlier RL of the same chemical
md_c2 <- left_join(md_c2, additions)
rm(additions)

#compress RLs into a single column, if applicable
md_c2 <- mutate(md_c2, ReportingLimit = ifelse(!is.na(RL), RL, 
                                        ifelse(!is.na(estRL), estRL, RESULTS)),
               RLmethod = ifelse(!is.na(RL), "reported", 
                          ifelse(!is.na(estRL), "estimated", "results")))
md_c2 <- left_join(md_c2, chem_list, by = c("CHEMICAL" = "chemVVL"))      #join with chemical MCLs
md_c2$chemNAME <- NULL

#calculate MCL index and track detections/non-detections, outliers
md_c2 <- md_c2 %>% mutate(detection = ifelse(newXMOD == "<", "non-detect", 
                                    ifelse(RESULTS <= 0, "non-detect",        #negative radioactive measurments are non-detects
                                    ifelse(RLmethod == "reported" & RESULTS <= RL, "non-detect", "detect")))) %>%
  mutate(RES = ifelse(detection == "non-detect", RLfun(ReportingLimit, CCV), RESULTS)) %>% #compress results to one number
  mutate(MCLindex = RES/CCV)                                            #calculate index
md_max <- summarize(group_by(md_c2 %>% dplyr::filter(detection == "detect"), CHEMICAL), avg = mean(RES), sd = sd(RES)) %>% 
  mutate(max = (avg + (10*sd))) %>% select(CHEMICAL, max)               #calculates outlier maximum
#write.csv(md_max, paste0(run, "outlierVALS_", Sys.Date(), ".csv"), row.names = F) 
md_c2 <- left_join(md_c2, md_max, by = "CHEMICAL") %>% 
  mutate(outlier_status = ifelse(is.na(max), "no outliers", 
                                 ifelse(RES > max & detection == "detect", "outlier", "standard")))

#clean table and output results
refcols <- c("WELL_ID", "SOURCE", "DATE", "RESULTS", "QUALIFER")
md_c <- md_c2[,c(refcols, setdiff(names(md_c2), refcols))]
write.csv(md_c, paste0("WQdata", Sys.Date(), ".csv"), row.names = F)

#domestic well wq results only
domwells <- allwells %>% filter(is.na(domesticdepth1) | domesticdepth1 == "yes")
domwq <- md_c %>% filter(WELL_ID %in% domwells$WELL_ID)
write.csv(domwq, paste0("domwellWQdata", Sys.Date(), ".csv"), row.names = F)
