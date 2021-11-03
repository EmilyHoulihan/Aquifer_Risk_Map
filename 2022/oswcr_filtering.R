library(tidyverse)
library(lubridate)
library(httr)
library(sf)
library(data.table)

st_drop_geometry <- function(x) {
  if(inherits(x,"sf")) {
    x <- st_set_geometry(x, NULL)
    class(x) <- 'data.frame'
  }
  return(x)
}

#load tabular oswcr data (takes several minutes)
#oswcr <- read.csv('https://data.cnra.ca.gov/dataset/647afc02-8954-426d-aabd-eff418d2652c/resource/8da7b93b-4e69-495d-9caa-335691a1896b/download/wellcompletionreports.csv', stringsAsFactors = F)
oswcr <- read.csv('C:/Users/EHoulihan/Documents/SB 200/ARM/ARM_project/OSWCR/wellcompletionreports.csv') #updated 9/24/2021

#isolate domestic use wells
#oswcr_dom <- oswcr %>% filter(str_detect(PlannedUseFormerUse, "omestic"))
oswcr_dom <- oswcr %>% filter(str_detect(PLANNEDUSEFORMERUSE, "omestic"), !str_detect(RECORDTYPE, "estruction"))

#read dates as dates
#oswcr_dom$DateWorkEnded <- lubridate::ymd_hms(oswcr_dom$DateWorkEnded)
oswcr_dom$DATEWORKENDED <- lubridate::mdy(oswcr_dom$DATEWORKENDED)

#append MTRS to each row
#oswcr_dom <- oswcr_dom %>% mutate(BM = ifelse(str_detect(BaselineMeridian, "Mount Diablo"), "M",
#                                        ifelse(str_detect(BaselineMeridian, "Humboldt"), "H",
#                                        ifelse(str_detect(BaselineMeridian, "San Bernardino"), "S", BaselineMeridian))),
#                                  MTRS = paste0(BM, Township, Range, Section))
oswcr_dom <- oswcr_dom %>% mutate(BM = ifelse(str_detect(BASELINEMERIDIAN, "Mount Diablo"), "M",
                                              ifelse(str_detect(BASELINEMERIDIAN, "Humboldt"), "H",
                                                     ifelse(str_detect(BASELINEMERIDIAN, "San Bernardino"), "S", BASELINEMERIDIAN))),
                                  MTRS = paste0(BM, TOWNSHIP, RANGE, SECTION))
#get statewide counts by section
url <- list(hostname = "gis.water.ca.gov/arcgis/rest/services",
            scheme = "https",
            path = "Environment/i07_WellCompletionReports/FeatureServer/1/query",
            query = list(where = "1=1",
                         outFields = "*",
                         returnGeometry = "true",
                         f = "geojson")) %>%
  setattr("class","url")
request <- build_url(url)
MTRS_sf <- st_read(request) %>% st_drop_geometry()
MTRS_sf_small <- MTRS_sf %>% select(MTRS, DomWellCount) %>% distinct()
MTRS_sf_small$DomWellCount[is.na(MTRS_sf_small$DomWellCount)] <- 0
write.csv(MTRS_sf_small, paste0("i07_WellCompletionReports", Sys.Date(), ".csv"), row.names = F)

#total number of domestic wells (only with real MTRS names)
oswcr_dom %>% filter(MTRS %in% MTRS_sf_small$MTRS) %>% nrow()

#total number of domestic wells post 1970
oswcr_dom %>% filter(MTRS %in% MTRS_sf_small$MTRS, DATEWORKENDED >= lubridate::ymd("1970-01-01")) %>% nrow()

#total number of domestic wells post 1970 including wells with no dates
oswcr_dom %>% filter(MTRS %in% MTRS_sf_small$MTRS) %>% 
  filter(DATEWORKENDED >= lubridate::ymd("1970-01-01") | is.na(DATEWORKENDED)) %>% nrow()

#exclude wells drilled prior to 1970
#keep wells without a date
oswcr_dom_date <- oswcr_dom %>%
  filter(is.na(DATEWORKENDED) | DATEWORKENDED >= lubridate::ymd("1970-01-01")) %>%
  group_by(MTRS) %>%
  summarize(dom_wells1970na = n())

#attach new (post 1970 counts) to original stats by section
oswcr_dom_final <- left_join(MTRS_sf_small, oswcr_dom_date)
head(oswcr_dom_final)
oswcr_dom_final <- oswcr_dom_final %>% distinct()
oswcr_dom_final$dom_wells1970na[is.na(oswcr_dom_final$dom_wells1970na)] <- 0
write.csv(oswcr_dom_final, paste0("Datasets/oswcr_dom_counts", Sys.Date(), ".csv"), row.names = F)
