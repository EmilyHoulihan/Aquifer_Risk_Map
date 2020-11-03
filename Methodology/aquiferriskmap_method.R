#Aquifer risk map methodology

#Map: https://gispublic.waterboards.ca.gov/portal/apps/webappviewer/index.html?id=b5488c0911314ef6a4796245dd584d0b
#Methodology: https://gispublic.waterboards.ca.gov/portal/home/item.html?id=5f94cf5ca2a04f66a3dcfab29cf6c764

#updated 10/21/2020 by Emily Houlihan, GAMA Unit (Emily.Houlihan@Waterboards.ca.gov)

#This script is provided as a documentation of methodology steps. If you would like to 
#re-run this script, you will need additional reference files. Please check
#the aquifer risk map metadata page, or contact Emily Houlihan.

#Components
#0. Load data
#1. Calculates well point averages and recent results
#2. Calculates square mile section averages and recent results
#3. Calculates census block group aggregated scores

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
if (!require('readxl')) install.packages('readxl'); library('readxl')

#0. Load water quality and census data####
welldata <- read.table("domwqdata_clean2020-10-29.txt", stringsAsFactors = F, header = T)
allwells <- read.delim("allwells2020-10-28.txt", stringsAsFactors = F)
domwells <- allwells %>% filter(domdepth == "yes")
#wellsmtrs <- read.dbf("Datasets/wells2020-08-07_MTRSjoin.dbf", as.is = T)
CBG <- read.dbf("Datasets/bg_mtrs_join.dbf", as.is = T)

#census information
cbg <- CBG %>% select(MTRS, DomWellCou, GEOID) %>% distinct()

#1. Calculate point risk scores (wells)####
welldata$DATE <- lubridate::ymd(welldata$DATE)
#wellsmtrs$MTRS[wellsmtrs$WELL_ID == "3110031-001"] <- "M15N16E24" #in Lake Tahoe, want to keep
domwelldata <- welldata %>% 
  filter(!outlier_status == "outlier", !CHEMICAL == "TCPR123_2", !is.na(MCLindex)) %>%
  select(WELL_ID, DATE, CHEMICAL, MCLindex) %>%
  mutate(year = year(DATE))
rm(welldata)

#CR6 adjustment - change to 10 ug/L (original has it at 20 ccv)
domwelldata <- domwelldata %>% mutate(MCLindex = ifelse(CHEMICAL == "CR6", MCLindex*2, MCLindex))

#de-cluster - calculate long-term average by year, then by well
domwellavg <- domwelldata %>% 
  group_by(CHEMICAL, WELL_ID, year) %>%
  summarize(yr_avg = mean(MCLindex))%>%
  group_by(CHEMICAL, WELL_ID) %>%
  summarize(well_avg = mean(yr_avg))

#de-cluster - calculate recent result groupings (last 2 years)
domwellre_long <- domwelldata %>% filter(DATE >= ymd("2018-01-01")) %>%
  mutate(re_status = ifelse(MCLindex > 1, 'over',
                            ifelse(MCLindex > 0.8, 'close', 'under'))) %>% 
  group_by(CHEMICAL, WELL_ID, re_status) %>% count()
domwellre_avg <- domwelldata %>% filter(DATE >= ymd("2018-01-01"), MCLindex > 1) %>% 
  group_by(WELL_ID, CHEMICAL) %>%
  summarize(re_avg_over = mean(MCLindex))
domwellre <- pivot_wider(domwellre_long, names_from = re_status, names_prefix = "re_", values_from = n) %>%
  left_join(., domwellre_avg)
domwellre[is.na(domwellre)] <- 0
rm(domwellre_avg, domwellre_long)

#joing all point data information (long-term average and recent results)
pointdata <- full_join(domwellavg, domwellre)
pointdata[is.na(pointdata)] <- 0
pointdata <- pointdata %>% mutate(avg_overMCL = ifelse(well_avg > 1 & re_avg_over > 0, (well_avg + re_avg_over)/2,
                                               ifelse(well_avg > 1 & re_avg_over == 0, well_avg,
                                                      ifelse(!well_avg > 1 & re_avg_over > 0, re_avg_over, 0))))
wellinfo <- allwells %>% select(WELL_ID, LATITUDE, LONGITUDE, Category, SOURCE)
pointdata <- left_join(pointdata, wellinfo)
#write.csv(pointdata, "pointdatadetailed_cr6adj.csv", row.names = F)
rm(domwelldata, welldata)

#determine risk for each point location
pdrisk <- pointdata %>% group_by(WELL_ID) %>% arrange(-avg_overMCL) %>%
  summarize(PRF1 = n_distinct(CHEMICAL[well_avg > 1 | re_over > 0]),
            PRF2 = n_distinct(CHEMICAL[well_avg <= 1 & well_avg > 0.8 | re_close > 0]),
            PRF3 = mean(avg_overMCL[avg_overMCL > 0]),
            PL1 = paste(CHEMICAL[well_avg > 1 | re_over > 0], collapse = "; "),
            PL2 = paste(CHEMICAL[well_avg <= 1 & well_avg > 0.8 | re_close > 0], collapse = "; "))
pdrisk <- pdrisk %>% mutate(risk = ifelse(PRF1 > 0, "high",
                                  ifelse(PRF1 == 0 & PRF2 > 0, "med", "low")))

pdrisk <- left_join(pdrisk, wellinfo)
write.csv(pdrisk, paste0("pointdata_risk", Sys.Date(), ".csv"), row.names = F)
rm(pointdata)

#2. Calculate section risk scores (square miles)####
#load file that contains neighbor reference data
neighborsec <- read.dbf("Datasets/MTRS_neighbor_all.dbf", as.is = T) %>% 
  select(MTRS, MTRS_2) %>% distinct() %>%
  filter(!MTRS == MTRS_2)

#de-cluster - calculate source section averages
domwellavg <- left_join(domwellavg, select(domwells, WELL_ID, MTRS)) %>%
  left_join(., neighborsec)
secavg <- domwellavg %>% group_by(CHEMICAL, MTRS) %>%
  summarize(SecDet = mean(well_avg))

#de-cluster - calculate neighbor section averages
neighbor_avg <- left_join(secavg, neighborsec) %>%
  group_by(CHEMICAL, MTRS_2) %>%
  summarize(nSecDet = mean(SecDet))
colnames(neighbor_avg)[colnames(neighbor_avg)== "MTRS_2"] <- "MTRS"
secavg <- full_join(secavg, neighbor_avg) %>% 
  mutate(method = ifelse(is.na(SecDet), "neighbor", "source"),
         SecDet = ifelse(is.na(SecDet), nSecDet, SecDet))
secavg$nSecDet <- NULL
rm(neighbor_avg)

#de-cluster - calculate source section recent results
domwellre <- left_join(domwellre, select(domwells, WELL_ID, MTRS))
secre <- domwellre %>% group_by(CHEMICAL, MTRS) %>%
  summarize(re_under = sum(re_under),
            re_close = sum(re_close),
            re_over = sum(re_over),
            re_avg_over = mean(re_avg_over[re_avg_over > 0]))

#de-cluster - calculate neighbor section recent results
neighbor_re <- left_join(secre, neighborsec) %>% 
  group_by(CHEMICAL, MTRS_2) %>%
  summarize(re_under = mean(re_under),
            re_close = mean(re_close),
            re_over = mean(re_over),
            re_avg_over = mean(re_avg_over, na.rm = T)) %>% 
  filter(!MTRS_2 %in% c(secre$MTRS))
colnames(neighbor_re)[colnames(neighbor_re) == "MTRS_2"] <- "MTRS"
secre <- rbind(secre, neighbor_re)
rm(neighbor_re)

#join average and recent section results
section_data <- left_join(secavg, secre)
section_data <- section_data %>% filter(!is.na(MTRS))
section_data[is.na(section_data)] <- 0

#calculate average magnitude
section_data <- section_data %>% mutate(avg_overMCL = ifelse(SecDet > 1 & re_avg_over > 0, (SecDet + re_avg_over)/2,
                                                 ifelse(SecDet > 1 & re_avg_over == 0, SecDet,
                                                 ifelse(!SecDet > 1 & re_avg_over > 0, re_avg_over, 0))))

#calculate section risk (summarize section_data)
sdrisk <- section_data %>% group_by(MTRS) %>% arrange(-SecDet) %>%
  summarize(SRF1 = n_distinct(CHEMICAL[SecDet > 1 | re_over > 0]),
            SRF2 = n_distinct(CHEMICAL[SecDet <= 1 & SecDet > 0.8 | re_close > 0]),
            SRF3 = mean(avg_overMCL[avg_overMCL > 0]),
            SL1 = paste(CHEMICAL[SecDet > 1 | re_over > 0], collapse = "; "),
            SL2 = paste(CHEMICAL[SecDet <= 1 & SecDet > 0.8 | re_close > 0], collapse = "; "))
#sdrisk[is.na(sdrisk)] <- 0
sdrisk <- sdrisk %>% mutate(WQriskbins = ifelse(SRF1 > 0, "high",
                                  ifelse(SRF2 > 0 & SRF1 == 0, "med", "low")))
write.csv(sdrisk, paste0("sectiondata_risk", Sys.Date(), ".csv"), row.names = F)
#rm(domwellavg, domwellre)

#3a. Census block group data - water quality risk scores####
#join full section data to census blocks
cdrisk <- left_join(section_data, cbg)
#calculate census risk factors 1-3
cdrisk <- cdrisk %>% group_by(GEOID) %>% arrange(-SecDet) %>% #arrange by magnitude
  summarize(CRF1 = n_distinct(CHEMICAL[SecDet > 1 | re_over > 0]),
            CRF2 = n_distinct(CHEMICAL[SecDet > 0.8 & SecDet <= 1 | re_close > 0]),
            CRF3 = mean(avg_overMCL[avg_overMCL > 0]),
            CL1 = paste(unique(CHEMICAL[SecDet > 1 | re_over > 0]), collapse = "; "),
            CL2 = paste(unique(CHEMICAL[SecDet > 0.8 & SecDet <= 1 | re_close > 0]), collapse = "; "))
cdrisk <- cdrisk %>% filter(!is.na(GEOID))
cdrisk[is.na(cdrisk)] <- 0

#join condensed section data to census blocks
cd1 <- left_join(sdrisk, cbg)
cd1 <- cd1 %>% group_by(GEOID) %>%
  summarize(highsec = n_distinct(MTRS[WQriskbins == "high"]),
            medsec = n_distinct(MTRS[WQriskbins == "med"]),
            datasec = n())
ct_area <- cbg %>% group_by(GEOID) %>%
  summarize(areasec = n())
cdrisk <- full_join(cdrisk, cd1) %>%
  full_join(., ct_area)
#calculate census risk factors 4-5 (here just displayed as "4") and overall wq score
cdrisk <- cdrisk %>% mutate(CRF4 = (highsec + medsec/2)/areasec*100,
                            wq_score = (CRF1 + CRF2/2 + CRF3/10)*(CRF4))

#convert water quality score to percentiles (remove 0 entries first)
census_risk_0 <- cdrisk %>% filter(is.na(wq_score) | wq_score == 0) %>%
  mutate(wq_per = 0)
census_risk_pos <- cdrisk %>% filter(wq_score > 0) %>%
  mutate(wq_per = ntile(wq_score, 100))

cdrisk_wq <- rbind(census_risk_0, census_risk_pos)

#edit table so it can be read in ArcGIS more easily
cdrisk_wq$datasec[is.na(cdrisk_wq$datasec)] <- 0
cdrisk_wq$CRF1[is.na(cdrisk_wq$CRF1)] <- 0
cdrisk_wq$CRF2[is.na(cdrisk_wq$CRF2)] <- 0
cdrisk_wq$CRF3[is.na(cdrisk_wq$CRF3)] <- 0
cdrisk_wq$highsec[is.na(cdrisk_wq$highsec)] <- 0
cdrisk_wq$medsec[is.na(cdrisk_wq$medsec)] <- 0
cdrisk_wq$datasec[is.na(cdrisk_wq$datasec)] <- 0
cdrisk_wq$CRF4[is.na(cdrisk_wq$CRF4)] <- 0
cdrisk_wq$wq_score[is.na(cdrisk_wq$wq_score)] <- 0
cdrisk_wq <- cdrisk_wq %>% 
  mutate(datacov = datasec/areasec*100)
cdrisk_wq$datacov <- round(cdrisk_wq$datacov, digits = 2)
cdrisk_wq$CRF3 <- round(cdrisk_wq$CRF3, digits = 2)
cdrisk_wq$wq_score <- round(cdrisk_wq$wq_score, digits = 2)
cdrisk_wq$CRF4 <- round(cdrisk_wq$CRF4, digits = 2)

#verify that census block groups are named correctly (opening a .csv in excel will remove the leading zeros)
cdrisk_wq <- cdrisk_wq %>% mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                                                             paste0("0", GEOID)))
write.csv(cdrisk_wq, paste0("censusdata_risk_wq", Sys.Date(), ".csv"), row.names = F)
rm(neighborsec, secavg, secre)

#3b. Census block group data - hazard + exposure (combined risk scores) ####
#load domestic well information
oswcr_dates <- read.csv("Datasets/oswcr_dom_date1970na.csv", stringsAsFactors = F)
cbg3 <- left_join(cbg, oswcr_dates)
cbg3[is.na(cbg3)] <- 0
cbg_dates <- cbg3 %>% group_by(GEOID) %>% summarize(oswcr_1970na = sum(dom_wells1970na))
cdrisk_com <- left_join(cdrisk_wq, cbg_dates)

#load state small system location information
statesmalls <- read.dbf("Datasets/SSWS_rcac_final_w_monterey_wqrisk.dbf", as.is = T)
ss_count <- statesmalls %>% group_by(GEOID) %>%
  summarize(statesmallcount = n())

#load demographic information
demog <- read.dbf("Datasets/2018bgmhi_dwr.dbf", as.is = T)
head(demog)

#add new risk factors to water quality scores
comprisk <- full_join(cdrisk_com, demog, by = c("GEOID" = "GEOID10")) %>%
  full_join(., ss_count)
comprisk$statesmallcount[is.na(comprisk$statesmallcount)] <- 0

#test count of domestic wells (sections duplicated across census block groups)
sum_test <- cbg %>% select(-GEOID) %>% distinct()
sum(sum_test$dom_wells1970na)

totaloswcr_1970na <- 287142
sumoswcr_1970na <- sum(comprisk$oswcr_1970na, na.rm = T)

#import block group info (counties)
bg_info <- read.dbf("Datasets/tl_2019_06_bg_WGS_geom_m.dbf", as.is = T)
bg_counties <- readxl::read_excel("Datasets/county_codes.xlsx")
bg_info <- left_join(bg_info, bg_counties, by = c("COUNTYFP" = "county_code"))
bg_info <- bg_info %>% mutate(area_sq_mi = area/2589988.1103)
bg_info_sm <- bg_info %>% select(GEOID, County, area_sq_mi)

comprisk <- left_join(comprisk, bg_info_sm)

comprisk <- comprisk %>% mutate(SDAC18 = ifelse(MHI18 < 42737 & !MHI18 == 0, "Y", "N"),
                                adj_dw = oswcr_1970na/sumoswcr_1970na*totaloswcr_1970na,
                                dwss_cou = statesmallcount + adj_dw,
                                dwss_den = dwss_cou/area_sq_mi)

comprisk_dwss_count_0 <- comprisk %>% filter(is.na(dwss_den) | dwss_den == 0) %>%
  mutate(dwss_per = 0)
comprisk_dwss_count_pos <- comprisk %>% filter(dwss_den > 0) %>%
  mutate(dwss_per = ntile(dwss_den, 100))
comprisk <- rbind(comprisk_dwss_count_0, comprisk_dwss_count_pos)

comprisk <- comprisk %>% mutate(prio1 = (wq_per + dwss_per))

comprisk$wq_per[is.na(comprisk$wq_per)] <- 0
comprisk$prio1[is.na(comprisk$prio1)] <- 0

comprisk_prio1_0 <- comprisk %>% filter(is.na(prio1) | prio1 == 0) %>%
  mutate(prio1_per = 0)
comprisk_prio1_pos <- comprisk %>% filter(prio1 > 0) %>%
  mutate(prio1_per = ntile(prio1, 100))

comprisk <- rbind(comprisk_prio1_0, comprisk_prio1_pos)

comprisk$prio1_per[is.na(comprisk$prio1_per)] <- 0

comprisk_write <- comprisk %>% select(GEOID, wq_per, prio1_per, wq_score, dwss_per, prio1,
                                      Pop18, DAC18, SDAC18, adj_dw, statesmallcount, 
                                      dwss_cou, dwss_den, datacov, County)

colnames(comprisk_write)[colnames(comprisk_write) == "wq_per"] <- "wq_risk"
colnames(comprisk_write)[colnames(comprisk_write) == "prio1_per"] <- "comb_risk"
colnames(comprisk_write)[colnames(comprisk_write) == "dwss_per"] <- "ex_risk"

comprisk_write <- comprisk_write %>% mutate(DAC_status = ifelse(DAC18 == "Data Not Available", "Data Not Available",
                                                                ifelse(SDAC18 == "Y", "SDAC",
                                                                       ifelse(DAC18 == "Y", "DAC", "None"))))
comprisk_write <- comprisk_write %>% select(-c("DAC18", "SDAC18", "prio1"))
head(comprisk_write)

comprisk_write <- comprisk_write %>% filter(!is.na(GEOID))
comprisk_write <- comprisk_write %>% mutate(DAC_status = ifelse(is.na(DAC_status), "Data Not Available", DAC_status))
comprisk_write[is.na(comprisk_write)] <- 0
comprisk_write$adj_dw <- round(comprisk_write$adj_dw, digits = 2)
comprisk_write$dwss_cou <- round(comprisk_write$dwss_cou, digits = 2)
comprisk_write$dwss_den <- round(comprisk_write$dwss_den, digits = 3)
write.csv(comprisk_write, paste0("comprisk_", Sys.Date(), ".csv"), row.names = F)

#Extra - make layers for individual chemicals ####
section_data[is.na(section_data)] <- 0
no3n <- sdata3 %>% filter(CHEMICAL == "NO3N")
as <- sdata3 %>% filter(CHEMICAL == "AS")
cr6 <- sdata3 %>% filter(CHEMICAL == "CR6")
tcp <- sdata3 %>% filter(CHEMICAL == "TCPR123")
u <- sdata3 %>% filter(CHEMICAL == "U")

chemical <- cr6 %>% mutate(WQriskbins = ifelse(SecDet > 1 | re_over > 0, "high",
                                               ifelse(SecDet > 0.8 & re_over == 0, "med", "low")))
write.csv(chemical, "section_risk_cr6adj.csv", row.names = F)

