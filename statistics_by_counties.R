#summary statistics

library(tidyverse)
library(readxl)
library(janitor)

bg_risk <- read.csv("combinedrisk_blockgroups.csv", stringsAsFactors = F)
bg_risk <- bg_risk %>% dplyr::mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                                                 paste0("0", GEOID)))
bg_risk <- bg_risk %>% dplyr::mutate(county_code_fips = substr(GEOID, 3, 5))

#get county names
fips_link <- 'https://www2.census.gov/programs-surveys/popest/geographies/2019/all-geocodes-v2019.xlsx'
temp <- tempfile()
download.file(url = fips_link, destfile = temp, mode = 'wb')
fips_codes <- read_xlsx(temp, skip = 4)
unlink(temp)
fips_codes <- clean_names(fips_codes)
fips_ca <- fips_codes %>% dplyr::filter(state_code_fips == "06") %>% 
  dplyr::select(county_code_fips, area_name_including_legal_statistical_area_description)

#join county names to aquifer risk map data
bg_risk_join <- dplyr::left_join(bg_risk, fips_ca, by = "county_code_fips")
colnames(bg_risk_join)[colnames(bg_risk_join) == "area_name_including_legal_statistical_area_description"] <- "county"

#bring in water quality detailed risk information
bg_risk_wq <- read.csv("waterquality_blockgroups.csv", stringsAsFactors = F)
bg_risk_wq <- bg_risk_wq %>% dplyr::mutate(GEOID = ifelse(grepl("^0...........", GEOID), GEOID,
                                                    paste0("0", GEOID)))
bg_risk_wq_small <- bg_risk_wq %>% dplyr::select(GEOID, datacov, CL1, CL2)

#join water quality metadata to combined risk data
bg_risk_join_wq <- dplyr::left_join(bg_risk_join, bg_risk_wq_small, by = "GEOID")

#summarize by county
county_sum <- bg_risk_join_wq %>% group_by(county) %>%
  summarize(avg_wq_risk = mean(wq_risk),
            avg_ex_risk = mean(ex_risk),
            avg_comb_risk = mean(comb_risk),
            count_comb_risk_over90 = sum(comb_risk >= 90),
            avg_datacov = mean(datacov, na.rm = T))

county_sum %>% arrange(-count_comb_risk_over90)
county_sum %>% arrange(avg_datacov)
county_sum %>% arrange(-avg_wq_risk)
county_sum %>% arrange(-avg_comb_risk)

