# load the packages we need
library(openair)
library(tidyverse)
library(openairmaps)

# get AURN meta data
aurn_meta <- importMeta(source = "aurn", all = TRUE)

# import data from MR
mary <- importAURN(site = "my1", year = 2021:2022)

# view map of AURN, AQE and locally managed data
networkMap(
  source = c("aurn", "aqe", "local"),
  control = "site_type",
  date = "2000",
  provider = c("OpenStreetMap", "CartoDB.Positron")
)

sunderland <- importAURN(site = c("sunr", "sun2"), year = 2019)

ne_scot_no2_sites <- 
  filter(aurn_meta,
         zone == "North East Scotland",
         end_date == "ongoing",
         variable == "NO2")

scot_codes <- ne_scot_no2_sites$code

scot_no2 <- importAURN(site = scot_codes, pollutant = "no2", year = 2022)

aurn_meta

networkMap("local")

london <- importKCL(year = 2021:2022, pollutant = "no2")

local_meta <- importMeta("local", all = TRUE)

london_meta <- filter(local_meta,
                      zone == "Greater London")

london_meta <- distinct(london_meta, site, .keep_all = TRUE)

local_london <- importLocal(site = london_meta$code, year = 2021)

filter(local_london, nox > 200)
 
mean(local_london$no2, na.rm = TRUE)

local_london_dropnano2 <- drop_na(local_london, no2)

select(local_meta, starts_with("a"))

filter(local_meta, !stringr::str_detect(string = site, pattern = "a"))

# access annual means

annual_2021 <- importAURN(data_type = "annual", year = 2021)

no2_annual <- select(annual_2021, code, site, no2, no2_capture)

o3_annual <- select(annual_2021, code, site, starts_with("o3"))
names(o3_annual)

aurn_2020 <- importAURN(data_type = "annual", year = 2020, meta = TRUE)
no2_aurn <- mutate(aurn_2020, network = "AURN")

aqe_2020 <- importAQE(data_type = "annual", year = 2020, meta = TRUE)
no2_aqe <- mutate(aqe_2020, network = "AQE")

local_2020 <- importLocal(data_type = "annual", year = 2020, meta = TRUE)
no2_local <- mutate(local_2020, network = "Local")

all_no2 <- bind_rows(no2_aurn, no2_aqe, no2_local)
all_no2 <- filter(all_no2, no2_capture >= 0.85)

summary_no2 <- group_by(all_no2, network) 
summary_no2 <- summarise(summary_no2, no2 = mean(no2))

summary_no2 <- group_by(all_no2, network, site_type) %>% 
  summarise(no2_mean = mean(no2), no2_median = median(no2))






