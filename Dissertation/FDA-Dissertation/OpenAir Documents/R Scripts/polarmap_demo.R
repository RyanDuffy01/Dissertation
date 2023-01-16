
library(openair)
library(openairmaps)

# Lots of sites around SE England -----------------------------------------

# choose some sites in the SE
select_sites <- c("NO12",  "cant",  "kc1", "brt3", "sdy", "hope")

# import AURN - remember meta = TRUE!
aq_meas <- importAURN(site = select_sites, year = 2018, meta = TRUE)

# polar plots
polarPlot(
  mydata = aq_meas,
  pollutant = "no2",
  type = "site"
)

# polar map (complicated example!)
aq_meas %>%
  # get season column to split by
  cutData("season") %>%
  # build a popup (site, site type, average NO2)
  buildPopup(c("site", "site_type", "no2"),
             control = "season",
             latitude = "latitude", longitude = "longitude",
             names = c("Site" = "site", "Site Type" = "site_type")) %>%
  # make polar map
  polarMap(
    latitude = "latitude", 
    longitude = "longitude",
    pollutant = "no2",
    control = "season", # layer control menu using season
    limits = c(0, 40), # all markers share the same limits
    cols = "inferno", # new colour scheme? See `?openair::openColours`
    popup = "popup" # use popup we just built!
  )

# windrose map? (can use six different polar-type plots as markers!)
aq_meas %>%
  windroseMap(
    latitude = "latitude",
    longitude = "longitude",
    label = "site",
    cols = "viridis",
    provider = "CartoDB.Positron",
    max.freq = 20 # use same radial scale for all
  )

# Look at different years (e.g., over lockdown) ---------------------------

# import multiple years
my1 <- importAURN("my1", year = 2018:2021, meta = TRUE)

my1 %>%
  # get "year" column
  cutData(type = "year") %>%
  # make map
  polarMap(
    pollutant = "nox",
    control = "year", # control by year
    limits = c(0, 400)
  )


