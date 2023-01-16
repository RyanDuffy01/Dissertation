
library(openair)
library(ggplot2)
library(ggmap)
library(leaflet)

# Import Data -------------------------------------------------------------

aurn <-
  importAURN(
    year = 2018,
    data_type = "annual",
    to_narrow = TRUE,
    meta = TRUE
  )

# Static maps -------------------------------------------------------------

## Vector/Polygon

# get polygon of the UK
mapdat <- map_data("world", "uk")

# filter data
statmap_data <- 
  aurn %>%
  filter(species == "no2") %>%
  arrange(value)

# static plot with polygon
ggplot(statmap_data) +
  geom_polygon(data = mapdat, 
               aes(x = long, y = lat, group = group),
               fill = "grey90") +
  geom_point(aes(x = longitude, y = latitude, color = value),
             alpha = .5) +
  coord_map() +
  scale_color_viridis_c() +
  theme_classic() +
  labs(color = quickText("NO2 (ug/m3)"))

## Raster/Basemap

# define bounding box (play with it a bit!)
left = min(aurn$longitude) * 1.1
right = max(aurn$longitude) * 1.5
down = min(aurn$latitude) * .95
up = max(aurn$latitude) * 1

# import basemap
basemap <-
  get_stamenmap(
    c(left, down, right, up),
    zoom = 6,
    maptype = "terrain-background",
    color = "bw"
  )

# use in ggplot
ggmap(basemap) +
  geom_point(
    data = filter(aurn, species == "no2") %>% arrange(value),
    aes(x = longitude, y = latitude, color = value)
  ) +
  scale_color_gradientn(colors = openColours()) +
  labs(color = quickText("NO2 (ug/m3)"))

## Another more local example

# just look at london
london_no2 <- 
  filter(aurn, str_detect(site, "^London"), species == "no2") %>% 
  arrange(value) %>%
  drop_na(value)

# define bounding box
left = min(london_no2$longitude)-.1
right = max(london_no2$longitude)+.1
down = min(london_no2$latitude)-.1
up = max(london_no2$latitude)+.1

# get basemap
london_basemap <-
  get_stamenmap(
    c(left, down, right, up),
    zoom = 11
  )

# ggplot
ggmap(london_basemap) +
  geom_point(
    data = london_no2,
    aes(x = longitude, y = latitude, fill = value),
    shape = 21,
    size = 5,
    color = "white",
    stroke = 2
  ) +
  ggrepel::geom_label_repel(
    data = london_no2,
    aes(x = longitude, y = latitude, label = site),
    size = 2.5
  ) +
  scale_fill_viridis_c(option = "inferno")

# Interactive Maps --------------------------------------------------------
# these all use leaflet

## markers
aurn %>%
  filter(species == "o3",
         !is.na(value)) %>%
  leaflet() %>%
  addTiles() %>%
  addMarkers(label = ~site,
             popup = ~paste("<b>Ozone:</b>", signif(value,4), "ug/m3"))

## circle markers (coloured)
# define map data
map_data <- 
  aurn %>%
  filter(species == "o3",
         !is.na(value)) %>%
  arrange(desc(value))

# create colour
pal <- colorNumeric(openColours(), map_data$value, reverse = TRUE)

# create map
map_data %>%
  leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    radius = 4,
    color = ~ pal(value),
    label = ~ site,
    popup = ~ paste("<b>Ozone:</b>", signif(value, 4), "ug/m3")
  ) %>%
  addLegend(pal = pal, values = map_data$value, 
            title = openairmaps::quickTextHTML("Ozone (O3)<br>(ug/m3)"))
