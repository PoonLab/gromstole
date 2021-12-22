library(sp)
library(rgdal)
library(OpenStreetMap)
library(raster)
library(here)

loc2 <- read.csv(here("data", "maps", "sampling-locations.csv"))
loc2 <- loc2[complete.cases(loc2),]
loc_merc <- projectMercator(loc2$lat, loc2$lon)

ontario = openmap(c(45.75,-84.25),c(42,-74), type="osm")

plot(ontario)
points(loc_merc, pch = 16, col = "darkorchid", cex = 1.25)
points(loc_merc, cex = 1.25)



