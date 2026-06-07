library(raster)
library(rgdal)
files <- list.files(pattern='.tif$', full.names=TRUE)
s <- stack(files)
points <- read.csv(file = "coords.csv", header = TRUE)

# if your points are lonlat
points <- SpatialPoints(points[, c('longitude', 'latitude')], 
                        proj4string=CRS('+proj=longlat +datum=WGS84'))
                        df <- extract(s, points)
                        write.csv(df, 'file.csv')
                  