###Calculate species area and niche breadth for each species
###Used for ecology driving diversification of Saussurea---Xu zhang 2020/10/28

#Remove all data
rm(list=ls())

###request packages
library(ENMTools)
library(rgbif)
library(dplyr)
library(sf)
library(spData)
library(tmap)
library(purrr)
library(readr)
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid
library(plyr)
library(rgeos)
library(PBSmapping)
library(geosphere)
library(RapidPolygonLookup)

sp<- "xx" ####define the species name
dat<- occ_search(scientificName = sp) ####search locations in gbif 
sp_location <- dat$data[c("decimalLongitude","decimalLatitude")] %>% 
filter(decimalLongitude != "NA" | decimalLatitude != "NA") %>%
data.frame() ###select longtitude and latitude and filter thoes without location information 
sp_location<- unique(sp_location) ### delete duplications
x<- nrow(sp_location)
if (x>10)
{
outlier_location <- sapply(sp_location,function(X){which(X%in%boxplot.stats(X)$out)})outlier_location <- sapply(sp_location,function(X){which(X%in%boxplot.stats(X)$out)})
todel <- (sort(unique(unlist(outlier_location))))
presence.points <- sp_location[-todel,]
} else
{
presence.points<- sp_location
}
names(presence.points)<- c("Longitude", "Latitude")
write.csv(presence.points,"xx loc.csv") ##write to .csv file 

number_loc<- nrow(sp_location)
number_loc_used<-nrow(presence.points)
xy.sp = SpatialPoints(presence.points)
xy.sp_buffer<- gBuffer(xy.sp, width = 5)
area<- areaPolygon(xy.sp_buffer)###calculate area after buffer 5km as polygons
#area<- geosphere::areaPolygon(xy) ###calculate area using presence points
env_world <- raster::getData('worldclim', var='bio', res=2.5)
x1<-floor(min(presence.points$Longitude))-2
x2<-ceiling(max(presence.points$Longitude))+2
y1<-floor(min(presence.points$Latitude))-2
y2<-ceiling(max(presence.points$Latitude))+2
env <- crop(env_world, extent(x1,x2,y1,y2))
plot(env[[1]])
Saussurea <- enmtools.species()
Saussurea <- enmtools.species(species.name = sp, 
                              presence.points = presence.points)
Saussurea$range <- background.raster.buffer(Saussurea$presence.points, 50000, mask = env)
Saussurea$background.points <- background.points.buffer(points = Saussurea$presence.points,
                                                        radius = 20000, n = 1000, mask = env[[1]])
names(Saussurea)
check.species(Saussurea)
#interactive.plot.enmtools.species(Saussurea)
#raster.cor.matrix(env)
#raster.cor.plot(env)
env <- env[[c("x", "x", "x","x","x")]]
##GLM##
Saussurea.glm <- enmtools.glm(species = Saussurea, env = env,  test.prop = 0.2)
Saussurea.glm
breadth_glm<- raster.breadth(Saussurea.glm)
Saussurea.mx <- enmtools.maxent(Saussurea, env, test.prop = 0.2,rts.reps = 5)
Saussurea.mx
breadth_mx<- raster.breadth(Saussurea.mx)
write.csv(data.frame(sp, number_loc,number_loc_used,area, breadth_glm$B1,breadth_glm$B2,breadth_mx$B1,breadth_mx$B2),
          "xx NB.csv")


####extract 19 bios for each species###
library(raster)

rm(list=ls())
location<-read.table("./loc/Saussurea alaschanica loc.csv",header=T,sep=",",row.names=1)

worldbio<-getData('worldclim', var = "bio", res =2.5 )#获取Bioclimatic variable
bio<-extract(worldbio, location)
biomean<-as.data.frame(apply(bio,2,mean))