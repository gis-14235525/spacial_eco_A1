# GEOG71922 Assessment 1
# Species Distribution Modelling for Meles meles

setwd("./14235525_A1")

install.packages(c("dismo","glmnet","maxnet","terra","randomForest","sf","rnaturalearth","dplyr","precrec"))

# Load the libraries that we will need straight away.
library(terra)  
library(sf)      
library(dismo)   

# read in the badger data
meles <- read.csv("Melesmeles.csv")

# view the first six rows of the data
head(meles)

#check the class of the loaded .csv file.
class(meles)

# remove records with missing latitude and longtitude values
meles <- meles[!is.na(meles$Latitude),]
meles <- meles[!is.na(meles$Longitude),]

# remove unconfirmed records
meles <- meles[meles$Identification.verification.status != "Unconfirmed",]

#remove all points with uncertainty > 1000m
meles <- meles[meles$Coordinate.uncertainty_m <= 1000,]

# make spatial points layer

# create coordinates object
meles.latlong=data.frame(x=meles$Longitude,y=meles$Latitude)

# use coordinates object to create our spatial points object
meles.sp=st_as_sf(meles.latlong,coords=c("x","y"),crs="epsg:4326")

# load in the study area polygon
scot=st_read("scotSamp.shp")

# load in the land cover map and clip to the study area 
LCM=rast("LCMUK.tif")

# plus a little more buffer(because we will lose a small amount of data in the next step)
LCM=crop(LCM,st_buffer(scot, dist=1000))

# aggregate LCM raster 
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")

# project badger data to the same CRS as the LCM
meles.sp=st_transform(meles.sp,crs(LCM))

# crop the badger points to the study area
melesFin=meles.sp[scot,]

# mask the LCM to the study area boundary
LCM=crop(LCM,scot,mask=TRUE)

# inspect
plot(LCM)
plot(melesFin$geometry,add=T)

# Re-classifying the raster

# access levels of the raster by treating them as categorical data
LCM = as.factor(LCM)

levels(LCM)

# create a vector object called reclass
reclass = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values
RCmatrix = cbind(levels(LCM)[[1]], reclass)

RCmatrix = RCmatrix[,2:3]

# apply function to make sure new columns are numeric
RCmatrix = apply(RCmatrix, 2, FUN = as.numeric)

# use the classify() function to assign new values
broadleaf = classify(LCM, RCmatrix)

# inspect the new data visually
plot(broadleaf)
plot(melesFin$geometry, add = TRUE)


# characteristic scale

# use previous badger points, but transform back to EPSG:4326
meles.scale <- st_transform(meles.sp, 4326)

# convert to terra vector 
meles.scale <- vect(meles.scale)

#First set the exetent to something workable
studyExtent<-c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min x, max x, min y, max y

#now crop our points to this area
C<-crop(meles.scale,studyExtent)

# read in LSM raster data
LCM_scale = rast("LCMUK.tif")

# create an object containing the crs information for LCM inside the project
melesFin<-project(C,crs(LCM_scale))

# crop the land cover data to the extent of our points data (plus 5km to allow space for the buffers we will create in subsequent steps)
melesCoords<-crds(melesFin)

x.min <- min(melesCoords[,1]) - 5000
x.max <- max(melesCoords[,1]) + 5000
y.min <- min(melesCoords[,2]) - 5000
y.max <- max(melesCoords[,2]) + 5000

extent.new <- ext(x.min, x.max, y.min, y.max)

LCM_scale <- crop(LCM_scale$LCMUK_1, extent.new)

plot(LCM_scale)
plot(melesFin,add = TRUE)
