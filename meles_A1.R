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

# use all cleaned badger points
meles.scale <- st_as_sf(meles.latlong, coords = c("x","y"), crs = "epsg:4326")

# convert to terra vector 
meles.scale <- vect(meles.scale)

# read in LCM raster data
LCM_scale <- rast("LCMUK.tif")

# project all cleaned badger points to the same CRS as the LCM
melesFin <- project(meles.scale, crs(LCM_scale))

# crop the land cover data to the extent of all badger points (plus 5km to allow space for the buffers we will create in subsequent steps)
melesCoords<-crds(melesFin)

x.min <- min(melesCoords[,1]) - 5000
x.max <- max(melesCoords[,1]) + 5000
y.min <- min(melesCoords[,2]) - 5000
y.max <- max(melesCoords[,2]) + 5000

extent.new <- ext(x.min, x.max, y.min, y.max)

LCM_scale <- crop(LCM_scale$LCMUK_1, extent.new)


# Generating (pseudo-)absence points
set.seed(11)
back.xy <- spatSample(LCM_scale, size=1000,as.points=TRUE) 

# create presence and absence data frames
Abs <- data.frame(crds(back.xy), Pres = 0)
Pres <- data.frame(crds(melesFin), Pres = 1)


# bind the two data frames by row
melesData <- rbind(Pres, Abs)

# convert to sf for buffer analysis
melesSF <- st_as_sf(melesData, coords = c("x","y"), crs = "EPSG:27700")


# reclassify original LCM to broadleaf woodland
LCM_scale <- as.factor(LCM_scale)

#create an vector object
reclass <- c(0, 1, rep(0, nrow(levels(LCM_scale)[[1]]) - 2))

# combine with the LCM categories into a matrix of old and new values
RCmatrix <- cbind(levels(LCM_scale)[[1]], reclass)
RCmatrix <- RCmatrix[,2:3]
RCmatrix <- apply(RCmatrix, 2, FUN = as.numeric)

broadleaf <- classify(LCM_scale, RCmatrix)

#Function for automating whole dataset

landBuffer <- function(speciesData, r){         
  
  #buffer each point
  melesBuffer <- st_buffer(speciesData, dist=r)                     
  
  #crop the woodland layer to the buffer extent
  bufferlandcover <- crop(broadleaf, melesBuffer)              
  
  # now extract the raster values (which should all be 1 for woodland and 0 for everything else) within each buffer and sum to get number of woodland cells inside the buffers.
  masklandcover <- extract(bufferlandcover, melesBuffer,fun="sum")   
  
  #get woodland area (625 is the area in metres of each cell of our 25m raster)
  landcoverArea <- masklandcover$LCMUK_1*625  
  
  # convert to precentage cover (we use the st_area() function from the sf package to get the area of our buffer) but convert to a numeric object (because sf applies units i.e. metres which then cant be entered into numeric calculations)
  percentcover <- landcoverArea/as.numeric(st_area(melesBuffer))*100 
  
  # return the result
  return(percentcover)                                       
}

# loop
radii<-seq(100,2000,by=100)

resList=list()

for(i in radii){
  res.i=landBuffer(speciesData=melesSF,r=i)
  res.i
  resList[[i/100]]=res.i
}

#collect all results together
resFin=do.call("cbind",resList)

#convert to data frame
glmData=data.frame(resFin)

#assign more intuitive column names
colnames(glmData)=paste0("radius",radii)

#add in the presences data
glmData$Pres<-melesData$Pres


#make an empty data frame to then add a sequence of results from glms (general linear models) that test the relationship between broadleaf cover and red squirrel presence for the different buffer sizes

#init empty data frame
glmRes=data.frame(radius=NA,loglikelihood=NA)

#for loop to iterate over radius values and run a general linear model with glm()
for(i in radii){
  n.i=paste0("Pres~","radius",i,sep ="")
  glm.i=glm(formula(n.i),family = "binomial",data = glmData)
  ll.i=as.numeric(logLik(glm.i))
  glmRes=rbind(glmRes,c(i,ll.i))
}

#remove the NAs in the first row
glmRes=glmRes[!is.na(glmRes),]

#use the which.max function to subset the dataframe to the just the row containing the max log likelihood value
opt<-glmRes[which.max(glmRes$loglikelihood),]



