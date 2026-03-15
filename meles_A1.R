# GEOG71922 Assessment 1
# Species Distribution Modelling for Meles meles

setwd("./14235525_A1")

#install.packages(c("dismo","glmnet","maxnet","terra","randomForest","sf","precrec","mlr"))

# Load the libraries that we will need straight away.
library(terra)  
library(sf)      
library(dismo) 
library(maxnet)
library(glmnet)
library(randomForest)
library(precrec)

# read in the badger data
meles <- read.csv("Melesmeles.csv")

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
meles_study=meles.sp[scot,]

# mask the LCM to the study area boundary
LCM=crop(LCM,scot,mask=TRUE)

plot(LCM)
plot(meles_study, add = TRUE)

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
plot(meles_study, add = TRUE)


# characteristic scale

# convert study-area badger points to terra vector
melesScale <- vect(meles_study)

# read LCM raster
LCM_scale <- rast("LCMUK.tif")

# crop LCM to study area(test a 100-2000m buffer around the study area)
LCM_scale <- crop(LCM_scale, st_buffer(scot, dist = 2000))
LCM_scale <- LCM_scale$LCMUK_1


# Generating (pseudo-)absence points
# generate pseudo-absence points inside the study area
sampleMask <- crop(LCM_scale, scot, mask = TRUE)

set.seed(11)
back.xy <- spatSample(sampleMask, size=2000,as.points=TRUE, na.rm = TRUE)

# create presence and absence data frames
Abs <- data.frame(crds(back.xy), Pres = 0)
Pres <- data.frame(crds(melesScale), Pres = 1)


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
plot(crop(broadleaf, scot, mask = TRUE))
plot(melesScale, add = TRUE)
plot(back.xy, add = TRUE, col = "red", pch = 16, cex = 0.4)

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



# set environmental covariates

# 1.create woodland covariate at optimum radius
# create an vector object called reclass
reclass = c(rep(0,20),1,1)

# combine with the LCM categories into a matrix of old and new values.
RCmatrix=cbind(levels(LCM)[[1]],reclass)
RCmatrix=RCmatrix[,2:3]

# apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)

# Use the classify() function to asssign new values to LCM with our reclassification matrix
broadleaf=classify(LCM, RCmatrix)

# broadleaf focal layer at optimum radius
nPix=round((opt$radius)/res(LCM)[1])
nPix=(nPix*2)+1

# buiild weights matrix
weightsMatrix=matrix(1:nPix^2,nrow=nPix,ncol=nPix)

#get focal cell 
x=ceiling(ncol(weightsMatrix)/2)
y=ceiling(nrow(weightsMatrix)/2)

focalCell=weightsMatrix[x,y]

indFocal=which(weightsMatrix==focalCell,arr.ind = TRUE)

#compute distances
distances=list()

for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=TRUE)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
}

#add distance values to the weights matrix
weightsMatrix[]=unlist(distances)

#set cells outside search radius to NA
weightsMatrix[weightsMatrix>(opt$radius)]=NA

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrix[!is.na(weightsMatrix)]=1/length(weightsMatrix[!is.na(weightsMatrix)])

#sum neighbourhood values from all surrounding cells
lcm_wood_opt=focal(broadleaf,w=weightsMatrix,fun="sum")

# 2.create urban covariate at 2300m as practice
reclassUrban <- c(rep(0, 20), 1, 1)

RCmatrixUrban <- cbind(levels(LCM)[[1]], reclassUrban)
RCmatrixUrban <- RCmatrixUrban[, 2:3]
RCmatrixUrban <- apply(RCmatrixUrban, 2, FUN = as.numeric)

urban <- classify(LCM, RCmatrixUrban)

nPixUrban <- round(2300 / res(LCM)[1])
nPixUrban <- (nPixUrban * 2) + 1

weightsMatrixUrban <- matrix(1:nPixUrban^2, nrow = nPixUrban, ncol = nPixUrban)

x <- ceiling(ncol(weightsMatrixUrban) / 2)
y <- ceiling(nrow(weightsMatrixUrban) / 2)

focalCell <- weightsMatrixUrban[x, y]
indFocal <- which(weightsMatrixUrban == focalCell, arr.ind = TRUE)

distancesUrban <- list()

for(i in 1:nPixUrban^2){
  ind.i <- which(weightsMatrixUrban == i, arr.ind = TRUE)
  diffX <- abs(ind.i[1,1] - indFocal[1,1]) * res(LCM)[1]
  diffY <- abs(ind.i[1,2] - indFocal[1,2]) * res(LCM)[1]
  dist.i <- sqrt(diffX^2 + diffY^2)
  distancesUrban[[i]] <- dist.i
}

weightsMatrixUrban[] <- unlist(distancesUrban)
weightsMatrixUrban[weightsMatrixUrban > 2300] <- NA
weightsMatrixUrban[!is.na(weightsMatrixUrban)] <- 1 / length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

lcm_urban_2300 <- focal(urban, w = weightsMatrixUrban, fun = "sum")

# elevation
demScot <- rast("demScotland (1).tif")
demScot <- terra::resample(demScot, lcm_wood_opt)

# stack the covariate layers together
allEnv=c(lcm_wood_opt,lcm_urban_2300,demScot)
names(allEnv)=c("broadleaf","urban","elev")


# extract covariates to points

# create background points
set.seed(11)

# sample background - one point for every cell (9775)
back = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back=back[!is.na(back$broadleaf),]
back=st_as_sf(back,crs="EPSG:27700")
# get environmental covariates at presence locations
eP=terra::extract(allEnv,meles_study)

# bind together the presence data using cbind() which binds together objects by column (i.e. with different columns but the same number of rows)
Pres.cov=st_as_sf(cbind(eP,meles_study))
Pres.cov$Pres=1

# remove the first column which is just an ID field.
Pres.cov=Pres.cov[,-1]

# get coordinates for spatial cross-validation later
coordsPres=st_coordinates(Pres.cov)

# drop geometry column using st_drop_geometry()
Back.cov=st_as_sf(data.frame(back,Pres=0))

# get coordinates of background points for cross validation later
coordsBack=st_coordinates(back)

# combine
coords=data.frame(rbind(coordsPres,coordsBack))

# assign coumn names
colnames(coords)=c("x","y")

# combine pres and background
all.cov=rbind(Pres.cov,Back.cov)

# add coordinates
all.cov=cbind(all.cov,coords)

# remove any NAs and geometry
all.cov=na.omit(all.cov)
all.cov <- st_drop_geometry(all.cov)


# GLM model
# specify the model
glm_badger=glm(Pres~broadleaf+urban+elev,binomial(link='logit'),
                data=all.cov)

# predict and inspect the output
prGLM=predict(allEnv,glm_badger,type="response")

# plot
plot(prGLM)

# build new data frame based on mean of elev and urban but varying values for broadleaf. 
glmNew=data.frame(broadleaf=seq(0,max(all.cov$broadleaf),length=1000),
                  elev=mean(all.cov$elev),
                  urban=mean(all.cov$urban))


# use type = "response" for probability-scale predictions and chose to return the standard error of the prediction (se.fit=TRUE)   
preds = predict(glm_badger, newdata = glmNew, type = "response", se.fit = TRUE)
glmNew$fit = preds$fit
glmNew$se = preds$se.fit

head(glmNew)


# maxnet evaluation(5folds)

#set number of folds to use
folds=5

# partition presence and background data and assign to folds using the kfold() function.
Pres.cov=all.cov[all.cov$Pres==1,]
Back.cov=all.cov[all.cov$Pres==0,]

kfold_pres = kfold(Pres.cov, folds)
kfold_back = kfold(Back.cov, folds)

eMax=list()

for (i in 1:folds) {
  train = Pres.cov[kfold_pres!= i,]
  test = Pres.cov[kfold_pres == i,]
  backTrain=Back.cov[kfold_back!=i,]
  backTest=Back.cov[kfold_back==i,]
  dataTrain=rbind(train,backTrain)
  dataTest=rbind(test,backTest)
  maxnetMod=maxnet(dataTrain$Pres, dataTrain[, c("broadleaf", "urban", "elev")],
                   maxnet.formula(dataTrain$Pres, dataTrain[, vars], classes = "lq"))
  
  pred_p <- predict(maxnetMod, dataTest[dataTest$Pres == 1, vars], type = "cloglog")
  pred_a <- predict(maxnetMod, dataTest[dataTest$Pres == 0, vars], type = "cloglog")
  
  eMax[[i]] <- evaluate(p = pred_p, a = pred_a)
}


#print the result
aucMax = sapply(eMax, function(x){slot(x, 'auc')} )
print(mean(aucMax))

Opt_Max = sapply(eMax, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptMax = mean(Opt_Max)
print(Mean_OptMax)


# maxnet model
maxnet_badger <- maxnet(
  st_drop_geometry(all.cov)$Pres,
  st_drop_geometry(all.cov)[, c("broadleaf", "urban", "elev")],
  maxnet.formula(
    st_drop_geometry(all.cov)$Pres,
    st_drop_geometry(all.cov)[, c("broadleaf", "urban", "elev")],
    classes = "lq"
  )
)

plot(maxnet_badger, type = "cloglog")

# maxnet prediction map
prMax <- terra::predict(allEnv, maxnet_badger, clamp = FALSE, type = "cloglog", na.rm = TRUE)

plot(prMax, main = "Maxnet cloglog prediction")
plot(meles_study, add = TRUE)

# maxnet presence or absence map
maxPA <- prMax > Mean_OptMax
par(mfrow = c(1,2))
plot(prMax, main = "Maxnet probability")
plot(maxPA, main = "Maxnet presence/absence")
plot(meles_study, add = TRUE)
par(mfrow = c(1,1))

