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

# 1.read in the badger data
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


# 2.load in the study area polygon
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


# 3.Re-classifying the raster
# access levels of the raster by treating them as categorical data
LCM = as.factor(LCM)
levels(LCM)
# broadleaf raster
reclass = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values
RCmatrix = cbind(levels(LCM)[[1]], reclass)

RCmatrix = RCmatrix[,2:3]

# apply function to make sure new columns are numeric
RCmatrix = apply(RCmatrix, 2, FUN = as.numeric)

# use the classify() function to assign new values
broadleaf = classify(LCM, RCmatrix)

# 4. characteristic scale
# optimum radius (2000 m) previously identified from badger characteristic-scale analysis; fixed here to avoid rerunning the computationally intensive optimisation loop
opt_radius <- 2000


# 5.set environmental covariates

# 5.1.create woodland covariate at optimum radius
nPix=round((opt_radius)/res(broadleaf)[1])
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
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(broadleaf)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(broadleaf)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
}

#add distance values to the weights matrix
weightsMatrix[]=unlist(distances)

#set cells outside search radius to NA
weightsMatrix[weightsMatrix>(opt_radius)]=NA

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrix[!is.na(weightsMatrix)]=1/length(weightsMatrix[!is.na(weightsMatrix)])

#sum neighbourhood values from all surrounding cells
lcm_wood_opt=focal(broadleaf,w=weightsMatrix,fun="sum")

# 5.2.create urban covariate at 2300m as practice
reclassUrban <- c(rep(0, 20), 1, 1)

RCmatrixUrban <- cbind(levels(LCM)[[1]], reclassUrban)
CmatrixUrban <- RCmatrixUrban[, 2:3]
RCmatrixUrban <- apply(RCmatrixUrban, 2, FUN = as.numeric)

urban <- classify(LCM, RCmatrixUrban)

nPixUrban <- round(2300/ res(urban)[1])
nPixUrban <- (nPixUrban * 2) + 1

weightsMatrixUrban <- matrix(1:nPixUrban^2, nrow = nPixUrban, ncol = nPixUrban)

x <- ceiling(ncol(weightsMatrixUrban) / 2)
y <- ceiling(nrow(weightsMatrixUrban) / 2)

focalCell <- weightsMatrixUrban[x, y]
indFocal <- which(weightsMatrixUrban == focalCell, arr.ind = TRUE)

distancesUrban <- list()

for(i in 1:nPixUrban^2){
  ind.i <- which(weightsMatrixUrban == i, arr.ind = TRUE)
  diffX <- abs(ind.i[1,1] - indFocal[1,1]) * res(urban)[1]
  diffY <- abs(ind.i[1,2] - indFocal[1,2]) * res(urban)[1]
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


# 6.extract covariates to points

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
Pres.cov=Pres.cov[, -1]

# create background data
Back.cov=st_as_sf(data.frame(back, Pres = 0))

# combine presence and background
all.cov=rbind(Pres.cov, Back.cov)

all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)


# 7.GLM model
# specify the model
glm_badger=glm(Pres~broadleaf+urban+elev,binomial(link='logit'),
                data=all.cov)

# predict and inspect the output
prGLM=predict(allEnv,glm_badger,type="response")

folds=5

# partition presence and background data and assign to folds using the kfold() function.
Pres.glm=all.cov[all.cov$Pres==1,]
Back.glm=all.cov[all.cov$Pres==0,]

kfold_pres_glm = kfold(Pres.cov, folds)
kfold_back_glm = kfold(Back.cov, folds)

eGLM=list()

vars <- c("broadleaf", "urban", "elev")

for(i in 1:folds){
  train=Pres.glm[kfold_pres_glm != i, ]
  test=Pres.glm[kfold_pres_glm == i, ]
  backTrain=Back.glm[kfold_back_glm != i, ]
  backTest=Back.glm[kfold_back_glm == i, ]
  dataTrain=rbind(train, backTrain)
  dataTest=rbind(test, backTest)
  glm_i=glm(Pres ~ broadleaf + urban + elev,
               binomial(link = "logit"),
               data = dataTrain)
  
  pred_p=predict(glm_i, newdata = dataTest[dataTest$Pres == 1, ], type = "response")
  pred_a=predict(glm_i, newdata = dataTest[dataTest$Pres == 0, ], type = "response")
  
  eGLM[[i]]=evaluate(p = pred_p, a = pred_a)
}


#print the result
aucGLM=sapply(eGLM, function(x){slot(x, 'auc')} )
Mean_AucGLM=mean(aucGLM)

Opt_GLM=sapply(eGLM, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptGLM=mean(Opt_GLM)

glmPA=prGLM > Mean_OptGLM

# glm ROC plots
par(mfrow = c(2,3))
for(i in 1:folds){
  plot(eRF[[i]], "ROC", col = "red", lwd = 2)
}
par(mfrow = c(1,1))


# 8.maxnet model
Pres.max=all.cov[all.cov$Pres == 1, ]
Back.max=all.cov[all.cov$Pres == 0, ]

kfold_pres=kfold(Pres.max, folds)
kfold_back=kfold(Back.max, folds)

eMax=list()

for(i in 1:folds){
  train=Pres.max[kfold_pres != i, ]
  test=Pres.max[kfold_pres == i, ]
  
  backTrain=Back.max[kfold_back != i, ]
  backTest=Back.max[kfold_back == i, ]
  
  dataTrain=rbind(train, backTrain)
  dataTest=rbind(test, backTest)
  
  maxnetMod=maxnet(
    dataTrain$Pres,
    dataTrain[, vars],
    maxnet.formula(dataTrain$Pres, dataTrain[, vars], classes = "lq")
  )
  
  eMax[[i]] = evaluate(
    p = dataTest[dataTest$Pres == 1, vars],
    a = dataTest[dataTest$Pres == 0, vars],
    model = maxnetMod
  )
}

aucMax=sapply(eMax, function(x){slot(x, "auc")})
Mean_AucMax=mean(aucMax)

Opt_Max=sapply(eMax, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptMax=mean(Opt_Max)

maxnet_badger=maxnet(
  all.cov$Pres,
  all.cov[, vars],
  maxnet.formula(all.cov$Pres, all.cov[, vars], classes = "lq")
)

prMax=predict(allEnv, maxnet_badger, clamp = FALSE, type = "cloglog", na.rm = TRUE)
maxPA=prMax > Mean_OptMax



# 9.random forest model

Pres.rf = all.cov[all.cov$Pres == 1, ]
Back.rf = all.cov[all.cov$Pres == 0, ]

kfold_pres_rf = kfold(Pres.rf, folds)
kfold_back_rf = kfold(Back.rf, folds)

eRF = list()

for(i in 1:folds){
  train = Pres.rf[kfold_pres_rf != i, ]
  test = Pres.rf[kfold_pres_rf == i, ]
  
  backTrain = Back.rf[kfold_back_rf != i, ]
  backTest = Back.rf[kfold_back_rf == i, ]
  
  dataTrain = rbind(train, backTrain)
  dataTest = rbind(test, backTest)
  
  dataTrain$Pres = as.factor(dataTrain$Pres)
  
  rfMod = randomForest(Pres ~ broadleaf + urban + elev,
                       data = dataTrain,
                       ntree = 500)
  
  pred_p = predict(rfMod,
                   newdata = dataTest[dataTest$Pres == 1, ],
                   type = "prob")[,2]
  
  pred_a = predict(rfMod,
                   newdata = dataTest[dataTest$Pres == 0, ],
                   type = "prob")[,2]
  
  eRF[[i]] = evaluate(p = pred_p, a = pred_a)
}

aucRF = sapply(eRF, function(x){slot(x, "auc")})
Mean_AucRF = mean(aucRF)

Opt_RF = sapply(eRF, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptRF = mean(Opt_RF)


all.cov.rf = all.cov
all.cov.rf$Pres = as.factor(all.cov.rf$Pres)

rf_badger = randomForest(Pres ~ broadleaf + urban + elev,
                         data = all.cov.rf,
                         ntree = 500)

# random forest prediction
prRF = terra::predict(allEnv, rf_badger, type = "prob", index = 2)

rfPA = prRF > Mean_OptRF

# rf ROC plots
par(mfrow = c(2,3))
for(i in 1:folds){
  plot(eRF[[i]], "ROC", col = "red", lwd = 2)
}
par(mfrow = c(1,1))

# final outputs
plot(prGLM, main = "GLM probability")
plot(prMax, main = "Maxnet probability")
plot(prRF, main = "Random Forest probability")
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

# 1.read in the badger data
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


# 2.load in the study area polygon
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


# 3.Re-classifying the raster
# access levels of the raster by treating them as categorical data
LCM = as.factor(LCM)
levels(LCM)
# broadleaf raster
reclass = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values
RCmatrix = cbind(levels(LCM)[[1]], reclass)

RCmatrix = RCmatrix[,2:3]

# apply function to make sure new columns are numeric
RCmatrix = apply(RCmatrix, 2, FUN = as.numeric)

# use the classify() function to assign new values
broadleaf = classify(LCM, RCmatrix)

# 4. characteristic scale
# optimum radius (2000 m) previously identified from badger characteristic-scale analysis; fixed here to avoid rerunning the computationally intensive optimisation loop
opt_radius <- 2000


# 5.set environmental covariates

# 5.1.create woodland covariate at optimum radius
nPix=round((opt_radius)/res(broadleaf)[1])
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
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(broadleaf)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(broadleaf)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
}

#add distance values to the weights matrix
weightsMatrix[]=unlist(distances)

#set cells outside search radius to NA
weightsMatrix[weightsMatrix>(opt_radius)]=NA

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrix[!is.na(weightsMatrix)]=1/length(weightsMatrix[!is.na(weightsMatrix)])

#sum neighbourhood values from all surrounding cells
lcm_wood_opt=focal(broadleaf,w=weightsMatrix,fun="sum")

# 5.2.create urban covariate at 2300m as practice
reclassUrban <- c(rep(0, 20), 1, 1)

RCmatrixUrban <- cbind(levels(LCM)[[1]], reclassUrban)
CmatrixUrban <- RCmatrixUrban[, 2:3]
RCmatrixUrban <- apply(RCmatrixUrban, 2, FUN = as.numeric)

urban <- classify(LCM, RCmatrixUrban)

nPixUrban <- round(2300/ res(urban)[1])
nPixUrban <- (nPixUrban * 2) + 1

weightsMatrixUrban <- matrix(1:nPixUrban^2, nrow = nPixUrban, ncol = nPixUrban)

x <- ceiling(ncol(weightsMatrixUrban) / 2)
y <- ceiling(nrow(weightsMatrixUrban) / 2)

focalCell <- weightsMatrixUrban[x, y]
indFocal <- which(weightsMatrixUrban == focalCell, arr.ind = TRUE)

distancesUrban <- list()

for(i in 1:nPixUrban^2){
  ind.i <- which(weightsMatrixUrban == i, arr.ind = TRUE)
  diffX <- abs(ind.i[1,1] - indFocal[1,1]) * res(urban)[1]
  diffY <- abs(ind.i[1,2] - indFocal[1,2]) * res(urban)[1]
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


# 6.extract covariates to points

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
Pres.cov=Pres.cov[, -1]

# create background data
Back.cov=st_as_sf(data.frame(back, Pres = 0))

# combine presence and background
all.cov=rbind(Pres.cov, Back.cov)

all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)


# 7.GLM model
# specify the model
glm_badger=glm(Pres~broadleaf+urban+elev,binomial(link='logit'),
               data=all.cov)

# predict and inspect the output
prGLM=predict(allEnv,glm_badger,type="response")

folds=5

# partition presence and background data and assign to folds using the kfold() function.
Pres.glm=all.cov[all.cov$Pres==1,]
Back.glm=all.cov[all.cov$Pres==0,]

kfold_pres_glm = kfold(Pres.cov, folds)
kfold_back_glm = kfold(Back.cov, folds)

eGLM=list()

vars <- c("broadleaf", "urban", "elev")

for(i in 1:folds){
  train=Pres.glm[kfold_pres_glm != i, ]
  test=Pres.glm[kfold_pres_glm == i, ]
  backTrain=Back.glm[kfold_back_glm != i, ]
  backTest=Back.glm[kfold_back_glm == i, ]
  dataTrain=rbind(train, backTrain)
  dataTest=rbind(test, backTest)
  glm_i=glm(Pres ~ broadleaf + urban + elev,
            binomial(link = "logit"),
            data = dataTrain)
  
  pred_p=predict(glm_i, newdata = dataTest[dataTest$Pres == 1, ], type = "response")
  pred_a=predict(glm_i, newdata = dataTest[dataTest$Pres == 0, ], type = "response")
  
  eGLM[[i]]=evaluate(p = pred_p, a = pred_a)
}


#print the result
aucGLM=sapply(eGLM, function(x){slot(x, 'auc')} )
Mean_AucGLM=mean(aucGLM)

Opt_GLM=sapply(eGLM, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptGLM=mean(Opt_GLM)

glmPA=prGLM > Mean_OptGLM

# ROC plots
par(mfrow = c(2,3))
for(i in 1:folds){
  plot(eRF[[i]], "ROC", col = "red", lwd = 2)
}
par(mfrow = c(1,1))


# 8.maxnet model
Pres.max=all.cov[all.cov$Pres == 1, ]
Back.max=all.cov[all.cov$Pres == 0, ]

kfold_pres=kfold(Pres.max, folds)
kfold_back=kfold(Back.max, folds)

eMax=list()

for(i in 1:folds){
  train=Pres.max[kfold_pres != i, ]
  test=Pres.max[kfold_pres == i, ]
  
  backTrain=Back.max[kfold_back != i, ]
  backTest=Back.max[kfold_back == i, ]
  
  dataTrain=rbind(train, backTrain)
  dataTest=rbind(test, backTest)
  
  maxnetMod=maxnet(
    dataTrain$Pres,
    dataTrain[, vars],
    maxnet.formula(dataTrain$Pres, dataTrain[, vars], classes = "lq")
  )
  
  eMax[[i]] = evaluate(
    p = dataTest[dataTest$Pres == 1, vars],
    a = dataTest[dataTest$Pres == 0, vars],
    model = maxnetMod
  )
}

aucMax=sapply(eMax, function(x){slot(x, "auc")})
Mean_AucMax=mean(aucMax)

Opt_Max=sapply(eMax, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptMax=mean(Opt_Max)

maxnet_badger=maxnet(
  all.cov$Pres,
  all.cov[, vars],
  maxnet.formula(all.cov$Pres, all.cov[, vars], classes = "lq")
)

prMax=predict(allEnv, maxnet_badger, clamp = FALSE, type = "cloglog", na.rm = TRUE)
maxPA=prMax > Mean_OptMax

# ROC plots
par(mfrow = c(2,3))
for(i in 1:folds){
  plot(eRF[[i]], "ROC", col = "red", lwd = 2)
}
par(mfrow = c(1,1))

# 9.random forest model

Pres.rf = all.cov[all.cov$Pres == 1, ]
Back.rf = all.cov[all.cov$Pres == 0, ]

kfold_pres_rf = kfold(Pres.rf, folds)
kfold_back_rf = kfold(Back.rf, folds)

eRF = list()

for(i in 1:folds){
  train = Pres.rf[kfold_pres_rf != i, ]
  test = Pres.rf[kfold_pres_rf == i, ]
  
  backTrain = Back.rf[kfold_back_rf != i, ]
  backTest = Back.rf[kfold_back_rf == i, ]
  
  dataTrain = rbind(train, backTrain)
  dataTest = rbind(test, backTest)
  
  dataTrain$Pres = as.factor(dataTrain$Pres)
  
  rfMod = randomForest(Pres ~ broadleaf + urban + elev,
                       data = dataTrain,
                       ntree = 500)
  
  pred_p = predict(rfMod,
                   newdata = dataTest[dataTest$Pres == 1, ],
                   type = "prob")[,2]
  
  pred_a = predict(rfMod,
                   newdata = dataTest[dataTest$Pres == 0, ],
                   type = "prob")[,2]
  
  eRF[[i]] = evaluate(p = pred_p, a = pred_a)
}

aucRF = sapply(eRF, function(x){slot(x, "auc")})
Mean_AucRF = mean(aucRF)

Opt_RF = sapply(eRF, function(x){ x@t[which.max(x@TPR + x@TNR)] })
Mean_OptRF = mean(Opt_RF)


all.cov.rf = all.cov
all.cov.rf$Pres = as.factor(all.cov.rf$Pres)

rf_badger = randomForest(Pres ~ broadleaf + urban + elev,
                         data = all.cov.rf,
                         ntree = 500)

# random forest prediction
prRF = terra::predict(allEnv, rf_badger, type = "prob", index = 2)

rfPA = prRF > Mean_OptRF

# ROC plots
par(mfrow = c(2,3))
for(i in 1:folds){
  plot(eRF[[i]], "ROC", col = "red", lwd = 2)
}
par(mfrow = c(1,1))

# final outputs
plot(prGLM, main = "GLM probability")
plot(prMax, main = "Maxnet probability")
plot(prRF, main = "Random Forest probability")
plot(rfPA, main = "Random Forest presence/absence")
Mean_AucGLM
Mean_OptGLM
Mean_AucMax
Mean_AucRF
Mean_OptRF
