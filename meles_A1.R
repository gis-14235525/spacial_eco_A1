# GEOG71922 Assessment 1
# Species Distribution Modelling for Meles meles

setwd("./GEOG71922/14235525_A1")

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

