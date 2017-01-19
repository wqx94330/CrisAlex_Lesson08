###Lesson8
#19/01/2017
##Team CrisAlex

library(sp)
library(raster)
source("R/rm_extrvalues.R")

#download the files
Git_L="https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/raw/gh-pages/data/"

download.file(paste0(Git_L,"GewataB1.rda"), destfile = "data/GewataB1.rda", method = "auto")
download.file(paste0(Git_L,"GewataB2.rda"), destfile = "data/GewataB2.rda", method = "auto")
download.file(paste0(Git_L,"GewataB3.rda"), destfile = "data/GewataB3.rda", method = "auto")
download.file(paste0(Git_L,"GewataB4.rda"), destfile = "data/GewataB4.rda", method = "auto")
download.file(paste0(Git_L,"GewataB5.rda"), destfile = "data/GewataB5.rda", method = "auto")
download.file(paste0(Git_L,"GewataB7.rda"), destfile = "data/GewataB7.rda", method = "auto")
download.file(paste0(Git_L,"vcfGewata.rda"), destfile = "data/vcfGewata.rda", method = "auto")
download.file(paste0(Git_L,"trainingPoly.rda"), destfile = "data/trainingPoly.rda", method = "auto")

#load all the data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")

## Build a brick containing all data
alldata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(alldata) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")


#remove the values higher than 100 in VCF
alldata$VCF[alldata$VCF > 100] <- NA


#Put together all the rasters and visualize theirs correlations
alldata_ext <- brick(band1, band2, band3, band4, band5, band7, alldata$VCF)
pairs(alldata_ext)

## Extract all data to a data.frame
df <- as.data.frame(getValues(addLayer(band1, band2, band3, band4, band5, band7, alldata$VCF)))

#Create the model with the function Lineal regression all data related with VCF
model <- lm(df$VCF ~ band1 + band2 + band3 + band4 + band5 + band7, data = df)
summary(model)

#Create a new model
model2 <- lm(df$VCF ~ band3 + band4 + band5, data = df)
summary(model2)

#now with the generated model we predict the VCF
predict_VCF <- predict(alldata, model2)
summary(predict_VCF)
predict_VCF[predict_VCF < 0] <- NA

#Visualize the predicted VCf and the original VCF
par(mfrow = c(1,2))
plot(predict_VCF, main = "Predicted VCF")
plot(alldata$VCF, main = "Original VCF")
par(mfrow = c(1,1))

#Compare these two raster
diff_raster <- (alldata$VCF - predict_VCF)
plot(diff_raster, main = "Difference between Predicted VCF and Original VCF")


#Calculate the RMSE
predict_VCF_df <- as.data.frame(predict_VCF)
ori_VFC_df <- as.data.frame(alldata$VCF)
RMSE <- sqrt( mean( (ori_VFC_df - predict_VCF_df)^2 , na.rm = TRUE ) )
print(paste("The RMSE between the predicted VCF and the Original VCF is", RMSE))


#####
####Differences of three classes RMSE in predicted VCF and the original VCF

load("data/trainingPoly.rda")
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
trainingPoly@data

#Rasterize the training areas
classes <- rasterize(trainingPoly, alldata, field='Code')

#take the mean of the classes in each VCF(predicted and original)
predictVCF_zonal <- zonal(predict_VCF, classes, fun = mean, na.rm = TRUE)
originalVCF_zonal <- zonal(alldata$VCF, classes, fun = mean, na.rm = TRUE)

predictVCF_zonal_df <- as.data.frame(predictVCF_zonal)
originalVCF_zonal_df <- as.data.frame(originalVCF_zonal)

#Calculate the RMSE for each class
RMSE_class1 <- sqrt( mean( (originalVCF_zonal_df[1,2] - predictVCF_zonal_df[1,2])^2 , na.rm = TRUE ) )
RMSE_class2 <- sqrt( mean( (originalVCF_zonal_df[2,2] - predictVCF_zonal_df[2,2])^2 , na.rm = TRUE ) )
RMSE_class3 <- sqrt( mean( (originalVCF_zonal_df[3,2] - predictVCF_zonal_df[3,2])^2 , na.rm = TRUE ) )

print(paste("The Cropland RMSE between the predicted VCF and the Original VCF is", RMSE_class1))
print(paste("The Forest RMSE between the predicted VCF and the Original VCF is", RMSE_class2))
print(paste("The Wetland RMSE between the predicted VCF and the Original VCF is", RMSE_class3))


#Create a matrix with the final results of the classes RMSE
RMSE_matrix <- matrix(c(1,2,3, round(RMSE_class1, digits = 2), round(RMSE_class2, digits = 2), round(RMSE_class3, digits = 2), "cropland", "forest", "wetland"), nrow = 3, ncol = 3)
colnames(RMSE_matrix) <- c("Code", "RMSE", "Class")
print(RMSE_matrix)














































































