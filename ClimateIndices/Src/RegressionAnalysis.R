library(ncdf4)
library(chron)
library(tidyverse)
library(kohonen) # fitting
library(aweSOM) # plotting
library(SOMbrero) # plotting
library(paletteer) #colors
library(PNWColors) #more colors
library(here) #navigating folders
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)       #basic mapping functions and some data

library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(PBSmapping)
set.seed(1234)
#here::i_am("som/sst_soms_spring.r")

#### Importing data ####

# start by loading NE Pacific SST/SLP - this may need to be moved into project if not locally stored
nc <- nc_open(here("som/copernicus_jul10.nc"))
#nc.wind <- nc_open(here("som/copernicus_jul10_wind.nc"))
#climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

#### SST Upwelling Analysis Index ####

# get times and 
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$spring_year <- dates$year
dates$spring_year[which(dates$month %in% 7:12)] <- NA
dates$spring_year[which(dates$month %in% 1:3)] <- NA # incremenent 

var_name = "sst"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
X <- X[,-to_drop]# drop cells with NAs 
spring_lat <- lat[-to_drop]
spring_lon <- lon[-to_drop]

# Add block for spring average calculations
spring_years <- unique(dates$spring_year)
spring_years <- sort(spring_years[-which(is.na(spring_years))])
X_spring <- matrix(NA, length(spring_years), ncol(X))
for(i in 1:length(spring_years)) {
  X_spring[i,] = colMeans(X[which(dates$spring_year == spring_years[i]),])
}

# remove annual mean for each year
for(i in 1:nrow(X_spring)) {
  X_spring[i,] = X_spring[i,] - mean(X_spring[i,], na.rm=T)
}
X_spring <- scale(X_spring, center=TRUE, scale=FALSE)


sst_anomaly <- matrix(NA, length(spring_years), ncol(X))
mean_sst <- apply(X_spring,2,mean)
#sd_sst <- apply(X_spring,2,sd)
for(i in 1:ncol(X)) {
  sst_anomaly[,i] = (X_spring[,i]-mean_sst[i])
}
#assigning names to time allows you to divide data by ERA

dimnames(sst_anomaly) <- list(as.character(spring_years),paste("N", spring_lat, "E", spring_lon, sep=""))
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))

SST1 <- sst_anomaly[rownames(sst_anomaly) %in% 1967:1988,] 
SST2 <- sst_anomaly[rownames(sst_anomaly) %in% 1989:2012,]
SST3 <- sst_anomaly[rownames(sst_anomaly) %in% 2013:2023,]


#### Climate SST####
PDO.FMA <- climate_dat%>% filter(season=='Spring'&region=='GoA')%>% 
  dplyr::select(Year_lag, seasonal_PDO) 
rownames(PDO.FMA)<- PDO.FMA$Year_lag

PDO1 <- PDO.FMA[rownames(PDO.FMA) %in% 1967:1988, 'seasonal_PDO'] 
PDO2 <- PDO.FMA[rownames(PDO.FMA) %in% 1989:2012, 'seasonal_PDO']
PDO3 <- PDO.FMA[rownames(PDO.FMA) %in% 2013:2022, 'seasonal_PDO']
PDO4 <- PDO.FMA[rownames(PDO.FMA) %in% 1967:2022, 'seasonal_PDO']

ONI.FMA <- climate_dat%>% filter(season=='Spring'&region=='GoA')%>% 
  dplyr::select(Year_lag, seasonal_ONI) 
rownames(ONI.FMA)<- ONI.FMA$Year_lag

ONI1 <- ONI.FMA[rownames(ONI.FMA) %in% 1967:1988, 'seasonal_ONI'] 
ONI2 <- ONI.FMA[rownames(ONI.FMA) %in% 1989:2012, 'seasonal_ONI']
ONI3 <- ONI.FMA[rownames(ONI.FMA) %in% 2013:2022, 'seasonal_ONI']
ONI4 <- ONI.FMA[rownames(ONI.FMA) %in% 1967:2022, 'seasonal_ONI']

NPGO.FMA <- climate_dat%>% filter(season=='Spring'&region=='GoA')%>% 
  dplyr::select(Year_lag, seasonal_NPGO) 
rownames(NPGO.FMA)<- NPGO.FMA$Year_lag

NPGO1 <- NPGO.FMA[rownames(NPGO.FMA) %in% 1967:1988, 'seasonal_NPGO'] 
NPGO2 <- NPGO.FMA[rownames(NPGO.FMA) %in% 1989:2012, 'seasonal_NPGO']
NPGO3 <- NPGO.FMA[rownames(NPGO.FMA) %in% 2013:2022, 'seasonal_NPGO']
NPGO4 <- NPGO.FMA[rownames(NPGO.FMA) %in% 1967:2022, 'seasonal_NPGO']

NPH.FMA <- climate_dat%>% filter(season=='Spring'&region=='GoA')%>% 
  dplyr::select(Year_lag, seasonal_NPH) 
rownames(NPH.FMA)<- NPH.FMA$Year_lag

NPH1 <- NPH.FMA[rownames(NPH.FMA) %in% 1967:1988, 'seasonal_NPH'] 
NPH2 <- NPH.FMA[rownames(NPH.FMA) %in% 1989:2012, 'seasonal_NPH']
NPH3 <- NPH.FMA[rownames(NPH.FMA) %in% 2013:2022, 'seasonal_NPH']
NPH4 <- NPH.FMA[rownames(NPH.FMA) %in% 1967:2022, 'seasonal_NPH']


# make objects to catch results
# make objects to catch results
PDO.regr1 <- PDO.regr2<- PDO.regr3<- NA
ONI.regr1 <- ONI.regr2<- ONI.regr3<- NA
NPH.regr1 <- NPH.regr2<- NPH.regr3<- NA
NPGO.regr1 <- NPGO.regr2<- NPGO.regr3<- NA

# now loop through each cell

for(i in 1:ncol(SST1)){
  #  i <- 1
  mod <- lm(PDO1~SST1[,i])
  PDO.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(PDO2~SST2[,i])
  PDO.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(PDO3~SST3[,i])
  PDO.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(ONI1~SST1[,i])
  ONI.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(ONI2~SST2[,i])
  ONI.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(ONI3~SST3[,i])
  ONI.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(NPGO1~SST1[,i])
  NPGO.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(NPGO2~SST2[,i])
  NPGO.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(NPGO3~SST3[,i])
  NPGO.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(NPH1~SST1[,i])
  NPH.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(NPH2~SST2[,i])
  NPH.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(NPH3~SST3[,i])
  NPH.regr3[i] <- summary(mod)$coef[2,1]
}
# calculate era differences for each cell

PDO.diff3 <- PDO.regr3 - PDO.regr2
PDO.diff2 <- PDO.regr3 - PDO.regr1
PDO.diff1 <- PDO.regr2 - PDO.regr1

NPGO.diff3 <- NPGO.regr3 - NPGO.regr2
NPGO.diff2 <- NPGO.regr3 - NPGO.regr1
NPGO.diff1 <- NPGO.regr2 - NPGO.regr1

NPH.diff3 <- NPH.regr3 - NPH.regr2
NPH.diff2 <- NPH.regr3 - NPH.regr1
NPH.diff1 <- NPH.regr2 - NPH.regr1

ONI.diff3 <- ONI.regr3 - ONI.regr2
ONI.diff2 <- ONI.regr3 - ONI.regr1
ONI.diff1 <- ONI.regr2 - ONI.regr1

X_PDO_SST<- as.data.frame(PDO.regr1)%>%
  cbind(as.data.frame(PDO.regr2))%>%
  cbind(as.data.frame(PDO.regr3))%>%
  cbind(as.data.frame(ONI.regr1))%>%
  cbind(as.data.frame(ONI.regr2))%>%
  cbind(as.data.frame(ONI.regr3))%>%
  cbind(as.data.frame(NPGO.regr1))%>%
  cbind(as.data.frame(NPGO.regr2))%>%
  cbind(as.data.frame(NPGO.regr3))%>%
  cbind(as.data.frame(NPH.regr1))%>%
  cbind(as.data.frame(NPH.regr2))%>%
  cbind(as.data.frame(NPH.regr3))%>%
  rename('A. PDO 1967 - 1988'= PDO.regr1,'B. PDO 1989 - 2012' = PDO.regr2, 'C. PDO 2013 - 2023' = PDO.regr3,
         'D. ONI 1967 - 1988'= ONI.regr1,'E. ONI 1989 - 2012' = ONI.regr2, 'F. ONI 2013 - 2023' = ONI.regr3,
         'G. NPGO 1967 - 1988'= NPGO.regr1,'H. NPGO 1989 - 2012' = NPGO.regr2, 'I. NPGO 2013 - 2023' = NPGO.regr3,
         'J. NPH 1967 - 1988'= NPH.regr1,'K. NPH 1989 - 2012' = NPH.regr2, 'L. NPH 2013 - 2023' = NPH.regr3)

X_Diff_SST<- as.data.frame(PDO.diff1)%>%
  cbind(as.data.frame(PDO.diff2))%>%
  cbind(as.data.frame(PDO.diff3))%>%
  cbind(as.data.frame(ONI.diff1))%>%
  cbind(as.data.frame(ONI.diff2))%>%
  cbind(as.data.frame(ONI.diff3))%>%
  cbind(as.data.frame(NPGO.diff1))%>%
  cbind(as.data.frame(NPGO.diff2))%>%
  cbind(as.data.frame(NPGO.diff3))%>%
  cbind(as.data.frame(NPH.diff1))%>%
  cbind(as.data.frame(NPH.diff2))%>%
  cbind(as.data.frame(NPH.diff3))%>%
  
  rename('A. PDO 1989:2012 - 1967:1988'= PDO.diff1,'B. PDO 2012:2023 - 1967:1988'= PDO.diff2,'C. PDO 2012:2023 - 1989:2012'= PDO.diff3,
         'D. ONI 1989:2012 - 1967:1988'= ONI.diff1,'E. ONI 2012:2023 - 1967:1988'= ONI.diff2,'F. ONI 2012:2023 - 1989:2012'= ONI.diff3,
         'G. NPGO 1989:2012 - 1967:1988'= NPGO.diff1,'H. NPGO 2012:2023 - 1967:1988'= NPGO.diff2,'I. NPGO 2012:2023 - 1989:2012'= NPGO.diff3,
         'J. NPH 1989:2012 - 1967:1988'= NPH.diff1,'K. NPH 2012:2023 - 1967:1988'= NPH.diff2,'L. NPH 2012:2023 - 1989:2012'= NPH.diff3)
X_PDO_SST$latitude <- spring_lat 
X_PDO_SST$longitude <- spring_lon+360


X_PDO_SST<- X_PDO_SST%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")


X_Diff_SST$latitude <- spring_lat 
X_Diff_SST$longitude <- spring_lon+360

X_Diff_SST<- X_Diff_SST%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")

world <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
#plot code
ggplot() + 
  geom_raster(data=X_Diff_SST, aes(x=longitude,y=latitude,fill = coefficient)) + 
  facet_wrap(~analysis, ncol = 3) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(120,135), ylim=c(0,60)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SST")+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 


#### SLP ####

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
# get times
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$winter_year <- dates$year
dates$winter_year[which(dates$month %in% 7:12)] <- NA
dates$winter_year[which(dates$month %in% 1:3)] <- NA

var_name = "msl"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
if(length(to_drop) > 0) {
  X <- X[,-to_drop]# drop cells with NAs 
  winter_lat <- lat[-to_drop]
  winter_lon <- lon[-to_drop]
} else {
  winter_lat <- lat
  winter_lon <- lon
}
# Add block for winter average calculations
winter_years <- unique(dates$winter_year)
winter_years <- sort(winter_years[-which(is.na(winter_years))])
X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,] = colMeans(X[which(dates$winter_year == winter_years[i]),])
}
dim(X_winter)
dimnames(X_winter) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))

# remove annual mean for each year
for(i in 1:nrow(X_winter)) {
  X_winter[i,] = X_winter[i,] - mean(X_winter[i,], na.rm=T)
}
X_winter <- scale(X_winter, center=TRUE, scale=FALSE)


#X_winter<- X_winter/100
slp_anomaly <- matrix(NA, length(winter_years), ncol(X))
mean_slp <- apply(X_winter,2,mean)
#sd_slp <- apply(X_winter,2,sd)
for(i in 1:ncol(X)) {
  slp_anomaly[,i] = (X_winter[,i]-mean_slp[i])
}

#slp_anomaly <-X_winter

dimnames(slp_anomaly) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))
dim(slp_anomaly)

SLP1 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:1988,] 
SLP2 <- slp_anomaly[rownames(slp_anomaly) %in% 1989:2012,]
SLP3 <- slp_anomaly[rownames(slp_anomaly) %in% 2013:2023,]
SLP4 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:2023,]
#### Climate SLP####

PDO.regr1 <- PDO.regr2<- PDO.regr3<- PDO.regr4<-NA
ONI.regr1 <- ONI.regr2<- ONI.regr3<-ONI.regr4<- NA
NPH.regr1 <- NPH.regr2<- NPH.regr3<-NPH.regr4<- NA
NPGO.regr1 <- NPGO.regr2<- NPGO.regr3<-NPGO.regr4<- NA

# now loop through each cell

for(i in 1:ncol(SLP1)){
  #  i <- 1
  mod <- lm(SLP1[,i] ~ PDO1)
  PDO.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP2[,i] ~ PDO2)
  PDO.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP3[,i] ~ PDO3)
  PDO.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP4[,i] ~ PDO4)
  PDO.regr4[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ ONI1)
  ONI.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ ONI2)
  ONI.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ ONI3)
  ONI.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP4[,i] ~ ONI4)
  ONI.regr4[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP1[,i] ~ NPGO1)
  NPGO.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NPGO2)
  NPGO.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ NPGO3)
  NPGO.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP4[,i] ~ NPGO4)
  NPGO.regr4[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ NPH1)
  NPH.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NPH2)
  NPH.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ NPH3)
  NPH.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP4[,i] ~ NPH4)
  NPH.regr4[i] <- summary(mod)$coef[2,1]
}
# calculate era differences for each cell

PDO.diff3 <- PDO.regr3 - PDO.regr4
PDO.diff2 <- PDO.regr2 - PDO.regr4
PDO.diff1 <- PDO.regr1 - PDO.regr4

NPGO.diff3 <- NPGO.regr3 - NPGO.regr4
NPGO.diff2 <- NPGO.regr2 - NPGO.regr4
NPGO.diff1 <- NPGO.regr1 - NPGO.regr4

NPH.diff3 <- NPH.regr3 - NPH.regr4
NPH.diff2 <- NPH.regr2 - NPH.regr4
NPH.diff1 <- NPH.regr1 - NPH.regr4

ONI.diff3 <- ONI.regr3 - ONI.regr4
ONI.diff2 <- ONI.regr2 - ONI.regr4
ONI.diff1 <- ONI.regr1 - ONI.regr4



X_PDO_SLP<- as.data.frame(PDO.regr1)%>%
  cbind(as.data.frame(PDO.regr2))%>%
  cbind(as.data.frame(PDO.regr3))%>%
  cbind(as.data.frame(PDO.regr4))%>%
  cbind(as.data.frame(ONI.regr1))%>%
  cbind(as.data.frame(ONI.regr2))%>%
  cbind(as.data.frame(ONI.regr3))%>%
  cbind(as.data.frame(ONI.regr4))%>%
  cbind(as.data.frame(NPGO.regr1))%>%
  cbind(as.data.frame(NPGO.regr2))%>%
  cbind(as.data.frame(NPGO.regr3))%>%
  cbind(as.data.frame(NPGO.regr4))%>%
  cbind(as.data.frame(NPH.regr1))%>%
  cbind(as.data.frame(NPH.regr2))%>%
  cbind(as.data.frame(NPH.regr3))%>%
  cbind(as.data.frame(NPH.regr4))%>%

  rename('A. PDO 1967 - 1988'= PDO.regr1,'B. PDO 1989 - 2012' = PDO.regr2, 
         'C. PDO 2013 - 2023' = PDO.regr3,'D. PDO 1967 - 2023' = PDO.regr4,
         
         'E. ONI 1967 - 1988'= ONI.regr1,'F. ONI 1989 - 2012' = ONI.regr2,
         'G. ONI 2013 - 2023' = ONI.regr3, 'H. ONI 1967 - 2023' = ONI.regr4,
         
         'I. NPGO 1967 - 1988'= NPGO.regr1,'J. NPGO 1989 - 2012' = NPGO.regr2,
         'K. NPGO 2013 - 2023' = NPGO.regr3,'L. NPGO 1967 - 2023' = NPGO.regr4,
         
         'M. NPH 1967 - 1988'= NPH.regr1,'N. NPH 1989 - 2012' = NPH.regr2, 
         'O. NPH 2013 - 2023' = NPH.regr3, 'P. NPH 1967 - 2023' = NPH.regr4)


X_PDO_SLP2<- as.data.frame(PDO.regr1)%>%
  cbind(as.data.frame(PDO.regr2))%>%
  cbind(as.data.frame(PDO.regr3))%>%
  cbind(as.data.frame(ONI.regr1))%>%
  cbind(as.data.frame(ONI.regr2))%>%
  cbind(as.data.frame(ONI.regr3))%>%
  cbind(as.data.frame(NPGO.regr1))%>%
  cbind(as.data.frame(NPGO.regr2))%>%
  cbind(as.data.frame(NPGO.regr3))%>%
  cbind(as.data.frame(NPH.regr1))%>%
  cbind(as.data.frame(NPH.regr2))%>%
  cbind(as.data.frame(NPH.regr3))%>%

  rename('A. PDO 1967 - 1988'= PDO.regr1,'B. PDO 1989 - 2012' = PDO.regr2, 
         'C. PDO 2013 - 2023' = PDO.regr3,
         
         'D. ONI 1967 - 1988'= ONI.regr1,'E. ONI 1989 - 2012' = ONI.regr2,
         'F. ONI 2013 - 2023' = ONI.regr3,
         
         'G. NPGO 1967 - 1988'= NPGO.regr1,'H. NPGO 1989 - 2012' = NPGO.regr2,
         'I. NPGO 2013 - 2023' = NPGO.regr3,
         
         'J. NPH 1967 - 1988'= NPH.regr1,'K. NPH 1989 - 2012' = NPH.regr2, 
         'L. NPH 2013 - 2023' = NPH.regr3)

X_Diff_SLP<- as.data.frame(PDO.diff1)%>%
  cbind(as.data.frame(PDO.diff2))%>%
  cbind(as.data.frame(PDO.diff3))%>%
  cbind(as.data.frame(ONI.diff1))%>%
  cbind(as.data.frame(ONI.diff2))%>%
  cbind(as.data.frame(ONI.diff3))%>%
  cbind(as.data.frame(NPGO.diff1))%>%
  cbind(as.data.frame(NPGO.diff2))%>%
  cbind(as.data.frame(NPGO.diff3))%>%
  cbind(as.data.frame(NPH.diff1))%>%
  cbind(as.data.frame(NPH.diff2))%>%
  cbind(as.data.frame(NPH.diff3))%>%
  
  rename('A. PDO 1967:1988 - 1967:2023'= PDO.diff1,'B. PDO 1989:2012 - 1967:2023'= PDO.diff2,'C. PDO 2013:2023 - 1967:2023'= PDO.diff3,
         'D. ONI 1967:1988 - 1967:2023'= ONI.diff1,'E. ONI 1989:2012 - 1967:2023'= ONI.diff2,'F. ONI 2013:2023 - 1967:2023'= ONI.diff3,
         'G. NPGO 1967:1988 - 1967:2023'= NPGO.diff1,'H. NPGO 1989:2012 - 1967:2023'= NPGO.diff2,'I. NPGO 2013:2023 - 1967:2023'= NPGO.diff3,
         'J. NPH 1967:1988 - 1967:2023'= NPH.diff1,'K. NPH 1989:2012 - 1967:2023'= NPH.diff2,'L. NPH 2013:2023 - 1967:2023'= NPH.diff3)

X_PDO_SLP$latitude <- winter_lat 
X_PDO_SLP$longitude <- winter_lon+360

X_PDO_SLP<- X_PDO_SLP%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")



X_Diff_SLP$latitude <- winter_lat 
X_Diff_SLP$longitude <- winter_lon+360

X_Diff_SLP<- X_Diff_SLP%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")

world <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
#plot code
ggplot() + 
  geom_raster(data=X_PDO_SLP, aes(x=longitude,y=latitude,fill = coefficient)) + 
  facet_wrap(~analysis, ncol = 4) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(220,250), ylim=c(20,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SLP")+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 



##### Data output for analysis summary qmd ####

cor_indices <- X_PDO_SST%>%
  mutate(var = 'SST')%>%
  rbind(X_PDO_SLP%>%mutate(var='SLP'))

#cor_cc <- X_cc%>%
#  mutate(var = 'SST')%>%
#  rbind(anomaly_cc_slp%>%mutate(var='SLP'))

cor_diff<- X_Diff_SLP%>%
 mutate(var='SLP')

saveRDS(cor_indices, file = here('data/physical/correlation_analysis_indices2.rds'))

#saveRDS(cor_cc, file = here('data/physical/correlation_analysis.rds'))

saveRDS(cor_diff, file = here('data/physical/correlation_analysis_diff.rds'))

ddd <-readRDS('data/physical/correlation_analysis.rds')



#### plotting SLP Spring means ####
# get lat/long
nc <- nc_open(here("som/copernicus_jul10.nc"))

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y)) 

# get times
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$winter_year <- dates$year
dates$winter_year[which(dates$month %in% 7:12)] <- NA
dates$winter_year[which(dates$month %in% 1:3)] <- NA
var_name = "msl"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
if(length(to_drop) > 0) {
  X <- X[,-to_drop]# drop cells with NAs 
  winter_lat <- lat[-to_drop]
  winter_lon <- lon[-to_drop]
} else {
  winter_lat <- lat
  winter_lon <- lon
}
# Add block for winter average calculations
winter_years <- unique(dates$winter_year)
winter_years <- sort(winter_years[-which(is.na(winter_years))])
X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,]= colMeans(X[which(dates$winter_year == winter_years[i]),])/100
}


#X_winter<- X_winter/100
slp_anomaly <- matrix(NA, length(winter_years), ncol(X))
mean_slp <- apply(X_winter,2,mean)
#sd_slp <- apply(X_winter,2,sd)
for(i in 1:ncol(X)) {
  slp_anomaly[,i] = (X_winter[,i]-mean_slp[i])
}
dimnames(slp_anomaly) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))
dim(slp_anomaly)


####plotting slp dat
SLP6 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:1988,] 
SLP7 <- slp_anomaly[rownames(slp_anomaly) %in% 1989:2012,]
SLP8 <- slp_anomaly[rownames(slp_anomaly) %in% 2013:2023,]

for(i in 1:ncol(SLP6)){
  #  i <- 1
  SLP6[i] <- mean(SLP6[,i])
  SLP7[i] <- mean(SLP7[,i])
  SLP8[i] <- mean(SLP8[,i])
}
dim(colMeans(SLP6))
SLP.dat<- as.data.frame(colMeans(SLP6))%>%
  cbind(as.data.frame(colMeans(SLP7)))%>%
  cbind(as.data.frame(colMeans(SLP8)))%>%
  rename(Era1='colMeans(SLP6)', Era2='colMeans(SLP7)',Era3='colMeans(SLP8)')

SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "Anomaly")


SLP.anom<- ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = Anomaly)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  # coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  scale_fill_gradientn(colours = c("blue","white","red")) + 
  ggtitle("SLP Anomalies (Spring)")+
  scale_x_continuous(breaks = c(190, 210, 230))+
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = Anomaly), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.anom
pdf("Output/SLPanom.pdf", 4,8) 
SLP.anom
dev.off()


#X_winter<- X_winter/100
slp_anomaly <- matrix(NA, length(winter_years), ncol(X))
mean_slp <- apply(X_winter,2,mean)
#sd_slp <- apply(X_winter,2,sd)
for(i in 1:ncol(X)) {
  slp_anomaly[,i] = (X_winter[,i])
}
dimnames(slp_anomaly) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))
dim(slp_anomaly)


####plotting slp dat
SLP6 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:1988,] 
SLP7 <- slp_anomaly[rownames(slp_anomaly) %in% 1989:2012,]
SLP8 <- slp_anomaly[rownames(slp_anomaly) %in% 2013:2023,]

for(i in 1:ncol(SLP6)){
  #  i <- 1
  SLP6[i] <- mean(SLP6[,i])
  SLP7[i] <- mean(SLP7[,i])
  SLP8[i] <- mean(SLP8[,i])
}
dim(colMeans(SLP6))
SLP.dat<- as.data.frame(colMeans(SLP6))%>%
  cbind(as.data.frame(colMeans(SLP7)))%>%
  cbind(as.data.frame(colMeans(SLP8)))%>%
  rename(Era1='colMeans(SLP6)', Era2='colMeans(SLP7)',Era3='colMeans(SLP8)')


SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "SLP")



SLP.plot<-ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = SLP)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  #coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  scale_fill_gradientn(colours = c("blue","white","red")) + 
  scale_x_continuous(breaks = c(190, 210, 230))+

  ggtitle("SLP (Spring)")+  
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SLP),col="lightgrey", lwd=0.5)+
  #geom_text_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SLP), stroke = 0.15)+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.plot


pdf("Output/SLPMeans.pdf", 4,8) 
SLP.plot
dev.off()



SLP.dat<- as.data.frame(apply(SLP6, 2, sd))%>%
  cbind(as.data.frame(apply(SLP7, 2, sd)))%>%
  cbind(as.data.frame(apply(SLP8, 2, sd)))%>%
  rename(Era1='apply(SLP6, 2, sd)', Era2='apply(SLP7, 2, sd)',Era3='apply(SLP8, 2, sd)')


SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "SD")



SLP.plot<-ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = SD)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  #coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  scale_fill_gradientn(colours = c("#704D9E", "#CF63A6", "#F7A086", "#F3E79A")) + 
  ggtitle("SLP SD (Spring)")+
  #geom_contour(colour='white')+
  scale_x_continuous(breaks = c(190, 210, 230))+
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SD), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.plot
pdf("Output/SLPSDs.pdf", 4,8) 
SLP.plot
dev.off()

####plotting slp dat for 2015 and 2022 ####


SLP6 <- slp_anomaly[rownames(slp_anomaly) %in% 2022,] 
SLP7 <- slp_anomaly[rownames(slp_anomaly) %in% 2015,]
SLP8 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:2023,]

for(i in 1:ncol(SLP8)){
  #  i <- 1
  SLP8[i] <- mean(SLP8[,i])
}

SLP.dat<- as.data.frame(SLP6)%>%
  cbind(as.data.frame(SLP7))%>%
  cbind(as.data.frame(colMeans(SLP8)))%>%
  rename('2022'='SLP6', '2015'='SLP7','1967-2013'='colMeans(SLP8)')


SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "SLP")



SLP.plot<-ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = SLP)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  #coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  scale_fill_gradientn(colours = c("blue","white","red")) + 
  scale_x_continuous(breaks = c(190, 210, 230))+

  ggtitle("SLP (Spring)")+  
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SLP),col="lightgrey", lwd=0.5)+
  #geom_text_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SLP), stroke = 0.15)+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.plot
#### plotting SLP Winter means ####
# get lat/long
nc <- nc_open(here("som/copernicus_jul10.nc"))

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y)) 

# get times
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$winter_year <- dates$year
dates$winter_year[which(dates$month %in% 4:10)] <- NA
#dates$winter_year[which(dates$month %in% 1:3)] <- NA
var_name = "msl"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
if(length(to_drop) > 0) {
  X <- X[,-to_drop]# drop cells with NAs 
  winter_lat <- lat[-to_drop]
  winter_lon <- lon[-to_drop]
} else {
  winter_lat <- lat
  winter_lon <- lon
}
# Add block for winter average calculations
winter_years <- unique(dates$winter_year)
winter_years <- sort(winter_years[-which(is.na(winter_years))])
X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,]= colMeans(X[which(dates$winter_year == winter_years[i]),])/100
}


#X_winter<- X_winter/100
slp_anomaly <- matrix(NA, length(winter_years), ncol(X))
mean_slp <- apply(X_winter,2,mean)
#sd_slp <- apply(X_winter,2,sd)
for(i in 1:ncol(X)) {
  slp_anomaly[,i] = (X_winter[,i]-mean_slp[i])
}
dimnames(slp_anomaly) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))
dim(slp_anomaly)


####plotting slp dat
SLP6 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:1988,] 
SLP7 <- slp_anomaly[rownames(slp_anomaly) %in% 1989:2012,]
SLP8 <- slp_anomaly[rownames(slp_anomaly) %in% 2013:2023,]

for(i in 1:ncol(SLP6)){
  #  i <- 1
  SLP6[i] <- mean(SLP6[,i])
  SLP7[i] <- mean(SLP7[,i])
  SLP8[i] <- mean(SLP8[,i])
}
dim(colMeans(SLP6))
SLP.dat<- as.data.frame(colMeans(SLP6))%>%
  cbind(as.data.frame(colMeans(SLP7)))%>%
  cbind(as.data.frame(colMeans(SLP8)))%>%
  rename(Era1='colMeans(SLP6)', Era2='colMeans(SLP7)',Era3='colMeans(SLP8)')

SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "Anomaly")


SLP.anom<- ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = Anomaly)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  # coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  scale_fill_gradientn(colours = c("blue","white","red")) + 
  ggtitle("SLP Anomalies (Winter)")+
  scale_x_continuous(breaks = c(190, 210, 230))+
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = Anomaly), col="lightgrey",lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.anom
pdf("Output/SLPanomWinter.pdf", 4,8) 
SLP.anom
dev.off()


#X_winter<- X_winter/100
slp_anomaly <- matrix(NA, length(winter_years), ncol(X))
mean_slp <- apply(X_winter,2,mean)
#sd_slp <- apply(X_winter,2,sd)
for(i in 1:ncol(X)) {
  slp_anomaly[,i] = (X_winter[,i])
}
dimnames(slp_anomaly) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))
dim(slp_anomaly)


####plotting slp dat
SLP6 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:1988,] 
SLP7 <- slp_anomaly[rownames(slp_anomaly) %in% 1989:2012,]
SLP8 <- slp_anomaly[rownames(slp_anomaly) %in% 2013:2023,]

for(i in 1:ncol(SLP6)){
  #  i <- 1
  SLP6[i] <- mean(SLP6[,i])
  SLP7[i] <- mean(SLP7[,i])
  SLP8[i] <- mean(SLP8[,i])
}
dim(colMeans(SLP6))
SLP.dat<- as.data.frame(colMeans(SLP6))%>%
  cbind(as.data.frame(colMeans(SLP7)))%>%
  cbind(as.data.frame(colMeans(SLP8)))%>%
  rename(Era1='colMeans(SLP6)', Era2='colMeans(SLP7)',Era3='colMeans(SLP8)')


SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "SLP")



SLP.plot<-ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = SLP)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  #coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  scale_fill_gradientn(colours = c("blue","white","red")) + 
  ggtitle("SLP (Winter)")+
  scale_x_continuous(breaks = c(190, 210, 230))+
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SLP),col="lightgrey", lwd=0.5)+
  #geom_text_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SLP), stroke = 0.15)+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.plot


pdf("Output/SLPMeansWinter.pdf", 4,8) 
SLP.plot
dev.off()



SLP.dat<- as.data.frame(apply(SLP6, 2, sd))%>%
  cbind(as.data.frame(apply(SLP7, 2, sd)))%>%
  cbind(as.data.frame(apply(SLP8, 2, sd)))%>%
  rename(Era1='apply(SLP6, 2, sd)', Era2='apply(SLP7, 2, sd)',Era3='apply(SLP8, 2, sd)')


SLP.dat$latitude <- winter_lat 
SLP.dat$longitude <- winter_lon+360
SLP.dat<- SLP.dat%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "SD")



SLP.plot<-ggplot() + 
  geom_raster(data=SLP.dat, aes(x=longitude,y=latitude,fill = SD)) + 
  facet_wrap(~analysis, ncol = 1) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  #coord_sf(xlim=c(140,240), ylim=c(0,60)) +
  scale_fill_gradientn(colours = c("#704D9E", "#CF63A6", "#F7A086", "#F3E79A")) + 
  ggtitle("SLP SD (Winter)")+
  #geom_contour(colour='white')+
  scale_x_continuous(breaks = c(190, 210, 230))+
  geom_contour(data=SLP.dat, aes(x=longitude,y=latitude,z = SD), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP.plot
pdf("Output/SLPSDsWinter.pdf", 4,8) 
SLP.plot
dev.off()


iceland <- ne_countries(scale = 10, country = "USA", returnclass = "sf")

#info of our spatial vector object
iceland

require(rnaturalearth)
world <- ne_countries(scale = "small", returnclass = "sf")
world_map <- world %>% 
  dplyr::group_by(continent) %>% 
  dplyr::select(continent) %>% 
  st_buffer(dist=.001) 

plot(ne_coastline())



#### SLP Detrended####

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
# get times
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$winter_year <- dates$year
dates$winter_year[which(dates$month %in% 7:12)] <- NA
dates$winter_year[which(dates$month %in% 1:3)] <- NA

var_name = "msl"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
if(length(to_drop) > 0) {
  X <- X[,-to_drop]# drop cells with NAs 
  winter_lat <- lat[-to_drop]
  winter_lon <- lon[-to_drop]
} else {
  winter_lat <- lat
  winter_lon <- lon
}

# remove annual mean for each year
for(i in 1:nrow(X_winter)) {
  X_winter[i,] = X_winter[i,] - mean(X_winter[i,], na.rm=T)
}
X_winter <- scale(X_winter)




#slp_anomaly <-X_winter

dimnames(slp_anomaly) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))
dim(slp_anomaly)

SLP1 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:1988,] 
SLP2 <- slp_anomaly[rownames(slp_anomaly) %in% 1989:2012,]
SLP3 <- slp_anomaly[rownames(slp_anomaly) %in% 2013:2023,]
SLP4 <- slp_anomaly[rownames(slp_anomaly) %in% 1967:2023,]
#### Climate SLP Detrended ####

PDO.regr1 <- PDO.regr2<- PDO.regr3<- PDO.regr4<-NA
ONI.regr1 <- ONI.regr2<- ONI.regr3<-ONI.regr4<- NA
NPH.regr1 <- NPH.regr2<- NPH.regr3<-NPH.regr4<- NA
NPGO.regr1 <- NPGO.regr2<- NPGO.regr3<-NPGO.regr4<- NA

# now loop through each cell

for(i in 1:ncol(SLP1)){
  #  i <- 1
  mod <- lm(SLP1[,i] ~ PDO1)
  PDO.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP2[,i] ~ PDO2)
  PDO.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP3[,i] ~ PDO3)
  PDO.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP4[,i] ~ PDO4)
  PDO.regr4[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ ONI1)
  ONI.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ ONI2)
  ONI.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ ONI3)
  ONI.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP4[,i] ~ ONI4)
  ONI.regr4[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP1[,i] ~ NPGO1)
  NPGO.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NPGO2)
  NPGO.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ NPGO3)
  NPGO.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP4[,i] ~ NPGO4)
  NPGO.regr4[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ NPH1)
  NPH.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NPH2)
  NPH.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ NPH3)
  NPH.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP4[,i] ~ NPH4)
  NPH.regr4[i] <- summary(mod)$coef[2,1]
}
# calculate era differences for each cell

PDO.diff3 <- PDO.regr3 - PDO.regr4
PDO.diff2 <- PDO.regr2 - PDO.regr4
PDO.diff1 <- PDO.regr1 - PDO.regr4

NPGO.diff3 <- NPGO.regr3 - NPGO.regr4
NPGO.diff2 <- NPGO.regr2 - NPGO.regr4
NPGO.diff1 <- NPGO.regr1 - NPGO.regr4

NPH.diff3 <- NPH.regr3 - NPH.regr4
NPH.diff2 <- NPH.regr2 - NPH.regr4
NPH.diff1 <- NPH.regr1 - NPH.regr4

ONI.diff3 <- ONI.regr3 - ONI.regr4
ONI.diff2 <- ONI.regr2 - ONI.regr4
ONI.diff1 <- ONI.regr1 - ONI.regr4



X_PDO_SLP<- as.data.frame(PDO.regr1)%>%
  cbind(as.data.frame(PDO.regr2))%>%
  cbind(as.data.frame(PDO.regr3))%>%
  cbind(as.data.frame(PDO.regr4))%>%
  cbind(as.data.frame(ONI.regr1))%>%
  cbind(as.data.frame(ONI.regr2))%>%
  cbind(as.data.frame(ONI.regr3))%>%
  cbind(as.data.frame(ONI.regr4))%>%
  cbind(as.data.frame(NPGO.regr1))%>%
  cbind(as.data.frame(NPGO.regr2))%>%
  cbind(as.data.frame(NPGO.regr3))%>%
  cbind(as.data.frame(NPGO.regr4))%>%
  cbind(as.data.frame(NPH.regr1))%>%
  cbind(as.data.frame(NPH.regr2))%>%
  cbind(as.data.frame(NPH.regr3))%>%
  cbind(as.data.frame(NPH.regr4))%>%

  rename('A. PDO 1967 - 1988'= PDO.regr1,'B. PDO 1989 - 2012' = PDO.regr2, 
         'C. PDO 2013 - 2023' = PDO.regr3,'D. PDO 1967 - 2023' = PDO.regr4,
         
         'E. ONI 1967 - 1988'= ONI.regr1,'F. ONI 1989 - 2012' = ONI.regr2,
         'G. ONI 2013 - 2023' = ONI.regr3, 'H. ONI 1967 - 2023' = ONI.regr4,
         
         'I. NPGO 1967 - 1988'= NPGO.regr1,'J. NPGO 1989 - 2012' = NPGO.regr2,
         'K. NPGO 2013 - 2023' = NPGO.regr3,'L. NPGO 1967 - 2023' = NPGO.regr4,
         
         'M. NPH 1967 - 1988'= NPH.regr1,'N. NPH 1989 - 2012' = NPH.regr2, 
         'O. NPH 2013 - 2023' = NPH.regr3, 'P. NPH 1967 - 2023' = NPH.regr4)


X_PDO_SLP2<- as.data.frame(PDO.regr1)%>%
  cbind(as.data.frame(PDO.regr2))%>%
  cbind(as.data.frame(PDO.regr3))%>%
  cbind(as.data.frame(ONI.regr1))%>%
  cbind(as.data.frame(ONI.regr2))%>%
  cbind(as.data.frame(ONI.regr3))%>%
  cbind(as.data.frame(NPGO.regr1))%>%
  cbind(as.data.frame(NPGO.regr2))%>%
  cbind(as.data.frame(NPGO.regr3))%>%
  cbind(as.data.frame(NPH.regr1))%>%
  cbind(as.data.frame(NPH.regr2))%>%
  cbind(as.data.frame(NPH.regr3))%>%

  rename('A. PDO 1967 - 1988'= PDO.regr1,'B. PDO 1989 - 2012' = PDO.regr2, 
         'C. PDO 2013 - 2023' = PDO.regr3,
         
         'D. ONI 1967 - 1988'= ONI.regr1,'E. ONI 1989 - 2012' = ONI.regr2,
         'F. ONI 2013 - 2023' = ONI.regr3,
         
         'G. NPGO 1967 - 1988'= NPGO.regr1,'H. NPGO 1989 - 2012' = NPGO.regr2,
         'I. NPGO 2013 - 2023' = NPGO.regr3,
         
         'J. NPH 1967 - 1988'= NPH.regr1,'K. NPH 1989 - 2012' = NPH.regr2, 
         'L. NPH 2013 - 2023' = NPH.regr3)

X_Diff_SLP<- as.data.frame(PDO.diff1)%>%
  cbind(as.data.frame(PDO.diff2))%>%
  cbind(as.data.frame(PDO.diff3))%>%
  cbind(as.data.frame(ONI.diff1))%>%
  cbind(as.data.frame(ONI.diff2))%>%
  cbind(as.data.frame(ONI.diff3))%>%
  cbind(as.data.frame(NPGO.diff1))%>%
  cbind(as.data.frame(NPGO.diff2))%>%
  cbind(as.data.frame(NPGO.diff3))%>%
  cbind(as.data.frame(NPH.diff1))%>%
  cbind(as.data.frame(NPH.diff2))%>%
  cbind(as.data.frame(NPH.diff3))%>%
  
  rename('A. PDO 1967:1988 - 1967:2023'= PDO.diff1,'B. PDO 1989:2012 - 1967:2023'= PDO.diff2,'C. PDO 2013:2023 - 1967:2023'= PDO.diff3,
         'D. ONI 1967:1988 - 1967:2023'= ONI.diff1,'E. ONI 1989:2012 - 1967:2023'= ONI.diff2,'F. ONI 2013:2023 - 1967:2023'= ONI.diff3,
         'G. NPGO 1967:1988 - 1967:2023'= NPGO.diff1,'H. NPGO 1989:2012 - 1967:2023'= NPGO.diff2,'I. NPGO 2013:2023 - 1967:2023'= NPGO.diff3,
         'J. NPH 1967:1988 - 1967:2023'= NPH.diff1,'K. NPH 1989:2012 - 1967:2023'= NPH.diff2,'L. NPH 2013:2023 - 1967:2023'= NPH.diff3)

X_PDO_SLP$latitude <- winter_lat 
X_PDO_SLP$longitude <- winter_lon+360

X_PDO_SLP<- X_PDO_SLP%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")



X_Diff_SLP$latitude <- winter_lat 
X_Diff_SLP$longitude <- winter_lon+360

X_Diff_SLP<- X_Diff_SLP%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")

world <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
#plot code
ggplot() + 
  geom_raster(data=X_PDO_SLP, aes(x=longitude,y=latitude,fill = coefficient)) + 
  facet_wrap(~analysis, ncol = 4) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(220,250), ylim=c(20,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SLP")+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 

global_trend
