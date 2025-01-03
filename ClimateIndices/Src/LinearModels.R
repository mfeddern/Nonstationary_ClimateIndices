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
library(rstan)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(PBSmapping)
library(bayestestR)
library(dplyr)
set.seed(1234)

##### Writing Bayes DFA model function ####
col<-pnw_palette("Sunset2",3,type="discrete")
linmod <- function(input_dat,yvar,columns,RESULTSDF, periods, coef) {
  for(i in 1:length(columns)){
    mod<-lm(input_dat[,which(colnames(input_dat) == yvar)]~input_dat[,columns[i]]* as.factor(input_dat$period),data=input_dat)
    results<- data.frame(coef(summary(mod)))%>%
      select(Estimate,Std..Error)%>%
      mutate(period=periods, coef=coef, Index=colnames(data[columns[i]]))%>%
      pivot_wider(names_from = coef, values_from = c(Estimate, Std..Error))%>%
      rename(se_alpha=Std..Error_alpha, se_beta=Std..Error_beta)%>%
      mutate(alpha=ifelse(period==1,Estimate_alpha,Estimate_alpha+Estimate_alpha[1]),
             beta=ifelse(period==1,Estimate_beta,Estimate_beta+Estimate_beta[1]))%>%
      select(period, Index, alpha, beta, se_alpha, se_beta)
    RESULTSDF<-rbind(results,RESULTSDF)
  }
  return(RESULTSDF)
}  


lev <- function(input_dat, yvar, columns,RESULTSDF, periods, coef) {
  for(i in 1:length(columns)){
    mod<-lm(input_dat[,which(colnames(input_dat) == yvar)]~input_dat[,columns[i]]* as.factor(data$period),data=input_dat)
    leverage<- as.data.frame(hatvalues(mod))%>%
      cbind(year=unique(input_dat$Year_lag))%>%
      mutate(Index=colnames(data[columns[i]]))%>%
      mutate(hat=hatvalues(mod))
    RESULTS_leverage<-rbind(leverage,RESULTS_leverage)
  }
  return(RESULTS_leverage)
}

linmodlev <- function(input_dat, yvar,columns,RESULTSDF, periodslev,coeflev,base) {
  for(i in 1:length(columns)){
    mod<-lm(input_dat[,which(colnames(input_dat) == yvar)]~input_dat[,columns[i]]* as.factor(input_dat$period),data=input_dat)
    leverage<- as.data.frame(hatvalues(mod))%>%
      cbind(year=unique(input_dat$Year_lag))%>%
      mutate(Index=colnames(input_dat[columns[i]]))%>%
      mutate(hat=hatvalues(mod))%>%
      filter(hat>3*((2+1)/length(input_dat$Year_lag)))
    
    datalev <- input_dat[ ! input_dat$Year_lag %in% leverage$year, ]
    mod2<-lm(datalev[,which(colnames(datalev) == yvar)]~datalev[,columns[i]]* as.factor(datalev$period),data=datalev)
    results<- data.frame(coef(summary(mod2)))%>%
      select(Estimate,Std..Error)%>%
      mutate(period=periodslev, coef=coeflev, Index=colnames(datalev[columns[i]]))%>%
      pivot_wider(names_from = coef, values_from = c(Estimate, Std..Error))%>%
      rename(se_alpha=Std..Error_alpha, se_beta=Std..Error_beta)%>%
      mutate(alpha=ifelse(period==1,Estimate_alpha,Estimate_alpha+Estimate_alpha[base]),
             beta=ifelse(period==1,Estimate_beta,Estimate_beta+Estimate_beta[base]))%>%
      select(period, Index, alpha, beta, se_alpha, se_beta)
    RESULTSDF<-rbind(results,RESULTSDF)
  }
  return(RESULTSDF)
}  


#### CALCOFI ####
######Spring #####
coef <- c("alpha", "beta", "alpha", "alpha", "beta", "beta")
periods<- c(1,1,2,3,2,3)
base<-1
data <- dat%>%filter(season=="Spring", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)



RESULTSDF <-NULL
RESULTS_leverage <-data.frame()
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))

CALCOFI<- linmod(data,"estimate",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="CALCOFI",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                      ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 0, leverage="Included")%>%
 bind_rows(linmodlev(data,"estimate",columns, RESULTSDF,periods, coef,base)%>%
  mutate(survey="CALCOFI",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 0, leverage="Excluded"))

leverage<-lev(data,"estimate",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="CALCOFI",
          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
         ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 0)

ggplot(RREASoffset, aes(x = beta,col = as.factor(leverage),y=as.factor(period)), group=as.factor(period)) +
  theme_bw() +
  facet_wrap(.~Index, ncol = 2, scales='free') +
  geom_point(alpha = 0.7) +
  geom_errorbarh(aes(xmin=beta-2*se_beta, xmax=beta+2*se_beta,))+
  #scale_colour_manual(values = c(col[1])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

#running models with 1-year offset

CALCOFIoffset<- linmod(data,"estimateoffset1",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="CALCOFI",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 1, leverage="Included")%>%
  bind_rows(linmodlev(data,"estimate",columns, RESULTSDF,periods, coef,base)%>%
              mutate(survey="CALCOFI",
                     Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                       ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
              mutate(region="SCC",Season="Spring", lag = 1, leverage="Excluded"))

#### RREAS ####
######Spring #####

data <- dat%>%filter(season=="Spring", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)
coeflev <- c("alpha", "beta", "alpha", "beta")
periodslev<- c(2,2,3,3)

RREAS<- linmod(data,"estimate",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="RREAS",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 0, leverage="Included")%>%
  bind_rows(linmodlev(data%>%filter(Year_lag>1988),"estimate",columns, RESULTSDF,periodslev,coeflev,1)%>%
              mutate(survey="RREAS",
                     Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                       ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
              mutate(region="SCC",Season="Spring", lag = 0, leverage="Excluded"))

#running models with 1-year offset

RREASoffset<- linmod(data,"estimateoffset1",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="RREAS",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 1, leverage="Included")%>%
  bind_rows(linmodlev(data%>%filter(Year_lag>1988),"estimate",columns, RESULTSDF,periodslev,coeflev,base)%>%
              mutate(survey="RREAS",
                     Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                       ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
              mutate(region="SCC",Season="Spring", lag = 1, leverage="Excluded"))


#### Copepod ####
###### Spring ####
climpdo <- climate_dat_cop%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)

data <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
                seasonal_copepod_southern)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)%>%
  mutate(full=1)

northern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_northern, data$period)
  northern <-rbind(northern,scaled.anom)
}
northern<-northern%>%mutate(survey="N. Copepod",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

northern_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_northern, data$full)
  northern_full <-rbind(northern_full,scaled.anom)
}

northern_full<-northern_full%>%mutate(survey="N. Copepod",
                                      Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                        ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

index.names <- unique(northern$Index)
period.names <- unique(northern$period)


northern<- bind_rows(northern,northern_full)%>%
  mutate(region="NCC",Season="Spring", lag = 0)
overlap.northern <- NA
i<-1
for(i in 1:4){
  temp <- northern%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
  ov2 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])
  ov1 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
  temp2<-rbind(ov3, ov2, ov1)
  overlap.northern <-rbind(temp2,overlap.northern)
}
overlap.northern<-overlap.northern%>%mutate(Survey="N. Copepods", region="NCC", season="Spring", offset=0) 


southern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_southern, data$period)
  southern <-rbind(southern,scaled.anom)
}
southern<-southern%>%mutate(survey="Southern Copepod (NCC)",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

southern_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_southern, data$full)
  southern_full <-rbind(southern_full,scaled.anom)
}

southern_full<-southern_full%>%mutate(survey="S. Copepod (NCC)",
                                      Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                        ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

index.names <- unique(southern$Index)
period.names <- unique(southern$period)
southern<- bind_rows(southern,southern_full)%>%
  mutate(region="NCC",Season="Spring", lag = 0)
overlap.southern <- NA

for(i in 1:4){
  temp <- southern%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
  ov2 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])
  ov1 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
  temp2<-rbind(ov3, ov2, ov1)
  overlap.southern <-rbind(temp2,overlap.southern)
}
overlap.southern<-overlap.southern%>%mutate(Survey="S. Copepods", region="NCC", season="Spring", offset=0) 

###### Winter ####
climpdo <- climate_dat_cop%>%filter(season=="Winter")%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)

data <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
                seasonal_copepod_southern)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)%>%
  mutate(full=1)

northernW <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_northern, data$period)
  northernW <-rbind(northernW,scaled.anom)
}
northernW<-northernW%>%mutate(survey="N. Copepod (NCC)",
                              Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

northernW_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_northern, data$full)
  northernW_full <-rbind(northernW_full,scaled.anom)
}

northernW_full<-northernW_full%>%mutate(survey="N. Copepod (NCC)",
                                        Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                          ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

index.names <- unique(northernW$Index)
period.names <- unique(northernW$period)


northernW<- bind_rows(northernW,northernW_full)%>%
  mutate(region="NCC",Season="Winter", lag = 0)
overlap.northernW <- NA
i<-1
for(i in 1:4){
  temp <- northernW%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
  ov2 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])
  ov1 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
  temp2<-rbind(ov3, ov2, ov1)
  overlap.northernW <-rbind(temp2,overlap.northernW)
}
overlap.northernW<-overlap.northernW%>%mutate(Survey="N. Copepods", region="NCC", season="Winter", offset=0) 


southernW <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_southern, data$period)
  southernW <-rbind(southernW,scaled.anom)
}
southernW<-southernW%>%mutate(survey="S. Copepod (NCC)",
                              Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

southernW_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$seasonal_copepod_southern, data$full)
  southernW_full <-rbind(southernW_full,scaled.anom)
}

southernW_full<-southernW_full%>%mutate(survey="S. Copepod (NCC)",
                                        Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                          ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

index.names <- unique(southernW$Index)
period.names <- unique(southernW$period)
southernW<- bind_rows(southernW,southernW_full)%>%
  mutate(region="NCC",Season="Winter", lag = 0)
overlap.southernW <- NA
for(i in 1:4){
  temp <- southernW%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
  ov2 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])
  ov1 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
  temp2<-rbind(ov3, ov2, ov1)
  overlap.southernW <-rbind(temp2,overlap.southernW)
}
overlap.southernW<-overlap.southernW%>%mutate(Survey="S. Copepods", region="NCC", season="Winter", offset=0) 






#### Upwelling Model Runs ####

###### spring #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Spring"
data<- climate_dat%>%filter(season=="Spring")%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()

upwelling <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_bakun_seasonally, data$era.region)
  upwelling <-rbind(upwelling,scaled.anom)
}

upwelling<-upwelling%>%mutate(survey="Upwelling",era.region=period,
                              Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region,period)%>%distinct())%>%
  mutate(Season='Spring', lag=0)


unique(upwelling$period)



upwelling_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_bakun_seasonally, data$region)
  upwelling_full <-rbind(upwelling_full,scaled.anom)
}

reg<-unique(data%>%dplyr::select(region))%>%
  rename(region2=region)%>%
  mutate(region=as.numeric(as.factor(region2)))

upwelling_full<-upwelling_full%>%mutate(survey="Upwelling",region=period,
                                        Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                          ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  left_join(reg)%>%mutate(period=4)%>%dplyr::select(!region)%>%rename(region=region2)%>%
  mutate(Season='Spring', lag=0)


upwelling <- upwelling%>%dplyr::select(!era.region)%>%mutate(period=as.numeric(as.factor(period)))%>%
  bind_rows(upwelling_full%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,region,period,Season,lag))

index.names <- unique(upwelling$Index)
region.names<-unique(upwelling$region)
period.names <- unique(upwelling$period)
overlap.up <- NA
for(j in 1:4){
  temp1 <- upwelling%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(1), period2=c(4), Index=index.names[i], region=region.names[j])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==2)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(2), period2=c(4), Index=index.names[i], region=region.names[j])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(3), period2=c(4), Index=index.names[i], region=region.names[j])
    
    temp2<-rbind(ov1,ov2,ov3,ov4,ov5,ov6)
    overlap.up <-rbind(temp2,overlap.up)
  }
  
}

overlap.up


#### Upwelling Model Runs ####

###### winter #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Winter"
data<- climate_dat%>%filter(season=="Winter")%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()

upwellingW <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_bakun_seasonally, data$era.region)
  upwellingW <-rbind(upwellingW,scaled.anom)
}

upwellingW<-upwellingW%>%mutate(survey="Upwelling",era.region=period,
                                Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region,period)%>%distinct())%>%
  mutate(Season='Winter', lag=0)


unique(upwellingW$period)



upwellingW_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_bakun_seasonally, data$region)
  upwellingW_full <-rbind(upwellingW_full,scaled.anom)
}

reg<-unique(data%>%dplyr::select(region))%>%
  rename(region2=region)%>%
  mutate(region=as.numeric(as.factor(region2)))

upwellingW_full<-upwellingW_full%>%mutate(survey="Upwelling",region=period,
                                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  left_join(reg)%>%mutate(period=4)%>%dplyr::select(!region)%>%rename(region=region2)%>%
  mutate(Season='Winter', lag=0)


upwellingW <- upwellingW%>%dplyr::select(!era.region)%>%mutate(period=as.numeric(as.factor(period)))%>%
  bind_rows(upwellingW_full%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,region,period,Season,lag))

index.names <- unique(upwellingW$Index)
region.names<-unique(upwellingW$region)
period.names <- unique(upwellingW$period)
overlap.upW <- NA
for(j in 1:4){
  temp1 <- upwellingW%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(1), period2=c(4), Index=index.names[i], region=region.names[j])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==2)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(2), period2=c(4), Index=index.names[i], region=region.names[j])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(3), period2=c(4), Index=index.names[i], region=region.names[j])
    
    temp2<-rbind(ov1,ov2,ov3,ov4,ov5,ov6)
    overlap.upW <-rbind(temp2,overlap.upW)
  }
  
}

overlap.upW


##### Phenology Model Runs #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))%>%
  mutate(region=fct_relevel(region,c("Northern CC","Central CC","Southern CC")))


season <- "Spring"
eras <- data.frame(era.region2=seq(1,9), era.region=seq(4,12))
data<- climate_dat%>%filter(season=="Spring"&region!='GoA'&Year_lag<2023&period!=4)%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()%>%
  left_join(eras)
index.names <-c("PDO", "NPGO", "NPH", "ONI")
region.names<-unique(climate_dat$region)
period.names <- unique(climate_dat$period)

data.phe.lm<-data%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, season, period,region,stand_tumi,stand_lusi,stand_sti,stand_bakun_seasonally,
                seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, period,region,stand_tumi,stand_lusi,stand_sti,stand_bakun_seasonally), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  pivot_longer(!c(Year_lag, season,period,region, Index_Value, Index_Name), 
               names_to = "Up_Name", values_to = "Up_Value")%>%
  distinct()

###### STI #######

sti_dat <-ggplot(data = data.phe.lm%>%filter(Up_Name=="stand_sti"&period!=4), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "STI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
sti_dat
pdf(file = "Output/Supplemental/FigureS10_STIlinearregression.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
sti_dat
dev.off()



STI<-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_sti, data$era.region2)
  STI <-rbind(STI,scaled.anom)
}
STI<-STI%>%mutate(survey="STI",era.region2=period,
                  Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                    ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region2,period)%>%distinct())%>%
  mutate(Season='Spring', lag=0, era.region=era.region2,period=as.numeric(as.factor(period)))


STI_beta <-ggplot(STI, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS11_STIbeta.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
STI_beta
dev.off()

STI_alpha <-ggplot(STI, aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS12_STIalpha.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
STI_alpha
dev.off()

index.names <-unique(STI$Index)
region.names<-unique(STI$region)
period.names <- unique(STI$period)

postplot(STI, STI$beta)
postplot(STI,STI$alpha)
as.numeric(STI$period)
overlap.sti <- NA
for(j in 2:4){
  temp1 <- STI%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    temp2<-rbind(ov1,ov2,ov3)
    overlap.sti <-rbind(temp2,overlap.sti)
  }
  
}

overlap.sti<-overlap.sti%>%mutate(survey="STI")


###### LUSI #######
lusi_dat <-ggplot(data = data.phe.lm%>%filter(Up_Name=="stand_lusi"&period!=4), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "LUSI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
lusi_dat
pdf(file = "Output/Supplemental/FigureS13_LUSIlinearregression.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
lusi_dat
dev.off()

LUSI<-NULL
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_lusi, data$era.region2)
  LUSI <-rbind(LUSI,scaled.anom)
}
LUSI<-LUSI%>%mutate(survey="LUSI",era.region2=period,
                    Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                      ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region2,period)%>%distinct())%>%
  mutate(Season='Spring', lag=0, era.region=era.region2,period=as.numeric(as.factor(period)))



LUSI_beta <-ggplot(LUSI, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS14_LUSIbeta.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
LUSI_beta
dev.off()

LUSI_alpha <-ggplot(LUSI, aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS15_LUSIalpha.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
LUSI_alpha
dev.off()

overlap.lusi <- NA
for(j in 2:4){
  temp1 <- LUSI%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    temp2<-rbind(ov1,ov2,ov3)
    overlap.lusi <-rbind(temp2,overlap.lusi)
  }
  
}

overlap.lusi
overlap.lusi<-overlap.lusi%>%mutate(survey="LUSI")

###### TUMI #######
tumi_dat <-ggplot(data = data.phe.lm%>%filter(Up_Name=="stand_tumi"&period!=4), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "TUMI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
tumi_dat
pdf(file = "Output/Supplemental/FigureS7_TUMIlinearregression.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
tumi_dat
dev.off()

TUMI<-NULL
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_tumi, data$era.region2)
  TUMI <-rbind(TUMI,scaled.anom)
}
TUMI<-TUMI%>%mutate(survey="TUMI",era.region2=period,
                    Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                      ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region2,period)%>%distinct())%>%
  mutate(Season='Spring', lag=0, era.region=era.region2,period=as.numeric(as.factor(period)))


TUMI_beta <-ggplot(TUMI, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS8_TUMIbeta.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
TUMI_beta
dev.off()

TUMI_alpha <-ggplot(TUMI, aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS9_TUMIalpha.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
TUMI_alpha
dev.off()


overlap.tumi <- NA
for(j in 2:4){
  temp1 <- TUMI%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    temp2<-rbind(ov1,ov2,ov3)
    overlap.tumi <-rbind(temp2,overlap.tumi)
  }
  
}

overlap.tumi<-overlap.tumi%>%mutate(survey="TUMI")

overlap.phe<-na.omit(overlap.tumi)%>%
  add_row(na.omit(overlap.lusi))%>%
  add_row(na.omit(overlap.sti))


overlap.phe%>%filter(period1==1,period2==3)%>%summarise(mean=mean(ov))
overlap.phe%>%filter(period1==1,period2==3)%>%summarise(sd=sd(ov))

overlap.phe%>%filter(period1==3,period2==2)%>%summarise(mean=mean(ov))
overlap.phe%>%filter(period1==3,period2==2)%>%summarise(sd=sd(ov))

overlap.phe%>%filter(region=="Northern CC",period1==3,period2==2)
overlap.phe%>%filter(region=="Southern CC",period1==3,period2==2)
overlap.phe%>%filter(region=="Central CC",period1==3,period2==2)

overlap.phe%>%filter(region=="Northern CC",period1==3,period2==2, survey=='TUMI'|survey=="STI", Index=="ONI"|Index=="NPGO"|Index=="PDO")%>%summarise(mean=mean(ov))
overlap.phe%>%filter(region=="Central CC",period1==3,period2==2, survey=='TUMI'|survey=="STI", Index=="ONI"|Index=="NPGO"|Index=="PDO")%>%summarise(mean=mean(ov))
overlap.phe%>%filter(region=="Southern CC",period1==3,period2==2, survey=='TUMI'|survey=="STI", Index=="ONI"|Index=="NPGO"|Index=="PDO")%>%summarise(mean=mean(ov))

overlap.phe%>%filter(region=="Northern CC",period1==3,period2==2, survey=='TUMI', Index=="ONI"|Index=="NPGO")%>%summarise(mean=mean(ov))
overlap.phe%>%filter(region=="Central CC",period1==3,period2==2, survey=='TUMI'|survey=="STI", Index=="ONI"|Index=="NPGO"|Index=="PDO")%>%summarise(mean=mean(ov))
overlap.phe%>%filter(region=="Southern CC",period1==3,period2==2, survey=='TUMI'|survey=="STI", Index=="ONI"|Index=="NPGO"|Index=="PDO")%>%summarise(mean=mean(ov))

overlap.phe%>%filter(region=="Central CC",period1==3,period2==2, survey=="STI", Index=="NPH")%>%summarise(mean=mean(ov))


overlap.phe%>%filter(ov>0.7&period1==3&period2==2)


overlap.phe%>%filter(period1==3&period2==2)%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe%>%filter(period1==3&period2==2)%>%
  group_by(region)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe%>%filter(period1==3&period2==2)%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe%>%filter(period1==3&period2==2)%>%
  group_by(Index, region)%>%
  summarise(mean=mean(ov), sd=sd(ov))
##### WINTER Phenology Model Runs #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Winter"
eras <- data.frame(era.region2=seq(1,9), era.region=seq(4,12))
data<- climate_dat%>%filter(season=="Winter"&region!='GoA'&Year_lag<2023)%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()%>%
  left_join(eras)

data.phe.lmW<-data%>%filter(season=="Winter")%>%
  dplyr::select(Year_lag, season, period,region,stand_tumi,stand_lusi,stand_sti,stand_bakun_seasonally,
                seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, period,region,stand_tumi,stand_lusi,stand_sti,stand_bakun_seasonally), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  pivot_longer(!c(Year_lag, season,period,region, Index_Value, Index_Name), 
               names_to = "Up_Name", values_to = "Up_Value")%>%
  distinct()
ggplot(data = data.phe.lm%>%filter(region=="Northern CC"), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~Up_Name, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
###### STI #######
STIW<-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_sti, data$era.region2)
  STIW <-rbind(STIW,scaled.anom)
}
STIW<-STIW%>%mutate(survey="STI",era.region2=period,
                    Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                      ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region2,period)%>%distinct())%>%
  mutate(Season='Winter', lag=0, era.region=era.region2)%>%
  mutate(region=fct_relevel(region,c("Northern CC","Central CC","Southern CC")))


ggplot(STIW, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_wrap(region~Index, ncol = 3, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

postplot(STIW, STIW$beta)
postplot(STIW,STIW$alpha)


###### LUSI #######
LUSIW<-NULL
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_lusi, data$era.region2)
  LUSIW <-rbind(LUSIW,scaled.anom)
}
LUSIW<-LUSIW%>%mutate(survey="LUSI",era.region2=period,
                      Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                        ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region2,period)%>%distinct())%>%
  mutate(Season='Winter', lag=0, era.region=era.region2)%>%
  mutate(region=fct_relevel(region,c("Northern CC","Central CC","Southern CC")))


ggplot(LUSIW, aes(x = beta, fill = as.factor(period), group=as.factor(era.region))) +
  theme_bw() +
  facet_wrap(region~Index, ncol = 4, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")


###### TUMI #######
TUMIW<-NULL
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_tumi, data$era.region2)
  TUMIW <-rbind(TUMIW,scaled.anom)
}
TUMIW<-TUMIW%>%mutate(survey="TUMI",era.region2=period,
                      Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                        ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region2,period)%>%distinct())%>%
  mutate(Season='Winter', lag=0, era.region=era.region2)%>%
  mutate(region=fct_relevel(region,c("Northern CC","Central CC","Southern CC")))


ggplot(TUMIW, aes(x = beta, fill = as.factor(period), group=as.factor(era.region))) +
  theme_bw() +
  facet_wrap(region~Index, ncol = 4, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

#####BIOLOGICAL VERSUS UPWELLING ####
###### Spring #####


upwelling_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))%>%
  filter(season=="Spring"&region!='GoA'&Year_lag<2023)%>%
  
  mutate(region=ifelse(region=='Southern CC','SCC',
                       ifelse(region=="Northern CC", 'NCC','CCC')))%>%
  dplyr::select(Year_lag,  period,region, season, stand_bakun_seasonally,
                stand_tumi,stand_lusi,stand_sti)

dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
  mutate(region=ifelse(trend=='CALCOFI','SCC',
                       ifelse(trend=="RREAS", 'CCC',0)))%>%
  dplyr::select(Year_lag, season,region, trend,  estimate)%>%
  filter(region!=0)%>%
  distinct()%>%
  filter(season=="Spring")%>%
  merge(upwelling_dat%>%filter(season=="Spring"))%>%
  distinct()

climate_dat_RREAS<-filter(dfa, trend=="RREAS", period==2|period==3)
climate_dat_CALCOFI<-filter(dfa, trend=="CALCOFI")


climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))%>%
  filter(season=='Summer')%>%
  dplyr::select(Year_lag, seasonal_copepod_northern,seasonal_copepod_southern)%>%
  pivot_longer(!Year_lag, names_to = "trend", values_to = "estimate")%>%
  mutate(region="NCC")%>%
  merge(upwelling_dat%>%filter(season=="Spring"))
# mutate(region="NCC", trend=ifelse(trend=='seasonal_copepod_northern', 'N. Copepod', 'S. Copepod'))
climate_dat_cop_northern<-filter(climate_dat_cop, trend=="seasonal_copepod_northern")
climate_dat_cop_southern<-filter(climate_dat_cop, trend=="seasonal_copepod_southern")

bioup_NCC_northern <-NULL
columns<-c(which(colnames(climate_dat_cop_northern) == "stand_bakun_seasonally"),
           which(colnames(climate_dat_cop_northern) == "stand_lusi"),
           which(colnames(climate_dat_cop_northern) == "stand_tumi"),
           which(colnames(climate_dat_cop_northern) == "stand_sti"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_northern, columns[i],  climate_dat_cop_northern$estimate, climate_dat_cop_northern$period)
  bioup_NCC_northern <-rbind(bioup_NCC_northern,scaled.anom)
}
bioup_NCC_northern_spring<-mutate(bioup_NCC_northern, period=ifelse(period==1,2,3),
                                  Survey="N. Copepod",region="NCC",
                                  Index=ifelse(Index=="stand_bakun_seasonally", "Upwelling", 
                                               ifelse(Index=="stand_sti", "STI",
                                                      ifelse(Index=="stand_tumi", "TUMI","LUSI"))))


#bioup_NCC_northern_full <-NULL
#columns<-c(which(colnames(climate_dat_cop_northern) == "stand_bakun_seasonally"))
#for(i in 1:length(columns)){
#  bayeslinmod(climate_dat_cop_northern, columns[i],  climate_dat_cop_northern$estimate, climate_dat_cop_northern$trend)
#  bioup_NCC_northern_full <-rbind(bioup_NCC_northern_full,scaled.anom)
#}
#bioup_NCC_northern_spring_full<-mutate(bioup_NCC_northern_full, period=4,
#                                       Survey="N. Copepod",region="NCC")


bioup_NCC_southern <-NULL
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_southern, columns[i],  climate_dat_cop_southern$estimate, climate_dat_cop_southern$period)
  bioup_NCC_southern <-rbind(bioup_NCC_southern,scaled.anom)
}
bioup_NCC_southern_spring<-mutate(bioup_NCC_southern, period=ifelse(period==1,2,3),
                                  Survey="S. Copepod",region="NCC",
                                  Index=ifelse(Index=="stand_bakun_seasonally", "Upwelling", 
                                               ifelse(Index=="stand_sti", "STI",
                                                      ifelse(Index=="stand_tumi", "TUMI","LUSI"))))

#bioup_NCC_southern_full <-NULL
#columns<-c(which(colnames(climate_dat_cop_southern) == "stand_bakun_seasonally"))
#for(i in 1:length(columns)){
#  bayeslinmod(climate_dat_cop_southern, columns[i],  climate_dat_cop_southern$estimate, climate_dat_cop_southern$trend)
#  bioup_NCC_southern_full <-rbind(bioup_NCC_southern_full,scaled.anom)
#}
#bioup_NCC_southern_spring_full<-mutate(bioup_NCC_southern_full, period=4,
#                           Survey="S. Copepod",region="NCC")

bioup_SCC <-NULL
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_CALCOFI, columns[i],  climate_dat_CALCOFI$estimate, climate_dat_CALCOFI$period)
  bioup_SCC <-rbind( bioup_SCC,scaled.anom)
}
bioup_SCC_spring<-mutate(bioup_SCC,
                         Survey="CALCOFI",region="SCC",
                         Index=ifelse(Index=="stand_bakun_seasonally", "Upwelling", 
                                      ifelse(Index=="stand_sti", "STI",
                                             ifelse(Index=="stand_tumi", "TUMI","LUSI"))))

#bioup_SCC <-NULL
#columns<-c(which(colnames(climate_dat_CALCOFI) == "stand_bakun_seasonally"))
#for(i in 1:length(columns)){
#  bayeslinmod(climate_dat_CALCOFI, columns[i],  climate_dat_CALCOFI$estimate, climate_dat_CALCOFI$trend)
#   bioup_SCC <-rbind( bioup_SCC,scaled.anom)
#}
#bioup_SCC_spring_full<-mutate(bioup_SCC,Survey="CALCOFI",region="SCC", period=4)

bioup_CCC <-NULL
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_RREAS, columns[i],  climate_dat_RREAS$estimate, climate_dat_RREAS$period)
  bioup_CCC <-rbind(bioup_CCC,scaled.anom)
}

bioup_CCC_spring<-mutate(bioup_CCC, period=ifelse(period==1,2,3),Survey="RREAS",region="CCC",
                         Index=ifelse(Index=="stand_bakun_seasonally", "Upwelling", 
                                      ifelse(Index=="stand_sti", "STI",
                                             ifelse(Index=="stand_tumi", "TUMI","LUSI"))))

#bioup_CCC_full <-NULL
#columns<-c(which(colnames(climate_dat_RREAS) == "stand_bakun_seasonally"))
#for(i in 1:length(columns)){
#  bayeslinmod(climate_dat_RREAS, columns[i],  climate_dat_RREAS$estimate, climate_dat_RREAS$trend)
#  bioup_CCC_full <-rbind(bioup_CCC_full,scaled.anom)
#}

#bioup_CCC_spring_full<-mutate(bioup_CCC_full, period=4,Survey="RREAS",region="CCC")


bioup_spring<- bind_rows(bioup_CCC_spring,bioup_SCC_spring,
                         bioup_NCC_southern_spring,
                         bioup_NCC_northern_spring)%>%
  mutate(season="Spring", lag=0)%>%
  mutate(Survey= ifelse(Survey=='RREAS', "RREAS (CCC)",
                        ifelse(Survey=="CALCOFI", "CALCOFI (SCC)",
                               ifelse(Survey=="N. Copepod", "N. Copepod (NCC)",
                                      ifelse(Survey=="S. Copepod","S. Copepod (NCC)",Survey)))))%>%
  distinct()%>%
  mutate(Survey= fct_relevel(Survey, "N. Copepod (NCC)", 
                             "S. Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))


BioUpbeta<- ggplot(bioup_spring%>%filter(Index!="Upwelling"), aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~Survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3], 'grey'),name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS20_BioUpbeta.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
BioUpbeta
dev.off()

BioUpalpha<-ggplot(bioup_spring%>%filter(Index!="Upwelling"), aes(x = alpha, fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~Survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3], 'grey'),name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")


pdf(file = "Output/Supplemental/FigureS21_BioUpalpha.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
BioUpalpha
dev.off()

#### Full Data Results ####
upwelling
STI<-STI%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,
                         region,period,Season,lag)%>%mutate(period=as.numeric(as.factor(period)))
LUSI<-LUSI%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,
                           region,period,Season,lag)%>%mutate(period=as.numeric(as.factor(period)))
TUMI<-TUMI%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,
                           region,period,Season,lag)%>%mutate(period=as.numeric(as.factor(period)))

STIW<-STIW%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,
                           region,period,Season,lag)%>%mutate(period=as.numeric(as.factor(period)))
LUSIW<-LUSIW%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,
                             region,period,Season,lag)%>%mutate(period=as.numeric(as.factor(period)))
TUMIW<-TUMIW%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,
                             region,period,Season,lag)%>%mutate(period=as.numeric(as.factor(period)))

Full_Results <- bind_rows(CALCOFI_Woffset,CALCOFI_W,CALCOFIoffset,CALCOFI,
                          RREAS_Woffset,RREAS_W,RREASoffset,RREAS,
                          #upwelling,
                          TUMI,STI,LUSI,TUMIW,STIW,LUSIW,
                          northernW%>%mutate(period=ifelse(period==1,2,ifelse(period==2,3,4))),
                          northern%>%mutate(period=ifelse(period==1,2,ifelse(period==2,3,4))), 
                          southernW%>%mutate(period=ifelse(period==1,2,ifelse(period==2,3,4))),
                          southern%>%mutate(period=ifelse(period==1,2,ifelse(period==2,3,4))))

overlap.upW<-overlap.upW%>%mutate(region=ifelse(region=='Southern CC','SCC',
                                                ifelse(region=="Northern CC", 'NCC','CCC')),
                                  season="Spring", offset=0, Survey="Upwelling")%>%
  dplyr::select(ov,period1,period2,Index,Survey,region,season,offset)

overlap.up<-overlap.up%>%mutate(region=ifelse(region=='Southern CC','SCC',
                                              ifelse(region=="Northern CC", 'NCC','CCC')),
                                season="Spring", offset=0, Survey="Upwelling")%>%
  dplyr::select(ov,period1,period2,Index,Survey,region,season,offset)

overlap_full<-bind_rows(na.omit(overlap.CALCOFI_Woffset),na.omit(overlap.CALCOFI_W),na.omit(overlap.CALCOFIoffset),na.omit(overlap.CALCOFI),
                        na.omit(overlap.RREAS_Woffset),na.omit(overlap.RREASoffset),na.omit(overlap.RREAS_W),na.omit(overlap.RREAS),
                        na.omit(overlap.northern),na.omit(overlap.northernW),
                        na.omit(overlap.southern),na.omit(overlap.southernW),
                        na.omit(overlap.up), na.omit(overlap.upW))


saveRDS(overlap_full, file = 'data/overlap_Results.rds')

saveRDS(Full_Results, file = here('data/Full_Results.rds'))

Violin_Data<-bind_rows(Full_Results,bioup_spring%>%rename(survey=Survey))%>%
  mutate(survey=ifelse(survey=="N. Copepod (NCC)","N. Copepod",
                       ifelse(survey=="S. Copepod (NCC)","S. Copepod",
                              ifelse(survey=="Southern Copepod (NCC)","S. Copepod",survey))))
saveRDS(Violin_Data, file = here('data/Violin_Data.rds'))

climate_dat_cop%>%filter(region=="NCC")%>%
  dplyr::select(seasonal_ONI, Year_lag)%>%
  distinct()

Violin_Data%>%filter(Season=="Winter"&survey=="TUMI")

###### supplment plots spring ####
dat_plot<-climate_dat_cop%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, season, period,seasonal_copepod_northern,seasonal_copepod_southern, seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI, 
         'S. Copepod (NCC)'=seasonal_copepod_southern,  'N. Copepod (NCC)'=seasonal_copepod_northern)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, 'N. Copepod (NCC)','S. Copepod (NCC)', period), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  pivot_longer(!c(Year_lag, season,period, Index_Value, Index_Name), 
               names_to = "trend", values_to = "estimate")%>%
  distinct()%>%
  bind_rows(dat%>%filter(season=="Spring")%>%
              dplyr::select(Year_lag, season, period,estimate, trend,seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
              rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
              distinct()%>%
              pivot_longer(!c(Year_lag, season, period,estimate, trend), 
                           names_to = "Index_Name", values_to = "Index_Value")%>%
              filter(trend=="RREAS"|trend=="CALCOFI")%>%
              mutate(trend= ifelse(trend=='RREAS', "RREAS (CCC)",ifelse(trend=="CALCOFI", "CALCOFI (SCC)", NA)))%>%
              distinct())%>%
  mutate(trend = fct_relevel(trend, "N. Copepod (NCC)", 
                             "S. Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))


bio_dat<-ggplot(data = dat_plot%>%filter(season=="Spring"&period!=4), 
                aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
bio_dat
pdf(file = "Output/Supplemental/FigureS16_BIOlinearregression.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio_dat
dev.off()


bio<-Full_Results%>%
  mutate(survey= ifelse(survey=='RREAS', "RREAS (CCC)",
                        ifelse(survey=="CALCOFI", "CALCOFI (SCC)",
                               ifelse(survey=="N. Copepod", "N. Copepod (NCC)",
                                      ifelse(survey=="Southern Copepod (NCC)","S. Copepod (NCC)",survey)))))%>%
  distinct()%>%
  mutate(survey = fct_relevel(survey, "N. Copepod (NCC)", 
                              "S. Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))
bio.beta<-ggplot(data = bio%>%
                   filter(Season=="Spring"&period!=4&lag==0)%>%
                   filter(survey=="RREAS (CCC)"|survey=="CALCOFI (SCC)"|survey=="N. Copepod (NCC)"|survey=="S. Copepod (NCC)"),
                 aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  #  xlim(c(-1,1))+
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS17_Biobeta.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio.beta
dev.off()

bio.alpha<-ggplot(data = bio%>%
                    filter(Season=="Spring"&period!=4&lag==0)%>%
                    filter(survey=="RREAS (CCC)"|survey=="CALCOFI (SCC)"|survey=="N. Copepod (NCC)"|survey=="S. Copepod (NCC)"),
                  aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  #  xlim(c(-1,1))+
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS18_Bioalpha.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio.alpha
dev.off()


###### supplment plots winter ####
dat_plot<-climate_dat_cop%>%filter(season=="Winter")%>%
  dplyr::select(Year_lag, season, period,seasonal_copepod_northern,seasonal_copepod_southern, seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI, 
         'S. Copepod (NCC)'=seasonal_copepod_southern,  'N. Copepod (NCC)'=seasonal_copepod_northern)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, 'N. Copepod (NCC)','S. Copepod (NCC)', period), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  pivot_longer(!c(Year_lag, season,period, Index_Value, Index_Name), 
               names_to = "trend", values_to = "estimate")%>%
  distinct()%>%
  bind_rows(dat%>%filter(season=="Winter")%>%
              dplyr::select(Year_lag, season, period,estimate, trend,seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
              rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
              distinct()%>%
              pivot_longer(!c(Year_lag, season, period,estimate, trend), 
                           names_to = "Index_Name", values_to = "Index_Value")%>%
              filter(trend=="RREAS"|trend=="CALCOFI")%>%
              mutate(trend= ifelse(trend=='RREAS', "RREAS (CCC)",ifelse(trend=="CALCOFI", "CALCOFI (SCC)", NA)))%>%
              distinct())%>%
  mutate(trend = fct_relevel(trend, "N. Copepod (NCC)", 
                             "S. Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))


bio_dat<-ggplot(data = dat_plot%>%filter(season=="Winter"&period!=4), 
                aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
bio_dat
pdf(file = "Output/Supplemental/FigureS31_BIOlinearregressionW.pdf",   # The directory you want to save the file in
    width =7, # The width of the plot in inches
    height = 5)
bio_dat
dev.off()


bio<-Full_Results%>%
  mutate(survey= ifelse(survey=='RREAS', "RREAS (CCC)",
                        ifelse(survey=="CALCOFI", "CALCOFI (SCC)",
                               ifelse(survey=="N. Copepod", "N. Copepod (NCC)",
                                      ifelse(survey=="Southern Copepod (NCC)","S. Copepod (NCC)",survey)))))%>%
  distinct()%>%
  mutate(survey = fct_relevel(survey, "N. Copepod (NCC)", 
                              "S. Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))
bio.beta<-ggplot(data = bio%>%
                   filter(Season=="Winter"&period!=4&lag==0)%>%
                   filter(survey=="RREAS (CCC)"|survey=="CALCOFI (SCC)"|survey=="N. Copepod (NCC)"|survey=="S. Copepod (NCC)"),
                 aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  #  xlim(c(-1,1))+
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS32_BiobetaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio.beta
dev.off()

bio.alpha<-ggplot(data = bio%>%
                    filter(Season=="Winter"&period!=4&lag==0)%>%
                    filter(survey=="RREAS (CCC)"|survey=="CALCOFI (SCC)"|survey=="N. Copepod (NCC)"|survey=="S. Copepod (NCC)"),
                  aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  #  xlim(c(-1,1))+
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS33_BioalphaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio.alpha
dev.off()

data<- climate_dat%>%filter(season=="Winter"&region!='GoA'&Year_lag<2023&period!=4)%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()%>%
  left_join(eras)%>%
  mutate(region=fct_relevel(region,c("Northern CC","Central CC","Southern CC")))
index.names <-c("PDO", "NPGO", "NPH", "ONI")
region.names<-unique(climate_dat$region)
period.names <- unique(climate_dat$period)

data.phe.lm<-data%>%filter(season=="Winter")%>%
  dplyr::select(Year_lag, season, period,region,stand_tumi,stand_lusi,stand_sti,stand_bakun_seasonally,
                seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, period,region,stand_tumi,stand_lusi,stand_sti,stand_bakun_seasonally), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  pivot_longer(!c(Year_lag, season,period,region, Index_Value, Index_Name), 
               names_to = "Up_Name", values_to = "Up_Value")%>%
  distinct()

sti_datW <-ggplot(data = data.phe.lm%>%filter(Up_Name=="stand_sti"&period!=4), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "STI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
sti_datW
pdf(file = "Output/Supplemental/FigureS25_STIlinearregressionW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
sti_datW
dev.off()


tumi_datW <-ggplot(data = data.phe.lm%>%filter(Up_Name=="stand_tumi"&period!=4), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "TUMI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
tumi_datW
pdf(file = "Output/Supplemental/FigureS22_TUMIlinearregressionW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
tumi_datW
dev.off()


lusi_datW <-ggplot(data = data.phe.lm%>%filter(Up_Name=="stand_lusi"&period!=4), aes(y = Up_Value, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "LUSI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
lusi_datW
pdf(file = "Output/Supplemental/FigureS28_lusilinearregressionW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
lusi_datW
dev.off()

STI_betaW <-ggplot(STIW, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS26_STIbetaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
STI_betaW
dev.off()

STI_alphaW <-ggplot(STIW, aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS27_STIalphaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
STI_alphaW
dev.off()

TUMI_betaW <-ggplot(TUMIW, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS23_TUMIbetaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
TUMI_betaW
dev.off()

TUMI_alphaW <-ggplot(TUMIW, aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS24_TUMIalphaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
TUMI_alphaW
dev.off()


LUSI_betaW <-ggplot(LUSIW, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS29_LUSIbetaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
LUSI_betaW
dev.off()

LUSI_alphaW <-ggplot(LUSIW, aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")

pdf(file = "Output/Supplemental/FigureS30_LUSIalphaW.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
LUSI_alphaW
dev.off()

#### Temporal Lag ####
dat_plot<-dat%>%filter(season=="Spring"&trend!="SEA")%>%
  dplyr::select(Year_lag, estimateoffset1,season, period,trend, seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag,estimateoffset1,trend, season, period), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  distinct()%>%
  mutate(trend= ifelse(trend=='RREAS', "RREAS (CCC)",ifelse(trend=="CALCOFI", "CALCOFI (SCC)", NA)))%>%
  distinct()

bio_dat<-ggplot(data =na.omit(dat_plot%>%filter(season=="Spring")), 
                aes(y = estimateoffset1, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period)), alpha = 0.5) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  #geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("1-year Lag")
bio_dat
pdf(file = "Output/Supplemental/FigureS33_BIOlinearregressionLAG.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio_dat
dev.off()


bio<-Full_Results%>%
  mutate(survey= ifelse(survey=='RREAS', "RREAS (CCC)",
                        ifelse(survey=="CALCOFI", "CALCOFI (SCC)",
                               ifelse(survey=="N. Copepod", "N. Copepod (NCC)",
                                      ifelse(survey=="Southern Copepod (NCC)","S. Copepod (NCC)",survey)))))%>%
  distinct()%>%
  mutate(survey = fct_relevel(survey, "N. Copepod (NCC)", 
                              "S. Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))
bio.beta<-ggplot(data = bio%>%
                   filter(Season=="Spring"&period!=4&lag==1)%>%
                   filter(survey=="RREAS (CCC)"|survey=="CALCOFI (SCC)"|survey=="N. Copepod (NCC)"|survey=="S. Copepod (NCC)"),
                 aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  #  xlim(c(-1,1))+
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS34_Biobeta.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio.beta
dev.off()

bio.alpha<-ggplot(data = bio%>%
                    filter(Season=="Spring"&period!=4&lag==1)%>%
                    filter(survey=="RREAS (CCC)"|survey=="CALCOFI (SCC)"|survey=="N. Copepod (NCC)"|survey=="S. Copepod (NCC)"),
                  aes(x = alpha, fill = as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  #  xlim(c(-1,1))+
  scale_fill_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")
pdf(file = "Output/Supplemental/FigureS34_Bioalpha.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
bio.alpha
dev.off()
