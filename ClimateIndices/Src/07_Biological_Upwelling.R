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
set.seed(1234)


##### Spring #####


upwelling_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))%>%
  mutate(region=ifelse(region=='Southern CC','SCC',
                ifelse(region=="Northern CC", 'NCC','CCC')))%>%
  dplyr::select(Year_lag,  period,region, season, stand_bakun_seasonally)

dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
    mutate(region=ifelse(trend=='CALCOFI','SCC',
                ifelse(trend=="RREAS", 'CCC',0)))%>%
  dplyr::select(Year_lag, region, trend,  estimate)%>%
  filter(region!=0)%>%
  distinct()%>%
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
columns<-c(which(colnames(climate_dat_cop_northern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_northern, columns[i],  climate_dat_cop_northern$estimate, climate_dat_cop_northern$period)
  bioup_NCC_northern <-rbind(bioup_NCC_northern,scaled.anom)
}
bioup_NCC_northern_spring<-mutate(bioup_NCC_northern, period=ifelse(period==1,2,3),
                                  Survey="N. Copepod",region="NCC")

bioup_NCC_northern_full <-NULL
columns<-c(which(colnames(climate_dat_cop_northern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_northern, columns[i],  climate_dat_cop_northern$estimate, climate_dat_cop_northern$trend)
  bioup_NCC_northern_full <-rbind(bioup_NCC_northern_full,scaled.anom)
}
bioup_NCC_northern_spring_full<-mutate(bioup_NCC_northern_full, period=4,
                                       Survey="N. Copepod",region="NCC")

 
bioup_NCC_southern <-NULL
columns<-c(which(colnames(climate_dat_cop_southern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_southern, columns[i],  climate_dat_cop_southern$estimate, climate_dat_cop_southern$period)
  bioup_NCC_southern <-rbind(bioup_NCC_southern,scaled.anom)
}
bioup_NCC_southern_spring<-mutate(bioup_NCC_southern, period=ifelse(period==1,2,3),
                           Survey="S. Copepod",region="NCC")

bioup_NCC_southern_full <-NULL
columns<-c(which(colnames(climate_dat_cop_southern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_southern, columns[i],  climate_dat_cop_southern$estimate, climate_dat_cop_southern$trend)
  bioup_NCC_southern_full <-rbind(bioup_NCC_southern_full,scaled.anom)
}
bioup_NCC_southern_spring_full<-mutate(bioup_NCC_southern_full, period=4,
                           Survey="S. Copepod",region="NCC")
 
bioup_SCC <-NULL
columns<-c(which(colnames(climate_dat_CALCOFI) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_CALCOFI, columns[i],  climate_dat_CALCOFI$estimate, climate_dat_CALCOFI$period)
   bioup_SCC <-rbind( bioup_SCC,scaled.anom)
}
bioup_SCC_spring<-mutate(bioup_SCC,
                           Survey="CALCOFI",region="SCC")

bioup_SCC <-NULL
columns<-c(which(colnames(climate_dat_CALCOFI) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_CALCOFI, columns[i],  climate_dat_CALCOFI$estimate, climate_dat_CALCOFI$trend)
   bioup_SCC <-rbind( bioup_SCC,scaled.anom)
}
bioup_SCC_spring_full<-mutate(bioup_SCC,
                           Survey="CALCOFI",region="SCC", period=4)

bioup_CCC <-NULL
columns<-c(which(colnames(climate_dat_RREAS) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_RREAS, columns[i],  climate_dat_RREAS$estimate, climate_dat_RREAS$period)
  bioup_CCC <-rbind(bioup_CCC,scaled.anom)
}

bioup_CCC_spring<-mutate(bioup_CCC, period=ifelse(period==1,2,3),Survey="RREAS",region="CCC")

bioup_CCC_full <-NULL
columns<-c(which(colnames(climate_dat_RREAS) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_RREAS, columns[i],  climate_dat_RREAS$estimate, climate_dat_RREAS$trend)
  bioup_CCC_full <-rbind(bioup_CCC_full,scaled.anom)
}

bioup_CCC_spring_full<-mutate(bioup_CCC_full, period=4,Survey="RREAS",region="CCC")


bioup_spring<- bind_rows(bioup_CCC_spring,bioup_SCC_spring,bioup_NCC_southern_spring,bioup_NCC_northern_spring,
                         bioup_CCC_spring_full,bioup_SCC_spring_full,
                         bioup_NCC_southern_spring_full,bioup_NCC_northern_spring_full)


 ggplot(bioup_spring, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(.~Survey, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1],col[2], col[3], 'grey')) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")

 ##### Winter #####


upwelling_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))%>%
  mutate(region=ifelse(region=='Southern CC','SCC',
                ifelse(region=="Northern CC", 'NCC','CCC')))%>%
  dplyr::select(Year_lag,  period,region, season, stand_bakun_seasonally)

dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
    mutate(region=ifelse(trend=='CALCOFI','SCC',
                ifelse(trend=="RREAS", 'CCC',0)))%>%
  dplyr::select(Year_lag, region, trend,  estimate)%>%
  filter(region!=0)%>%
  distinct()%>%
  merge(upwelling_dat%>%filter(season=="Winter"))%>%
  distinct()
climate_dat_RREAS<-filter(dfa, trend=="RREAS", period==2|period==3)
climate_dat_CALCOFI<-filter(dfa, trend=="CALCOFI")


bioup_NCC_northern <-NULL
columns<-c(which(colnames(climate_dat_cop_northern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_northern, columns[i],  climate_dat_cop_northern$estimate, climate_dat_cop_northern$period)
  bioup_NCC_northern <-rbind(bioup_NCC_northern,scaled.anom)
}
bioup_NCC_northern_winter<-mutate(bioup_NCC_northern, period=ifelse(period==1,2,3),
                                  Survey="N. Copepod",region="NCC")

bioup_NCC_northern_full <-NULL
columns<-c(which(colnames(climate_dat_cop_northern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_northern, columns[i],  climate_dat_cop_northern$estimate, climate_dat_cop_northern$trend)
  bioup_NCC_northern_full <-rbind(bioup_NCC_northern_full,scaled.anom)
}
bioup_NCC_northern_winter_full<-mutate(bioup_NCC_northern_full, period=4,
                                       Survey="N. Copepod",region="NCC")

 
bioup_NCC_southern <-NULL
columns<-c(which(colnames(climate_dat_cop_southern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_southern, columns[i],  climate_dat_cop_southern$estimate, climate_dat_cop_southern$period)
  bioup_NCC_southern <-rbind(bioup_NCC_southern,scaled.anom)
}
bioup_NCC_southern_winter<-mutate(bioup_NCC_southern, period=ifelse(period==1,2,3),
                           Survey="S. Copepod",region="NCC")

bioup_NCC_southern_full <-NULL
columns<-c(which(colnames(climate_dat_cop_southern) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_cop_southern, columns[i],  climate_dat_cop_southern$estimate, climate_dat_cop_southern$trend)
  bioup_NCC_southern_full <-rbind(bioup_NCC_southern_full,scaled.anom)
}
bioup_NCC_southern_winter_full<-mutate(bioup_NCC_southern_full, period=4,
                           Survey="S. Copepod",region="NCC")
 
bioup_SCC <-NULL
columns<-c(which(colnames(climate_dat_CALCOFI) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_CALCOFI, columns[i],  climate_dat_CALCOFI$estimate, climate_dat_CALCOFI$period)
   bioup_SCC <-rbind( bioup_SCC,scaled.anom)
}
bioup_SCC_winter<-mutate(bioup_SCC,
                           Survey="CALCOFI",region="SCC")

bioup_SCC <-NULL
columns<-c(which(colnames(climate_dat_CALCOFI) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_CALCOFI, columns[i],  climate_dat_CALCOFI$estimate, climate_dat_CALCOFI$trend)
   bioup_SCC <-rbind( bioup_SCC,scaled.anom)
}
bioup_SCC_winter_full<-mutate(bioup_SCC,
                           Survey="CALCOFI",region="SCC", period=4)

bioup_CCC <-NULL
columns<-c(which(colnames(climate_dat_RREAS) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_RREAS, columns[i],  climate_dat_RREAS$estimate, climate_dat_RREAS$period)
  bioup_CCC <-rbind(bioup_CCC,scaled.anom)
}

bioup_CCC_winter<-mutate(bioup_CCC, period=ifelse(period==1,2,3),Survey="RREAS",region="CCC")

bioup_CCC_full <-NULL
columns<-c(which(colnames(climate_dat_RREAS) == "stand_bakun_seasonally"))
for(i in 1:length(columns)){
  bayeslinmod(climate_dat_RREAS, columns[i],  climate_dat_RREAS$estimate, climate_dat_RREAS$trend)
  bioup_CCC_full <-rbind(bioup_CCC_full,scaled.anom)
}

bioup_CCC_winter_full<-mutate(bioup_CCC_full, period=4,Survey="RREAS",region="CCC")


bioup_winter<- bind_rows(bioup_CCC_winter,bioup_SCC_winter,bioup_NCC_southern_winter,bioup_NCC_northern_winter,
                         bioup_CCC_winter_full,bioup_SCC_winter_full,
                         bioup_NCC_southern_winter_full,bioup_NCC_northern_winter_full)


 ggplot(bioup_winter, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(.~Survey, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1],col[2], col[3], 'grey')) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")
 
 bioup<- bind_rows(bioup_winter%>%mutate(Season='Winter'),
                   bioup_spring%>%mutate(Season='Spring'))%>%
   mutate(lag=0,Index="Upwelling")%>%
   rename(survey=Survey)%>%
   dplyr::select(alpha,beta,period,Index,yfirst,ylast, survey,region,Season,lag)


 ggplot(bioup, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(Season~survey, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1],col[2], col[3], 'grey')) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")
 
 ggplot(bioup%>%dplyr::filter(period!=4), aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(Season~survey, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1],col[2], col[3], 'grey')) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")
 
 
 