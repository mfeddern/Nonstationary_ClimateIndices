library(dplyr)
library(nord)
library(tidyr)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(MCMCvis)
library(HDInterval)
library(reshape2)
library(tidyverse)
library(dplyr)
library(rstan)
library(here)
library(PNWColors)
library(corrplot)
library(maps)       #basic mapping functions and some data
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(rgdal)
library(colorspace)
library(PBSmapping) #powerful mapping functions developed by Pacific Biological Station
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#### Importing Upwelling ####

here::i_am("Output/BayesianLinearModels.Rmd")


bakundat <-read.csv(here('data/physical/Bakun/erdUI246hr_d68d_e898_8529.csv'))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI276hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI306hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI336hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI366hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI396hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI426hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI456hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI486hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI516hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI546hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI576hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI606hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI616hr_d68d_e898_8529.csv')))

#converting dates
bakun <- bakundat%>%
  add_column('Year'=as.numeric(format(as.Date(bakundat$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(bakundat$time),"%m")))%>%
  add_column("Day"=as.numeric(format(as.Date(bakundat$time),"%d")))%>%
  add_column("YearDay"=as.numeric(yday(format(as.Date(bakundat$time)))))

#selecting stations and assigning regions
bakun_region <- bakun%>%
  filter(station_id=='36N'|station_id=='39N'|station_id=='42N')%>%
  mutate(region="Central CC")%>%
  bind_rows(bakun%>%
              filter(station_id=='33N'|station_id=='30N'|station_id=='27N'|station_id=='24N')%>%
              mutate(region="Southern CC"))%>%
  bind_rows(bakun%>%
              filter(station_id=='45N'|station_id=='48N')%>%
              mutate(region="Northern CC"))%>%
  bind_rows(bakun%>%
              filter(station_id=='51N'|station_id=='54N'|station_id=='57N'|station_id=='60N')%>%
              mutate(region="GoA"))

#selecting months to assign a season and adding lag for nov/dec for winter
bakun_season <- bakun_region%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(bakun_region%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(bakun_region%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

#creating an annual standardization
bakun_summ_annual <- bakun_season%>%
  group_by(region, Year_lag) %>%
  summarise(annual_mean = mean(na.omit(upwelling_index)))%>%
  ungroup()%>%
  group_by(region)%>%
  mutate(mean = mean(annual_mean), sd=sd(annual_mean))%>%
  ungroup()%>%
  mutate(stand_bakun_annual = (annual_mean-mean)/sd)%>%
  select(Year_lag, region,stand_bakun_annual)

#creating a seasonal standardization
bakun_summ_seasonal <- bakun_season%>%
  group_by(region, Year_lag,  season) %>%
  summarise(season_mean = mean(na.omit(upwelling_index)))%>%
  ungroup()%>%
  group_by(region, season)%>%
  mutate(mean = mean(season_mean), sd=sd(season_mean))%>%
  ungroup()%>%
  mutate(stand_bakun_seasonally = (season_mean-mean)/sd)%>%
  select(Year_lag, season, region,stand_bakun_seasonally)

#creating a monthly standardization (note this uses YearLag)
bakun_summ<- bakun_season%>%
  group_by(Month, region, Year_lag) %>%
  summarise(monthly_mean = mean(na.omit(upwelling_index)))%>%
  mutate(stand_bakun_monthly = (monthly_mean-mean(na.omit(monthly_mean)))/sd(na.omit(monthly_mean)))%>%
  mutate(season=if_else(Month == 11|Month ==12|Month ==1|Month ==2|Month ==3, "Winter",
                        if_else(Month ==4|Month ==5|Month ==6, "Spring",
                                if_else(Month ==7|Month ==8, "Summer", "Autumn"))))%>%
  left_join(bakun_summ_seasonal)%>%
  left_join(bakun_summ_annual)


bakun_time <-bakun_summ%>%
  filter(Year_lag>1963 & Year_lag<1989)%>%
  mutate(period='1')%>%
  bind_rows(bakun_summ%>%
              filter(Year_lag>=1989 & Year_lag<2014)%>%
              mutate(period='2'))%>%
  bind_rows(bakun_summ%>%
              filter(Year_lag>2013)%>%
              mutate(period='3'))%>%
  mutate(era.region = if_else(region == "GoA"&period == 1, 1, 
                        if_else(region == "GoA"&period == 2, 2, 
                           if_else(region == "GoA"&period == 3, 3,
                              if_else(region == "Northern CC"&period == 1, 4,
                                  if_else(region == "Northern CC"&period == 2, 5,
                                      if_else(region == "Northern CC"&period == 3, 6,
                                           if_else(region == "Central CC"&period == 1, 7,
                                                if_else(region == "Central CC"&period == 2, 8,
                                                     if_else(region == "Central CC"&period == 3, 9, 
                                                          if_else(region == "Southern CC"&period == 1, 10,
                                                              if_else(region == "Southern CC"&period == 2, 11,
                                                                   if_else(region == "Southern CC"&period == 3, 12,
                                                                        13)))))))))))))

bakun_monthly <- bakun_time%>%
  select(Month, region, Year_lag, monthly_mean, stand_bakun_monthly)
bakun_time%>%filter(Year_lag==2023&season=='Spring')


##### Generating Cumulative Upwelling Plots #####
bakun_daily <- bakun_region%>%
  group_by( station_id,Year, YearDay,region)%>%
  summarise(upwelling_index_sum = mean(na.omit(upwelling_index)))


bakun_cum <- bakun_daily%>%
  group_by(station_id,Year,region)%>%
  reframe(upwelling_index_cum = cumsum(upwelling_index_sum))%>%
  add_column(YearDay=bakun_daily$YearDay)
period2=data.frame(period2=c('1967 - 1988', '1989 - 2013', '2014 - 2022'), period=c('1','2','3'))

bakun_time <-bakun_cum%>%
  filter(Year>1963 & Year<1989)%>%
  mutate(period='1')%>%
  bind_rows(bakun_cum%>%
              filter(Year>1989 & Year<2014)%>%
              mutate(period='2'))%>%
  bind_rows(bakun_cum%>%
              filter(Year>2013)%>%
              mutate(period='3'))%>%
  left_join(period2)


ggplot(data = bakun_time, aes(x = YearDay, y = upwelling_index_cum, col=Year)) +
  facet_wrap(.~station_id, ncol = 3, scales='free') +
  geom_line(size=0.75,aes(group=as.numeric(Year))) +
  # geom_smooth(col='red',aes(group=period)) +
  #scale_x_continuous(name = "Day") +
  #scale_color_manual(values=col3)+
  theme_bw()

bakun_time[which.min(bakun_time$upwelling_index_cum),]

CUM<-ggplot(data = bakun_time%>%filter(station_id!='60N'&
                                         station_id!='24N'&
                                         station_id!='27N'&
                                         station_id!='30N'&
                                         station_id!='51N'&
                                         station_id!='54N'&
                                         station_id!='57N'), aes(x = YearDay, y = upwelling_index_cum)) +
  facet_wrap(station_id~region, scales='free', ncol=3) +
  #  ggtitle(paste(bakun_time$station_id, "/n", bakun_time$region))+
  geom_line(size=0.75,aes(group=as.numeric(Year)),col='grey') +
  geom_smooth(aes(group=period2, colour=period2)) +
  #scale_x_continuous(name = "Day") +
  scale_color_manual(values=col[1:3], name="Period")+
  xlab("Julian Day")+
  ylab(expression('CUI ' ~ m^3 ~ '/ s / 100m'))+
  theme_bw()
CUM

pdf(file = "Output/Figures/CumulativeUpwelling.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
CUM
#thetaplot
dev.off()

#### Importing Phenology Data ####

TUMIdat <-read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_TUMI_48N.csv'))%>%
  mutate(station_id='48N')%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_TUMI_45N.csv'))%>%mutate(station_id='45N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_TUMI_42N.csv'))%>%mutate(station_id='42N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_TUMI_39N.csv'))%>%mutate(station_id='39N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_TUMI_36N.csv'))%>%mutate(station_id='36N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_TUMI_33N.csv'))%>%mutate(station_id='33N'))

TUMIdat <- TUMIdat%>%
  add_column('Year'=as.numeric(format(as.Date(TUMIdat$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(TUMIdat$time),"%m")))%>%
  add_column("Day"=as.numeric(format(as.Date(TUMIdat$time),"%d")))

corrplot.mixed(cor(pivot_wider(TUMIdat%>%select(station_id,tumi, Year),names_from = station_id, values_from = tumi)%>%
               select(!Year)),
         order="alphabet", title="TUMI", mar=c(0,0,3,0))

STIdat <-read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_STI_48N.csv'))%>%
  mutate(station_id='48N')%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_STI_45N.csv'))%>%mutate(station_id='45N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_STI_42N.csv'))%>%mutate(station_id='42N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_STI_39N.csv'))%>%mutate(station_id='39N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_STI_36N.csv'))%>%mutate(station_id='36N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_STI_33N.csv'))%>%mutate(station_id='33N'))

STIdat <- STIdat%>%
  add_column('Year'=as.numeric(format(as.Date(STIdat$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(STIdat$time),"%m")))%>%
  add_column("Day"=as.numeric(format(as.Date(STIdat$time),"%d")))
STIdat%>%filter(station_id=='48N')

corrplot.mixed(cor(na.omit(scale(pivot_wider(STIdat%>%select(station_id,sti, Year),names_from = station_id, values_from = sti)%>%
                     select(!Year), center=TRUE, scale=TRUE))),
               order="alphabet", title="STI", mar=c(0,0,3,0))



LUSIdat <-read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_LUSI_48N.csv'))%>%
  mutate(station_id='48N')%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_LUSI_45N.csv'))%>%mutate(station_id='45N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_LUSI_42N.csv'))%>%mutate(station_id='42N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_LUSI_39N.csv'))%>%mutate(station_id='39N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_LUSI_36N.csv'))%>%mutate(station_id='36N'))%>%
  bind_rows(read.csv(here('data/physical/Upwelling_Phenology/cciea_OC_LUSI_33N.csv'))%>%mutate(station_id='33N'))

LUSIdat <- LUSIdat%>%
  add_column('Year'=as.numeric(format(as.Date(LUSIdat$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(LUSIdat$time),"%m")))%>%
  add_column("Day"=as.numeric(format(as.Date(LUSIdat$time),"%d")))

corrplot.mixed(cor(pivot_wider(LUSIdat%>%select(station_id,lusi, Year),names_from = station_id, values_from = lusi)%>%
               select(!Year)),
               order="alphabet", title="LUSI", mar=c(0,0,3,0))

phendat <- LUSIdat%>%left_join(STIdat)%>%
  left_join(TUMIdat)

phen_region <- phendat%>%
  filter(station_id=='36N'|station_id=='39N')%>%
  mutate(region="Central CC")%>%
  bind_rows(phendat%>%
              filter(station_id=='33N'|station_id=='30N'|station_id=='27N'|station_id=='24N')%>%
              mutate(region="Southern CC"))%>%
  bind_rows(phendat%>%
              filter(station_id=='42N'|station_id=='45N'|station_id=='48N')%>%
              mutate(region="Northern CC"))%>%
  bind_rows(phendat%>%
              filter(station_id=='51N'|station_id=='54N'|station_id=='57N'|station_id=='60N')%>%
              mutate(region="GoA"))

cor(phendat)


phen_stand <- phen_region%>%
  group_by(region, Year) %>%
  summarise(lusi_mean = mean(lusi), sti_mean=mean(sti), tumi_mean=mean(tumi))%>%
  ungroup()%>%
  group_by(region)%>%
  mutate(lusi_annual_mean = mean(lusi_mean), lusi_annual_sd=sd(lusi_mean),
         sti_annual_mean = mean(sti_mean), sti_annual_sd=sd(sti_mean),
         tumi_annual_mean = mean(tumi_mean), tumi_annual_sd=sd(tumi_mean))%>%
  ungroup()%>%
  mutate(stand_lusi = (lusi_mean-lusi_annual_mean)/lusi_annual_sd,
         stand_sti = (sti_mean-sti_annual_mean)/sti_annual_sd,
        stand_tumi = (tumi_mean-tumi_annual_mean)/tumi_annual_sd)%>%
  select(Year, region,stand_tumi,stand_sti,stand_lusi)%>%
  rename(Year_lag=Year)

##### PDO ##### 
#### Climate Indices ####

PDO <- read.csv(here('data/physical/PDO.csv'))%>%
  pivot_longer(!Year, names_to = "Month", values_to = "PDO")%>%
  mutate(Month = as.integer(factor(Month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

PDO_seasonal <-  PDO%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(PDO%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(PDO%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_PDO = mean(PDO))


PDO_annual <-  PDO%>%
  group_by(Year_lag)%>%
  summarise(annual_PDO = mean(PDO))

##### NPGO #####
NPGO <- read.csv(here('data/physical/NPGO.csv'))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

NPGO_seasonal <-  NPGO%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(NPGO%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(NPGO%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_NPGO = mean(NPGO))


NPGO_annual <-  NPGO%>%
  group_by(Year_lag)%>%
  summarise(annual_NPGO = mean(NPGO))

##### North Pacific High #####
NPH <- read.csv(here('data/physical/NPH.csv'))
NPH<-NPH%>%
  add_column('Year'=as.numeric(format(as.Date(NPH$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(NPH$time),"%m")))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

NPH_seasonal <-  NPH%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(NPH%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(NPH%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_NPH = mean(nph_area))

  
NPH_seasonal<-NPH_seasonal%>%
  group_by(season)%>%
  mutate(seasonal_NPH=(seasonal_NPH-mean(seasonal_NPH))/sd(seasonal_NPH))


  ggplot(data = NPH_seasonal%>%filter(season=="Spring"), aes(x = Year_lag, y = NPH_stand)) +
    geom_line(size=0.75)+
    geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
    theme_bw()
  
  
  
  NPH_annual <-  NPH%>%
  group_by(Year_lag)%>%
  summarise(annual_NPH = mean(nph_area))

  
##### ENSO #####

ONI <- read.csv(here('data/physical/ONI.csv'))%>%
  pivot_longer(!Year, names_to = "Month", values_to = "ONI")%>%
  mutate(Month = as.integer(factor(Month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

ONI_seasonal <-  ONI%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(ONI%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(ONI%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_ONI = mean(na.omit(ONI)))

ONI_seasonal%>%group_by(season)%>%
  summarise(mean = mean(na.omit(seasonal_ONI)), sd=sd(na.omit(seasonal_ONI)))


ONI_annual <-  ONI%>%
  group_by(Year_lag)%>%
  summarise(annual_ONI = mean(na.omit(ONI)))

##### NPI #####

NPI <- read.csv(here('data/physical/NPI_monthly.csv'))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))
NPI_seasonal<-NPI%>%
  filter(Year>1941)%>%
filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(NPI%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(NPI%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_NPI = mean(NPI))

NPI_seasonal<-NPI_seasonal%>%
  group_by(season)%>%
  mutate(NPI_stand=(seasonal_NPI-mean(seasonal_NPI))/sd(seasonal_NPI))


### Copepod Data ####
cop <- read.csv(here('data/biological/copepod_biomass_anomaly.csv'))

copepod<-cop%>%
  add_column('Year'=as.numeric(format(as.Date(cop$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(cop$time),"%m")))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))
 
copepod_annual<-copepod%>%
  group_by(Year_lag)%>%
  summarise(annual_copepod_northern = mean(northern_biomass_anomaly), 
            annual_copepod_southern = mean(southern_biomass_anomaly))

copepod_seasonal<-copepod%>%  
filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(copepod%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(copepod%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_copepod_northern = mean(northern_biomass_anomaly), 
            seasonal_copepod_southern = mean(southern_biomass_anomaly))
  
copepod_seasonal<-copepod_seasonal%>%  
  mutate(period=ifelse(Year_lag>2012,3,2))

### Compile  with Upwelling Indices to Dataframe ####
unique(bakun_time$region)
climate_dat <- bakun_time%>%
  left_join(phen_stand, by=c('region', 'Year_lag'))%>%
  merge(PDO, by=c('Month', 'Year_lag'))%>%
  # merge(PDO_annual, by=c('Year_lag'))%>%
  merge(PDO_seasonal, by=c('season', 'Year_lag'))%>%
  merge(NPH, by=c('Month', 'Year_lag'))%>%
  merge(NPH_annual, by=c('Year_lag'))%>%
  merge(NPH_seasonal, by=c('season', 'Year_lag'))%>%
  merge(ONI, by=c('Month', 'Year_lag'))%>%
  merge(ONI_annual, by=c('Year_lag'))%>%
  merge(ONI_seasonal, by=c('season', 'Year_lag'))%>%
# merge(copepod_annual, by=c('Year_lag'))%>%
#  merge(copepod_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPGO, by=c('Month', 'Year_lag'))%>%
  left_join(NPGO_annual, by=c('Year_lag'))%>%
  left_join(NPGO_seasonal, by=c('season', 'Year_lag'))

colnames(climate_dat)
climate_dat <- climate_dat%>%
  select(Year_lag, region, season, stand_bakun_seasonally,period,stand_tumi,           
         stand_sti, stand_lusi,
         era.region, seasonal_PDO, annual_NPGO, seasonal_NPGO, annual_ONI, seasonal_ONI,
         annual_NPH, seasonal_NPH)%>%
  distinct()

ggplot(climate_dat, aes(x=seasonal_PDO, y=seasonal_NPI))+
  geom_point()

saveRDS(climate_dat, file = here('data/physical/climate_dat_upwelling.rds'))


### Compile  with Upwelling Indices and Copepod into Dataframe ####
dfa.trend<-readRDS("data/biological/dfa.trends.rds")%>%
  rename(Year_lag=time)

climate_dat <- PDO_seasonal%>%
  merge(PDO_annual, by=c('Year_lag'))%>%
  #merge(PDO, by=c('Year_lag'))%>%
  #merge(NPH, by=c('Month', 'Year_lag'))%>%
  merge(NPH_annual, by=c('Year_lag'))%>%
  merge(NPH_seasonal, by=c('season', 'Year_lag'))%>%
 # merge(ONI, by=c('Month', 'Year_lag'))%>%
  merge(ONI_annual, by=c('Year_lag'))%>%
  merge(ONI_seasonal, by=c('season', 'Year_lag'))%>%
  # merge(copepod%>%select(Year_lag, period), by=c('Year_lag'))%>%
  merge(dfa.trend, by=c('Year_lag'))%>%
  #left_join(NPGO, by=c('Month', 'Year_lag'))%>%
  left_join(NPGO_annual, by=c('Year_lag'))%>%
  left_join(NPI_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPGO_seasonal, by=c('season', 'Year_lag'))%>%
  #add_column(region='NCC')%>%
  distinct()


climate_dat <- climate_dat%>%
  select(Year_lag, season,estimate, lower, upper, trend, era,
         annual_PDO, seasonal_PDO, annual_NPGO, seasonal_NPGO, annual_ONI, seasonal_ONI,
         annual_NPH, seasonal_NPH,NPI_stand)%>%
  distinct()

ggplot(climate_dat, aes(x=seasonal_PDO, y=estimate))+
  geom_point()

saveRDS(climate_dat, file = here('data/physical/climate_dat_dfa.rds'))



### Compile  with Upwelling Indices and Copepod into Dataframe ####

climate_dat <- PDO_seasonal%>%
  merge(PDO_annual, by=c('Year_lag'))%>%
  #merge(PDO, by=c('Year_lag'))%>%
  #merge(NPH, by=c('Month', 'Year_lag'))%>%
  merge(NPH_annual, by=c('Year_lag'))%>%
  merge(NPH_seasonal, by=c('season', 'Year_lag'))%>%
  # merge(ONI, by=c('Month', 'Year_lag'))%>%
  merge(ONI_annual, by=c('Year_lag'))%>%
  merge(ONI_seasonal, by=c('season', 'Year_lag'))%>%
  # merge(copepod%>%select(Year_lag, period), by=c('Year_lag'))%>%
  #merge(copepod_annual, by=c('Year_lag'))%>%
  merge(copepod_seasonal, by=c('season', 'Year_lag'))%>%
  #left_join(NPGO, by=c('Month', 'Year_lag'))%>%
  left_join(NPGO_annual, by=c('Year_lag'))%>%
  left_join(NPI_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPGO_seasonal, by=c('season', 'Year_lag'))%>%
  add_column(region='NCC')%>%
  distinct()


climate_dat <- climate_dat%>%
  select(Year_lag, region, period, season,
         annual_PDO, seasonal_PDO, annual_NPGO, seasonal_NPGO, annual_ONI, seasonal_ONI,
         annual_NPH, seasonal_NPH,NPI_stand, 
         seasonal_copepod_northern, seasonal_copepod_southern)%>%
  distinct()

ggplot(climate_dat, aes(x=seasonal_PDO, y=seasonal_NPH))+
  geom_point()


saveRDS(climate_dat, file = here('data/physical/climate_dat_cop.rds'))





### Compile  Indices to Dataframe ####

climate_dat <-PDO_seasonal%>%
  merge(ONI_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPGO_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPH_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPI, by=c('Year_lag'))

climate_dat <-climate_dat%>%  
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))

climate_dat <-climate_dat%>%  
  filter(Year_lag<1979)%>%
  mutate(period2='1')%>%
  bind_rows(climate_dat%>%
              filter(Year_lag>=1979)%>%
              mutate(period2='2'))%>%
  select(Year_lag, season, period, period2,seasonal_PDO, seasonal_NPGO, seasonal_ONI)%>%
  distinct()

ggplot(climate_dat, aes(x=seasonal_PDO, y=seasonal_ONI))+
  geom_point()

saveRDS(climate_dat, file = here('data/physical/climate_dat.rds'))


#### creating dataframe for monthly data ####


climate_dat_monthly <-PDO%>%
  merge(ONI, by=c('Month', 'Year_lag', 'Year'))%>%
  merge(NPH, by=c('Month', 'Year_lag', 'Year'))%>%
  merge(NPGO, by=c('Month', 'Year_lag', 'Year'))%>%
  rename(NPH=nph_area)%>%
  merge(bakun_monthly, by=c('Month', 'Year_lag'))

climate_dat_monthly2 <-climate_dat_monthly%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))

  
saveRDS(climate_dat_monthly2, file = here('data/physical/climate_dat_monthly.rds'))

