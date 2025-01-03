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

climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))
dat <- dfa%>%dplyr::select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))
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

leverage_CALCOFI<-lev(data,"estimate",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="CALCOFI",
          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
         ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 0)



#running models with 1-year offset
#### RREAS ####
######Spring #####

data <- dat%>%filter(season=="Spring", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)
coeflev <- c("alpha", "beta", "alpha", "beta")
periodslev<- c(2,2,3,3)
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))

RREAS<- linmod(data,"estimate",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="RREAS",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="SCC",Season="Spring", lag = 0, leverage="Included")%>%
  bind_rows(linmodlev(data%>%filter(Year_lag>1988),"estimate",columns, RESULTSDF,periodslev,coeflev,1)%>%
              mutate(survey="RREAS",
                     Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                       ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
              mutate(region="CCC",Season="Spring", lag = 0, leverage="Excluded"))

#running models with 1-year offset

leverage_RREAS<-lev(data,"estimate",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="RREAS",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="CCC",Season="Spring", lag = 0)

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
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))

NCOP<- linmod(data,"seasonal_copepod_northern",columns, RESULTSDF,periodslev, coeflev)%>%
  mutate(survey="N. Copepods",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="NCC",Season="Spring", lag = 0, leverage="Included")%>%
  bind_rows(linmodlev(data%>%filter(Year_lag>1988),"seasonal_copepod_northern",columns, RESULTSDF,periodslev,coeflev,1)%>%
              mutate(survey="N. Copepods",
                     Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                       ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
              mutate(region="NCC",Season="Spring", lag = 0, leverage="Excluded"))

leverage_NCOP<-lev(data,"seasonal_copepod_northern",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="N. Copepods",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="NCC",Season="Spring", lag = 0)


SCOP<- linmod(data,"seasonal_copepod_southern",columns, RESULTSDF,periodslev, coeflev)%>%
  mutate(survey="S. Copepods",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="NCC",Season="Spring", lag = 0, leverage="Included")%>%
  bind_rows(linmodlev(data%>%filter(Year_lag>1988),"seasonal_copepod_northern",columns, RESULTSDF,periodslev,coeflev,1)%>%
              mutate(survey="S. Copepods",
                     Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                       ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
              mutate(region="NCC",Season="Spring", lag = 0, leverage="Excluded"))

leverage_NCOP<-lev(data,"seasonal_copepod_southern",columns, RESULTSDF,periods, coef)%>%
  mutate(survey="S. Copepods",
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(region="NCC",Season="Spring", lag = 0)

BIO<- CALCOFI%>%
  bind_rows(RREAS)%>%
  bind_rows(NCOP)%>%
  bind_rows(SCOP)
##### Phenology Model Runs #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))%>%
  mutate(region=fct_relevel(region,c("Northern CC","Central CC","Southern CC")))
season <- "Spring"
data<- climate_dat%>%filter(season=="Spring"&region!='GoA'&Year_lag<2023&period!=4)%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
regions <- unique(data$region)
responses<-c("stand_tumi","stand_sti","stand_lusi")
survey<- c("TUMI", "STI", "LUSI")
PHE<- data.frame()
regionslab<- c("SCC", "NCC", "CCC")
for(k in 1:length(responses)){
  temp<-data.frame()
  for(j in 1:length(regions)){
    results<- linmod(data%>%filter(region==regions[j]),responses[k],columns, RESULTSDF,periods, coef)%>%
      mutate(survey=survey[k],
             Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                               ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
      mutate(region=regionslab[j],Season="Spring", lag = 0, leverage="Included")%>%
      bind_rows(linmodlev(data%>%filter(region==regions[j]),"stand_tumi",columns, RESULTSDF,periods,coef,1)%>%
                  mutate(survey=survey[k],
                         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
                  mutate(region=regionslab[j],Season="Spring", lag = 0, leverage="Excluded"))
    temp<-rbind(results,temp)
  } 
  PHE<-rbind(PHE,temp)
}

##### Plots ###### 

###### re leveling #######
PHE$region <- factor(PHE$region, levels = c("NCC", "CCC", "SCC"))
BIO$region <- factor(BIO$region, levels = c("NCC", "CCC", "SCC"))

PHE$survey <- factor(PHE$survey, levels = c("STI", "TUMI", "LUSI"))
BIO$survey <- factor(BIO$survey, levels = c("N. Copepods", "S. Copepods", "RREAS", "CALCOFI"))

###### Intercept plots #######
dodge<-1.25
PHE_intercept <-ggplot(PHE, aes(x = alpha, group=as.factor(period),lty=as.factor(leverage),col = as.factor(leverage),y=Index, shape=as.factor(period))) +
  theme_bw() +
  facet_grid(region~survey) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey"))+
  geom_errorbarh(aes(xmin=alpha-2*se_alpha, xmax=alpha+2*se_alpha),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "",
       y = "Climate Index")+
  theme(legend.position = "none")

BIO_intercept<-ggplot(BIO, aes(x = alpha, group=as.factor(period),lty=as.factor(leverage),col = as.factor(leverage),y=Index, shape=as.factor(period))) +
  theme_bw() +
  geom_errorbarh(aes(xmin=alpha-2*se_alpha, xmax=alpha+2*se_alpha),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  facet_wrap(~survey,ncol=1) +
  scale_colour_manual(name="Leverage Years", values=c("palevioletred", "darkgrey"))+
  scale_linetype_manual(name="Leverage Years",values=c(2,1))+
  scale_shape_manual(name="Period", values=c(15,16,17),labels=c("1967 - 1988", "1989 - 2012", "2013 - 2022"))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "",
       y = "")

pdf("Output/LeverageIntercepts.pdf", 11,6.5) 
intfig<-ggarrange(PHE_intercept,BIO_intercept, ncol = 2, labels = c("A", "B"), 
          widths=c(4,3), heights=c(2,2))
annotate_figure(intfig, 
                bottom = "Intercept")
dev.off()
##### Slope Plot #####
PHE_slope<-ggplot(PHE, aes(x = beta, group=as.factor(period),lty=as.factor(leverage),col = as.factor(leverage),y=Index, shape=as.factor(period))) +
  theme_bw() +
  facet_grid(region~survey) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey"))+
  geom_errorbarh(aes(xmin=beta-2*se_beta, xmax=beta+2*se_beta),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "",
       y = "Climate Index")+
  theme(legend.position = "none")

BIO_slope<-ggplot(BIO, aes(x = beta, group=as.factor(period),lty=as.factor(leverage),col = as.factor(leverage),y=Index, shape=as.factor(period))) +
  theme_bw() +
  geom_errorbarh(aes(xmin=beta-2*se_beta, xmax=beta+2*se_beta),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  geom_point(alpha = 0.79,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  facet_wrap(~survey,ncol=1) +
  scale_colour_manual(name="Leverage Years", values=c("palevioletred", "darkgrey"))+
  scale_linetype_manual(name="Leverage Years",values=c(2,1))+
  scale_shape_manual(name="Period", values=c(15,16,17),labels=c("1967 - 1988", "1989 - 2012", "2013 - 2022"))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "",
       y = "")

pdf("Output/LeverageSlopes.pdf", 11,6.5) 
slopefig<-ggarrange(PHE_slope,BIO_slope, ncol = 2, labels = c("A", "B"), 
          widths=c(4,3), heights=c(2,2))
annotate_figure(slopefig, 
                bottom = "Slope")
dev.off()
