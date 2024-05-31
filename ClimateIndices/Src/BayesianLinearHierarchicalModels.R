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
warmups <- 1000
total_iterations <- 10000
max_treedepth <-  10
n_chains <-  3
n_cores <- 4
adapt_delta <- 0.95
col<-pnw_palette("Sunset2",3,type="discrete")

bayeslinmod <- function(dat, ind,xyz, period) {
  N <- length(period)
  NP <- as.numeric(length(unique(period)))
  P<- as.numeric(as.factor(period))
  K <- 1
  y <-xyz
  x<-dat[ind]
  data <- list(N = N, #total number of opservations
               NP=NP, #total number of time periods
               P = P, #time period pointer vector
               K = K,#number of covariates, starting with one but can add for final model structure
               x=x, 
               y=y #response variable
  )
  
  bhfit <- stan(
    file = here::here("Src/BayesianLinearHierarchicalModels.stan"),
    data = data,
    chains = n_chains,
    warmup = warmups,
    iter = total_iterations,
    cores = n_cores,
    refresh = 250,
    control = list(max_treedepth = max_treedepth,
                   adapt_delta = adapt_delta))
  posterior<-data.frame(summary(bhfit, prob=c(0.025, 0.25,0.75, 0.975, 0.1, 0.9))$summary)
  scaled<-NULL
  alphaP<-posterior[1:NP,]%>%
    add_column(
      period=rep(seq(1,NP),1))
  
  betaP<-posterior[seq(NP+1,NP*2),]%>%
    add_column(
      period=rep(seq(1,NP),1))
  n<- 10000
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled <-rbind(scaled,cbind(tempalpha, tempbeta))
  }
  scaled.anom<<- scaled%>%
    bind_cols(period=rep(rep(seq(1,NP), each = n),1), Index=rep(colnames(dat[ind]), n*NP),
              yfirst=rep(min(dat$Year_lag), n*NP),ylast=rep(max(dat$Year_lag), n*NP))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}

postplot <- function(dat, param){
  ggplot(dat, aes(x = param, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1], col[2], col[3])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")
}

##### DFA Trend Model Runs ####
climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))
dat <- dfa%>%dplyr::select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))

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
  dplyr::select(Year_lag, seasonal_copepod_northern,seasonal_copepod_southern,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  pivot_longer(!c(Year_lag,  NPH, NPGO,PDO,ONI), 
               names_to = "trend", values_to = "estimate")%>%
  mutate(season="Spring")%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")



dat.long<- dat%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, season, estimate,trend,estimateoffset1,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

dat.lm<-data%>%
  bind_rows(dat.long)%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
  filter(trend!="SEA")%>%
  mutate(trend=ifelse(trend=="CALCOFI","CALCOFI (SCC)",
                      ifelse(trend=="RREAS", "RREAS (CCC)", 
                          ifelse(trend=="seasonal_copepod_southern","Southern Copepod (NCC)","Northern Copepod (NCC)"))))%>%
  mutate(trend = fct_relevel(trend, "Northern Copepod (NCC)", 
                              "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))

survey.lm<-ggplot(data = dat.lm%>%filter(period!=4), aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
   geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
 # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
survey.lm



pdf(file = "Output/Figures/biological.lm.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
survey.lm
dev.off()
#### CALCOFI ####
######Spring #####

data <- dat%>%filter(season=="Spring", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

CALCOFI <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  CALCOFI <-rbind(CALCOFI,scaled.anom)
}

CALCOFI<-CALCOFI%>%mutate(survey="CALCOFI",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(CALCOFI, CALCOFI$beta)
postplot(CALCOFI, CALCOFI$alpha)

CALCOFI_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFI_full <-rbind(CALCOFI_full,scaled.anom)
}

CALCOFI_full<-CALCOFI_full%>%mutate(survey="CALCOFI",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFI_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFI<- bind_rows(CALCOFI,CALCOFI_full)%>%
  mutate(region="SCC",Season="Spring", lag = 0)%>%
  filter(period!=4)

CALCOFI_beta<-ggplot(CALCOFI, aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")

index.names <- unique(CALCOFI$Index)
period.names <- unique(CALCOFI$period)

overlap.CALCOFI <- NA
  for(i in 1:4){
    temp <- CALCOFI%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
 #   ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
 ##   ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
  #  ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3)
    overlap.CALCOFI <-rbind(temp2,overlap.CALCOFI)
  }
overlap.CALCOFI<-overlap.CALCOFI%>%mutate(Survey="CALCOFI", region="SCC", season="Spring", offset=0) 

#running models with 1-year offset

CALCOFIoffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  CALCOFIoffset <-rbind(CALCOFIoffset,scaled.anom)
}
CALCOFIoffset<-CALCOFIoffset%>%mutate(survey="CALCOFI",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(CALCOFIoffset, CALCOFIoffset$beta)
postplot(CALCOFIoffset, CALCOFIoffset$alpha)

CALCOFIoffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFIoffset_full <-rbind(CALCOFIoffset_full,scaled.anom)
}

CALCOFIoffset_full<-CALCOFIoffset_full%>%mutate(survey="CALCOFI",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFIoffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFIoffset<- bind_rows(CALCOFIoffset,CALCOFIoffset_full)%>%
  mutate(region="SCC",Season="Spring", lag = 1)

index.names <- unique(CALCOFIoffset$Index)
period.names <- unique(CALCOFIoffset$period)
overlap.CALCOFIoffset <- NA
  for(i in 1:4){
    temp <- CALCOFIoffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFIoffset <-rbind(temp2,overlap.CALCOFIoffset)
  }
overlap.CALCOFIoffset<-overlap.CALCOFIoffset%>%mutate(Survey="CALCOFI", region="SCC", season="Spring", offset=1) 

###### Winter #####

data <- dat%>%filter(season=="Winter", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

CALCOFI_W <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  CALCOFI_W <-rbind(CALCOFI_W,scaled.anom)
}

CALCOFI_W<-CALCOFI_W%>%mutate(survey="CALCOFI",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(CALCOFI_W, CALCOFI_W$beta)
postplot(CALCOFI_W, CALCOFI_W$alpha)

CALCOFI_W_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFI_W_full <-rbind(CALCOFI_W_full,scaled.anom)
}

CALCOFI_W_full<-CALCOFI_W_full%>%mutate(survey="CALCOFI",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFI_W_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFI_W<- bind_rows(CALCOFI_W,CALCOFI_W_full)%>%
  mutate(region="SCC",Season="Winter", lag = 0)

index.names <- unique(CALCOFI_W$Index)
period.names <- unique(CALCOFI_W$period)
overlap.CALCOFI_W <- NA
  for(i in 1:4){
    temp <- CALCOFI_W%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFI_W <-rbind(temp2,overlap.CALCOFI_W)
  }
overlap.CALCOFI_W<-overlap.CALCOFI_W%>%mutate(Survey="CALCOFI", region="SCC", season="Winter", offset=0) 


#running models with 1-year offset

CALCOFI_Woffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  CALCOFI_Woffset <-rbind(CALCOFI_Woffset,scaled.anom)
}
CALCOFI_Woffset<-CALCOFI_Woffset%>%mutate(survey="CALCOFI",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(CALCOFI_Woffset, CALCOFI_Woffset$beta)
postplot(CALCOFI_Woffset, CALCOFI_Woffset$alpha)

CALCOFI_Woffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFI_Woffset_full <-rbind(CALCOFI_Woffset_full,scaled.anom)
}

CALCOFI_Woffset_full<-CALCOFI_Woffset_full%>%mutate(survey="CALCOFI",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFI_Woffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFI_Woffset<- bind_rows(CALCOFI_Woffset,CALCOFI_Woffset_full)%>%
  mutate(region="SCC",Season="Winter", lag = 1)

index.names <- unique(CALCOFI_Woffset$Index)
period.names <- unique(CALCOFI_Woffset$period)
overlap.CALCOFI_Woffset <- NA
 for(i in 1:4){
    temp <- CALCOFI_Woffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFI_Woffset <-rbind(temp2,overlap.CALCOFI_Woffset)
  }
overlap.CALCOFI_Woffset<-overlap.CALCOFI_Woffset%>%mutate(Survey="CALCOFI", region="SCC", season="Winter", offset=1)



#### RREAS ####
######Spring #####

data <- dat%>%filter(season=="Spring", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

RREAS <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  RREAS <-rbind(RREAS,scaled.anom)
}

RREAS<-RREAS%>%mutate(survey="RREAS",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(RREAS, RREAS$beta)
postplot(RREAS, RREAS$alpha)

RREAS_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREAS_full <-rbind(RREAS_full,scaled.anom)
}

RREAS_full<-RREAS_full%>%mutate(survey="RREAS",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREAS, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREAS<- bind_rows(RREAS,RREAS_full)%>%
  mutate(region="CCC",Season="Spring", lag = 0)

index.names <- unique(RREAS$Index)
period.names <- unique(RREAS$period)
overlap.RREAS <- NA
for(i in 1:4){
    temp <- RREAS%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREAS <-rbind(temp2,overlap.RREAS)
  }
overlap.RREAS<-overlap.RREAS%>%mutate(Survey="RREAS", region="CCC", season="Spring", offset=0) 

#running models with 1-year offset

RREASoffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  RREASoffset <-rbind(RREASoffset,scaled.anom)
}
RREASoffset<-RREASoffset%>%mutate(survey="RREAS",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(RREASoffset, RREASoffset$beta)
postplot(RREASoffset, RREASoffset$alpha)

RREASoffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREASoffset_full <-rbind(RREASoffset_full,scaled.anom)
}

RREASoffset_full<-RREASoffset_full%>%mutate(survey="RREAS",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREASoffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREASoffset<- bind_rows(RREASoffset,RREASoffset_full)%>%
  mutate(region="CCC",Season="Spring", lag = 1)

index.names <- unique(RREASoffset$Index)
period.names <- unique(RREASoffset$period)
overlap.RREASoffset <- NA
  for(i in 1:4){
    temp <- RREASoffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREASoffset <-rbind(temp2,overlap.RREASoffset)
  }
overlap.RREASoffset<-overlap.RREASoffset%>%mutate(Survey="RREAS", region="CCC", season="Spring", offset=1) 

###### Winter #####

data <- dat%>%filter(season=="Winter", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

RREAS_W <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  RREAS_W <-rbind(RREAS_W,scaled.anom)
}

RREAS_W<-RREAS_W%>%mutate(survey="RREAS",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(RREAS_W, RREAS_W$beta)
postplot(RREAS_W, RREAS_W$alpha)

RREAS_W_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREAS_W_full <-rbind(RREAS_W_full,scaled.anom)
}

RREAS_W_full<-RREAS_W_full%>%mutate(survey="RREAS",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREAS_W_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREAS_W<- bind_rows(RREAS_W,RREAS_W_full)%>%
  mutate(region="CCC",Season="Winter", lag = 0)

index.names <- unique(RREAS_W$Index)
period.names <- unique(RREAS_W$period)
overlap.RREAS_W <- NA
 for(i in 1:4){
    temp <- RREAS_W%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREAS_W <-rbind(temp2,overlap.RREAS_W)
  }
overlap.RREAS_W<-overlap.RREAS_W%>%mutate(Survey="RREAS", region="CCC", season="Winter", offset=0) 


#running models with 1-year offset

RREAS_Woffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  RREAS_Woffset <-rbind(RREAS_Woffset,scaled.anom)
}
RREAS_Woffset<-RREAS_Woffset%>%mutate(survey="RREAS",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(RREAS_Woffset, RREAS_Woffset$beta)
postplot(RREAS_Woffset, RREAS_Woffset$alpha)

RREAS_Woffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREAS_Woffset_full <-rbind(RREAS_Woffset_full,scaled.anom)
}

RREAS_Woffset_full<-RREAS_Woffset_full%>%mutate(survey="RREAS",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREAS_Woffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREAS_Woffset<- bind_rows(RREAS_Woffset,RREAS_Woffset_full)%>%
  mutate(region="CCC",Season="Winter", lag = 1)

index.names <- unique(RREAS_Woffset$Index)
period.names <- unique(RREAS_Woffset$period)
overlap.RREAS_Woffset <- NA
  for(i in 1:4){
    temp <- RREAS_Woffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREAS_Woffset <-rbind(temp2,overlap.RREAS_Woffset)
  }
overlap.RREAS_Woffset<-overlap.RREAS_Woffset%>%mutate(Survey="RREAS", region="CCC", season="Winter", offset=1) 


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
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Spring"
eras <- data.frame(era.region2=seq(1,9), era.region=seq(4,12))
data<- climate_dat%>%filter(season=="Spring"&region!='GoA'&Year_lag<2023)%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()%>%
  left_join(eras)

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
  ggtitle("Spring")
###### STI #######
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


ggplot(STI, aes(x = beta, fill = as.factor(period))) +
  theme_bw() +
  facet_wrap(region~Index, ncol = 3, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

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


ggplot(LUSI, aes(x = beta, fill = as.factor(period), group=as.factor(era.region))) +
  theme_bw() +
  facet_wrap(region~Index, ncol = 4, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

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


ggplot(TUMI, aes(x = beta, fill = as.factor(period), group=as.factor(era.region))) +
  theme_bw() +
  facet_wrap(region~Index, ncol = 4, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")


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
  mutate(Season='Winter', lag=0, era.region=era.region2)


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
  mutate(Season='Winter', lag=0, era.region=era.region2)


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
  mutate(Season='Winter', lag=0, era.region=era.region2)


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
  mutate(season="Spring", lag=0)


 ggplot(bioup_spring, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(Index~Survey, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1],col[2], col[3], 'grey')) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")

 
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
          upwelling,
          upwelling,TUMI,STI,LUSI,TUMIW,STIW,LUSIW,
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

