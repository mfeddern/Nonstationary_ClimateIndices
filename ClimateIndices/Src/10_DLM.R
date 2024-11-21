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
library(dplyr)
library(MARSS)
set.seed(1234)

##### Data ####
climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))

daty <- dfa%>%dplyr::select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  #rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa%>%select(-estimate, -lower, -upper), by=c('trend', 'Year_lag','season'))%>%
  select(Year_lag, trend, season, estimate, seasonal_PDO, seasonal_ONI, seasonal_NPH, seasonal_NPGO)%>%
  mutate(offset=1)

datx <- dfa%>%
  select(Year_lag, trend, season, estimate, seasonal_PDO, seasonal_ONI, seasonal_NPH, seasonal_NPGO)%>%
  mutate(offset=0)%>%
  add_row(daty)

Index_Names<- c("PDO","NPGO","NPH","ONI")
Phen_Names<- c('TUMI','LUSI','STI')
Param_name<-c('alpha', 'beta')
region<-c("Northern CC","Southern CC","Central CC")
season<-c("Spring","Winter")
Trend_Name<- c("CALCOFI", "RREAS")
Cop_Name<- c("Northern Copepods", "Southern Copepods")
offset<-c(0,1)


#### Biology - Climate Models: RREAS/CalCOFI ####

#assigning names to create and object

Bio_Clim_DLM1<-NULL
datx%>%
  pivot_wider(names_from = trend, values_from = estimate)

for(l in 1:length(season)){
  dd<- datx%>%
      filter(season==season[l])%>%
    distinct()
  for(k in 1:2){
    dd2<-dd%>%
      filter(offset==offset[k])
      for(j in 1:2){
        data<-dd2%>%
          filter(trend==Trend_Name[j])
        for(i in 1:4){
        ## get time indices
        years <- data[, 'Year_lag']
        ## number of years of data
        TT <- length(years)
        ## get response variable: logit(survival)
        dat <- matrix(c(data[, 'estimate'], data[, 'Year_lag']), 
                      nrow = 2, byrow=TRUE)
        #get predictor
        Index<- matrix(c(data[, 'seasonal_PDO'],data[, 'seasonal_NPGO'],data[, 'seasonal_NPH'],data[, 'seasonal_ONI']), 
                       nrow = 4, ncol=TT)
        U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
        Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
        diag(Q) <- c(0.01,0.01)  ## 2x2; diag = (0.1,0.1)
        #diag(Q) <- c("q.alpha","q.beta")  ## 2x2; diag = (0.1,0.1)
        
        ## for observation eqn
        Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
        Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
        Z[1, 2, ] <- scale(Index[i, ]) ## Nx1; predictor variable
        
        ## list of model matrices & vectors
        mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
        
        ## fit univariate DLM
        fit_2 <- MARSS(dat[1,], inits = inits_list, model = mod_list)
        ## get estimates of alpha
        alpha_hat <- fit_2$states[1,]
        ## get estimates of eta
        beta_hat <- fit_2$states[2,]
        ###get alpha se
        alpha_se <- fit_2$states.se[1,]
        ###get alpha se
        beta_se <- fit_2$states.se[2,]
        
        temp<-cbind(as.numeric(dat[2,]),as.numeric(beta_hat),as.numeric(alpha_hat),as.numeric(beta_se),as.numeric(alpha_se),
                    Trend_Name[j],Index_Names[i],season[l],as.numeric(offset[k]))
        
        Bio_Clim_DLM1<-rbind(Bio_Clim_DLM1,temp)
        ## plot the estimated level and drift
        par(mfrow = c(2,1), mai = c(0.8, 0.8, 0.2, 0.2), omi = c(0, 0, 0, 0))
        ## plot alpha
        plot.ts(alpha_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(alpha[t]))
        ## plot eta
        plot.ts(beta_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(beta[t]))
        
      }  
    }
    
  }
}



colnames(Bio_Clim_DLM1)<-c("Year", "beta", "alpha","beta_se", "alpha_se","Trend","Climate_Index", "season","offset")
Bio_Clim_DLM1<-data.frame(Bio_Clim_DLM1)
                                  
#### Biology - Climate Models: Cop ####

climpdo <- climate_dat_cop%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH, season)%>%
  distinct()%>%
  filter(Year_lag<2023)

dd1 <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
                seasonal_copepod_southern)%>%
  left_join(climpdo%>%filter(season=="Spring"))%>%
  distinct()

dd2 <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
                seasonal_copepod_southern)%>%
  left_join(climpdo%>%filter(season=="Winter"))%>%
  distinct()

dd<-dd1%>%bind_rows(dd2)
Bio_Clim_DLM2<-NULL
for(l in 1:length(season)){
  data<- dd%>%
    filter(season==season[l])%>%
    distinct()
    for(j in 1:2){
      for(i in 1:4){
        ## get time indices
        years <- data[, 'Year_lag']
        ## number of years of data
        TT <- length(years)
        ## get response variable: logit(survival)
        dat <- matrix(c(data[, 'seasonal_copepod_northern'], 
                        data[, 'seasonal_copepod_southern'],
                        data[, 'Year_lag']), 
                      nrow = 3, byrow=TRUE)
        #get predictor
        Index<- matrix(c(data[, 'seasonal_PDO'],data[, 'seasonal_NPGO'],data[, 'seasonal_NPH'],data[, 'seasonal_ONI']), 
                       nrow = 4, ncol=TT)
        U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
        Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
        diag(Q) <- c(0.01,0.01)  ## 2x2; diag = (0.1,0.1)
        #diag(Q) <- c("q.alpha","q.beta")  ## 2x2; diag = (0.1,0.1)
        
        ## for observation eqn
        Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
        Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
        Z[1, 2, ] <- scale(Index[i, ]) ## Nx1; predictor variable
        
        ## list of model matrices & vectors
        mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
        
        ## fit univariate DLM
        fit_2 <- MARSS(dat[j,], inits = inits_list, model = mod_list)
        ## get estimates of alpha
        alpha_hat <- fit_2$states[1,]
        ## get estimates of eta
        beta_hat <- fit_2$states[2,]
        ###get alpha se
        alpha_se <- fit_2$states.se[1,]
        ###get alpha se
        beta_se <- fit_2$states.se[2,]
        
        temp<-cbind(as.numeric(dat[3,]),as.numeric(beta_hat),as.numeric(alpha_hat),as.numeric(beta_se),as.numeric(alpha_se),
                    Cop_Name[j],Index_Names[i],season[l])
        
        Bio_Clim_DLM2<-rbind(Bio_Clim_DLM2,temp)
        ## plot the estimated level and drift
        par(mfrow = c(2,1), mai = c(0.8, 0.8, 0.2, 0.2), omi = c(0, 0, 0, 0))
        ## plot alpha
        plot.ts(alpha_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(alpha[t]))
        ## plot eta
        plot.ts(beta_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(beta[t]))
        
      }  
    }
    
  }


unique(factor(Bio_Clim_DLM3$Trend))

colnames(Bio_Clim_DLM2)<-c("Year", "beta", "alpha","beta_se", "alpha_se","Trend","Climate_Index", "season")
Bio_Clim_DLM2<-data.frame(Bio_Clim_DLM2)
Bio_Clim_DLM2<-Bio_Clim_DLM2%>%mutate(offset='0')
Bio_Clim_DLM3<-Bio_Clim_DLM1%>%bind_rows(Bio_Clim_DLM2)
Bio_Clim_DLM3$Trend <- factor(Bio_Clim_DLM3$Trend,
                                      levels = c("Southern Copepods",
                                                 "Northern Copepods",
                                                 "RREAS", "CALCOFI"))



#####Biology - Phenology Models  ####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))

ddat<- left_join(climate_dat%>%filter(season=="Spring"&region!='GoA'&Year_lag<2023)%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                       ifelse(region=="Southern CC", "SCC", "CCC")))%>%
  select(Year_lag, region, stand_tumi, stand_sti, stand_lusi,),
  bind_rows(dd%>%filter(season=="Spring")%>%
  select(Year_lag, region,seasonal_copepod_northern, seasonal_copepod_southern)%>%
  pivot_longer(!c(Year_lag,region),names_to = 'trend', values_to = 'estimate'),
datx%>%filter(season=="Spring",offset==0, trend!="SEA")%>%select(Year_lag, trend, estimate)%>%
  mutate(region=ifelse(trend=='CALCOFI',"SCC","CCC"))))

ddat <- ddat[complete.cases(ddat), ]

Trend_Name2<- c("CALCOFI", "RREAS","seasonal_copepod_southern", "seasonal_copepod_northern")
region2<-c("SCC", "CCC", "NCC","NCC")

Phen_Bio_DLM<-NULL
  for(k in 1:4){
    data<-ddat%>%
      filter(trend==Trend_Name2[k])
    
    ## get time indices
    years <- data[, 'Year_lag']
    ## number of years of data
    TT <- length(years)
    ## get response variable: logit(survival)
    dat <- matrix(c(data[, 'estimate'], 
                    data[, 'Year_lag']), 
                  nrow = 2, byrow=TRUE)
    #get predictor
    Index<- matrix(c(data[, 'stand_tumi'],data[, 'stand_sti'],data[, 'stand_lusi']), 
                   nrow = 3, ncol=TT)
      for(j in 1:3){
        U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
        Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
        diag(Q) <- c(0.01,0.01)  ## 2x2; diag = (0.1,0.1)
        #diag(Q) <- c("q.alpha","q.beta")  ## 2x2; diag = (0.1,0.1)
        
        ## for observation eqn
        Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
        Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
        Z[1, 2, ] <- scale(Index[j, ]) ## Nx1; predictor variable
        
        ## list of model matrices & vectors
        mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
        
        ## fit univariate DLM
        fit_2 <- MARSS(dat[1,], inits = inits_list, model = mod_list)
        ## get estimates of alpha
        alpha_hat <- fit_2$states[1,]
        ## get estimates of eta
        beta_hat <- fit_2$states[2,]
        ###get alpha se
        alpha_se <- fit_2$states.se[1,]
        ###get alpha se
        beta_se <- fit_2$states.se[2,]
        
        temp<-cbind(as.numeric(dat[2,]),as.numeric(beta_hat),as.numeric(alpha_hat),as.numeric(beta_se),as.numeric(alpha_se),
                    Phen_Names[j],Trend_Name2[k], region2[k])
        
        Phen_Bio_DLM<-rbind(Phen_Bio_DLM,temp)
        ## plot the estimated level and drift
        par(mfrow = c(2,1), mai = c(0.8, 0.8, 0.2, 0.2), omi = c(0, 0, 0, 0))
        ## plot alpha
        plot.ts(alpha_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(alpha[t]))
        ## plot eta
        plot.ts(beta_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(beta[t]))
    }
}

colnames(Phen_Bio_DLM)<-c("Year", "beta", "alpha","beta_se", "alpha_se","Upwelling_Index","Trend", "region")
Phen_Bio_DLM<-data.frame(Phen_Bio_DLM)

Phen_Bio_DLM<-Phen_Bio_DLM%>%
  mutate(Trend=ifelse(Trend=="CALCOFI","CALCOFI",
                       ifelse(Trend=="RREAS", "RREAS", 
                        ifelse(Trend=="seasonal_copepod_northern", "N. Cop.",
                        "S. Cop."))))
Phen_Bio_DLM$Trend <- factor(Phen_Bio_DLM$Trend,
                                      levels = c("N. Cop.","S. Cop.",
                                                 "RREAS", "CALCOFI"))


##### Phenology-Climate Model Runs #####

climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
Phen_Clim_DLM<-NULL
Phen_Clim_region<-NULL


for(l in 1: length(season)){
  dd<- climate_dat%>%filter(season==season[l]&region!='GoA'&Year_lag<2023)%>%
    #  filter(Year_lag!=2022&Year_lag!=2015)%>%
    distinct()
  for(k in 1:length(region)){
    data<-dd%>%
      filter(region==region[k])
    
    ## get time indices
    years <- data[, 'Year_lag']
    ## number of years of data
    TT <- length(years)
    ## get response variable: logit(survival)
    dat <- matrix(c(data[, 'stand_tumi'], 
                    data[, 'stand_lusi'],
                    data[, 'stand_sti'],
                    data[, 'Year_lag']), 
                  nrow = 4, byrow=TRUE)
    #get predictor
    Index<- matrix(c(data[, 'seasonal_PDO'],data[, 'seasonal_NPGO'],data[, 'seasonal_NPH'],data[, 'seasonal_ONI']), 
                   nrow = 4, ncol=TT)
    for(i in 1:4){
      for(j in 1:3){
        
        U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
        Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
        diag(Q) <- c(0.01,0.01)  ## 2x2; diag = (0.1,0.1)
        #diag(Q) <- c("q.alpha","q.beta")  ## 2x2; diag = (0.1,0.1)
        
        ## for observation eqn
        Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
        Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
        Z[1, 2, ] <- scale(Index[i, ]) ## Nx1; predictor variable
        
        ## list of model matrices & vectors
        mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and equal")
        
        ## fit univariate DLM
        fit_2 <- MARSS(dat[j,], inits = inits_list, model = mod_list)
        ## get estimates of alpha
        alpha_hat <- fit_2$states[1,]
        ## get estimates of eta
        beta_hat <- fit_2$states[2,]
        ###get alpha se
        alpha_se <- fit_2$states.se[1,]
        ###get alpha se
        beta_se <- fit_2$states.se[2,]
        
        temp<-cbind(as.numeric(dat[4,]),as.numeric(beta_hat),as.numeric(alpha_hat),as.numeric(beta_se),as.numeric(alpha_se),
                    Phen_Names[j],Index_Names[i],region[k],season[l])
        
        Phen_Clim_DLM<-rbind(Phen_Clim_DLM,temp)
        ## plot the estimated level and drift
        par(mfrow = c(2,1), mai = c(0.8, 0.8, 0.2, 0.2), omi = c(0, 0, 0, 0))
        ## plot alpha
        plot.ts(alpha_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(alpha[t]))
        ## plot eta
        plot.ts(beta_hat, las = 1, lwd = 2, col = "blue",
                ylab = expression(beta[t]))
        
      }  
    }
    
  }
}



colnames(Phen_Clim_DLM)<-c("Year", "beta", "alpha","beta_se", "alpha_se","Upwelling_Index","Climate_Index", "region","season")
Phen_Clim_DLM_Region<-data.frame(Phen_Clim_DLM)
Phen_Clim_DLM_Region<-Phen_Clim_DLM_Region%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                       ifelse(region=="Southern CC", "SCC", "CCC")))
Phen_Clim_DLM_Region$region <- factor(Phen_Clim_DLM_Region$region,
                                      levels = c("NCC","CCC","SCC"))




#### Making Plots BETAS ####

###### Spring ####

Phen_Bio_DLM_line<-ggplot(data=Phen_Bio_DLM, aes(y = as.numeric(beta), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Upwelling_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
 scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="")+
  theme(legend.position="none")
Phen_Bio_DLM_line

Bio_Clim_DLM_line<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Spring"&offset==0), aes(y = as.numeric(beta), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Slope")+
  theme(legend.position="bottom")
Bio_Clim_DLM_line


Phen_Clim_DLM_line<-ggplot(data=Phen_Clim_DLM_Region%>%filter(season=="Spring"), aes(y = as.numeric(beta), fill=region,x= as.numeric(Year),col = as.factor(region))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_grid(Climate_Index~Upwelling_Index,  scales='free') +
  geom_line() +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="")+
  labs(x = "Year",
       y = "Slope")+
  theme(legend.position="bottom")
Phen_Clim_DLM_line

z.plot <- ggplot()+theme_void()

pdf("Output/DLMbetaspring.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_line,Bio_Clim_DLM_line, ncol = 3, labels = c("A", "B", "C"),
          ggarrange(Phen_Bio_DLM_line,z.plot,nrow=2,labels = c("",""), 
                    heights = c(6.25,1)), widths=c(3,1.5,1.5))

dev.off()



###### Winter ####


Bio_Clim_DLM_lineW<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Winter"&offset==0), aes(y = as.numeric(beta), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Slope")+
  theme(legend.position="right")
Bio_Clim_DLM_lineW


Phen_Clim_DLM_lineW<-ggplot(data=Phen_Clim_DLM_Region%>%filter(season=="Winter"), aes(y = as.numeric(beta), fill=region,x= as.numeric(Year),col = as.factor(region))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_grid(Climate_Index~Upwelling_Index,  scales='free') +
  geom_line() +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="")+
  labs(x = "Year",
       y = "Slope")+
  theme(legend.position="right")
Phen_Clim_DLM_lineW

pdf("Output/DLMbetawinter.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_lineW,Bio_Clim_DLM_lineW, ncol = 2, labels = c("A", "B"),
                    heights = c(6.25,1), widths=c(3,3))

dev.off()



###### Offset ####
Bio_Clim_DLM_lineOff<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Spring"&offset==1), aes(y = as.numeric(beta), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[2], col[3])) +
  scale_linetype_manual(values = c(1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Slope")+
  theme(legend.position="right")
Bio_Clim_DLM_lineOff
pdf("Output/DLMbetaoff.pdf", 4,5) 
Bio_Clim_DLM_lineOff
dev.off()


###### Smoothed Spring ####

Phen_Bio_DLM_line<-ggplot(data=Phen_Bio_DLM, aes(y = as.numeric(beta), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Upwelling_Index,  scales='free', nrow=4) +
  geom_smooth(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="")+
  theme(legend.position="none")
Phen_Bio_DLM_line

Bio_Clim_DLM_line<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Spring"&offset==0), aes(y = as.numeric(beta), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_smooth(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Slope")+
  theme(legend.position="bottom")
Bio_Clim_DLM_line


Phen_Clim_DLM_line<-ggplot(data=Phen_Clim_DLM_Region%>%filter(season=="Spring"), aes(y = as.numeric(beta), fill=region,x= as.numeric(Year),col = as.factor(region))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_grid(Climate_Index~Upwelling_Index,  scales='free') +
  geom_smooth() +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(beta)-as.numeric(beta_se), ymax = as.numeric(beta)+as.numeric(beta_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="")+
  labs(x = "Year",
       y = "Slope")+
  theme(legend.position="bottom")
Phen_Clim_DLM_line


pdf("Output/DLMbetaspringsmooth.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_line,Bio_Clim_DLM_line, ncol = 3, labels = c("A", "B", "C"),
          ggarrange(Phen_Bio_DLM_line,z.plot,nrow=2,labels = c("",""), 
                    heights = c(6.25,1)), widths=c(3,1.5,1.5))

dev.off()



#### Making Plots Alphas ####

###### Spring ####

Phen_Bio_DLM_lineA<-ggplot(data=Phen_Bio_DLM, aes(y = as.numeric(alpha), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Upwelling_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(alpha)-as.numeric(alpha_se), ymax = as.numeric(alpha)+as.numeric(alpha_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="")+
  theme(legend.position="none")
Phen_Bio_DLM_lineA

Bio_Clim_DLM_lineA<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Spring"&offset==0), aes(y = as.numeric(alpha), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(alpha)-as.numeric(alpha_se), ymax = as.numeric(alpha)+as.numeric(alpha_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Intercept")+
  theme(legend.position="bottom")
Bio_Clim_DLM_lineA


Phen_Clim_DLM_lineA<-ggplot(data=Phen_Clim_DLM_Region%>%filter(season=="Spring"), aes(y = as.numeric(alpha), fill=region,x= as.numeric(Year),col = as.factor(region))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_grid(Climate_Index~Upwelling_Index,  scales='free') +
  geom_line() +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(alpha)-as.numeric(alpha_se), ymax = as.numeric(alpha)+as.numeric(alpha_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="")+
  labs(x = "Year",
       y = "Intercept")+
  theme(legend.position="bottom")
Phen_Clim_DLM_lineA


pdf("Output/DLMalphaspringAlpha.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_lineA,Bio_Clim_DLM_lineA, ncol = 3, labels = c("A", "B", "C"),
          ggarrange(Phen_Bio_DLM_lineA,z.plot,nrow=2,labels = c("",""), 
                    heights = c(6.25,1)), widths=c(3,1.5,1.5))

dev.off()



###### Winter ####


Bio_Clim_DLM_lineW_A<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Winter"&offset==0), aes(y = as.numeric(alpha), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[1],col[1], col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(alpha)-as.numeric(alpha_se), ymax = as.numeric(alpha)+as.numeric(alpha_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[1], col[2], col[3])) +
  scale_linetype_manual(values = c(2,1,1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Intercept")+
  theme(legend.position="right")
Bio_Clim_DLM_lineW_A


Phen_Clim_DLM_lineW_A<-ggplot(data=Phen_Clim_DLM_Region%>%filter(season=="Winter"), aes(y = as.numeric(alpha), fill=region,x= as.numeric(Year),col = as.factor(region))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_grid(Climate_Index~Upwelling_Index,  scales='free') +
  geom_line() +
  scale_fill_manual(values = c(col[1],col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(alpha)-as.numeric(alpha_se), ymax = as.numeric(alpha)+as.numeric(alpha_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[1],col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="")+
  labs(x = "Year",
       y = "Intercept")+
  theme(legend.position="right")
Phen_Clim_DLM_lineW_A

pdf("Output/DLMalphawinter.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_lineW_A,Bio_Clim_DLM_lineW_A, ncol = 2, labels = c("A", "B"),
          heights = c(6.25,1), widths=c(3,3))

dev.off()



###### Offset ####
Bio_Clim_DLM_lineOffA<-ggplot(data=Bio_Clim_DLM3%>%filter(season=="Spring"&offset==1), aes(y = as.numeric(alpha), fill=Trend,x= as.numeric(Year),col = as.factor(Trend))) +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  geom_vline(xintercept = 1988, lty = 3) +
  geom_vline(xintercept = 2012, lty = 3) +
  theme_bw() +
  facet_wrap(~Climate_Index,  scales='free', nrow=4) +
  geom_line(aes(linetype=Trend)) +
  scale_fill_manual(values = c(col[2], col[3])) +
  geom_ribbon(aes(ymin = as.numeric(alpha)-as.numeric(alpha_se), ymax = as.numeric(alpha)+as.numeric(alpha_se)), alpha=0.1,colour = NA)+
  scale_colour_manual(values = c(col[2], col[3])) +
  scale_linetype_manual(values = c(1,1))+
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  guides(fill = FALSE)+
  labs(colour="", linetype="", x="Year", y="Intercept")+
  theme(legend.position="right")
Bio_Clim_DLM_lineOffA
pdf("Output/DLMalphaoff.pdf", 4,5) 
Bio_Clim_DLM_lineOffA
dev.off()


#### Summary Stats ####
Phen_Bio_DLM<-Phen_Bio_DLM%>%mutate(era=ifelse(Year<1988,1, 
                      ifelse(Year>2012,3,2)))

Phen_Bio_DLM2<-Phen_Bio_DLM%>% group_by(Trend,Upwelling_Index,era)%>%
  summarise(mean_beta=mean(as.numeric(beta)), mean_alpha=mean(as.numeric(alpha)),
            sd_beta=sd(beta), sd_alpha=sd(alpha))
  
Phen_Bio_DLM_point <- ggplot(data=Phen_Bio_DLM2, aes(x=mean_beta, y=Trend, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_beta-sd_beta, xmax=mean_beta+sd_beta, col=as.factor(era)), width=.2) +
  facet_wrap(.~Upwelling_Index, ncol=1) +
  xlab("Beta") +
  ylab("") +
  theme(legend.position="bottom")

Bio_Clim_DLM3<-Bio_Clim_DLM3%>%mutate(era=ifelse(Year<1988,1, 
                                               ifelse(Year>2012,3,2)))

Bio_Clim_DLM32<-Bio_Clim_DLM3%>% group_by(Climate_Index,Trend,era, offset,season)%>%
  summarise(mean_beta=mean(as.numeric(beta)), mean_alpha=mean(as.numeric(alpha)),
            sd_beta=sd(beta), sd_alpha=sd(alpha))

Phen_Bio_DLM_point<-ggplot(data=Bio_Clim_DLM32%>%filter(season=="Spring"&offset==0), aes(x=mean_beta, y=Trend, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_beta-sd_beta, xmax=mean_beta+sd_beta, col=as.factor(era)), width=.2) +
  facet_wrap(.~Climate_Index, ncol=1) +
  xlab("Beta") +
  ylab("") +
  theme(legend.position="bottom")

Bio_Clim_DLM32<-Bio_Clim_DLM3%>% group_by(Climate_Index,Trend,era, offset,season)%>%
  summarise(mean_beta=mean(as.numeric(beta)), mean_alpha=mean(as.numeric(alpha)),
            sd_beta=sd(beta), sd_alpha=sd(alpha))

Bio_Clim_DLM_point<-ggplot(data=Bio_Clim_DLM32%>%filter(season=="Spring"&offset==0), aes(x=mean_beta, y=Trend, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_beta-sd_beta, xmax=mean_beta+sd_beta, col=as.factor(era)), width=.2) +
  facet_wrap(.~Climate_Index, ncol=1) +
  xlab("Beta") +
  ylab("") +
  theme(legend.position="bottom")

Phen_Clim_DLM_Region<-Phen_Clim_DLM_Region%>%mutate(era=ifelse(Year<1988,1, 
                                               ifelse(Year>2012,3,2)))

Phen_Clim_DLM_Region2<-Phen_Clim_DLM_Region%>% group_by(Climate_Index,Upwelling_Index,region,era, season)%>%
  summarise(mean_beta=mean(as.numeric(beta)), mean_alpha=mean(as.numeric(alpha)),
            sd_beta=sd(beta), sd_alpha=sd(alpha))


Phen_Clim_DLM_point<-ggplot(data=Phen_Clim_DLM_Region2%>%filter(season=="Spring"), aes(x=mean_beta, y=region, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_beta-sd_beta, xmax=mean_beta+sd_beta, col=as.factor(era)), width=.2) +
  facet_grid(Climate_Index~Upwelling_Index) +
  xlab("Beta") +
  ylab("") +
  theme(legend.position="bottom")

pdf("Output/DLMbetaspringPOINT.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_point,Bio_Clim_DLM_point, ncol = 3, labels = c("A", "B", "C"),
          ggarrange(Phen_Bio_DLM_point,z.plot,nrow=2,labels = c("",""), 
                    heights = c(6.25,1)), widths=c(3,1.5,1.5))

dev.off()


Phen_Bio_DLM_pointA<-ggplot(data=Bio_Clim_DLM32%>%filter(season=="Spring"&offset==0), aes(x=mean_alpha, y=Trend, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_alpha-sd_alpha, xmax=mean_alpha+sd_alpha, col=as.factor(era)), width=.2) +
  facet_wrap(.~Climate_Index, ncol=1) +
  xlab("alpha") +
  ylab("") +
  theme(legend.position="bottom")

Bio_Clim_DLM_pointA<-ggplot(data=Bio_Clim_DLM32%>%filter(season=="Spring"&offset==0), aes(x=mean_alpha, y=Trend, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_alpha-sd_alpha, xmax=mean_alpha+sd_alpha, col=as.factor(era)), width=.2) +
  facet_wrap(.~Climate_Index, ncol=1) +
  xlab("alpha") +
  ylab("") +
  theme(legend.position="bottom")



Phen_Clim_DLM_pointA<-ggplot(data=Phen_Clim_DLM_Region2%>%filter(season=="Spring"), aes(x=mean_alpha, y=region, coulour=as.factor(era))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 1, col='grey') +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Era")+
  geom_point(aes(col=as.factor(era)))+
  geom_errorbar(aes(xmin=mean_alpha-sd_alpha, xmax=mean_alpha+sd_alpha, col=as.factor(era)), width=.2) +
  facet_grid(Climate_Index~Upwelling_Index) +
  xlab("alpha") +
  ylab("") +
  theme(legend.position="bottom")

pdf("Output/DLMalphaspringPOINT.pdf", 11,6) 
ggarrange(Phen_Clim_DLM_pointA,Bio_Clim_DLM_pointA, ncol = 3, labels = c("A", "B", "C"),
          ggarrange(Phen_Bio_DLM_pointA,z.plot,nrow=2,labels = c("",""), 
                    heights = c(6.25,1)), widths=c(3,1.5,1.5))

dev.off()
