library(dplyr)
library(nord)
library(tidyr)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
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
library(rgdal)
library(colorspace)
library(PBSmapping) #powerful mapping functions developed by Pacific Biological Station
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(bayestestR)
library(ggh4x)
library(ggridges)

#### Data Reorg ####
Violin_Data <-readRDS('data/Violin_Data.rds')

Violin_indices <- filter(Violin_Data, survey!="Upwelling",survey!="TUMI",
                         survey!="LUSI",survey!="STI") 
Violin_index<- filter(Violin_Data, Index=="Upwelling"|Index=="TUMI"|
                         Index=="LUSI"|Index=="STI") 
Violin_upwelling <- filter(Violin_Data, survey=="Upwelling"|survey=="TUMI"|
                             survey=="LUSI"|survey=="STI") 


ratio_biological<-rbind(Violin_indices%>%filter(period==2)%>%group_by(Index,survey, lag, Season)%>%
                       add_column((Violin_indices%>%filter(period==3)%>%
                                     group_by(Index,survey, lag, Season))$beta-(Violin_indices%>%
                                            filter(period==2)%>%
                                            group_by(Index,survey, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%
  mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2")%>%
#  bind_rows(rbind(Violin_indices%>%filter(period==3)%>%group_by(Index,region, lag, Season)%>%
#                       add_column((Violin_indices%>%filter(period==3)%>%
#                                     group_by(Index,region, lag, Season))$beta-(Violin_indices%>%
#                                            filter(period==4)%>%
#                                            group_by(Index,region, lag, Season))$beta))%>%
#                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - Full TS", Difference2="Era 3 - Full TS"))%>%
  bind_rows(rbind(Violin_indices%>%filter(period==2)%>%group_by(Index,survey, lag, Season)%>%
                       add_column((Violin_indices%>%filter(period==3)%>%
                                     group_by(Index,survey, lag, Season))$beta-(Violin_indices%>%
                                            filter(period==2)%>%
                                            group_by(Index,survey, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%
  mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2"))

ratio_index<-rbind(Violin_index%>%filter(period==2)%>%group_by(Index,survey, lag, Season)%>%
                          add_column((Violin_index%>%filter(period==3)%>%
                                        group_by(Index,survey, lag, Season))$beta-(Violin_index%>%
                                                                                     filter(period==2)%>%
                                                                                     group_by(Index,survey, lag, Season))$beta))%>%
  rename('beta_diff'=`... - ...`)%>%
  mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2")%>%
  #  bind_rows(rbind(Violin_index%>%filter(period==3)%>%group_by(Index,region, lag, Season)%>%
  #                       add_column((Violin_index%>%filter(period==3)%>%
  #                                     group_by(Index,region, lag, Season))$beta-(Violin_index%>%
  #                                            filter(period==4)%>%
  #                                            group_by(Index,region, lag, Season))$beta))%>%
  #                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - Full TS", Difference2="Era 3 - Full TS"))%>%
  bind_rows(rbind(Violin_index%>%filter(period==2)%>%group_by(Index,survey, lag, Season)%>%
                    add_column((Violin_index%>%filter(period==3)%>%
                                  group_by(Index,survey, lag, Season))$beta-(Violin_index%>%
                                                                               filter(period==2)%>%
                                                                               group_by(Index,survey, lag, Season))$beta))%>%
              rename('beta_diff'=`... - ...`)%>%
              mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2"))


ratio_upwelling<-rbind(Violin_upwelling%>%filter(period==2)%>%group_by(Index,region, survey,lag, Season)%>%
                       add_column((Violin_upwelling%>%filter(period==3)%>%
                                     group_by(Index,region, lag, Season))$beta-(Violin_upwelling%>%
                                            filter(period==2)%>%
                                            group_by(Index,region,survey, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2")%>%
#  bind_rows(rbind(Violin_upwelling%>%filter(period==3)%>%group_by(Index,region, lag, Season)%>%
#                       add_column((Violin_upwelling%>%filter(period==3)%>%
#                                     group_by(Index,region, lag, Season))$beta-(Violin_upwelling%>%
#                                            filter(period==4)%>%
#                                            group_by(Index,region, lag, Season))$beta))%>%
#                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - Full TS", Difference2="Era 3 - Full TS"))%>%
  bind_rows(rbind(Violin_upwelling%>%filter(period==3)%>%group_by(Index,region,survey, lag, Season)%>%
                    add_column((Violin_upwelling%>%filter(period==3)%>%
                                  group_by(Index,region,survey, lag, Season))$beta-(Violin_upwelling%>%
                                                                               filter(period==1)%>%
                                                                               group_by(Index,region, lag, Season))$beta))%>%
              rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - 1967:1988", Difference2="Era 3 - Era 1"))%>%
  bind_rows(rbind(Violin_upwelling%>%filter(period==3)%>%group_by(Index,region,survey, lag, Season)%>%
                    add_column((Violin_upwelling%>%filter(period==2)%>%
                                  group_by(Index,region, survey,lag, Season))$beta-(Violin_upwelling%>%
                                                                               filter(period==1)%>%
                                                                               group_by(Index,region, lag, Season))$beta))%>%
              rename('beta_diff'=`... - ...`)%>%mutate(Difference="1989:2012 - 1967:1988", Difference2="Era 2 - Era 1"))




col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ratio.up<-ratio_upwelling%>%filter(region!="GoA", season=="Spring"|Season=="Spring", lag==0)%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                      ifelse(region=="Southern CC", "SCC", "CCC")))
ratio.up$region <- factor(ratio.up$region,
                levels = c("SCC","CCC","NCC"))

ratio.upw<-ratio_upwelling%>%filter(region!="GoA", season=="Winter"|Season=="Winter", lag==0)%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                       ifelse(region=="Southern CC", "SCC", "CCC")))
ratio.upw$region <- factor(ratio.upw$region,
                          levels = c("SCC","CCC","NCC"))
unique(ratio.bio$Index)
ratio.bio<-ratio_biological%>%filter(Season=="Spring", lag==0)%>%
         mutate(survey=ifelse(survey=="CALCOFI (SCC)","CALCOFI",
                      ifelse(survey=="RREAS (CCC)", "RREAS",
                        ifelse(survey=="N. Copepod (NCC)", "N. Copepod", 
                               ifelse(survey=="S. Copepod (NCC)","S. Copepod",
                                      ifelse(survey=="Southern Copepod (NCC)","S. Copepod",
                                             ifelse(survey=="Northern Copepod (NCC)","N. Copepod",
                                                    ifelse(survey=="RREAS (SCC)","RREAS",
                                             survey))))))))
#ratio.bio<-ratio.bio%>%filter(Index=='PDO'|Index=='ONI'|Index=='NPH'|Index=='NPGO'|
#                                survey=='CALCOFI'|survey=='RREAS'|survey=='N. Copepod'|survey=="S. Copepod")
ratio.bio$region <- factor(ratio.bio$region,
                levels = c("SCC","CCC","NCC"))
ratio.bio$survey <- factor(ratio.bio$survey,
                levels = c("CALCOFI","RREAS","N. Copepod","S. Copepod", 
                           "Upwelling","STI","TUMI","LUSI"))
ratio_index$region <- factor(ratio_index$region,
                           levels = c("SCC","CCC","NCC"))
ratio_index$survey <- factor(ratio_index$survey,
                           levels = c("CALCOFI","RREAS","N. Copepod","S. Copepod", 
                                      "Upwelling","STI","TUMI","LUSI"))

#### Alternative Plot: Posteriors with Densities ####

Violin_Data <-readRDS('data/Violin_Data.rds')
col2<-pnw_palette("Sunset2",3,type="discrete")
Violin_Data<-Violin_Data%>%mutate(region=ifelse(region=="Northern CC","NCC",
                                            ifelse(region=="Southern CC", "SCC", "CCC")))
Violin_Data$region <- factor(Violin_Data$region,
                           levels = c("NCC","CCC","SCC"))
Violin_Data$survey <- factor(Violin_Data$survey,
                             levels = c("N. Copepod","S. Copepod","RREAS",
                                        "CALCOFI","TUMI", "STI", "LUSI", "Upwelling"))
density_phe<-ggplot(Violin_Data%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
         filter(Season=="Spring"&period!=4&lag==0), aes(x = beta,y=Index, fill=as.factor(period))) +
  theme_bw() +
  # geom_errorbar(aes(xmin=median_beta-sd_beta, xmax=median_beta+sd_beta),width = 0)+
  geom_density_ridges2(alpha = 0.4, scale=1, quantiles = c(0.1, 0.9))+
  ggh4x::facet_grid2(region~survey,scales = 'free') +
  scale_fill_manual(values = c(col2[1],col2[2], col2[3]))+
  xlab("Slope") +
  geom_vline(xintercept = 0, lty = 2) +
  ylab("Climate Index")+
  guides(fill=guide_legend(title="Period"))+
  theme(legend.position="bottom")
density_phe

density_bio<-ggplot(Violin_Data%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Season=="Spring"&period!=1&period!=4&lag==0), aes(x = beta,y=Index, fill=as.factor(period))) +
  theme_bw() +
  # geom_errorbar(aes(xmin=median_beta-sd_beta, xmax=median_beta+sd_beta),width = 0)+
  geom_density_ridges(alpha = 0.4, rel_min_height = 0.01,scale=1,rel_min_height = 0.01)+
  facet_wrap(~survey, ncol=1, scales='free') +
  xlab("Slope") +
  #xlim(c(-1.25,1.25))+
  scale_fill_manual(values = c(col2[2], col2[3]))+
  geom_vline(xintercept = 0, lty = 2) +
  ylab("") +
 # geom_text(data=overlap_data,aes(label=ov))+
  theme(legend.position="none")
density_bio
density_bio2<-ggplot(Violin_Data%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Index=="TUMI"|Index=="LUSI"|Index=="STI")%>%
         filter(period!=1&lag==0&period!=4&season=="Spring"), aes(x = beta,y=Index, fill=as.factor(period))) +
  theme_bw() +
  # geom_errorbar(aes(xmin=median_beta-sd_beta, xmax=median_beta+sd_beta),width = 0)+
  geom_density_ridges(alpha = 0.4, rel_min_height = 0.01,scale=1,rel_min_height = 0.01)+
  facet_wrap(~survey, ncol=1, scales='free') +
  xlab("Slope") +
  #xlim(c(-1.25,1.25))+
  scale_fill_manual(values = c(col2[2], col2[3]))+
  geom_vline(xintercept = 0, lty = 2) +
  ylab("Upwelling Index") +
  theme(legend.position="none")
density_bio2

pdf("Output/Density.pdf", 13,6) 
ggarrange(density_phe,density_bio,density_bio2, ncol = 3, labels = c("A", "B", "C"), 
                   widths=c(4,2,2))
dev.off()

col<-pnw_palette("Sunset2",3,type="discrete")
#### Posterior Means Spring ####
dodge<-0.6
mean_beta<-Violin_Data%>%group_by(region,period,survey, Season,Index,lag)%>%
  summarise(median_beta=median(beta),sd_beta_80=ci(beta,ci = 0.95, method = "HDI"),
            sd_beta_50=ci(beta,ci = 0.5, method = "HDI"))
mean_beta$region <- factor(mean_beta$region,
                           levels = c("NCC","CCC","SCC"))
mean_beta$survey <- factor(mean_beta$survey,
                           levels = c("N. Copepod","S. Copepod","RREAS", "CALCOFI",
                                      "Upwelling","STI","TUMI","LUSI"))
mean_beta$region <- factor(mean_beta$region,
                             levels = c("NCC","CCC","SCC"))
mean_beta$survey <- factor(mean_beta$survey,
                             levels = c("N. Copepod","S. Copepod","RREAS","CALCOFI", 
                                        "Upwelling","STI","TUMI","LUSI"))

dist_phe<-ggplot(mean_beta%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                  # filter(region=="SCC"|region=="CCC"|region=="NCC")%>%
         filter(Season=="Spring")%>%filter(lag==0), aes(x = median_beta, y=Index,col = as.factor(period))) +
  theme_bw() +
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3, position =ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8, position =ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75, position = ggstance::position_dodgev(height=dodge))+
  ggh4x::facet_grid2(region~survey) +
  ylab("Climate Index") +
  guides(col=guide_legend(title="Period"))+
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="bottom") 
dist_phe
dist_bio<-ggplot(mean_beta%>%
                   filter(survey=="CALCOFI"&period!=4)%>%
                   bind_rows(mean_beta%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
                   filter(Season=="Spring"&lag==0),
       aes(x =median_beta,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +

  ylab("") +
  scale_colour_manual(values = col)+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="none") 
dist_bio

dist_bio2<-ggplot(mean_beta%>%filter(survey=="CALCOFI")%>%
                    bind_rows(mean_beta%>%filter(survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod"&period!=1))%>%
         filter(Index=="TUMI"|Index=="LUSI"|Index=="STI")%>%
         filter(lag==0), 
       aes(x =median_beta,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.7,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  ylab("Upwelling Index") +
  scale_colour_manual(values = col)+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="none") 
dist_bio2 

pdf("Output/Distribution.pdf", 11,5) 
ggarrange(dist_phe,dist_bio,dist_bio2, ncol = 3, labels = c("A", "B", "C"), 
          widths=c(4,2,2), heights=c(2,2,1.75))
dev.off()


#### Posterior Means Intercepts ####
dodge<-0.6
mean_alpha<-Violin_Data%>%group_by(region,period,survey, Season,Index,lag)%>%
  summarise(median_alpha=median(alpha),sd_alpha_80=ci(alpha,ci = 0.95, method = "HDI"),
            sd_alpha_50=ci(alpha,ci = 0.5, method = "HDI"))
mean_alpha$region <- factor(mean_alpha$region,
                           levels = c("NCC","CCC","SCC"))
mean_alpha$survey <- factor(mean_alpha$survey,
                           levels = c("N. Copepod","S. Copepod","RREAS", "CALCOFI",
                                      "Upwelling","STI","TUMI","LUSI"))
mean_alpha$region <- factor(mean_alpha$region,
                           levels = c("NCC","CCC","SCC"))
mean_alpha$survey <- factor(mean_alpha$survey,
                           levels = c("N. Copepod","S. Copepod","RREAS","CALCOFI", 
                                      "Upwelling","STI","TUMI","LUSI"))

dist_phe<-ggplot(mean_alpha%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                   filter(Season=="Spring")%>%filter(lag==0), aes(x = median_alpha, y=Index,col = as.factor(period))) +
  theme_bw() +
  geom_errorbar(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),width = 0, lwd=0.3, position =ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_alpha_50$CI_low, xmax=sd_alpha_50$CI_high),width = 0, lwd=0.8, position =ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75, position = ggstance::position_dodgev(height=dodge))+
  ggh4x::facet_grid2(region~survey) +
  ylab("Climate Index") +
  guides(col=guide_legend(title="Period"))+
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="bottom") 
dist_phe
dist_bio<-ggplot(mean_alpha%>%
                   filter(survey=="CALCOFI"&period!=4)%>%
                   bind_rows(mean_alpha%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
                   filter(Season=="Spring"&lag==0),
                 aes(x =median_alpha,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_alpha_50$CI_low, xmax=sd_alpha_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  
  ylab("") +
  scale_colour_manual(values = col)+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="none") 
dist_bio
dist_bio2<-ggplot(mean_alpha%>%filter(survey=="CALCOFI")%>%
                    bind_rows(mean_alpha%>%filter(survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod"&period!=1))%>%
                    filter(Index=="TUMI"|Index=="LUSI"|Index=="STI")%>%
                    filter(lag==0), 
                  aes(x =median_alpha,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_alpha_50$CI_low, xmax=sd_alpha_50$CI_high),width = 0, lwd=0.7,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  ylab("Upwelling Index") +
  scale_colour_manual(values = col)+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="none") 
dist_bio2 

pdf("Output/DistributionIntercept.pdf", 11,5) 
ggarrange(dist_phe,dist_bio,dist_bio2, ncol = 3, labels = c("A", "B", "C"), 
          widths=c(4,2,2), heights=c(2,2,1.75))
dev.off()


overlap <- readRDS(here('data/Full_Results.rds'))
overlap%>%filter(Survey=="RREAS"&Index=="TUMI")

unique(overlap$Survey)
mean_alpha<-Violin_Data%>%group_by(region,period,survey, Season,Index,lag)%>%
  summarise(median_alpha=median(alpha),sd_alpha_80=ci(alpha,ci = 0.8, method = "HDI"),
            sd_alpha_50=ci(alpha,ci = 0.5, method = "HDI"))

dist_phe<-ggplot(mean_beta%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                   filter(Season=="Spring")%>%filter(period!=1&lag==0), aes(x = median_beta, y=Index, col = as.factor(period))) +
  theme_bw() +
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0.25, lwd=0.75)+
  #geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),alpha=0.7,width = 0)+
  geom_point(aes(col = as.factor(period)), cex=2.5)+
  ggh4x::facet_grid2(region~survey) +
  ylab("Climate Index") +
  scale_colour_manual(values = c(col2[2], col2[3]))+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="none") 

dist_bio<-ggplot(mean_beta%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
                   filter(Season=="Spring"&period!=1&period!=4&lag==0), 
                 aes(x =median_beta,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0.25, lwd=0.75)+
  #geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),alpha=0.7,width = 0)+
  geom_point(aes(col = as.factor(period)), cex=2.5)+
  facet_wrap(~survey,ncol=1) +
  ylab("") +
  scale_colour_manual(values = c(col2[2], col2[3]))+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="none") 


dist_bio2<-ggplot(mean_beta%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
                    filter(Index=="TUMI"|Index=="LUSI"|Index=="STI")%>%
                    filter(period!=1&lag==0), 
                  aes(x =median_beta,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0.25, lwd=0.75)+
  #geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),alpha=0.7,width = 0)+
  geom_point(aes(col = as.factor(period)), cex=2.5)+
  facet_wrap(~survey,ncol=1) +
  ylab("Upwelling Index") +
  scale_colour_manual(values = c(col2[2], col2[3]))+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="bottom") 


pdf("Output/Distribution.pdf", 11,6) 
ggarrange(dist_phe,dist_bio,dist_bio2, ncol = 3, labels = c("A", "B", "C"), 
          widths=c(4,2,2))
dev.off()


#### Posterior Means Winter ####

dist_phe<-ggplot(mean_beta%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                   filter(Season=="Winter")%>%filter(lag==0), aes(x = median_beta, y=Index,col = as.factor(period))) +
  theme_bw() +
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3, position =ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8, position =ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75, position = ggstance::position_dodgev(height=dodge))+
  ggh4x::facet_grid2(region~survey) +
  ylab("Climate Index") +
  guides(col=guide_legend(title="Period"))+
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="bottom") 
dist_phe
dist_bio<-ggplot(mean_beta%>%
                   filter(survey=="CALCOFI"&period!=4)%>%
                   bind_rows(mean_beta%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
                   filter(Season=="Winter"&lag==0),
                 aes(x =median_beta,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  
  ylab("") +
  scale_colour_manual(values = col)+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope") +
  theme(legend.position="none") 
dist_bio


pdf("Output/DistributionWinter.pdf", 8,5) 
ggarrange(dist_phe,dist_bio, ncol = 2, labels = c("A", "B"), 
          widths=c(4,2), heights=c(2,2))
dev.off()


#### Posterior Means Intercepts Winter ####

dist_phe<-ggplot(mean_alpha%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                   filter(Season=="Winter")%>%filter(lag==0), aes(x = median_alpha, y=Index,col = as.factor(period))) +
  theme_bw() +
  geom_errorbar(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),width = 0, lwd=0.3, position =ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_alpha_50$CI_low, xmax=sd_alpha_50$CI_high),width = 0, lwd=0.8, position =ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75, position = ggstance::position_dodgev(height=dodge))+
  ggh4x::facet_grid2(region~survey) +
  ylab("Climate Index") +
  guides(col=guide_legend(title="Period"))+
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="bottom") 
dist_phe
dist_bio<-ggplot(mean_alpha%>%
                   filter(survey=="CALCOFI"&period!=4)%>%
                   bind_rows(mean_alpha%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
                   filter(Season=="Winter"&lag==0),
                 aes(x =median_alpha,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_alpha_50$CI_low, xmax=sd_alpha_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  
  ylab("") +
  scale_colour_manual(values = col)+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="none") 
dist_bio
 

pdf("Output/DistributionInterceptWinter.pdf", 8,5) 
ggarrange(dist_phe,dist_bio,dist_bio2, ncol = 2, labels = c("A", "B"), 
          widths=c(4,2), heights=c(2,2))
dev.off()

#### Posterior Means Lag ####

dist_bio<-ggplot(mean_beta%>%
                   filter(survey=="CALCOFI"&period!=4)%>%
                   bind_rows(mean_beta%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
                   filter(Season=="Spring"&lag==1),
                 aes(x =median_beta,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  
  ylab("") +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="bottom") 
dist_bio


pdf("Output/DistributionLag.pdf", 6,5) 
ggarrange(dist_bio, ncol = 1)
dev.off()


#### Posterior Means Intercepts ####

dist_bio<-ggplot(mean_alpha%>%
                   filter(survey=="CALCOFI"&period!=4)%>%
                   bind_rows(mean_alpha%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
                   filter(Season=="Spring"&lag==1),
                 aes(x =median_alpha,y=Index, col=as.factor(period))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_alpha_50$CI_low, xmax=sd_alpha_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(aes(col = as.factor(period)), cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  
  ylab("") +
  scale_colour_manual(values = c(col[1],col[2], col[3]), name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022')) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Intercept") +
  theme(legend.position="bottom") 
dist_bio


pdf("Output/DistributionLagIntercept.pdf", 6,5) 
ggarrange(dist_bio, ncol = 1)
dev.off()


#### Difference Means ####
col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
mean_beta<-ratio_upwelling%>%group_by(region,survey, Season,Index,lag, Difference)%>%
  summarise(median_beta=median(beta_diff),sd_beta_80=ci(beta_diff,ci = 0.8, method = "HDI"),
            sd_beta_50=ci(beta_diff,ci = 0.5, method = "HDI"))%>%
mutate(region=ifelse(region=="Northern CC","NCC",
                     ifelse(region=="Southern CC", "SCC", "CCC")))
#ratio.up<-ratio.up%>%left_join(mean_beta)
mean_beta$region <- factor(mean_beta$region,
                          levels = c("NCC","CCC","SCC"))


diff_phe<-ggplot(mean_beta%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                   filter(Index=="NPGO"|Index=="PDO"|Index=="ONI"|Index=="NPH")%>%
                   filter(Season=="Spring")%>%filter(lag==0), aes(x = median_beta, y=Index, col=Difference)) +
  theme_bw() +
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3, position =ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8, position =ggstance::position_dodgev(height=dodge))+
  geom_point(cex=1.75, position = ggstance::position_dodgev(height=dodge))+
  ggh4x::facet_grid2(region~survey) +
  ylab("Climate Index") +
  guides(col=guide_legend(title="Period"))+
  scale_colour_manual(values=c(col4[4],col4[3],col4[1]), name="Period")+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope Difference") +
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position="bottom") 
diff_phe


mean_beta<-ratio_biological%>%group_by(region,survey, Season,Index,lag, Difference)%>%
  summarise(median_beta=median(beta_diff),sd_beta_80=ci(beta_diff,ci = 0.8, method = "HDI"),
            sd_beta_50=ci(beta_diff,ci = 0.5, method = "HDI"))
mean_beta$survey <- factor(mean_beta$survey,
                           levels = c("N. Copepod","S. Copepod","RREAS", "CALCOFI",
                                      "Upwelling","STI","TUMI","LUSI"))
mean_beta$region <- factor(mean_beta$region,
                           levels = c("NCC","CCC","SCC"))
diff_bio<-ggplot(mean_beta%>%filter(Season=="Spring"&lag==0)%>%
                   filter(Index=="PDO"|Index=="ONI"|Index=="NPGO"|Index=='NPH')%>%
                   filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod"),
                 aes(x =median_beta,y=Index, col=as.factor(Difference))) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.8,position = ggstance::position_dodgev(height=dodge))+
  geom_point(cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  ylab("") +
  scale_colour_manual(values = c(col4[1],col4[3],col4[4]))+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope Difference") +
  theme(legend.position="none") 
diff_bio



mean_beta<-ratio_index%>%group_by(region,survey, Season,Index,lag, Difference)%>%
  summarise(median_beta=median(beta_diff),sd_beta_80=ci(beta_diff,ci = 0.8, method = "HDI"),
            sd_beta_50=ci(beta_diff,ci = 0.5, method = "HDI"))

mean_beta$survey <- factor(mean_beta$survey,
                           levels = c("N. Copepod","S. Copepod", "RREAS","CALCOFI",
                                      "Upwelling","STI","TUMI","LUSI"))
mean_beta$region <- factor(mean_beta$region,
                           levels = c("NCC","CCC","SCC"))


diff_bio2<-ggplot(mean_beta%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
                    filter(Index=="TUMI"|Index=="LUSI"|Index=="STI")%>%
                    filter(lag==0), 
                  aes(x =median_beta,y=Index, col=Difference)) +
  theme_bw()+
  geom_errorbar(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),width = 0, lwd=0.3,position = ggstance::position_dodgev(height=dodge))+
  geom_errorbar(aes(xmin=sd_beta_50$CI_low, xmax=sd_beta_50$CI_high),width = 0, lwd=0.7,position = ggstance::position_dodgev(height=dodge))+
  geom_point(cex=1.75,position = ggstance::position_dodgev(height=dodge))+
  facet_wrap(~survey,ncol=1) +
  ylab("Upwelling Index") +
  scale_colour_manual(values = c(col4[1],col4[3],col4[4]))+
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Slope Difference") +
  theme(legend.position="none") 
diff_bio2 

pdf("Output/Difference_Distribution.pdf", 10.5,5) 
ggarrange(diff_phe,diff_bio,diff_bio2, ncol = 3, labels = c("A", "B", "C"), 
          widths=c(4,2,2), heights=c(2,2,1.75))
dev.off()

##### Posterior Mean Winter/Spring #####

dat<-mean_beta%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                  filter(lag==0)
dat$sd_beta_80.CI_high
spw_slope_A<-ggplot(dat, aes(x = median_beta, y=Index,lty=Season,col=Season, group=as.factor(period),shape=as.factor(period))) +
         theme_bw() +
  facet_grid(region~survey) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey")) +
  geom_errorbarh(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  scale_shape_manual(values=c(16,17,15),labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Climate Index",shape="Period")+
theme(legend.position = "none")

biodat<-mean_beta%>%
  filter(survey=="CALCOFI"&period!=4)%>%
  bind_rows(mean_beta%>%filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
  filter(lag==0)%>%
  filter(Index=="NPGO"|Index=="PDO"|Index=="ONI"|Index=="NPH")

biodatint<-mean_alpha%>%
  filter(survey=="CALCOFI"&period!=4)%>%
  bind_rows(mean_alpha%>%
              filter(survey=="RREAS"&period!=1&period!=4|survey=="N. Copepod"&period!=1&period!=4|survey=="S. Copepod"&period!=1&period!=4))%>%
  filter(lag==0)%>%
  filter(Index=="NPGO"|Index=="PDO"|Index=="ONI"|Index=="NPH")

spw_slope_B<-ggplot(biodat, aes(x = median_beta, y=Index,lty=Season,col=Season, group=as.factor(period),shape=as.factor(period))) +
  theme_bw() +
  facet_wrap(~survey,ncol=1) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey")) +
  geom_errorbarh(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  scale_shape_manual(values=c(16,17,15),labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Climate Index",shape="Period")

spw_slope_B

spw_int_A<-ggplot(mean_alpha%>%filter(survey=="TUMI"|survey=="STI"|survey=="LUSI")%>%
                    filter(lag==0), aes(x = median_alpha, y=Index,lty=Season,col=Season, group=as.factor(period),shape=as.factor(period))) +
  theme_bw() +
  facet_grid(region~survey) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey")) +
  geom_errorbarh(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  scale_shape_manual(values=c(16,17,15),labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Climate Index",shape="Period")+
  theme(legend.position = "none")
spw_int_A

spw_int_B<-ggplot(biodatint, aes(x = median_alpha, y=Index,lty=Season,col=Season, group=as.factor(period),shape=as.factor(period))) +
  theme_bw() +
  facet_wrap(~survey,ncol=1) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey")) +
  geom_errorbarh(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  scale_shape_manual(values=c(16,17,15),labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Climate Index",shape="Period")

spw_int_B


pdf("Output/WinterV2Slope.pdf", 11,5) 
ggarrange(spw_slope_A,spw_slope_B, ncol = 2, labels = c("A", "B"), 
          widths=c(4,3), heights=c(2,1.75))
dev.off()

pdf("Output/WinterV2Int.pdf", 11,5) 
ggarrange(spw_int_A,spw_int_B, ncol = 2, labels = c("A", "B"), 
          widths=c(4,3), heights=c(2,1.75))
dev.off()


#### lag alt ####

biodat2<-mean_beta%>%
  filter(survey=="CALCOFI"&period!=4)%>%
  bind_rows(mean_beta%>%filter(survey=="RREAS"&period!=1&period!=4))%>%
  filter(Season=="Spring")%>%
  filter(Index=="NPGO"|Index=="PDO"|Index=="ONI"|Index=="NPH")

biodatint2<-mean_alpha%>%
  filter(survey=="CALCOFI"&period!=4)%>%
  bind_rows(mean_alpha%>%
              filter(survey=="RREAS"&period!=1&period!=4))%>%
  filter(Season=="Spring")%>%
  filter(Index=="NPGO"|Index=="PDO"|Index=="ONI"|Index=="NPH")


lag_int<-ggplot(biodatint2, aes(x = median_alpha, y=Index,lty=as.factor(lag),col=as.factor(lag), group=as.factor(period),shape=as.factor(period))) +
  theme_bw() +
  facet_wrap(~survey,ncol=1) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey")) +
  geom_errorbarh(aes(xmin=sd_alpha_80$CI_low, xmax=sd_alpha_80$CI_high),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  scale_shape_manual(values=c(16,17,15),labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Climate Index",shape="Period",col="Lag",lty="Lag")

lag_int

lag_slope<-ggplot(biodat2, aes(x = median_beta, y=Index,lty=as.factor(lag),col=as.factor(lag), group=as.factor(period),shape=as.factor(period))) +
  theme_bw() +
  facet_wrap(~survey,ncol=1) +
  geom_point(alpha = 0.9,position = ggstance::position_dodgev(height=dodge),cex=2.5) +
  scale_colour_manual(values=c("palevioletred", "darkgrey")) +
  geom_errorbarh(aes(xmin=sd_beta_80$CI_low, xmax=sd_beta_80$CI_high),position = ggstance::position_dodgev(height=dodge),height=0.5, alpha=0.75)+
  scale_linetype_manual(values=c(2,1))+
  scale_shape_manual(values=c(16,17,15),labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Climate Index",shape="Period",col="Lag",lty="Lag")

lag_slope

pdf("Output/Lag_Slopw.pdf", 8,5) 
lag_slope
dev.off()

pdf("Output/Lag_Int.pdf", 8,5) 
lag_int
dev.off()
