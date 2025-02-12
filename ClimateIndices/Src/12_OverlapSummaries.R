library(gt)
overlap<-readRDS('data/overlap_Results_beta.rds')
overlapalpha<-readRDS('data/overlap_Results_alpha.rds')

overlapTab<-read.csv('climateov.csv')

#### Slope Overlap ####
overlap%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0)
overlap.phe%>%filter(ov<0.1)%>%
  group_by(Index)

#### Overlap Table ####
climateov<-bind_rows(overlapalpha%>%filter(ov<0.05, season=="Spring", period1!=4&offset==0)%>%
  mutate("Time Period 1" = ifelse(period1==1, "1967 - 1988", ifelse(period1==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Time Period 2"= ifelse(period2==1, "1967 - 1988", ifelse(period2==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Posterior Overlap"=round(ov,2))%>%
  rename("Driver"=Index, "Response"=Survey)%>%
  dplyr::select("Driver","region", "Response", "Time Period 1", "Time Period 2", "Posterior Overlap"),
overlap.phe.alpha%>%filter(ov<0.05)%>%
  mutate("Time Period 1" = ifelse(period1==1, "1967 - 1988", ifelse(period1==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Time Period 2"= ifelse(period2==1, "1967 - 1988", ifelse(period2==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Posterior Overlap"=round(ov,2))%>%
  rename("Driver"=Index, "Response"=survey)%>%
  dplyr::select("Driver", "region","Response", "Time Period 1", "Time Period 2", "Posterior Overlap"))

ovup<-overlap.bioup.alpha%>%filter(ov<0.05)%>%
  mutate("Time Period 1" = ifelse(period1==1, "1967 - 1988", ifelse(period1==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Time Period 2"= ifelse(period2==1, "1967 - 1988", ifelse(period2==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Posterior Overlap"=round(ov,2))%>%
  rename("Driver"=Index, "Response"=Survey)%>%
  dplyr::select("Driver", "Response", "Time Period 1", "Time Period 2", "Posterior Overlap")

write.csv(ovup,"ovup.csv")
write.csv(climateov,"climateov.csv")


slopeov<-bind_rows(overlap%>%filter(ov<0.05, season=="Spring", period1!=4&offset==0)%>%
  mutate("Time Period 1" = ifelse(period1==1, "1967 - 1988", ifelse(period1==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Time Period 2"= ifelse(period2==1, "1967 - 1988", ifelse(period2==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Posterior Overlap"=round(ov,2))%>%
  rename("Driver"=Index, "Response"=Survey)%>%
  dplyr::select("Driver", "Response", "Time Period 1", "Time Period 2", "Posterior Overlap"),
overlap.phe%>%filter(ov<0.05)%>%
  mutate("Time Period 1" = ifelse(period1==1, "1967 - 1988", ifelse(period1==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Time Period 2"= ifelse(period2==1, "1967 - 1988", ifelse(period2==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Posterior Overlap"=round(ov,2))%>%
  rename("Driver"=Index, "Response"=survey)%>%
  dplyr::select("Driver", "Response", "Time Period 1", "Time Period 2", "Posterior Overlap"))
write.csv(slopeov,"slopeov.csv")
overlap.bioup%>%filter(ov<0.05)%>%
  mutate("Time Period 1" = ifelse(period1==1, "1967 - 1988", ifelse(period1==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Time Period 2"= ifelse(period2==1, "1967 - 1988", ifelse(period2==2, "1989 - 2012", "2013 - 2022")))%>%
  mutate("Posterior Overlap"=round(ov,2))%>%
  rename("Driver"=Index, "Response"=Survey)%>%
  dplyr::select("Driver", "Response", "Time Period 1", "Time Period 2", "Posterior Overlap")





overlapalpha%>%filter(ov>0.1, season=="Spring", period1!=4&offset==0&Survey == "CALCOFI")
overlap%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0)

overlapalpha%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0&Survey == "RREAS"&period1!=3)
overlap%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0&Survey == "RREAS"&period1!=3)


overlapalpha%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0&Survey == "N. Copepods")
overlapalpha%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0&Survey == "S. Copepods")

### TUMI Intercepts ###
overlap.phe.alpha%>%filter(period1==3&period2==2&survey=="TUMI"&region=="Southern CC")%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe.alpha%>%filter(period1==1&period2==3&survey=="TUMI"&region=="Southern CC")%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe.alpha%>%filter(period1==1&period2==3&survey=="TUMI"&region=="Central CC")%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe.alpha%>%filter(period1==3&period2==2&survey=="TUMI"&region=="Central CC")%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe%>%filter(ov<0.1)




overlap.phe%>%filter(period1==3&period2==2)%>%
  group_by(Index)%>%
  summarise(mean=mean(ov), sd=sd(ov))

overlap.phe%>%filter(period1==3&period2==2)%>%
  group_by(Index, region)%>%
  summarise(mean=mean(ov), sd=sd(ov))