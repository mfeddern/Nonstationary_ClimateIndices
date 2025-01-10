overlap<-readRDS('data/overlap_Results_beta.rds')
overlapalpha<-readRDS('data/overlap_Results_alpha.rds')

overlapalpha%>%filter(ov>0.1, season=="Spring", period1!=4&offset==0&Survey == "CALCOFI")
overlap%>%filter(ov<0.1, season=="Spring", period1!=4&offset==0&Survey == "CALCOFI")

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