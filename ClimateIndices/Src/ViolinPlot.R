
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




##### Upwelling Violin ####
# plotting functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }

col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ratio.up<-ratio_upwelling%>%filter(region!="GoA", season=="Spring"|Season=="Spring", lag==0)%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                      ifelse(region=="Southern CC", "SCC", "CCC")))
ratio.up$region <- factor(ratio.up$region,
                levels = c("SCC","CCC","NCC"))
violin.up <-ggplot(ratio.up%>%filter(survey=="Upwelling"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("Upwelling Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.up

#### Violin Phenology ####

violin.sti <-ggplot(ratio.up%>%filter(survey=="STI"&Season=="Spring"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("STI Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.sti

violin.tumi <-ggplot(ratio.up%>%filter(survey=="TUMI"&Season=="Spring"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("TUMI Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.tumi


violin.lusi <-ggplot(ratio.up%>%filter(survey=="LUSI"&Season=="Spring"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("LUSI Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.lusi

colnames(ratio.up)
ratio.up$Difference2
violin.phe <-ggplot(ratio.up%>%filter(survey!="Upwelling"&Season=="Spring"&Difference2=="Era 3 - Era 2"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~survey)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~survey) +
  ylab("Phyical Conditions Posterior Difference") +
  xlab("") +
  theme(legend.position="none")+
  geom_hline(aes(yintercept=0), size=0.3) 
violin.phe


#### Violin Winter Phenology ####
ratio.upw<-ratio_upwelling%>%filter(region!="GoA", season=="Winter"|Season=="Winter", lag==0)%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                       ifelse(region=="Southern CC", "SCC", "CCC")))
ratio.upw$region <- factor(ratio.upw$region,
                          levels = c("SCC","CCC","NCC"))
violin.stiw <-ggplot(ratio.upw%>%filter(survey=="STI"&Season=="Winter"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("STI Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.stiw

violin.tumiw <-ggplot(ratio.upw%>%filter(survey=="TUMI"&Season=="Winter"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("TUMI Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.tumiw


violin.lusiw <-ggplot(ratio.upw%>%filter(survey=="LUSI"&Season=="Winter"), aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  facet_wrap(~Difference)+
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("LUSI Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.lusiw
#### Violin Biological 2 ####
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
violin.bio <-ggplot(ratio.bio%>%filter(Difference=="2013:2023 - 1989:2012"), 
                    aes(x=survey, y=beta_diff, fill=region)) +
  theme_bw() +
  scale_fill_manual(values=c(col4[3:1],col4[1]), name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("Biological Posterior \n Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="none")
violin.bio

#### Violin Biological ####
unique(ratio_index$Index)

#ratio.bio<-ratio.bio%>%filter(Index=='PDO'|Index=='ONI'|Index=='NPH'|Index=='NPGO'|
#                                survey=='CALCOFI'|survey=='RREAS'|survey=='N. Copepod'|survey=="S. Copepod")
ratio_index$region <- factor(ratio_index$region,
                           levels = c("SCC","CCC","NCC"))
ratio_index$survey <- factor(ratio_index$survey,
                           levels = c("CALCOFI","RREAS","N. Copepod","S. Copepod", 
                                      "Upwelling","STI","TUMI","LUSI"))
violin.bio2 <-ggplot(ratio_index%>%filter(Difference=="2013:2023 - 1989:2012" & Index!="Upwelling"), 
                    aes(x=survey, y=beta_diff, fill=region)) +
  theme_bw() +
  scale_fill_manual(values=c(col4[3:1],col4[1]), name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("Biological Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom",axis.text.y=element_blank())
violin.bio2


pdf("Output/ViolinV3.pdf", 11,6) 
ggarrange(violin.phe,violin.bio,violin.bio2,labels = c("A", "B", "C"),  
          font.label = list(size = 12, face="plain"), ncol=3,widths=c(3,2,2))
dev.off()


pdf("Output/ViolinV4.pdf", 11,6) 
ggarrange(violin.phe,violin.bio, ncol = 3, labels = c("A", "B", "C"),
          ggarrange(violin.bio2,z.plot,nrow=2,labels = c("",""), 
                    heights = c(6.25,1)), widths=c(3,2,1.5))
dev.off()


violin.bio
pdf(file = "Output/Figures/violin.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
violin.up
dev.off()
