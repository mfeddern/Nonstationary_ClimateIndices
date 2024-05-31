library(ggrepel)
library(strucchange)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(PNWColors)
library(ggspatial)

# reading in the data 
schroeder.nph <-read.csv('data/physical/year_mon_area_max_x_y_lon_lat_2023.csv')%>%
  mutate(era= ifelse(Year<=1988,1, ifelse(Year>2012,3,2)))%>% #assigning eras
  mutate(Year_win = if_else(Month == 11|Month ==12, Year+1, Year))%>% #creating a year offset fro winter
  mutate(dec.yr =as.numeric(as.character(Year)) + (as.numeric(Month)-0.5)/12, #creating a decimal year
         Area=Area/1000000)%>%  #converting area to smaller units
  mutate(area.anom = Area-mean(Area)) #creating an area anomaly

plot(schroeder.nph$dec.yr, schroeder.nph$Area, type="l") #check plot to visualize area through time
plot(schroeder.nph$dec.yr, schroeder.nph$Max, type="l") #check plot to visualize intensity through time


ggplot(data=schroeder.nph, 
       aes(Month,Max, group=Year))+
  facet_wrap(.~era, ncol = 3) +
  #geom_line()+
  geom_smooth(se=F, col='grey')+
  geom_smooth(aes(group=era))+
  theme_bw()

ggplot(data=schroeder.nph, 
       aes(Month,Area, group=Year))+
  facet_wrap(.~era, ncol = 3) +
  #geom_line()+
  geom_smooth(se=F, col='grey')+
  geom_smooth(aes(group=era))+
  theme_bw()

ggplot(data=schroeder.nph%>%filter(era==3), 
       aes(x,y, col=Max))+
  facet_wrap(.~Month, ncol = 3) +
  geom_point()+
  scale_colour_gradientn(colours = colorspace::diverge_hcl(7), limits=c(1010, 1030))+
  #geom_text()+
  stat_ellipse(level = 0.9) +
  ggtitle("2013 - 2023")+
  theme_bw()

ggplot(data=schroeder.nph%>%filter(era==2), 
       aes(x,y, col=Max))+
  facet_wrap(.~Month, ncol = 3) +
  geom_point()+
  scale_colour_gradientn(colours = colorspace::diverge_hcl(7), limits=c(1010, 1030))+
  #geom_text()+
  ggtitle("1989 - 2012")+
  stat_ellipse(level = 0.9) +
  theme_bw()
#visulaizing location with ellipses
ggplot(data=schroeder.nph%>%filter(Month==4|Month==5|Month==6), 
       aes(x,y, label = Year,col=as.factor(era)))+
  facet_wrap(.~Month, ncol = 3, scales='free') +
  geom_point()+
  geom_text()+
  stat_ellipse(level = 0.9) +
  theme_bw()
#examining monthrs through time
ggplot(data=schroeder.nph,
       aes(dec.yr,Max))+
  facet_wrap(.~Month, ncol = 3, scales='free') +
  geom_line()+
  theme_bw()

#creating a spring dataset
spring.schroeder <- schroeder.nph%>%
  filter(Month==4|Month==5|Month==6)%>%
  group_by(Year)%>%
  summarise(mean.x=mean(x), mean.y=mean(y),mean.max=mean(Max),mean.area=mean(Area),
            mean.area.anom=mean(area.anom))%>%
  mutate(era= ifelse(Year<=1988,1, ifelse(Year>2012,3,2)))%>%
  mutate(xy=sqrt(mean.x^2+mean.y^2), era=as.factor(era),
         era.lab = ifelse(era==1, '1967 - 1988', ifelse(era==2, "1989 - 2012", "2013 - 2023")))

#creating a winter dataset
winter.schroeder <- schroeder.nph%>%
  filter(Month==11|Month==12|Month==1|Month==3|Month==2)%>%
  group_by(Year_win)%>%
  summarise(mean.x=mean(x), mean.y=mean(y),mean.max=mean(Max),mean.area=mean(Area),
            mean.area.anom=mean(area.anom))%>%
  mutate(era= ifelse(Year_win<=1988,1, ifelse(Year_win>2012,3,2)))%>%
  mutate(xy=sqrt(mean.x^2+mean.y^2), era=as.factor(era),
         era.lab = ifelse(era==1, '1967 - 1988', ifelse(era==2, "1989 - 2012", "2013 - 2023")))

#calculating mean and SD values for location
means <- spring.schroeder%>%
  group_by(era.lab)%>%
  summarise(x=mean(mean.x),y=mean(mean.y),sd.x=sd(mean.x), sd.y=sd(mean.y),
            area=mean(mean.max),intensity=mean(mean.area),
            sd.area=sd(mean.max), sd.intensity=sd(mean.area),
            Year=0)%>%
  rename(mean.x=x, mean.y=y,mean.area=intensity, mean.max=area )
theme_set(theme_classic())


means2 <- winter.schroeder%>%
  group_by(era.lab)%>%
  summarise(x=mean(mean.x),y=mean(mean.y),sd.x=sd(mean.x), sd.y=sd(mean.y),
            area=mean(mean.max),intensity=mean(mean.area),sd.area=sd(mean.max), sd.intensity=sd(mean.area),
            Year=0)%>%
  rename(mean.x=x, mean.y=y,mean.area=intensity, mean.max=area )
theme_set(theme_classic())
#spring.schroeder<-winter.schroeder
#plotting location
col<-pnw_palette("Sunset2",3,type="discrete")
a.plot <-ggplot(data=spring.schroeder,aes(abs(mean.x-360),mean.y, label=Year,group=era.lab,col=era.lab))+
  geom_point(alpha=0.4)+
  #geom_text(col='grey')+
  #ggtitle("Center of North Pacific High") +
  geom_point(data=means)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means,aes(ymin = mean.y-sd.y, ymax=  mean.y+sd.y), width=0.5) +
  geom_errorbar(data=means,aes(xmin = abs(mean.x-360)-sd.x, xmax=  abs(mean.x-360)+sd.x), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  geom_hline(yintercept=mean(spring.schroeder$mean.y),lty=2, col='grey')+
  geom_vline(xintercept=mean(spring.schroeder$mean.x), lty=2, col='grey')+
  scale_x_reverse(lim=c(147,135))+
  ylab('Latitude (ºN)')+
  theme_bw() +
  xlab('Longitude (ºW)')+
  theme(legend.position = "none")
a.plot
means2<-means2%>%mutate(Year_win=Year)
a.plotw <-ggplot(data=winter.schroeder,aes(abs(mean.x-360),mean.y, label=Year_win,group=era.lab,col=era.lab))+
  geom_point(alpha=0.4)+
  #geom_text(col='grey')+
  #ggtitle("Center of North Pacific High") +
  geom_point(data=means2)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means2,aes(ymin = mean.y-sd.y, ymax=  mean.y+sd.y), width=0.5) +
  geom_errorbar(data=means2,aes(xmin = abs(mean.x-360)-sd.x, xmax=  abs(mean.x-360)+sd.x), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  geom_hline(yintercept=mean(winter.schroeder$mean.y),lty=2, col='grey')+
  geom_vline(xintercept=mean(winter.schroeder$mean.x), lty=2, col='grey')+
  scale_x_reverse(lim=c(147,125))+
  ylab('Latitude (ºN)')+
  theme_bw() +
  xlab('Longitude (ºW)')+
  theme(legend.position = "none")
a.plotw

k.plot <-ggplot(data=spring.schroeder,aes(abs(mean.x-360),mean.y, label=Year,group=era.lab,col=era.lab))+
  geom_point(alpha=0.4)+
  ggtitle("Center of North Pacific High") +
  geom_point(data=means)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means,aes(ymin = mean.y-sd.y, ymax=  mean.y+sd.y), width=0.5) +
  geom_errorbar(data=means,aes(xmin = abs(mean.x-360)-sd.x, xmax=  abs(mean.x-360)+sd.x), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  geom_hline(yintercept=mean(spring.schroeder$mean.y),lty=2, col='grey')+
  geom_vline(xintercept=mean(spring.schroeder$mean.x), lty=2, col='grey')+
  scale_x_reverse(lim=c(147,135))+
  ylab('Latitude (ºN)')+
  theme_bw() +
  xlab('Longitude (ºW)')
k.plot

j.plot <-ggplot(data=spring.schroeder,aes(y=mean.max,x=mean.area, label=Year,group=era.lab,col=era.lab))+
  geom_point(alpha=0.4)+
 # ggtitle("North Pacific High\n Areal Extent and Intensity") +
  geom_point(data=means)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means,aes(xmin = mean.area-sd.area, xmax= mean.area+sd.area), width=0.5) +
  geom_errorbar(data=means,aes(ymin = mean.max-sd.intensity, ymax= mean.max+sd.intensity), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  #scale_x_reverse(lim=c(147,135))+
  ylab('North Pacific High \n Intensity (hPa)')+
  theme_bw()+
  geom_hline(yintercept=mean(spring.schroeder$mean.max),lty=2, col='grey')+
  geom_vline(xintercept=mean(spring.schroeder$mean.area), lty=2, col='grey')+
  xlab(expression("North Pacific High Area "~(10^6 ~km^2)))
j.plot

j.plotw <-ggplot(data=winter.schroeder,aes(y=mean.max,x=mean.area, label=Year_win,group=era.lab,col=era.lab))+
  geom_point(alpha=0.4)+
  # ggtitle("North Pacific High\n Areal Extent and Intensity") +
  geom_point(data=means2)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means2,aes(xmin = mean.area-sd.area, xmax= mean.area+sd.area), width=0.5) +
  geom_errorbar(data=means2,aes(ymin = mean.max-sd.intensity, ymax= mean.max+sd.intensity), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  #scale_x_reverse(lim=c(147,135))+
  ylab('North Pacific High \n Intensity (hPa)')+
  theme_bw()+
  geom_hline(yintercept=mean(winter.schroeder$mean.max),lty=2, col='grey')+
  geom_vline(xintercept=mean(winter.schroeder$mean.area), lty=2, col='grey')+
  xlab(expression("North Pacific High Area "~(10^6 ~km^2)))
j.plotw

h.plot <-ggplot(data=spring.schroeder,aes(y=mean.max,x=mean.area, label=Year,group=era.lab,col=era.lab))+
  geom_point()+
  ggtitle("") +
  #  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(era.lab))) +
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  #geom_errorbar(data=means,aes(ymin = mean.area-sd.area, ymax= mean.area+sd.area), width=0.5) +
 # geom_errorbar(data=means,aes(xmin = mean.max-sd.intensity, xmax=  mean.max+sd.intensity), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  ylab('North Pacific High \n Intensity (hPa)')+
  theme_bw() +
  geom_text_repel(data=subset(spring.schroeder, era==3),
            aes(y=mean.max,x=mean.area,label=Year),col='black', max.overlaps = Inf, position = position_jitter(seed = 5))+
  geom_hline(yintercept=mean(spring.schroeder$mean.max),lty=2, col='grey')+
  geom_vline(xintercept=mean(spring.schroeder$mean.area), lty=2, col='grey')+
  xlab(expression("North Pacific High Area "~(10^6 ~km^2)))
h.plot

g.plot <-ggplot(data=spring.schroeder,aes(x=mean.max,y=mean.area, group=era.lab,col=era.lab))+
  #ggtitle("Center of North Pacific High") +
  geom_point(col='grey',data=subset(spring.schroeder, era!=3))+
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(era.lab))) +
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
 # geom_errorbar(data=means,aes(ymin = mean.y-sd.y, ymax=  mean.y+sd.y), width=0.5) +
 # geom_errorbar(data=means,aes(xmin = mean.x-sd.x, xmax=  mean.x+sd.x), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  xlab('North Pacific High \n Intensity (hPa)')+
  theme_bw() +
  geom_text(data=subset(spring.schroeder, era==3),
            aes(x=mean.max,y=mean.area,label=Year),col='black')+
  ylab(expression("North Pacific High Area "~(10^6 ~km^2)))
g.plot
#check TS to look at variables through time for spring
plot(spring.schroeder$Year, spring.schroeder$xy, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.y, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.x, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.max, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.area.anom, type="l")

#making a TS for break point analysis Y
y.ts <- ts(data=spring.schroeder%>%select(mean.y), 1967, 2023, frequency=1)
# fit breakpoint model
bp.y <- breakpoints(y.ts ~ 1)
summary(bp.y) 

#making a TS for break point analysis X
x.ts <- ts(data=spring.schroeder%>%select(mean.x), 1967, 2023, frequency=1)
# fit breakpoint model
bp.x <- breakpoints(x.ts ~ 1)
summary(bp.x) 

#making a TS for break point analysis Area
area.ts <- ts(data=spring.schroeder%>%select(mean.area.anom), 1967, 2023, frequency=1)
# fit breakpoint model
bp.area <- breakpoints(area.ts ~ 1)
summary(bp.area) 

#making a TS for break point analysis Intensity
max.ts <- ts(data=spring.schroeder%>%select(mean.max), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max) 

#creating a model to plot y breakpoints
mod <- lm(mean.y ~ as.factor(era), data=spring.schroeder)
# and get predicted values of the model to plot
pred <- predict(mod, se=T, newdata = spring.schroeder)
spring.schroeder$mean <- pred$fit

# now save for a combined plot

# set the colors to use - colorblind pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_classic())

b.plot <- ggplot(data=spring.schroeder, aes(Year, mean.y)) +
  geom_line(size=0.2) +
  xlab("")+
  geom_line(aes(Year, mean), color=cb[6], size=0.4) + 
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  ylab("y (ºN)") + ggtitle("Center of North Pacific High (ºN)") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  theme_bw() +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1966,2024)
b.plot


# and calculate standard deviation over 11-year rolling windows for intensity
NPH.max.sd <- rollapply(spring.schroeder%>%dplyr::select(mean.max), 11, sd, fill=NA)
plot(1967:2023, NPH.max.sd , type="l") #check plot

# and calculate standard deviation over 11-year rolling windows for area
NPH.area.sd <- rollapply(spring.schroeder%>%dplyr::select(mean.area.anom), 11, sd, fill=NA)
plot(1967:2023, NPH.area.sd, type="l") #check plot
# now fit a non-parametric regression

# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.max.sd))
# fit the model
mod <- gam(mean.max ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

c.plot <- ggplot(plot.dat, aes(year, mean.max)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4)  +  theme_bw() +
  ylab("Standard deviation \n (hPa)") +
  xlab("")+
 # ggtitle("North Pacific High Intensity Variability (Spring)") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
c.plot

# and calculate standard deviation over 11-year rolling windows for intensity
NPH.max.sd <- rollapply(winter.schroeder%>%dplyr::select(mean.max), 11, sd, fill=NA)
plot(1967:2024, NPH.max.sd , type="l") #check plot

# and calculate standard deviation over 11-year rolling windows for area
NPH.area.sd <- rollapply(winter.schroeder%>%dplyr::select(mean.area.anom), 11, sd, fill=NA)
plot(1967:2024, NPH.area.sd, type="l") #check plot
# now fit a non-parametric regression

# first, make a data frame
plot.dat <- data.frame(year=1972:2019, sd=na.omit(NPH.max.sd))
# fit the model
mod <- gam(mean.max ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

c.plotw <- ggplot(plot.dat, aes(year, mean.max)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4)  +  theme_bw() +
  ylab("Standard deviation \n (hPa)") +
  xlab("")+
  # ggtitle("North Pacific High Intensity Variability (Spring)") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
c.plotw

# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.area.sd))
# fit the model
mod <- gam(mean.area.anom ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

e.plot <- ggplot(plot.dat, aes(year, mean.area.anom)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4) +   theme_bw() +
  ylab(expression("Standard deviation "~(10^6 ~km^2))) +
  xlab("Year")+
  ggtitle("North Pacific High Area Variability") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
e.plot

#TS for summer data
max.ts <- ts(data=summer.schroeder%>%select(mean.max), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max)

NPH.max.sd <- rollapply(summer.schroeder%>%select(mean.max), 11, sd, fill=NA)
plot(1967:2023, NPH.max.sd , type="l")

# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.max.sd))
#plot.dat <- data.frame(year=1969:2021, sd=na.omit(NPH.max.sd))

# fit the model
mod <- gam(mean.max ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

d.plot <- ggplot(plot.dat, aes(year, mean.max)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4) +   theme_bw() +
  ylab("Standard deviation (hPa)") +
  xlab("Year")+
  ggtitle("North Pacific High Intensity Variability (Summer)") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
d.plot

png("Output/Fig 1.png", 6, 11, units="cm", res=300) 
ggarrange(a.plot, b.plot, c.plot, d.plot,labels = c("A", "B", "C", "D"),  
          font.label = list(size = 10, face="plain"), nrow=4)
dev.off()

pdf("Output/Fig 1.pdf", 6,11) 
ggarrange(a.plot, g.plot, b.plot, c.plot,labels = c("A", "B", "C", "D"),  
          font.label = list(size = 12, face="plain"), nrow=4)
dev.off()

pdf("Output/Fig 1v1.pdf", 6,11) 
ggarrange(a.plot, b.plot, labels = c("A", "B"),  
          font.label = list(size = 12, face="plain"), nrow=2)
dev.off()

pdf("Output/Fig 1v2.pdf", 6,11) 
ggarrange( g.plot, c.plot,d.plot,labels = c("A", "B", "C"),  
          font.label = list(size = 12, face="plain"), nrow=3)
dev.off()

pdf("Output/Fig 1v3.pdf", 5,7) 
ggarrange( a.plot, h.plot,c.plot,labels = c("A", "B", "C"),  
           font.label = list(size = 12, face="plain"), nrow=3)
dev.off()
#TS for winter data
max.ts <- ts(data=winter.schroeder%>%select(mean.max), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max)

NPH.max.sd <- rollapply(winter.schroeder%>%select(mean.max), 11, sd, fill=NA)
plot(1967:2023, NPH.max.sd , type="l")

# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.max.sd))
#plot.dat <- data.frame(year=1969:2021, sd=na.omit(NPH.max.sd))

# fit the model
mod <- gam(mean.max ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

f.plot <- ggplot(plot.dat, aes(year, mean.max)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4) +   theme_bw() +
  ylab("Standard deviation (hPa)") +
  xlab("Year")+
  ggtitle("North Pacific High Intensity Variability (Summer)") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
f.plot

z.plot<-ggplot()
  
  z.plot
##### Jumbo Plot ####

#### Map Data #####
col2<-pnw_palette("Sunset2",4,type="discrete")
col<-pnw_palette("Sunset2",3,type="discrete")
col3<-pnw_palette("Sunset2",8,type="continuous")

col<-pnw_palette("Sunset2",3,type="discrete")
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
climate_dat_cop <-readRDS(here('data/physical/climate_dat_cop.rds'))
bakunsites <- read.csv(here('data/physical/Bakun/MapLocations2.csv'))%>%
  mutate(longitude=longitude)
sites <- st_as_sf(data.frame(bakunsites[,1:2]), coords = c("longitude","latitude"), crs = 4326, 
                  agr = "constant")

#### Making the Map #####
world <- st_as_sf(map('world', plot=F, fill=T)) #base layer for land masses

map<-ggplot() +
  geom_polygon(aes(x=c(-105, -113, -127,-105,-105),
                   y=c(22.1, 22.1,34.4486,34.4486,20)),
               fill='white', col='black',alpha=0.6)+
  geom_polygon(aes(x=c(-110, -127,-130,-110,-110), 
                  # y=c(34.4486,34.4486,40.4401,40.4401,34.4486)),
               y=c(34.4486,34.4486,43,43,34.4486)),

               fill='white', col='black',alpha=0.6)+
  geom_polygon(aes(x=c(-110, -130,-130,-110,-110), 
                  # y=c(40.4401,40.4401,49.5,49.5,40.4401)),
               y=c(43,43,49.5,49.5,43)),
               fill='white', col='black',alpha=0.6)+
  geom_sf(data = world)+
  annotate("rect", xmin= -121.5, xmax = -109, ymin = 42, ymax = 48.8, 
           fill = 'white', col='black',size = 0.8, lwd=0.2) +
  geom_sf(fill='grey95') +
  geom_sf(data = sites, size = c(rep(2,68+35+12), rep(3,2)), 
          shape = c(rep(24,68), rep(21,35),rep(23,12),rep(22,2)), 
          col = c(rep('black',68+35+14)), 
          fill = c(rep(col3[3],68), rep(col3[7],35),rep(col3[5],12),rep(col3[8],2))) +
  coord_sf(xlim = c(-132, -108), ylim = c(26, 50), expand = FALSE)+
  ylab(" ")+
  xlab(" ")+
  annotation_scale()+
  annotation_north_arrow(which_north = "true",pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"))+
  annotate(geom = "text", x = c(-127.5,-118.5,-129), y = c(39,28,46.5), 
           label = str_wrap(c("Central", "Southern","Northern"), width = 20),
           fontface = "italic", color = "grey22", size = 3.75, angle=c('285', '311','270')) +
  annotate(geom = "text", x = c(-114,-114,-114,-114,-114), y = c(48,46.5,45, 44,43), 
           label = str_wrap(c("Upwelling Data","CC Regions", "Newport Line","CalCOFI", "RREAS"), width = 22),
           color = "grey22", size =3.5) +
  annotate(geom = "text", x = c(-120.5,-117), y = c(41,35), 
           label = str_wrap(c("Cape Mendocino","Point Conception"), width = 20),
           fontface = "italic", color = "grey22", size = 3) +
  annotate("rect", xmin= -121, xmax = -119, ymin = 46, ymax = 47, 
           fill = 'white', col='black',size = 0.8, lwd=0.5) +
  annotate("line", x= c(-124.1, -124.65), y = c(44.652, 44.652),col=col2[1],size = 0.8, lwd=1) +
  annotate("line", x= c(-120.5, -119.5), y = c(45, 45),col=col2[1],size = 0.8, lwd=1) +
  
  theme(panel.background = element_rect(fill = "lightsteelblue2"),
        panel.border = element_rect(fill = NA),panel.grid.major = element_line(colour = "transparent"))

map 

#### SLP Plots ####


ddd <-readRDS('data/physical/correlation_analysis_indices2.rds')%>%
  mutate(index = case_when(grepl("PDO", analysis) ~ "PDO",
                           grepl("NPGO", analysis) ~ "NPGO",
                           grepl("ONI", analysis) ~ "ONI",
                           grepl("NPH", analysis) ~ "NPH"))%>%
  mutate(era = case_when(grepl("1967 - 1988", analysis) ~ "Era 1",
                         grepl("1989 - 2012", analysis) ~ "Era 2",
                         grepl("2013 - 2023", analysis) ~ "Era 3"))%>%
  mutate(period = case_when(grepl("1967 - 1988", analysis) ~ "1967 - 1988",
                            grepl("1989 - 2012", analysis) ~ "1989 - 2012",
                            grepl("2013 - 2023", analysis) ~ "2013 - 2023"))


world <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
world3<-world%>%filter(ID=="Canada"|ID=="USA"|ID=="Mexico") #pulling ou USA,M and C
world3<-st_union(world3) #making them one multipoly object
buff<-buffer(vect(world3),width=130000) #adding a buffer
buff2<- st_as_sf(buff) #making it a sf object

buffer<- data.frame(st_coordinates(buff2))%>% #creating the buffer polygon
  filter(Y>30&Y<50)%>%
  mutate(X=X+360)%>%
  filter((X>220&X<250))%>%
 # filter(L1==1,L2==1,L3==1020)%>%
 # select(X,Y)%>%
  mutate(x=X,y=Y)
unique(coast$L1)
unique(coast$L2)
unique(coast$L3)
coast2<-st_coordinates(st_as_sf(buffer(vect(world3),width=1)))
coast<- data.frame(coast2)%>%
  filter(Y>30&Y<55)%>%
  mutate(X=X+360)%>%
  filter((X>220&X<250))%>%
  filter((L2==113))%>%
 # filter(L1==1,L2==1,L3==1020)%>%|L2==151 OR56
 # select(X,Y)%>%
  purrr::map_df(rev)%>%
  mutate(x=X,y=Y)

#no 31
buffer<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>31&y<48.5)
bufferNCC<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>40.4401&y<48.5)
bufferCCC<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>34.4486&y<40.5401)
bufferSCC<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>31&y<34.5486)
analysis2<-data.frame(analysis=unique(X_PDO_SLP$analysis), analysis2=c("1967 - 1988","1989 - 2012","2013 - 2023","X",
                                                            "1967 - 1988","1989 - 2012","2013 - 2023", "X",
                                                            "1967 - 1988", "1989 - 2012","2013 - 2023","X",
                                                            "1967 - 1988","1989 - 2012","2013 - 2023","X"),
                      Index=c(rep('PDO', 4), rep('ONI',4), rep('NPGO',4), rep('NPH',4))) 
analysis2<-merge(X_PDO_SLP,analysis2)%>%
  filter(analysis2!='X')

SLP_coast2<- ggplot() + 
 geom_raster(data=na.omit(ddd%>%filter(var=="SLP")), aes(x=longitude,y=latitude,fill = coefficient)) + 
   facet_grid(index~period) +
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(230,245), ylim=c(30,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("Spring SLP Anomalies (Pa) vs. Climate Indices")+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks = c(230,240))+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(strip.background=element_rect(colour="black",
                                    fill="white"),panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA))
SLP_coast2


world2 <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
#plot code
SLP.plot <-ggplot() + 
  geom_raster(data=na.omit(ddd%>%filter(var=="SLP")), aes(x=longitude,y=latitude,fill = coefficient)) + 
  facet_grid(index~period) + 
  geom_sf(data=world2, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(220,250), ylim=c(20,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SLP")+
  scale_x_continuous(breaks = c(225, 245))+
  scale_y_continuous(breaks = c(25, 35,45))+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"),panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA))
SLP.plot
#### JUMBO Plot #####
pdf("Output/Fig JumboV2.pdf", 11,8) 
ggarrange(ggarrange(map,ggarrange(a.plot,j.plot, nrow = 2, labels = c("B", "C")),
                    ncol=2,labels = c("A", "")), c.plot,labels = c("", "D"),nrow=2,heights=c(2,0.75))
dev.off()

pdf("Output/MapSpace.pdf", 8,5) 
ggarrange(map, z.plot,  ncol=2,labels = c("A", "B"))
dev.off()

pdf("Output/Fig SLP and NPH.pdf", 11,6.5) 
ggarrange(ggarrange(a.plot,j.plot,c.plot, nrow = 3, labels = c("A", "B", "C")),
          SLP_coast2,
                    ncol=2,labels = c("","D"))
dev.off()

pdf("Output/Fig NPH Winter.pdf", 6,6.5) 
ggarrange(a.plotw,j.plotw,c.plotw, nrow = 3, labels = c("A", "B", "C"))
dev.off()
