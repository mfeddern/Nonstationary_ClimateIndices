library(ggrepel)
library(strucchange)
library(ncdf4)
library(tidyverse)
library(ggpp)
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
library(sf)
library(terra)
library(raster)

# reading in the data 
schroeder.nph <-read.csv('data/physical/year_mon_area_max_x_y_lon_lat_2023.csv')%>%
  mutate(era= ifelse(Year<=1988,1, ifelse(Year>2012,3,2)))%>% #assigning eras
  mutate(Year_win = if_else(Month == 11|Month ==12, Year+1, Year))%>% #creating a year offset fro winter
  mutate(dec.yr =as.numeric(as.character(Year)) + (as.numeric(Month)-0.5)/12, #creating a decimal year
         Area=Area/1000000)%>%  #converting area to smaller units
  mutate(area.anom = Area-mean(Area)) #creating an area anomaly

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

spring.schroeder.Z<- spring.schroeder%>%
  mutate(Z.mean.area=scale(mean.area), 
         Z.mean.y=scale(mean.y), 
         Z.mean.max=scale(mean.max), 
         Z.mean.x=scale(mean.x))
z<-spring.schroeder.Z%>%filter(abs(Z.mean.area)>3|
                                 abs(Z.mean.y)>3|
                                 abs(Z.mean.x)>3|
                                 abs(Z.mean.max)>3)

Zplot<-spring.schroeder.Z%>%rename(Extent=Z.mean.area,
       Intensity = Z.mean.max,
       Latitude=Z.mean.y,
       Longitude = Z.mean.x)%>%
  dplyr::select(Year, Extent, Intensity, Latitude, Longitude)%>%
  pivot_longer(-Year,names_to = 'Index', values_to = 'Z_Score')
nino<-Zplot%>%filter(Year==1997|Year==1998|Year==2018|Year==2019)
ggplot(data=Zplot,aes(x=Year,y=Z_Score))+
  facet_wrap(~Index)+
  geom_line()+
  geom_point()+
  geom_point(data=nino,aes(x=Year,y=Z_Score),col="red")+
  geom_line(data=Zplot%>%filter(Year==1997|Year==1998),aes(x=Year,y=Z_Score),col="red")+
  geom_line(data=Zplot%>%filter(Year==2018|Year==2019),aes(x=Year,y=Z_Score),col="red")+
  #geom_text(col='grey')+
  ggtitle("North Pacific High") +
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7))

#spring.schroeder.Z%>%filter(abs(Z.mean.y)>2)

means <- spring.schroeder%>%filter(Year!=1992&Year!=1997&Year!=1998)%>%
  group_by(era.lab)%>%
  summarise(x=mean(mean.x),y=mean(mean.y),sd.x=sd(mean.x), sd.y=sd(mean.y),
            area=mean(mean.max),intensity=mean(mean.area),
            sd.area=sd(mean.max), sd.intensity=sd(mean.area),
            Year=0)%>%
  rename(mean.x=x, mean.y=y,mean.area=intensity, mean.max=area )


#spring.schroeder<-winter.schroeder
#plotting location
col<-pnw_palette("Sunset2",3,type="discrete")
a.plot <-ggplot(data=spring.schroeder%>%filter(Year!=1992&Year!=1997),aes(abs(mean.x-360),mean.y, label=Year,group=era.lab,col=era.lab))+
  geom_point(alpha=0.7,aes(shape=era.lab))+
  #geom_text(col='grey')+
  #ggtitle("Center of North Pacific High") +
  geom_point(data=means)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means,aes(ymin = mean.y-sd.y, ymax=  mean.y+sd.y), width=0.5) +
  geom_errorbar(data=means,aes(xmin = abs(mean.x-360)-sd.x, xmax=  abs(mean.x-360)+sd.x), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  geom_hline(yintercept=mean(spring.schroeder$mean.y),lty=2, col='grey')+
  geom_vline(xintercept=abs(mean(spring.schroeder$mean.x-360)), lty=2, col='grey')+
  scale_x_reverse(lim=c(147,139))+
  ylab('Latitude (ºN)')+
  theme_bw() +
  xlab('Longitude (ºW)')+
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.position = "none")
a.plot

j.plot <-ggplot(data=spring.schroeder%>%
                  filter(Year!=1992&Year!=1998),aes(shape=era.lab,y=mean.max,x=mean.area, label=Year,group=era.lab,col=era.lab))+
  geom_point(alpha=0.7)+
 # ggtitle("North Pacific High\n Areal Extent and Intensity") +
  geom_point(data=means)+
  geom_errorbar(data=means,aes(xmin = mean.area-sd.area, xmax= mean.area+sd.area), width=0.5) +
  geom_errorbar(data=means,aes(ymin = mean.max-sd.intensity, ymax= mean.max+sd.intensity), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  #scale_x_reverse(lim=c(147,135))+
  ylab('North Pacific High \n Intensity (hPa)')+
  theme_bw()+
  geom_hline(yintercept=mean(spring.schroeder$mean.max),lty=2, col='grey')+
  geom_vline(xintercept=mean(spring.schroeder$mean.area), lty=2, col='grey')+
  xlab(expression("North Pacific High Area "~(10^6 ~km^2)))+
  #geom_text_s(data=spring.schroeder%>%filter(Year>2020),colour='black',
  #            #colour.target = "segment",
   #          arrow = arrow(length = grid::unit(1.5, "mm")),
   #          point.padding = 0.4,
   #          # angle = 45,
    #          nudge_y = c(0.5,0.1,-0.2),
    #         nudge_x = c(0.3,0.5,0.5),
   #           show.legend = FALSE) +
  labs(col  = "", shape = "")+
   theme( panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.7), 
        legend.key.size = unit(3.5, "mm"),legend.text = element_text(size = 10))
j.plot


#making a TS for break point analysis Y
y.ts <- ts(data=spring.schroeder%>%dplyr::select(mean.y), 1967, 2023, frequency=1)
# fit breakpoint model
bp.y <- breakpoints(y.ts ~ 1)
summary(bp.y) 

#making a TS for break point analysis X
x.ts <- ts(data=spring.schroeder%>%dplyr::select(mean.x), 1967, 2023, frequency=1)
# fit breakpoint model
bp.x <- breakpoints(x.ts ~ 1)
summary(bp.x) 

#making a TS for break point analysis Area
area.ts <- ts(data=spring.schroeder%>%dplyr::select(mean.area.anom), 1967, 2023, frequency=1)
# fit breakpoint model
bp.area <- breakpoints(area.ts ~ 1)
summary(bp.area) 

#making a TS for break point analysis Intensity
max.ts <- ts(data=spring.schroeder%>%dplyr::select(mean.max), 1967, 2023, frequency=1)
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


# and calculate standard deviation over 11-year rolling windows for intensity
NPH.max.sd <- rollapply(spring.schroeder%>%
                          #filter(Year!=1998)%>%
                          dplyr::select(mean.max), 11, sd, fill=NA,na.rm=T)
plot(1967:2023, NPH.max.sd , type="l") #check plot
NPH.max.sd2 <- rollapply(spring.schroeder%>%dplyr::filter(Year!=1998)%>%dplyr::select(mean.max), 11, sd, fill=NA,na.rm=T)

# and calculate standard deviation over 11-year rolling windows for area
NPH.area.sd <- rollapply(spring.schroeder%>%dplyr::select(mean.area.anom), 11, sd, fill=NA)
plot(1967:2023, NPH.area.sd, type="l") #check plot
# now fit a non-parametric regression

# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.max.sd))%>%mutate(Year=year)
plot.dat3 <- data.frame(year=c(1972:1997,1999:2018), sd=na.omit(NPH.max.sd2))%>%
  mutate(Year=year)%>%
  merge(plot.dat%>%dplyr::select(year),all.y=TRUE)
plot.dat3[27,2] = (1.635+1.4681698)/2 
plot.dat3[27,3] = 1998 
# fit the model

mod <- gam(mean.max ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit 

merge(plot.dat3,plot.dat)

mod <- gam(mean.max ~ s(year), data=plot.dat3)
pred <- predict(mod, se=T, newdata = plot.dat3)
plot.dat3$mean <- pred$fit

plot.dat3<-plot.dat3%>%rename(mean.max.out=mean.max, mean.out=mean)
plot.dat2 <-spring.schroeder%>%
  dplyr::select(Year, mean.x, mean.y,mean.max, mean.area,era,era.lab)%>%
  rename(mean.max.full=mean.max)%>%
  left_join(plot.dat%>%dplyr::select(Year, mean, mean.max))%>%
  left_join(plot.dat3%>%dplyr::select(Year, mean.out, mean.max.out))

plot.dat2[is.na(plot.dat2)] <- 1.01
plot.dat2 <-plot.dat2%>%left_join(plot.dat%>%dplyr::select(Year,year))
max_first  <- 2.4#max(plot.dat2$mean.max)   # Specify max of first y axis
max_second <- max(plot.dat2$mean.max.full) # Specify max of second y axis
min_first  <- 0.5#min(plot.dat2$mean.max)   # Specify min of first y axis
min_second <- min(plot.dat2$mean.max.full) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x, scale, shift){
  return ((x)*scale - shift)
}

# Function to scale secondary variable values
inv_scale_function <- function(x, scale, shift){
  return ((x + shift)/scale)
}

summary(lm(mean.max.full~Year, data = plot.dat2))
i.plot<-ggplot(data=plot.dat2,aes(x=Year,y=mean.max))+
  geom_vline(xintercept=2012.5, lty=2, size=0.3)+
  geom_vline(xintercept=1988.5, lty=2, size=0.3)+
 # geom_smooth(aes(x=Year, y=inv_scale_function(mean.max.full, scale, shift)), 
 #             method='lm',col='gray', alpha=.5, se=F) +
  geom_line(size=0.2,aes(x=year))+
  geom_point(aes(x=Year, y=inv_scale_function(mean.max.full, scale, shift)), col='gray', alpha=.5) +
  geom_point(data=plot.dat2%>%filter(Year==1992|Year==1997|Year==1998),aes(x=Year, y=inv_scale_function(mean.max.full, scale, shift)), col='red', alpha=.5) +
  geom_line(aes(x=year, y=mean), color=cb[2], size=0.6) +
  geom_line(aes(x=year, y=mean.out), color=cb[6], size=0.6) +
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  scale_y_continuous(limits = c(min_first, max_first), sec.axis = sec_axis(~scale_function(., scale, shift), name='North Pacific High \n Intensity (hPa)')) +
  ylab('Standard deviation \n (hPa)')+
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
i.plot

##### Jumbo Plot ####

#### Map Data #####
col2<-pnw_palette("Sunset2",4,type="discrete")
col<-pnw_palette("Sunset2",3,type="discrete")
col3<-pnw_palette("Sunset2",8,type="continuous")

col<-pnw_palette("Sunset2",3,type="discrete")
climate_dat <-readRDS('data/physical/climate_dat_upwelling.rds')
climate_dat_cop <-readRDS('data/physical/climate_dat_cop.rds')
bakunsites <- read.csv('data/physical/Bakun/MapLocations3.csv')%>%
  mutate(longitude=longitude)
sites <- st_as_sf(data.frame(bakunsites[,1:2]), coords = c("longitude","latitude"), crs = 4326, 
                  agr = "constant")

#### Making the Map #####
world <- st_as_sf(map('world', plot=F, fill=T)) #base layer for land masses

map_og<-ggplot() +
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
  annotate("rect", xmin= -121.5, xmax = -109, ymin = 43, ymax = 48.8, 
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
  annotate(geom = "text", x = c(-114,-114,-114,-114), y = c(48,46.5,45, 44), 
           label = str_wrap(c("Upwelling Data","Newport Line","CalCOFI", "RREAS"), width = 22),
           color = "grey22", size =3.5) +
  
  annotate(geom = "text", x = c(-120.5,-117), y = c(41,35), 
           label = str_wrap(c("Cape Mendocino","Point Conception"), width = 20),
           fontface = "italic", color = "grey22", size = 3) +
  annotate("line", x= c(-124.1, -124.65), y = c(44.652, 44.652),col=col2[1],size = 0.8, lwd=1) +
  annotate("line", x= c(-120.5, -119.5), y = c(46.5, 46.5),col=col2[1],size = 0.8, lwd=1) +
  
  theme(panel.background = element_rect(fill = "lightsteelblue2"),
        panel.border = element_rect(fill = NA),panel.grid.major = element_line(colour = "transparent"))

map 

#### JUMBO Plot #####


mapPlot<-ggarrange(map, ncol=1,labels = c("A"))
NPHPlot<-ggarrange(a.plot,j.plot,i.plot, nrow = 3, labels = c("B", "C", "D"))

pdf("Output/Fig Map and NPH.pdf", 9.75,8) 
ggarrange(mapPlot,NPHPlot,ncol=2,widths=c(1,1.25))
dev.off()
