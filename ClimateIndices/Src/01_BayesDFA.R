library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(ggplot2)
library(stringr)
library(ggpubr)
#Organize SCC biology data
dat<- read.csv("data/biologydata_south.central_2023.csv")
ax<- 20
ti<-24
wid <- 28
plot_trends2 <- function(rotated_modelfit,
                         years = NULL,
                         highlight_outliers = FALSE,
                         threshold = 0.01) {
  rotated <- rotated_modelfit
  df <- dfa_trends(rotated, years = years)
  
  # make faceted ribbon plot of trends
  p1 <- ggplot(df, aes_string(x = "time", y = "estimate")) +
    geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), alpha = 0.4) +
    geom_line() +
    facet_wrap("trend_number") +
    xlab("Time") +
    ylab("")+
    
    theme(axis.text=element_text(size=ax),
          axis.title=element_text(size=ti,face="bold"))+
    theme_bw()
  
  if (highlight_outliers) {
    swans <- find_swans(rotated, threshold = threshold)
    df$outliers <- swans$below_threshold
    p1 <- p1 + geom_point(data = df[which(df$outliers), ], color = "red")
  }
  
  p1
}


plot_loadings2 <- function(rotated_modelfit,
                           names = NULL,
                           facet = TRUE,
                           violin = TRUE,
                           conf_level = 0.95,
                           threshold = NULL) {
  v <- dfa_loadings(rotated_modelfit,
                    summary = FALSE,
                    names = names,
                    conf_level = conf_level
  )
  df <- dfa_loadings(rotated_modelfit,
                     summary = TRUE,
                     names = names,
                     conf_level = conf_level
  )
  
  # filter values below threshold
  if (!is.null(threshold)) {
    df <- df[df$prob_diff0 >= threshold, ]
    v <- v[v$prob_diff0 >= threshold, ]
  }
  
  if (!violin) {
    p1 <- ggplot(df, aes_string(
      x = "name", y = "median", col = "trend",
      alpha = "prob_diff0"
    )) +
      geom_point(size = 3, position = position_dodge(0.3)) +
      geom_errorbar(aes_string(ymin = "lower", ymax = "upper"),
                    position = position_dodge(0.3), width = 0
      ) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() +
      xlab("Time Series") +
      ylab("Loading")+
      guides(fill="none", alpha='none')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = wid))+
      theme(legend.position="none")+
      theme_bw()
  }
  
  if (violin) {
    p1 <- ggplot(v, aes_string(
      x = "name", y = "loading", fill = "trend",
      alpha = "prob_diff0"
    )) +
      geom_violin(color = NA) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() +
      xlab("Time Series") +
      ylab("Loading")+
      theme(axis.text=element_text(size=ax),
            axis.title=element_text(size=ti,face="bold"))+
      guides(fill="none", alpha='none')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = wid))+
      theme_bw()
  }
  
  if (facet) {
    p1 <- p1 + facet_wrap(~trend, scales = "free_x")
  }
  
  p1
}



l.calcofi <- plot_loadings2(r.calcofi,names=namescalcofi$namesshort)
l.calcofi
n1 <- names(dat)[grepl('calcofi.',names(dat))]
n2 <- names(dat)[grepl('rreas.',names(dat))]
ids <- c(n1,n2,"ZALOPHUS.PUPCT","ZALOPHUS.PUPWT")
n3 <- names(dat)[grepl('SBRD.',names(dat))]
n4 <- names(dat)[grepl('ZALOPHUS.',names(dat))]

for(i in 1:ncol(dat)){
  if(names(dat)[i] %in% ids){
    dat[,i] <- log(dat[,i])
  }
}

#dat<-dat%>%select(c(year,n4,n3,n1))


#### CALCOFI####
y1calcofi <- 1951
y2calcofi<- 2022
dat.calcofi<-dat%>%select(c(year,n1)) #just setting an order so we know how to set the variance index for each survey
remelt = melt(dat.calcofi,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(FALSE),
                       var_index = c("survey"), num_trends = 1,
                       elpd_loo = NA, se_elpd_loo=NA)

varIndx = c(rep(1,length(n1)))
fit.mod.calcofi = fit_dfa(y = Y,
                  num_trends = 1,
                  iter=n_iter,
                  varIndx = varIndx,
                  chains=n_chains, estimate_nu=model_df$est_nu[1],
                  estimate_trend_ma = model_df$estimate_trend_ma[1],
                  estimate_trend_ar = model_df$estimate_trend_ar[1],
                  estimate_process_sigma = model_df$estimate_process_sigma[1],
                  seed=123)

namescalcofi<-data.frame(names)%>%
  mutate(namesshort= c("Snubnose smelt", "Slender blacksmelt", "Dogtooth lampfish", 
                       "Northern anchovy", "Cal. smoothtongue", "Eared blacksmelt", 
                       "North Pacific hake", "Cal. flashlightfish", "Sardinops",
                       "Northern lampfish", "Bigfin laternfish", "Blue laternfish",
                       "Mexican lampfish", "Lightfishes"))
pars = rstan::extract(fit.mod.calcofi$model)
r.calcofi <- rotate_trends(fit.mod.calcofi)
p.calcofi <- plot_trends2(r.calcofi,years =dat.calcofi$year)
p.calcofi
l.calcofi <- plot_loadings2(r.calcofi,names=namescalcofi$namesshort)
l.calcofi
is_converged(fit.mod.calcofi)
summary(fit.mod.calcofi)

#### RREAS ####
y1rreas <- 1983
y2rreas<- 2022
dat.RREAS<-dat%>%select(c(year,n2))%>%
  filter(year>=y1rreas&year<=y2rreas)#just setting an order so we know how to set the variance index for each survey
remelt = melt(dat.RREAS,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear


varIndx = c(rep(1,length(n2)))
fit.mod.RREAS = fit_dfa(y = Y,
                          num_trends = 1,
                          iter=n_iter,
                          varIndx = varIndx,
                          chains=n_chains, estimate_nu=model_df$est_nu[1],
                          estimate_trend_ma = model_df$estimate_trend_ma[1],
                          estimate_trend_ar = model_df$estimate_trend_ar[1],
                          estimate_process_sigma = model_df$estimate_process_sigma[1],
                          seed=123)

namesRREAs<-data.frame(names)%>%
  mutate(namesshort= c("Adult anchovy", "Adult sardine", "Krill", "Market squid",
                       "Myctophids", "Juv. anchovy", "Juv. boccacio rock.",
                       "Juv. chilipepper rock.", "Juv. Pacific hake", 
                       "Juv. halfbanded rock.", "Juv Pacific sand.",
                       "Juv. shortbelly rock.","Juv. speckled sand.", "Juv. widow rock."))
pars = rstan::extract(fit.mod.RREAS$model)
r.RREAS <- rotate_trends(fit.mod.RREAS)
p.RREAS <- plot_trends2(r.RREAS,years =dat.RREAS$year)
p.RREAS
l.RREAS <- plot_loadings2(r.RREAS,names=namesRREAs$namesshort)
l.RREAS
is_converged(fit.mod.RREAS)
summary(fit.mod.RREAS)


#### Seabird####
y1sea <- 1971
y2sea<- 2022

dat.SEA<-dat%>%select(c(year,n3))%>%
  filter(year>=y1sea&year<=y2sea) #just setting an order so we know how to set the variance index for each survey
remelt = melt(dat.SEA,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear


varIndx = c(rep(1,length(n3)))
fit.mod.SEA = fit_dfa(y = Y,
                        num_trends = 1,
                        iter=n_iter,
                        varIndx = varIndx,
                        chains=n_chains, estimate_nu=model_df$est_nu[1],
                        estimate_trend_ma = model_df$estimate_trend_ma[1],
                        estimate_trend_ar = model_df$estimate_trend_ar[1],
                        estimate_process_sigma = model_df$estimate_process_sigma[1],
                        seed=123)
namessea<-data.frame(names)%>%
  mutate(namesshort= c("Ashy-storm petrel", "Brandt's cormorant", "Cassin's auklet",
                       "Common murre", "Pelagic Cormorant", "Pigeon guillemot",
                       "Rhinoceros auklet","Western Gull"))

pars = rstan::extract(fit.mod.SEA$model)
r.SEA <- rotate_trends(fit.mod.SEA)
p.SEA <- plot_trends2(r.SEA, years =dat.SEA$year)
p.SEA
l.SEA <- plot_loadings2(r.SEA,names=namessea$namesshort)
l.SEA
is_converged(fit.mod.SEA, threshold = 1.05, parameters = c("sigma", "x", "Z"))
summary(fit.mod.SEA)

trend.sea <-dfa_trends(r.SEA, years =dat.SEA$year)%>%
  select(time, estimate, lower, upper)%>%
 # rename(sea.estimate = estimate, sea.lower=lower,sea.upper=upper)%>%
  mutate(trend = "SEA")%>%
  mutate(era=ifelse(time<=1988, 1, ifelse(time>1988&time<=2012, 2,3)))
trend.rreas <-dfa_trends(r.RREAS, years =dat.RREAS$year)%>%
  select(time, estimate, lower, upper)%>%
#  rename(rreas.estimate = estimate, rreas.lower=lower,rreas.upper=upper)%>%
  mutate(trend = "RREAS")%>%
  mutate(era=ifelse(time<=1988, 1, ifelse(time>1988&time<=2012, 2,3)))
trend.calcofi <-dfa_trends(r.calcofi, years =dat.calcofi$year)%>%
  select(time, estimate, lower, upper)%>%
 # rename(calcofi.estimate = estimate, calcofi.lower=lower,calcofi.upper=upper)%>%
  mutate(trend = "CALCOFI")%>%
  mutate(era=ifelse(time<=1988, 1, ifelse(time>1988&time<=2012, 2,3)))

dfa.trends <- trend.calcofi%>%
  add_row(trend.rreas)
saveRDS(dfa.trends, file = "dfa.trends.rds")

arranged <- ggarrange(p.calcofi, l.calcofi,p.RREAS,l.RREAS,ncol = 2, nrow = 2,labels = c("A", "B", "C", "D"),
                      heights=c(1,1.25,1.5))

pdf(file = "DFAtrendsloadings.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7)
arranged
dev.off()
