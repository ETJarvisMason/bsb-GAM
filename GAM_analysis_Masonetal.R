# Jarvis Mason, ET
# 26 FEB 2024


# load libraries ----------------------------------------------------------

library(tidyverse)
library(funtimes)
library(mgcv)
# source("http://www.sthda.com/upload/rquery_cormat.r")
library(corrplot)
library(patchwork)

# load theme pub function -------------------------------------------------

theme_Publication <- function(base_size=14, base_family="Helvetica") {

  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),#element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",# bottom
            legend.direction = "vertical",
            legend.key.size= unit(0.75, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_blank(),#element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


# load data ---------------------------------------------------------------

# 1) Barred Sand Bass CPFV harvest data (Bightwide, reconstructed 1947-1974)
# 2) King Harbor, CA SCUBA density survey data for adults and young-of-the-year juvenile recruits 1974-2022
# 3) Pt. Dume, CA shore station sea surface temperature data 1957-2022 (annual mean after taking summer monthly means, Jun-Aug)
# 4) Oceanic Nino Index 1957-2022 (annual mean of monthly values)

load(file="BSB_dat.Rdata")

# cross-correlation to identify any lags among time series (use package funtimes t --------
# uses bootstrap approach to account for potential autocorrelation
# NOTE:  NAs not allowed!

set.seed(022224)

# ---------are sst and adult densities correlated?

d <- dat %>% 
  dplyr::select(sst,ad) %>% 
  drop_na()

# 3 - 4 year lag
ccf_01_sst.ad <- ccf_boot(as.numeric(d$sst),
                      as.numeric(d$ad),lag.max=10)

# ----------are sst and harvest correlated?

d <- dat %>% 
  dplyr::select(sst,lands) %>% 
  drop_na()

# 4 - 10 year lag 
ccf_02_sst.lands <- ccf_boot(as.numeric(d$sst),
                      as.numeric(d$lands),lag.max=10)


# ----------are adult density data correlated with harvest data?
d <- dat %>% 
  dplyr::select(lands,ad) %>% 
  drop_na()

ccf_03_ad.lands <- ccf_boot(as.numeric(d$ad),
                            as.numeric(d$lands),lag.max=10)

ccf_03_ad.lands # adult densities lead harvest data from 0-5 y


# ----------are sst and YOY juvenile recruitment correlated?

d <- dat %>% 
  dplyr::select(sst,juv) %>% 
  drop_na()

# 0 lag is highest, but not outside blue region
ccf_04_sst.juv <- ccf_boot(as.numeric(d$sst),
                           as.numeric(d$juv),lag.max=3)

# ---------are oni and YOY juvenile recruitment correlated?

d <- dat %>% 
  dplyr::select(oni,juv) %>% 
  drop_na()

# 1 and 2 y lags highest, but not outside blue region
ccf_05_oni.juv <- ccf_boot(as.numeric(d$oni),
                      as.numeric(d$juv),lag.max=3)

# ---------are landings and YOY juvenile recruitment correlated?

d <- dat %>% 
  dplyr::select(lands,juv) %>% 
  drop_na()

# 0 and 1 y lags highest
ccf_06_lands.juv <- ccf_boot(as.numeric(d$lands),
                           as.numeric(d$juv),lag.max=10)

# ---------are adult densities and YOY juvenile recruitment correlated?

d <- dat %>% 
  dplyr::select(ad,juv) %>% 
  drop_na()

# 2-4 y lags highest??
ccf_07_ad.juv <- ccf_boot(as.numeric(d$ad),
                             as.numeric(d$juv),lag.max=10)



# fancy plots
ccf.plot01 <- ggplot(ccf_01_sst.ad,aes(Lag,r_P))+
  geom_ribbon(aes(x=Lag,ymin=lower_P,ymax=upper_P,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_P), size = 1.5))+
  scale_color_gradient(low = "yellow", high = "darkred") +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_P)) +
  geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype="dotted")+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,1))+
  scale_y_continuous(limits = c(-0.6,0.8))+
  xlab("Lag")+
  ylab("Cross-correlation\n coefficient")+
  labs(subtitle = "SST (t + Lag)\nand Adult Density (t)")+
  theme_Publication()+
  theme(legend.position = "")
ccf.plot01

ccf.plot02 <- ggplot(ccf_02_sst.lands,aes(Lag,r_P))+
  geom_ribbon(aes(x=Lag,ymin=lower_P,ymax=upper_P,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_P), size = 1.5))+
  scale_color_gradient(low = "yellow", high = "darkred") +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_P)) +
  geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype="dotted")+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,1))+
  scale_y_continuous(limits = c(-0.6,0.8))+
  xlab("Lag")+
  ylab("Cross-correlation\n coefficient")+
  labs(subtitle = "SST (t + Lag)\nand Harvest (t)")+
  theme_Publication()+
  theme(legend.position = "")
ccf.plot02

ccf.plot03 <- ggplot(ccf_03_ad.lands,aes(Lag,r_P))+
  geom_ribbon(aes(x=Lag,ymin=lower_P,ymax=upper_P,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_P), size = 1.5))+
  scale_color_gradient(low = "yellow", high = "darkred") +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_P)) +
  geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype="dotted")+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,1))+
  scale_y_continuous(limits = c(-0.6,0.8))+
  xlab("Lag")+
  ylab("Cross-correlation\n coefficient")+
  labs(subtitle = "Adult Density (t + Lag)\nand Harvest (t)")+
  theme_Publication()+
  theme(legend.position = "")
ccf.plot03

ccf.plot04 <- ggplot(ccf_04_sst.juv,aes(Lag,r_P))+
  geom_ribbon(aes(x=Lag,ymin=lower_P,ymax=upper_P,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_P), size = 1.5))+
  scale_color_gradient(low = "yellow", high = "darkred") +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_P)) +
  geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype="dotted")+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,1))+
  scale_y_continuous(limits = c(-0.4,0.6),breaks = seq(-0.4,0.6,0.2))+
  xlab("Lag")+
  ylab("Cross-correlation\n coefficient")+
  ggtitle("SST (t + Lag) and Recruits (t)")+
  # ggtitle("ONI (t + Lag) and Recruits (t)")+
  theme_Publication()+
  theme(legend.position = "")
ccf.plot04

ccf.plot05 <- ggplot(ccf_05_oni.juv,aes(Lag,r_P))+
  geom_ribbon(aes(x=Lag,ymin=lower_P,ymax=upper_P,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_P), size = 1.5))+
  scale_color_gradient(low = "yellow", high = "darkred") +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_P)) +
  geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype="dotted")+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,1))+
  scale_y_continuous(limits = c(-0.4,0.6),breaks = seq(-0.4,0.6,0.2))+
  xlab("Lag")+
  ylab("Cross-correlation\n coefficient")+
  # ggtitle("SST (t + Lag) and Recruits (t)")+
  ggtitle("ONI (t + Lag) and Recruits (t)")+
  theme_Publication()+
  theme(legend.position = "")
ccf.plot05

# Generalized Additive ModelS --------------------------------------------

gam.datwlags <- dat %>% 
  mutate(lands.l = lag(lands),
         ad.l = lag(ad))  # lag by one year
  
# Note there are missing data for some variables in some years; however, 
# we are leaving as is because the the default for running a gam   
# is simply to use only the ‘complete cases’.

# check out correlation matrix
justvars <- gam.datwlags %>%
  dplyr::select(-Year) %>% 
  drop_na()

C <- cor(justvars)
corrplot(C, method = "number")


# model based on top lags from cross correlation above (sst lag of zero, oni lag of one)
gam.m1 <- gam(juv ~ 
                   s(sst) +
                   s(oni.l),
                 data = gam.datwlags,
                 family = tw(),
                 method = "REML",
                 select = TRUE)
m <- gam.m1
gam.check(m)
summary(m)
m$aic

# check for concurvity
concurvity(m)
round(concurvity(m, full = FALSE)$estimate, 2) # pairwise

# Plot the partial responses
plot(m, pages = 1, scale = 0)
predsm1 <- visreg::visreg(m,
                          line = list(col = "red"),
                          points = list(cex = 1, pch = 1),
                          scale = "response",
                          type = "conditional",
                          gg=TRUE
)

predsm1
library(patchwork)
e<-predsm1[[1]] + theme_Publication()
f<-predsm1[[2]] + theme_Publication()


e+f + plot_layout(ncol=2)


# include harvest to account for changes in recruit densities due to fishing (no lags) -- **best model
gam.m2 <- gam(juv ~ 
                s(sst) +
                s(oni.l) +
                s(lands),
              data = gam.datwlags,
              family = tw(),
              method = "REML", 
              select = TRUE)
m <- gam.m2
gam.check(m)
summary(m)  
m$aic

# check for concurvity
concurvity(m)
round(concurvity(m, full = FALSE)$estimate, 2) # pairwise


# Plot the partial responses
plot(m, pages = 1, scale = 0)
predsm2 <- visreg::visreg(m,
                          line = list(col = "red"),
                          points = list(cex = 1, pch = 1),
                          scale = "response",
                          type = "conditional",
                          gg=TRUE
)

library(patchwork)
a<-predsm2[[1]] + ylab("Density") + theme_Publication()
b<-predsm2[[2]] + ylab("Density") + theme_Publication()
c<-predsm2[[3]] + ylab("Density") + theme_Publication()

fig7_top <- a+b+c + plot_layout(ncol=2)
fig7_top
  

gam.m2b <- gam(juv ~ 
                s(sst) +
                s(oni.l) +
                s(lands.l),
              data = gam.datwlags,
              family = tw(),
              method = "REML", 
              select = TRUE)
m <- gam.m2b
gam.check(m)
summary(m)  
m$aic

# check for concurvity
concurvity(m)
round(concurvity(m, full = FALSE)$estimate, 2) # pairwise


# add in adult densities
gam.m3 <- gam(juv ~ 
                s(sst) +
                s(ad) +
                s(lands) +
                s(oni.l),
              data = gam.datwlags,
              family = tw(),
              method = "REML",
              select = TRUE)
m <- gam.m3
gam.check(m)
summary(m) # adult densities do not add anything
m$aic
# check for concurvity
concurvity(m)
round(concurvity(m, full = FALSE)$estimate, 2) # concurvity high for s(ad) 



# Now lag adult densities by one year
gam.m4 <- gam(juv ~ 
                s(sst) +
                s(ad.l) +
                s(lands) +
                s(oni.l),
              data = gam.datwlags,
              family = tw(),
              method = "REML",     
              select = TRUE)
m <- gam.m4
gam.check(m) # edf not close to k'
summary(m) # no big improvement
m$aic

# check for concurvity
concurvity(m)
round(concurvity(m, full = FALSE)$estimate, 2) # concurvity moderate for s(ad), but not found to be important anyway

# Plot the partial responses

plot.gam(m, residuals = FALSE, pch =1, cex = 1, shade = 
           TRUE, shade.col = "lightblue", seWithMean = TRUE, 
         pages = 1, all.terms = TRUE, scheme = 2)
par(mfrow=c(1,4))



# no oni or ads
gam.m9 <- gam(juv ~ 
                s(sst) +
                s(lands),
              data = gam.datwlags,
              family = tw(),
              method = "REML",     
              select = TRUE)
m <- gam.m9
gam.check(m)
summary(m) 
m$aic
# check for concurvity
concurvity(m)
round(concurvity(m, full = FALSE)$estimate, 2)

plot(m)
## Predicted recruit density All years---------------------------------------------------

# reduced model (excludes adult densities)
m <- gam.m2
ggdf2 <- ggpredict(m, terms = c("sst","oni.l [-0.85,0,0.85]"))
plot(ggdf2)
nina <- expression(paste("La Ni","u+00F1","a"))
cols <- c("-0.85" ="skyblue", "0" = "lightgray", "0.85" = "red")
oni.labs <- c("-0.85" ="La Nina","0" ="Neutral","0.85" ="El Nino")
pred.plot.red.mod <- ggplot(ggdf2,aes(x, predicted,colour = group,color = group,fill=group))+
  geom_line(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low,ymax = conf.high,
                  fill=group), colour = NA,alpha = 0.3)+
  facet_wrap(~ group, labeller = labeller(group=oni.labs))+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  scale_y_continuous("Density",limits=c(0, 1.0),breaks = seq(0,1.0,0.5))+
  guides(fill=guide_legend(title="ONI (t-1)"),
         color=FALSE)+
  theme_Publication() +
  xlab("Sea surface temperature (\u00b0C)")+
  ggtitle("Predicted recruit density")
pred.plot.red.mod

# temperature only plot (includes lag on oni)
m <- gam.m1
ggdf1 <- ggpredict(m, terms = c("sst","oni.l [-0.85,0,0.85]"))
plot(ggdf1)

cols <- c("-0.85" ="skyblue", "0" = "lightgray", "0.85" = "red")
oni.labs <- c("-0.85" ="La Nina","0" ="Neutral","0.85" ="El Nino")
pred.plot.tempOnly.mod <- ggplot(ggdf1,aes(x, predicted,colour = group,color = group,fill=group))+
  geom_line(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low,ymax = conf.high,
                  fill=group), colour = NA,alpha = 0.3)+
  facet_wrap(~ group, labeller = labeller(group=oni.labs))+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  guides(fill=guide_legend(title="ONI (t-1)"),
         color=FALSE)+
  theme_Publication() +
  xlab("Sea surface temperature (\u00b0C)")+
  scale_y_continuous("Density",limits=c(0, 1.5),breaks = seq(0,1.5,0.5))+
  ggtitle("Predicted recruit density")
pred.plot.tempOnly.mod


