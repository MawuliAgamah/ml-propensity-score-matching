library("ggplot2") # visualization
library("MatchIt") # Matching
library("cobalt") #  covariate balance
library("gridExtra")
library('ggpubr')
library('haven')
library ('gmodels')
library ('MASS')
library('dplyr')
library('MatchItSE')
library('survey')
library('Matching')
library('rgenoud')
options(scipen=999)
set.seed(1234)
# KEY 
# 1 - nsw treated + cps control  (lalonde's original sample)
# 2 - nsw treated + psid control (lalonde's original sample)
# 3 - nsw treated + cps control  (dehejia & wahba sub-sample)
# 4 - nsw treated + psid control (dehejia & wahba sub-sample)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Load Data and match
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

# LOGIT

# un-adjusted logit data set's taken from python
# When matching the default estimand for the match-it function is the ATT , which we use. 


trimming.funct <- function(dataset){
  
  treatedUnits <- dataset[dataset$treat==1,]
  max_ps <- max(treatedUnits$propensity_score)
  min_ps <- min(treatedUnits$propensity_score)
  criteria_met_booleon <- between(dataset$propensity_score, min_ps, max_ps)
  trimmed_df <- dataset[criteria_met_booleon,]
  dropped_count <- nrow(dataset)- nrow(trimmed_df)
  print (c("dropped:", dropped_count))
  return(trimmed_df)
  
}

logitUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswCps_lalonde_ps_unmatched_LOGIT_FS1.csv')
logitUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswPsid_lalonde_ps_unmatched_LOGIT_FS1.csv')
logitUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswCps_dehWab_ps_unmatched_LOGIT_FS1.csv')
logitUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswPsid_dehWab_ps_unmatched_LOGIT_FS1.csv')

logitUndajusted1 <- trimming.funct(logitUndajusted1)
logitUndajusted2 <- trimming.funct(logitUndajusted2)
logitUndajusted3 <- trimming.funct(logitUndajusted3)
logitUndajusted4 <- trimming.funct(logitUndajusted4)

caliper1 = sd(logitUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(logitUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(logitUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(logitUndajusted4$propensity_score, na.rm = FALSE)*0.25

forumla1 = treat ~ age + education. + black + hispanic + married + nodegree + re75 + propensity_score
forumla2 = treat ~ age + education. + black + hispanic + married + nodegree + re74+ re75 + propensity_score

m_out_logit1 <- matchit(formula = forumla1, data = logitUndajusted1,method = "genetic",distance = logitUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_logit2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "genetic", distance = logitUndajusted2$propensity_score,caliper = caliper2, replace =  TRUE,pop.size = 50)
m_out_logit3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "genetic", distance = logitUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_logit4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "genetic", distance = logitUndajusted4$propensity_score,caliper = caliper4, replace =  TRUE,pop.size = 50)

# CART 
cartUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswCps_lalonde_ps_unmatched_CART_FS1.csv')
cartUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswPsid_lalonde_ps_unmatched_CART_FS1.csv')
cartUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswCps_dehWab_ps_unmatched_CART_FS1.csv')
cartUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswPsid_dehWab_ps_unmatched_CART_FS1.csv')


cartUndajusted1 <- trimming.funct(cartUndajusted1)
cartUndajusted2 <- trimming.funct(cartUndajusted2)
cartUndajusted3 <- trimming.funct(cartUndajusted3)
cartUndajusted4 <- trimming.funct(cartUndajusted4)

caliper1 = sd(cartUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(cartUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(cartUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(cartUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_cart1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "genetic",distance = cartUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_cart2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "genetic", distance = cartUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_cart3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "genetic", distance = cartUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_cart4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "genetic", distance = cartUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

# Forest
forestUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswCps_lalonde_ps_unmatched_FOREST_FS1.csv')
forestUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswPsid_lalonde_ps_unmatched_FOREST_FS1.csv')
forestUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswCps_dehWab_ps_unmatched_FOREST_FS1.csv')
forestUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswPsid_dehWab_ps_unmatched_FOREST_FS1.csv')


forestUndajusted1 <- trimming.funct(forestUndajusted1)
forestUndajusted2 <- trimming.funct(forestUndajusted2)
forestUndajusted3 <- trimming.funct(forestUndajusted3)
forestUndajusted4 <- trimming.funct(forestUndajusted4)

caliper1 = sd(forestUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(forestUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(forestUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(forestUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_forest1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "genetic",distance = forestUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_forest2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "genetic", distance = forestUndajusted2$propensity_score,caliper =caliper2,replace =  TRUE,pop.size = 50)
m_out_forest3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "genetic", distance = forestUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_forest4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "genetic", distance = forestUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

# Boost 

boostUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswCps_lalonde_ps_unmatched_BOOST_FS1.csv')
boostUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswPsid_lalonde_ps_unmatched_BOOST_FS1.csv')
boostUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswCps_dehWab_ps_unmatched_BOOST_FS1.csv')
boostUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswPsid_dehWab_ps_unmatched_BOOST_FS1.csv')

boostUndajusted1 <- trimming.funct(boostUndajusted1)
boostUndajusted2 <- trimming.funct(boostUndajusted2)
boostUndajusted3 <- trimming.funct(boostUndajusted3)
boostUndajusted4 <- trimming.funct(boostUndajusted4)

caliper1 = sd(boostUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(boostUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(boostUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(boostUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_boost1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "genetic", distance = boostUndajusted1$propensity_score,caliper =caliper1,replace =  TRUE,pop.size = 50)
m_out_boost2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "genetic", distance = boostUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_boost3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "genetic", distance = boostUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_boost4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "genetic", distance = boostUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

# ANN
# Load unadjusted dataet 
annUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswCps_lalonde_ps_unmatched_ANN_FS1.csv')
annUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswPsid_lalonde_ps_unmatched_ANN_FS1.csv')
annUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswCps_dehWab_ps_unmatched_ANN_FS1.csv')
annUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswPsid_dehWab_ps_unmatched_ANN_FS1.csv')

annUndajusted1 <- trimming.funct(annUndajusted1)
annUndajusted2 <- trimming.funct(annUndajusted2)
annUndajusted3 <- trimming.funct(annUndajusted3)
annUndajusted4 <- trimming.funct(annUndajusted4)

caliper1 = sd(annUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(annUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(annUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(annUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_ann1 <- matchit(formula = forumla1, data = annUndajusted1, method = "genetic", distance = annUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_ann2 <- matchit(formula = forumla1, data = annUndajusted2, method = "genetic", distance = annUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_ann3 <- matchit(formula = forumla2, data = annUndajusted3, method = "genetic", distance = annUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_ann4 <- matchit(formula = forumla2, data = annUndajusted4, method = "genetic", distance = annUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Matching basic summary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

print(summary(m_out_logit1,un = FALSE)) 
print(summary(m_out_logit2,un = FALSE)) 
print(summary(m_out_logit3,un = FALSE)) 
print(summary(m_out_logit4,un = FALSE)) 

print(summary(m_out_cart1,un = FALSE)) 
print(summary(m_out_cart2,un = FALSE)) 
print(summary(m_out_cart3,un = FALSE)) 
print(summary(m_out_cart4,un = FALSE)) 


print(summary(m_out_forest1,un = FALSE)) 
print(summary(m_out_forest2,un = FALSE)) 
print(summary(m_out_forest3,un = FALSE)) 
print(summary(m_out_forest4,un = FALSE))

print(summary(m_out_boost1,un = FALSE)) 
print(summary(m_out_boost2,un = FALSE)) 
print(summary(m_out_boost3,un = FALSE)) 
print(summary(m_out_boost4,un = FALSE)) 

print(summary(m_out_ann1,un = FALSE)) 
print(summary(m_out_ann2,un = FALSE)) 
print(summary(m_out_ann3,un = FALSE)) 
print(summary(m_out_ann4,un = FALSE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Covariate balance plots 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

# LOGIT - balance plots
balance_plot_logit1 <- love.plot(m_out_logit1,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
        xlab("Lalonde's sample with CPS control's")+
        theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
         panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
        ggtitle("")


balance_plot_logit2 <- love.plot(m_out_logit2,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line = TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")

balance_plot_logit3 <- love.plot(m_out_logit3,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with CPS control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_logit4 <- love.plot(m_out_logit4,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line = TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


logit_lalonde <-ggarrange(balance_plot_logit1,balance_plot_logit2, 
                          balance_plot_logit3,balance_plot_logit4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
logit_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/LOGIT_balance_plots_4grid.png', plot = last_plot(),
       dpi = 300)

# CART 
balance_plot_cart1 <- love.plot(m_out_cart1,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with CPS control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_cart2 <- love.plot(m_out_cart2,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line = TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")

balance_plot_cart3 <- love.plot(m_out_cart3,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_cart4 <- love.plot(m_out_cart4,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line = TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .1),
                                colors = c("#003366","#E31B23"),
                                sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


CART_lalonde <-ggarrange(balance_plot_cart1,balance_plot_cart2, 
                          balance_plot_cart3,balance_plot_cart4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
CART_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/CART_balance_plots_4grid.png',dpi = 300, plot = last_plot())

# Random Forest 
# forest balance plots
balance_plot_forest1 <- love.plot(m_out_forest1,
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line =TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled"),
                                thresholds = c(m = .25),
                                colors = c("#003366","#E31B23"),
                                sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with CPS control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_forest2 <- love.plot(m_out_forest2,
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line = TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled"),
                                thresholds = c(m = .25),
                                colors = c("#003366","#E31B23"),
                                sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")

balance_plot_forest3 <- love.plot(m_out_forest3,
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line =TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled"),
                                thresholds = c(m = .25),
                                colors = c("#003366","#E31B23"),
                                sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_forest4 <- love.plot(m_out_forest4,
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line = TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled"),
                                thresholds = c(m = .25),
                                colors = c("#003366","#E31B23"),
                                sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


forest_lalonde <-ggarrange(balance_plot_forest1,balance_plot_forest2, 
                         balance_plot_forest3,balance_plot_forest4,
                         ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
forest_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/FOREST_balance_plots_4grid.png',dpi = 300, plot = last_plot())

# Boosted trees 

# XGboost balance plots
balance_plot_boost1 <- love.plot(m_out_boost1,
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line =TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled"),
                                  thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with CPS control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_boost2 <- love.plot(m_out_boost2,
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line = TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled"),
                                  thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")

balance_plot_boost3 <- love.plot(m_out_boost3,
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line =TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled"),
                                  thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_boost4 <- love.plot(m_out_boost4,
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line = TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled"),
                                  thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


boost_lalonde <-ggarrange(balance_plot_boost1,balance_plot_boost2, 
                           balance_plot_boost3,balance_plot_boost4,
                           ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
boost_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/BOOST_balance_plots_4grid.png',dpi = 300, plot = last_plot())

# ANN
# ANN balance plots
balance_plot_ann1 <- love.plot(m_out_ann1,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with CPS control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_ann2 <- love.plot(m_out_ann2,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line = TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Lalonde's sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")

balance_plot_ann3 <- love.plot(m_out_ann3,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")


balance_plot_ann4 <- love.plot(m_out_ann4,
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line = TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled"),
                                 thresholds = c(m = .25),
                                 colors = c("#003366","#E31B23"),
                                 sample.names = c("unadjusted", "adjusted"))+
  xlab("Dehejia - Wahba sample with PSID control's")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=1),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


ann_lalonde <-ggarrange(balance_plot_ann1,balance_plot_ann2, 
                          balance_plot_ann3,balance_plot_ann4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
ann_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/ann_balance_plots_4grid.png',dpi = 300, plot = last_plot())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
#EQQ plot's
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
plot(m_out_ann3, type = "qq", interactive = FALSE,
     which.xs = c("propensity_score","education.","re75"))

plot(m_out_ann1, type = "qq", interactive = FALSE,
     which.xs = c("hispanic","married","re75"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Common support plot's
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

cs_logit_plt1 <- bal.plot(m_out_logit1, 
                              var.name = "propensity_score", 
                              which = "both",
                              type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                              mirror = TRUE)+
  xlab("Matched CPS control unit's")+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  ggtitle("")

cs_logit_plt2 <- bal.plot(m_out_logit2, 
                              var.name = "propensity_score", 
                              which = "both",
                              type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                              mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18)) +
  xlab("Matched PSID control unit's")+
  ylab("")+
  ggtitle("")
cs_logit_plt1
cs_cart_plt1 <- bal.plot(m_out_cart1, 
                              var.name = "propensity_score", 
                              which = "both",
                              type = "histogram", 
                         colors = c("#E31B23", "#003366"),
                              mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched CPS control unit's")+
  ggtitle("")

cs_cart_plt2 <- bal.plot(m_out_cart2, 
                              var.name = "propensity_score", 
                              which = "both",
                              type = "histogram", 
                         colors = c("#E31B23", "#003366"),
                              mirror = TRUE)+
       theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18)) +
  xlab("Matched PSID control unit's")+
  ylab("")+
  ggtitle("")
cs_cart_plt2
cs_forest_plt1 <- bal.plot(m_out_forest1, 
                         var.name = "propensity_score", 
                         which = "both",
                         type = "histogram", 
                         colors = c("#E31B23", "#003366"),
                         mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched PSID control unit's")+
  ggtitle("")
cs_forest_plt1
cs_forest_plt2 <- bal.plot(m_out_forest2, 
                           var.name = "propensity_score", 
                           which = "both",
                           type = "histogram", 
                           colors = c("#E31B23", "#003366"),
                           mirror = TRUE)+
        theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched PSID control unit's")+
  ggtitle("")

cs_forest_plt2
cs_boost_plt1 <- bal.plot(m_out_boost1, 
                           var.name = "propensity_score", 
                           which = "both",
                           type = "histogram",
                          colors = c("#E31B23", "#003366"),
                           mirror = TRUE)+
        theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched CPS control unit's")+
  ggtitle("")
cs_boost_plt1
cs_boost_plt2 <- bal.plot(m_out_boost2, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
        theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched psid control unit's")+
  ggtitle("")

cs_boost_plt2
cs_ann_plt1 <- bal.plot(m_out_ann1, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram",
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched CPS control unit's")+
  ggtitle("")
cs_ann_plt1
cs_ann_plt2 <- bal.plot(m_out_ann2, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched psid control unit's")+
  ylab("")+
  ggtitle("")
cs_ann_plt2
ggarrange(cs_logit_plt1,cs_logit_plt2,
          cs_cart_plt1,cs_cart_plt2,
          cs_forest_plt1,cs_forest_plt2,
          cs_boost_plt1,cs_boost_plt2,
          cs_ann_plt1,cs_ann_plt2,
          nrow = 5 , ncol = 2,
          common.legend = TRUE, 
          labels = c("Logit","","CART","","XG-BOOST","","Forest","","ANN"),
          font.label = list(size = 18, family = "Times New Roman"),
          legend="right")

ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/balance_plots_4model_lalonde.png',
       dpi = 300, 
       width = 20, height = 18,
       plot = last_plot())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Dehejia & Wahba sub sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cs_logit_plt3 <- bal.plot(m_out_logit3, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  xlab("Matched CPS control unit's")+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  ggtitle("")

cs_logit_plt4 <- bal.plot(m_out_logit4, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18)) +
  xlab("Matched PSID control unit's")+
  ylab("")+
  ggtitle("")

cs_cart_plt3 <- bal.plot(m_out_cart3, 
                         var.name = "propensity_score", 
                         which = "both",
                         type = "histogram", 
                         colors = c("#E31B23", "#003366"),
                         mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched CPS control unit's")+
  ggtitle("")

cs_cart_plt4 <- bal.plot(m_out_cart4, 
                         var.name = "propensity_score", 
                         which = "both",
                         type = "histogram", 
                         colors = c("#E31B23", "#003366"),
                         mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18)) +
  xlab("Matched PSID control unit's")+
  ylab("")+
  ggtitle("")

cs_forest_plt3 <- bal.plot(m_out_forest4, 
                           var.name = "propensity_score", 
                           which = "both",
                           type = "histogram", 
                           colors = c("#E31B23", "#003366"),
                           mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched PSID control unit's")+
  ggtitle("")

cs_forest_plt4 <- bal.plot(m_out_forest4, 
                           var.name = "propensity_score", 
                           which = "both",
                           type = "histogram", 
                           colors = c("#E31B23", "#003366"),
                           mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched PSID control unit's")+
  ylab("")+
  ggtitle("")


cs_boost_plt3 <- bal.plot(m_out_boost3, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram",
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched CPS control unit's")+
  ggtitle("")

cs_boost_plt4 <- bal.plot(m_out_boost4, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched psid control unit's")+
  ylab("")+
  ggtitle("")

cs_ann_plt3 <- bal.plot(m_out_ann3, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram",
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched CPS control unit's")+
  ylab("")+
  ggtitle("")

cs_ann_plt4 <- bal.plot(m_out_ann4, 
                          var.name = "propensity_score", 
                          which = "both",
                          type = "histogram", 
                          colors = c("#E31B23", "#003366"),
                          mirror = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(family = "Times New Roman",size = 18))+
  xlab("Matched psid control unit's")+
  ggtitle("")


ggarrange(cs_logit_plt2,cs_logit_plt4,
          cs_cart_plt2,cs_cart_plt3,
          cs_forest_plt2,cs_forest_plt3,
          cs_boost_plt2,cs_boost_plt3,
          cs_ann_plt2,cs_ann_plt3,
          nrow = 5 , ncol = 2,
          common.legend = TRUE, 
          labels = c("Logit","","CART","","XG-BOOST","","Forest","","ANN"),
          font.label = list(size = 18, family = "Times New Roman"),
          legend="right")

ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/balance_plots_4model_dehejiaWahba.png',
       dpi = 300, 
       width = 20, height = 18,
       plot = last_plot())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
# Get matching counts 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

# Logit 






#m_out_logit1 
#m_out_logit2 
#m_out_logit3 
#m_out_logit4 

forumla1


# CART 
#m_out_cart1 
#m_out_cart2 
#m_out_cart3 
#m_out_cart4 

# Forest
#m_out_forest1 
#m_out_forest2 
#m_out_forest3 
#m_out_forest4 

# Boost 
#m_out_boost1 
#m_out_boost2 
#m_out_boost3 
#m_out_boost4 

# ANN
#m_out_ann1 
#m_out_ann2
#m_out_ann3
#m_out_ann4 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
# WIP - Guido Imbens (2014) stratification_algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

# Create strata based on propensity score 
library("kableExtra")
library(dplyr)
library(gmodels)


x <- match.data(m_out_logit1)
#x <- match.data(m)
y$distance

x$propensity_score


# calculate log odds 
x$log_odds = log(x$propensity_score/(1-x$propensity_score))
#perform t-test of means 
t_test <- t.test(log_odds~treat,data=x)
t_stat <- t_test$statistic # get test statistic 

n <- 1 # initial number of strata 

x$strata <- ifelse(x$log_odds > median(x$log_odds), 1,2)
split1 <- x[x$strata==1,]

t_stat_s1 <- t.test(log_odds~treat,data=split1)
abs(t_stat_s1$statistic)

x$strata <- ifelse(x$log_odds > median(x$log_odds), 3,4)

split1 <- x[x$strata==3,]
t_stat_s1 <- t.test(log_odds~treat,data=x[x$strata==3,])
x$strata
x[x$strata==2,]

if (abs(t_stat) > n){ # check if block is inequality balanced , inadequate balance if t>n
  
  # split block and run t tests again
  x$strata <- ifelse(x$log_odds > median(x$log_odds), 1,2)
  
  split1 <- x[x$strata==1,]
  t.test(log_odds~treat,data=split1)
  x[x$strata==2,]
}else{
  print('no')
}


x$strata <- ifelse(x$log_odds > median(x$log_odds), 1,2)
x$strata
data = data %>% mutate(quantile = ntile(propensity_score, strata))


imbens_rubin_stratification_algorithm <- function(matchit_object,strata){
  
  # Begin with 1 strata 
  
    data <- match.data(matchit_object) # calculate log odds 
  
}

# t-statistic 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
#  Sub classification
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

# Recursive stratification to ensure there is always at least 1 treated unit in each strata.
# This could be improved upon to better test for common support (mean diff or statistically)
stratification_function <- function(matchit_object,strata){
  
  if(is.data.frame(matchit_object)){
    dataout1 <- matchit_object %>% mutate(quantile = ntile(propensity_score, strata)) # stratify
    x <- CrossTable(dataout1$treat, dataout1$quantile) # summary of strata
    if (0 %in% x$t){return(stratification_function(dataout1,strata-1))}
    else{return(dataout1)}
  }else{
    datain <- match.data(matchit_object) # match it object to data frame 
    dataout2 <- datain %>% mutate(quantile = ntile(propensity_score, strata)) # stratify
    x2 <- CrossTable(dataout2$treat, dataout2$quantile) # summary of strata
    if (0 %in% x2$t){return(stratification_function(dataout2,strata-1))}# check common support
    else{return(dataout2)}
  }}





# Run function over mathched data set's
# logit

stratifiedMatch_logit1<- stratification_function(m_out_logit1,5)
stratifiedMatch_logit2 <- stratification_function(m_out_logit2,5)
stratifiedMatch_logit3 <- stratification_function(m_out_logit3,5)
stratifiedMatch_logit4 <- stratification_function(m_out_logit4,5)
# cart
stratifiedMatch_cart1<- stratification_function(m_out_cart1,5)
stratifiedMatch_cart2 <- stratification_function(m_out_cart2,5)
stratifiedMatch_cart3 <- stratification_function(m_out_cart3,5)
stratifiedMatch_cart4 <- stratification_function(m_out_cart4,5)
# forest
stratifiedMatch_forest1<- stratification_function(m_out_forest1,5)
stratifiedMatch_forest2 <- stratification_function(m_out_forest2,5)
stratifiedMatch_forest3 <- stratification_function(m_out_forest3,5)
stratifiedMatch_forest4 <- stratification_function(m_out_forest4,5)
# boost
stratifiedMatch_boost1<- stratification_function(m_out_boost1,5)
stratifiedMatch_boost2 <- stratification_function(m_out_boost2,5)
stratifiedMatch_boost3 <- stratification_function(m_out_boost3,5)
stratifiedMatch_boost4 <- stratification_function(m_out_boost4,5)
# ann
stratifiedMatch_ann1<- stratification_function(m_out_ann1,5)
stratifiedMatch_ann2 <- stratification_function(m_out_ann2,5)
stratifiedMatch_ann3 <- stratification_function(m_out_ann3,5)
stratifiedMatch_ann4 <- stratification_function(m_out_ann4,5)


# Function to estimate ATT pooled across strata 
simple_att_pooled_estimator <-function(stratified_data){
  estimates <- list()
  for(i in unique(stratified_data$quantile)){
    block <- stratified_data[stratified_data$quantile==i,]
    Y_t <- block$re78[block$treat==1]
    Y_c <- block$re78[block$treat==0]
    n_t <- table(Y_t) %>% sum()
    n_c <- table(Y_c) %>% sum()
    tor <- (sum(Y_t)-sum(Y_c))/n_t
    estimates[[i]]<-tor
  } 
  return(Reduce("+",estimates))
  
}

simple_att_pooled_estimator(stratifiedMatch_logit2)
stratifiedMatch_logit1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
# Within strata balance tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
MatchBalance(forumla1, data=stratifiedMatch_logit1[stratifiedMatch_logit1$quantile == 3,],nboots=500)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
# Treatment effect estimation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----


# Benchmark experimental data 
benchmark.data <- read_dta('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/nsw.dta')
benchmark.data2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/nsw_dehejia_wahba.csv')
benchmark.experimental.data <- subset(benchmark.data,select = -c(data_id))
unmatched.data.psid.dehwab <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/unmatched data/Quasi_NswPsid_dehWab.csv')
benchmark.experimental.data$agesq = benchmark.experimental.data$age*benchmark.experimental.data$age # age squared 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # # Linear regressions - no controls
specification_1 = re78 ~ treat


# Experimental benchmark 
summary(lm(specification_1 ,benchmark.experimental.data))
summary(lm(specification_1 ,benchmark.data2))

# matching logit benchmark 

# Logit
summary(lm(re78 ~ treat ,match.data(m_out_logit1)))
summary(lm(re78 ~ treat ,match.data(m_out_logit2)))
summary(lm(re78 ~ treat ,match.data(m_out_logit3)))
summary(lm(re78 ~ treat ,match.data(m_out_logit4)))

# Cart
summary(lm(re78 ~ treat ,match.data(m_out_cart1)))
summary(lm(re78 ~ treat ,match.data(m_out_cart2)))
summary(lm(re78 ~ treat ,match.data(m_out_cart3)))
summary(lm(re78 ~ treat ,match.data(m_out_cart4)))
# Forest
summary(lm(re78 ~ treat ,match.data(m_out_forest1)))
summary(lm(re78 ~ treat ,match.data(m_out_forest2)))
summary(lm(re78 ~ treat ,match.data(m_out_forest3)))
summary(lm(re78 ~ treat ,match.data(m_out_forest4)))
# Boost
summary(lm(re78 ~ treat ,match.data(m_out_boost1)))
summary(lm(re78 ~ treat ,match.data(m_out_boost2)))
summary(lm(re78 ~ treat ,match.data(m_out_boost3)))
summary(lm(re78 ~ treat ,match.data(m_out_boost4)))

# ANN
summary(lm(re78 ~ treat ,match.data(m_out_ann1)))
summary(lm(re78 ~ treat ,match.data(m_out_ann2)))
summary(lm(re78 ~ treat ,match.data(m_out_ann3)))
summary(lm(re78 ~ treat ,match.data(m_out_ann4)))

# Linear regression with controls 
specification_2 = re78 ~ treat + age + agesq + nodegree+black+hispanic + re75
specification_3 = re78 ~ treat + age + agesq + nodegree+black+hispanic + re74 + re74
# Experimental benchmark 
summary(lm(specification_2 ,benchmark.experimental.data))
summary(lm(specification_3 ,benchmark.data2))

regression_controls <- function(data,ref_specification){
  data = match.data(data)
  data$agesq = data$age*data$age
  return(summary(lm(ref_specification , data)))
  #return(lm(ref_specification , data))
  
}


# Logit
regression_controls(m_out_logit1,specification_2)
regression_controls(m_out_logit2,specification_2)
regression_controls(m_out_logit3,specification_3)
regression_controls(m_out_logit4,specification_3)

# Cart
regression_controls(m_out_cart1,specification_2)
regression_controls(m_out_cart2,specification_2)
regression_controls(m_out_cart3,specification_3)
regression_controls(m_out_cart4,specification_3)
# Forest
regression_controls(m_out_forest1,specification_2)
regression_controls(m_out_forest2,specification_2)
regression_controls(m_out_forest3,specification_3)
regression_controls(m_out_forest4,specification_3)
# Boost
regression_controls(m_out_boost1,specification_2)
regression_controls(m_out_boost2,specification_2)
regression_controls(m_out_boost3,specification_3)
regression_controls(m_out_boost4,specification_3)

# ANN
regression_controls(m_out_ann1,specification_2)
regression_controls(m_out_ann2,specification_2)
regression_controls(m_out_ann3,specification_3)
regression_controls(m_out_ann4,specification_3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Addide and Imebns methods 


abadie_imbens_estimator <- function(object,dataset){
  
  x <- att(obj = object, Y = dataset$re78)
  y <- bootstrap.se(object, Y = dataset$re78, max.iter = 1000)
  y <- y[[2]]
  return(cat("ate:",x,"se:",y))
  
}

abadie_imbens_estimator(m_out_logit1,logitUndajusted1)
abadie_imbens_estimator(m_out_logit2,logitUndajusted2)
abadie_imbens_estimator(m_out_logit3,logitUndajusted3)
abadie_imbens_estimator(m_out_logit4,logitUndajusted4)

abadie_imbens_estimator(m_out_cart1,cartUndajusted1)
abadie_imbens_estimator(m_out_cart2,cartUndajusted2)
abadie_imbens_estimator(m_out_cart3,cartUndajusted3)
abadie_imbens_estimator(m_out_cart4,cartUndajusted4)

abadie_imbens_estimator(m_out_forest1,forestUndajusted1)
abadie_imbens_estimator(m_out_forest2,forestUndajusted2)
abadie_imbens_estimator(m_out_forest3,forestUndajusted3)
abadie_imbens_estimator(m_out_forest4,forestUndajusted4)

abadie_imbens_estimator(m_out_boost1,boostUndajusted1)
abadie_imbens_estimator(m_out_boost2,boostUndajusted2)
abadie_imbens_estimator(m_out_boost3,boostUndajusted3)
abadie_imbens_estimator(m_out_boost4,boostUndajusted4)


abadie_imbens_estimator(m_out_ann1, annUndajusted1)
abadie_imbens_estimator(m_out_ann2,annUndajusted2)
abadie_imbens_estimator(m_out_ann3,annUndajusted3)
abadie_imbens_estimator(m_out_ann4,annUndajusted4)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

# Horowitz-Thompson Weighted ATT estimator


#  STRATIEFIED Weighted ATT

#stratification_function <- function(matchit_object,strata){

#    if(is.data.frame(matchit_object)){
#      dataout1 <- datain %>% mutate(quantile = ntile(propensity_score, strata)) # stratify
#      x <- CrossTable(dataout1$treat, dataout1$quantile) # summary of strata
#        if (0 %in% x$t){return(stratification_function(dataout1,strata-1))}
#          else{return(dataout1)}
#    }else{
#      datain <- match.data(matchit_object) # match it object to data frame 
#      dataout2 <- datain %>% mutate(quantile = ntile(propensity_score, strata)) # stratify
#      x2 <- CrossTable(dataout2$treat, dataout2$quantile) # summary of strata
#      if (0 %in% x$t){return(stratification_function(dataout2,strata-1))}# check common support
#         else{return(dataout2)}
#}}

# weighted glm estimator 
weighted_regression_estimator <- function(matchit_object){
  options(survey.lonely.psu = 'adjust')
  stratified_data <- stratification_function(matchit_object,5)
  stratified_data$ID <- seq.int(nrow(stratified_data)) # create id column
  surveyDesign1 <- svydesign(ids =~ID,strata=~quantile,data = stratified_data,nest=F) # create survey design
  # replicate weight's for bootstrapped standard errors
  surveyDesign1.bootstrap <- as.svrepdesign(surveyDesign1,type=c("bootstrap"),replicates=100)

  outcomeModel2006Boot <- svyglm(re78~treat,surveyDesign1.bootstrap)
  summary(outcomeModel2006Boot)
  
}
  
# Logit
weighted_regression_estimator(m_out_logit1)
weighted_regression_estimator(m_out_logit2)
weighted_regression_estimator(m_out_logit3)
weighted_regression_estimator(m_out_logit4)

# Cart
weighted_regression_estimator(m_out_cart1)
weighted_regression_estimator(m_out_cart2)
weighted_regression_estimator(m_out_cart3)
weighted_regression_estimator(m_out_cart4)
# Forest
weighted_regression_estimator(m_out_forest1)
weighted_regression_estimator(m_out_forest2)
weighted_regression_estimator(m_out_forest3)
weighted_regression_estimator(m_out_forest4)
# Boost
weighted_regression_estimator(m_out_boost1)
weighted_regression_estimator(m_out_boost2)
weighted_regression_estimator(m_out_boost3)
weighted_regression_estimator(m_out_boost4)

# ANN
weighted_regression_estimator(m_out_ann1)
weighted_regression_estimator(m_out_ann2)
weighted_regression_estimator(m_out_ann3)
weighted_regression_estimator(m_out_ann4)

x <- match.data(m_out_forest1)
treated <- x[x$treat==1,]
control <- x[x$treat==0,]

mean(treated$re78) - mean(control$re78)
summary(control)
summary(treated)


by(x, x$treat, summary)


