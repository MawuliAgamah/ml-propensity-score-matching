
library("ggplot2") # visualization
library("MatchIt") # Matching
library("cobalt") #  covariate balance
library("gridExtra")
library('ggpubr')
library('haven')

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

library('Matching')
library('rgenoud')

logitUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/unmatched/nswCps_lalonde_ps_unmatched_LOGIT.csv')
logitUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/unmatched/nswPsid_lalonde_ps_unmatched_LOGIT.csv')
logitUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/unmatched/nswCps_dehWab_ps_unmatched_LOGIT.csv')
logitUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/unmatched/nswPsid_dehWab_ps_unmatched_LOGIT.csv')

logitUndajusted1$comparison_group <- "cps"
logitUndajusted2$comparison_group <- "psid"
logitUndajusted3$comparison_group <- "cps"
logitUndajusted4$comparison_group <- "psid"

logitUndajusted1$sample <- "lalonde"
logitUndajusted2$sample <- "lalonde"
logitUndajusted3$sample <- "dehejia_wahba"
logitUndajusted4$sample <- "dehejia_wahba"

forumla1 = treat ~ age + education. + black + hispanic + married + nodegree + re75 + re78 + propensity_score
forumla2 = treat ~ age + education. + black + hispanic + married + nodegree + re74+ re75 + re78 + propensity_score 

caliper1 = sd(logitUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(logitUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(logitUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(logitUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_logit1 <- matchit(formula = forumla1, data = logitUndajusted1, method = "genetic",distance = logitUndajusted1$propensity_score,caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_logit2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "genetic", distance = logitUndajusted2$propensity_score,caliper = 0.25, replace =  TRUE,pop.size = 50)
m_out_logit3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "genetic", distance = logitUndajusted3$propensity_score,caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_logit4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "genetic", distance = logitUndajusted4$propensity_score,caliper = 0.25, replace =  TRUE,pop.size = 50)

# CART 
cartUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/unmatched/nswCps_lalonde_ps_unmatched_CART.csv')
cartUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/unmatched/nswPsid_lalonde_ps_unmatched_CART.csv')
cartUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/unmatched/nswCps_dehWab_ps_unmatched_CART.csv')
cartUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/unmatched/nswPsid_dehWab_ps_unmatched_CART.csv')

cartUndajusted1$comparison_group <- "cps"
cartUndajusted2$comparison_group <- "psid"
cartUndajusted3$comparison_group <- "cps"
cartUndajusted4$comparison_group <- "psid"

cartUndajusted1$sample <- "lalonde"
cartUndajusted2$sample <- "lalonde"
cartUndajusted3$sample <- "dehejia_wahba"
cartUndajusted4$sample <- "dehejia_wahba"

m_out_cart1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "genetic", distance = cartUndajusted1$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_cart2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "genetic", distance = cartUndajusted2$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_cart3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "genetic", distance = cartUndajusted3$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_cart4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "genetic", distance = cartUndajusted4$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)

# Forest

forestUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/unmatched/nswCps_lalonde_ps_unmatched_FOREST.csv')
forestUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/unmatched/nswPsid_lalonde_ps_unmatched_FOREST.csv')
forestUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/unmatched/nswCps_dehWab_ps_unmatched_FOREST.csv')
forestUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/unmatched/nswPsid_dehWab_ps_unmatched_FOREST.csv')

forestUndajusted1$comparison_group <- "cps"
forestUndajusted2$comparison_group <- "psid"
forestUndajusted3$comparison_group <- "cps"
forestUndajusted4$comparison_group <- "psid"

forestUndajusted1$sample <- "lalonde"
forestUndajusted2$sample <- "lalonde"
forestUndajusted3$sample <- "dehejia_wahba"
forestUndajusted4$sample <- "dehejia_wahba"


caliper1 = sd(forestUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(forestUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(forestUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(forestUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_forest1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "genetic", distance = forestUndajusted1$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_forest2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "genetic", distance = forestUndajusted2$propensity_score,caliper ,caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_forest3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "genetic", distance = forestUndajusted3$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_forest4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "genetic", distance = forestUndajusted4$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)

# Boost 

boostUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/unmatched/nswCps_lalonde_ps_unmatched_BOOST.csv')
boostUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/unmatched/nswPsid_lalonde_ps_unmatched_BOOST.csv')
boostUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/unmatched/nswCps_dehWab_ps_unmatched_BOOST.csv')
boostUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/unmatched/nswPsid_dehWab_ps_unmatched_BOOST.csv')

boostUndajusted1$comparison_group <- "cps"
boostUndajusted2$comparison_group <- "psid"
boostUndajusted3$comparison_group <- "cps"
boostUndajusted4$comparison_group <- "psid"

boostUndajusted1$sample <- "lalonde"
boostUndajusted2$sample <- "lalonde"
boostUndajusted3$sample <- "dehejia_wahba"
boostUndajusted4$sample <- "dehejia_wahba"

caliper1 = sd(boostUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(boostUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(boostUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(boostUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_boost1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "genetic", distance = boostUndajusted1$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_boost2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "genetic", distance = boostUndajusted2$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_boost3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "genetic", distance = boostUndajusted3$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_boost4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "genetic", distance = boostUndajusted4$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)


# ANN
# Load unadjusted dataset 
annUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/ann/unmatched/nswCps_lalonde_ps_unmatched_ANN.csv')
annUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/ann/unmatched/nswPsid_lalonde_ps_unmatched_ANN.csv')
annUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/ann/unmatched/nswCps_dehWab_ps_unmatched_ANN.csv')
annUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/ann/unmatched/nswPsid_dehWab_ps_unmatched_ANN.csv')

annUndajusted1$comparison_group <- "cps"
annUndajusted2$comparison_group <- "psid"
annUndajusted3$comparison_group <- "cps"
annUndajusted4$comparison_group <- "psid"

annUndajusted1$sample <- "lalonde"
annUndajusted2$sample <- "lalonde"
annUndajusted3$sample <- "dehejia_wahba"
annUndajusted4$sample <- "dehejia_wahba"

caliper1 = sd(annUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(annUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(annUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(annUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_ann1 <- matchit(formula = forumla1, data = annUndajusted1, method = "genetic", distance = annUndajusted1$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_ann2 <- matchit(formula = forumla1, data = annUndajusted2, method = "genetic", distance = annUndajusted2$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_ann3 <- matchit(formula = forumla2, data = annUndajusted3, method = "genetic", distance = annUndajusted3$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)
m_out_ann4 <- matchit(formula = forumla2, data = annUndajusted4, method = "genetic", distance = annUndajusted4$propensity_score,caliper , caliper = 0.25,replace =  TRUE,pop.size = 50)

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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                abs = TRUE,
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
                                abs = TRUE,
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
                                abs = TRUE,
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
                                abs = TRUE,
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
                                  abs = TRUE,
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
                                  abs = TRUE,
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
                                  abs = TRUE,
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
                                  abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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
                                 abs = TRUE,
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

balance_plot_ann4

ann_lalonde <-ggarrange(balance_plot_ann1,balance_plot_ann2, 
                          balance_plot_ann3,balance_plot_ann4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
ann_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/ann_balance_plots_4grid.png',dpi = 300, plot = last_plot())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
#EQQ plot's
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
plot(m_out_ann2, type = "qq", interactive = FALSE,
     which.xs = c("age","education.","re75"))



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
# WIP - Guido Imben's (2014) stratification_algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

# Create strata based on propensity score 
library("kableExtra")
library(dplyr)
library(gmodels)


# Recursive stratification to ensure there no untreated unit's in each strata
# This could be improved upon to also statistically test for common support

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
# Sub classification
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

stratification_function <- function(matchit_object,strata){
  
  # Begin with 5 strata
  if(strata == 5){
      data <- match.data(matchit_object)
      data = data %>% mutate(quantile = ntile(propensity_score, strata))
      x <- CrossTable(data$treat, data$quantile)
          if (0 %in% x$t){ # check common support
            return(stratification_function(data,strata-1)) # if common support not met in a strata, run function again
          }else{
            return(data) # return stratified data 
          }
    }else{
      data = matchit_object %>% mutate(quantile = ntile(propensity_score, strata))
      x <- CrossTable(data$treat, data$quantile)
          if (0 %in% x$t){ # check common support
              return(stratification_function(data,strata-1)) # if common support not met in a strata, run function again
          }else{
            return(data) # return stratified data 
            }
    }
}

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


#logit

CrossTable(stratifiedMatch_logit1$treat, stratifiedMatch_logit1$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_logit2$treat, stratifiedMatch_logit2$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_logit3$treat, stratifiedMatch_logit3$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_logit4$treat, stratifiedMatch_logit4$quantile)  # Treated and control counts
# Cart
CrossTable(stratifiedMatch_cart1$treat,stratifiedMatch_cart1$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_cart2$treat,stratifiedMatch_cart2$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_cart3$treat,stratifiedMatch_cart3$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_cart4$treat,stratifiedMatch_cart4$quantile)  # Treated and control counts
# FOREST
CrossTable(stratifiedMatch_forest1$treat,stratifiedMatch_forest1$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_forest2$treat,stratifiedMatch_forest2$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_forest3$treat,stratifiedMatch_forest3$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_forest4$treat,stratifiedMatch_forest4$quantile)  # Treated and control counts
# BOOST
CrossTable(stratifiedMatch_boost1$treat,stratifiedMatch_boost1$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_boost2$treat,stratifiedMatch_boost2$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_boost3$treat,stratifiedMatch_boost3$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_boost4$treat,stratifiedMatch_boost4$quantile)  # Treated and control counts
# ANN
CrossTable(stratifiedMatch_ann1$treat,stratifiedMatch_ann1$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_ann2$treat,stratifiedMatch_ann2$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_ann3$treat,stratifiedMatch_ann3$quantile)  # Treated and control counts
CrossTable(stratifiedMatch_ann4$treat,stratifiedMatch_ann4$quantile)  # Treated and control counts

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
# Within strata balance tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
MatchBalance(forumla1, data=stratifiedMatch_logit1[stratifiedMatch_logit1$quantile == 3,],nboots=500)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
# Treatment effect estimation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----

# Benchmark experimental data 
benchmark.data <- read_dta('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/nsw.dta')
benchmark.experimental.data <- subset(benchmark.data,select = -c(data_id))
unmatched.data.psid.dehwab <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/unmatched data/Quasi_NswPsid_dehWab.csv')

# Benchmark for matching with match It package , nearest neighbor matching and logit propensity scores 
nn_logit_benchmark <- matchit(treat~
                              age + agesq + 
                              education. + educsq +  black + 
                              hispanic + married + nodegree +
                              re75 + re75 + u74+ u75, data =  unmatched.data.cps.lalonde,
                              caliper = 0.05, method = "nearest",distance = "logit")


# Function to estimate ATE pooled across strata 
`simple_ate_pooled_estimator <-function(stratified_data){
  estimates <- list()
  for(i in unique(stratified_data$quantile)){
    block <- stratified_data[stratified_data$quantile==i,]
    Y_t <- block$re78[block$treat==1]
    Y_c <- block$re78[block$treat==0]
    n_t <- table(Y_t) %>% sum()
    n_c <- table(Y_c) %>% sum()
    tor <- sum(Y_t/n_t)-sum(Y_c/n_c)
    estimates[[i]]<-tor
    } 
  return(Reduce("+",estimates))

}

simple_ate_pooled_estimator(stratifiedMatch_forest1)


summary(lm(re78 ~ treat + age + agesq + education. + nodegree + black + hispanic+re74+re75 ,match.data(m)))

benchmark.experimental.data$agesq = benchmark.experimental.data$age*benchmark.experimental.data$age
# benchmark 
summary(lm(re78 ~ treat + 
             age + agesq +
             education + nodegree + 
             black + hispanic+re75 ,
            benchmark.experimental.data))

stratifiedMatch_cart2$agesq = stratifiedMatch_cart2$age*stratifiedMatch_cart2$age

summary(lm(re78 ~ treat ,logitstrata1))





