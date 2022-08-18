
library("ggplot2") # visualization
library("MatchIt") # Matching
library("cobalt") #  covariate balance
library("gridExtra")
library('ggpubr')


# KEY 
# 1 - nsw treated + cps control  (lalonde's original sample)
# 2 - nsw treated + psid control (lalonde's original sample)
# 3 - nsw treated + cps control  (dehejia & wahba sub-sample)
# 4 - nsw treated + psid control (dehejia & wahba sub-sample)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Load Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

# LOGIT

# un-adjusted logit data set's taken from python

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

m_out_logit1 <- matchit(formula = forumla1, data = logitUndajusted1, method = "nearest", distance = logitUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_logit2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "nearest", distance = logitUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_logit3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "nearest", distance = logitUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_logit4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "nearest", distance = logitUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)

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

caliper1 = sd(cartUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(cartUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(cartUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(cartUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_cart1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "nearest", distance = cartUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_cart2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "nearest", distance = cartUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_cart3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "nearest", distance = cartUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_cart4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "nearest", distance = cartUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)

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

m_out_forest1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "nearest", distance = forestUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_forest2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "nearest", distance = forestUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_forest3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "nearest", distance = forestUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_forest4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "nearest", distance = forestUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)

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

m_out_boost1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "nearest", distance = boostUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_boost2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "nearest", distance = boostUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_boost3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "nearest", distance = boostUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_boost4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "nearest", distance = boostUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)


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

m_out_ann1 <- matchit(formula = forumla1, data = annUndajusted1, method = "nearest", distance = annUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_ann2 <- matchit(formula = forumla1, data = annUndajusted2, method = "nearest", distance = annUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_ann3 <- matchit(formula = forumla2, data = annUndajusted3, method = "nearest", distance = annUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_ann4 <- matchit(formula = forumla2, data = annUndajusted4, method = "nearest", distance = annUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)



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
        xlab("Lalonde's subsample with CPS controls")+
        theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
         panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Lalonde's sample with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Dehejia & Wahba's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Dehejia & Wahba's subsample with PSID control")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


logit_lalonde <-ggarrange(balance_plot_logit1,balance_plot_logit2, 
                          balance_plot_logit3,balance_plot_logit4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
logit_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/LOGIT_balance_plots_4grid.png', plot = last_plot(),dpi = 300)

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
  xlab("Lalonde's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Lalonde's sample with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Dehejia & Wahba's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Dehejia & Wahba's subsample with PSID control")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")


CART_lalonde <-ggarrange(balance_plot_cart1,balance_plot_cart2, 
                          balance_plot_cart3,balance_plot_cart4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")

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
  xlab("Lalonde's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Lalonde's sample with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Dehejia & Wahba's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Dehejia & Wahba's subsample with PSID control")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Lalonde's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Lalonde's sample with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Dehejia & Wahba's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Dehejia & Wahba's subsample with PSID control")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Lalonde's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Lalonde's sample with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
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
  xlab("Dehejia & Wahba's subsample with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="right",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5))+
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
  xlab("Dehejia & Wahba's subsample with PSID control")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman"),
        legend.position="none",
        panel.border=element_rect(size=2),
        panel.grid.minor.y = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        panel.grid.minor.x= element_line(colour="grey", size=0.5),
        panel.grid.major.x= element_line(colour="grey", size=0.5),
        axis.ticks.length=unit(-0.2, "cm"))+
  ggtitle("")

balance_plot_ann4

ann_lalonde <-ggarrange(balance_plot_ann1,balance_plot_ann2, 
                          balance_plot_ann3,balance_plot_ann4,
                          ncol = 2,nrow = 2, common.legend = TRUE, legend="right")
ann_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/ann_balance_plots_4grid.png',dpi = 300, plot = last_plot())


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

cs_ANN_plt3 <- bal.plot(m_out_ANN3, 
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

cs_ANN_plt4 <- bal.plot(m_out_ANN4, 
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

ggarrange(cs_logit_plt3,cs_logit_plt4,
          cs_cart_plt3,cs_cart_plt4,
          cs_forest_plt3,cs_forest_plt4,
          cs_boost_plt3,cs_boost_plt4,
          cs_ANN_plt3,cs_ANN_plt4,
          nrow = 4 , ncol = 2,
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

m_out_logit1 <- matchit(formula = forumla1, data = logitUndajusted1, method = "nearest", distance = logitUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_logit2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "nearest", distance = logitUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_logit3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "nearest", distance = logitUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_logit4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "nearest", distance = logitUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)

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

caliper1 = sd(cartUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(cartUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(cartUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(cartUndajusted4$propensity_score, na.rm = FALSE)*0.25

m_out_cart1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "nearest", distance = cartUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_cart2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "nearest", distance = cartUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_cart3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "nearest", distance = cartUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_cart4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "nearest", distance = cartUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)

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

m_out_forest1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "nearest", distance = forestUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_forest2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "nearest", distance = forestUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_forest3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "nearest", distance = forestUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_forest4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "nearest", distance = forestUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)

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

m_out_boost1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "nearest", distance = boostUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_boost2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "nearest", distance = boostUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_boost3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "nearest", distance = boostUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_boost4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "nearest", distance = boostUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)


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

m_out_ann1 <- matchit(formula = forumla1, data = annUndajusted1, method = "nearest", distance = annUndajusted1$propensity_score,caliper = caliper1,replace = FALSE)
m_out_ann2 <- matchit(formula = forumla1, data = annUndajusted2, method = "nearest", distance = annUndajusted2$propensity_score,caliper = caliper2, replace = FALSE)
m_out_ann3 <- matchit(formula = forumla2, data = annUndajusted3, method = "nearest", distance = annUndajusted3$propensity_score,caliper = caliper3,replace = FALSE)
m_out_ann4 <- matchit(formula = forumla2, data = annUndajusted4, method = "nearest", distance = annUndajusted4$propensity_score,caliper = caliper4, replace = FALSE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#---
# Treatment Effect Estimation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Divide propensity scores into five strata using quantiles 

quantiles1 <- quantile(LogitMatched1$propensity_score,prob = seq(0,1,1/5))
LogitMatched1$strata <- cut(LogitMatched1$propensity_score, breaks = quantiles1,include.lowest = TRUE)

levels(LogitMatched1$strata )

levels(LogitMatched1$strata) <- 1:length(levels(LogitMatched1$strata )) # rename strata labels 

# Examine common support (cs) - no strata should have 0 treated or 0 control units 
xtabs(~LogitMatched1$treat+LogitMatched1$strata)  # cs met 


# Check covariate balance with MatchIt

forumla1 = treat ~ age + education. + black + hispanic + married + nodegree + re75 + re78 + propensity_score
forumla2 = treat ~ age + education. + black + hispanic + married + nodegree + re74+ re75 + re78 + propensity_score 

covariates <- c("treat","age","education.","black","hispanic","married","nodegree","re75","re78")
covariates_dehWab <- c("treat","age","education.","black","hispanic","married","nodegree","re74","re75","re78")

LogitMatched1[c("propensity_score","strata","treat")]

balance_formula1 <- paste(covariates , collapse  =  "+")
balance_formula2 <- paste(covariates_dehWab , collapse = "+")
stratification <- matchit()













