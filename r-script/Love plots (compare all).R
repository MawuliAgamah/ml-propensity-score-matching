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
# 1 - nsw treated + CPS control  ( Lalonde's original sample)
# 2 - nsw treated + PSID control ( Lalonde's original sample)
# 3 - nsw treated + CPS control  ( Dehejia & Wahba sub-sample)
# 4 - nsw treated + PSID control ( Dehejia & Wahba sub-sample)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Load Data and match
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# LOGIT

# un-adjusted logit data set's taken from python
# When matching the default estimated for the match-it function is the ATT , which we use. 

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

forumla1 = treat ~ age + education. + black + hispanic + married + nodegree + re75 + propensity_score
forumla2 = treat ~ age + education. + black + hispanic + married + nodegree + re74+ re75 + propensity_score

# logit 

# Get comparison selected data

logitUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswCps_lalonde_ps_unmatched_LOGIT.csv')
logitUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswPsid_lalonde_ps_unmatched_LOGIT.csv')
logitUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswCps_dehWab_ps_unmatched_LOGIT.csv')
logitUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswPsid_dehWab_ps_unmatched_LOGIT.csv')

logitUndajusted1 <- trimming.funct(logitUndajusted1)
logitUndajusted2 <- trimming.funct(logitUndajusted2)
logitUndajusted3 <- trimming.funct(logitUndajusted3)
logitUndajusted4 <- trimming.funct(logitUndajusted4)

caliper1 = sd(logitUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(logitUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(logitUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(logitUndajusted4$propensity_score, na.rm = FALSE)*0.25

# nearest neighbor 

m_out_logit_nn1 <- matchit(formula = forumla1, data = logitUndajusted1, method = "nearest",ratio=1,distance = logitUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE)
m_out_logit_nn2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "nearest",ratio=1, distance = logitUndajusted2$propensity_score,caliper = caliper2, replace =  TRUE)
m_out_logit_nn3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "nearest",ratio=1, distance = logitUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_logit_nn4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "nearest",ratio=1, distance = logitUndajusted4$propensity_score,caliper = caliper4, replace =  TRUE)

# genetic matching 

m_out_logit_gen1 <- matchit(formula = forumla1, data = logitUndajusted1,method = "genetic",distance = logitUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_logit_gen2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "genetic", distance = logitUndajusted2$propensity_score,caliper = caliper2, replace =  TRUE,pop.size = 50)
m_out_logit_gen3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "genetic", distance = logitUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_logit_gen4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "genetic", distance = logitUndajusted4$propensity_score,caliper = caliper4, replace =  TRUE,pop.size = 50)# 

# Get non-featured selected data
logitUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswCps_lalonde_ps_unmatched_LOGIT_FS1.csv')
logitUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswPsid_lalonde_ps_unmatched_LOGIT_FS1.csv')
logitUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswCps_dehWab_ps_unmatched_LOGIT_FS1.csv')
logitUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/logit/unmatched/nswPsid_dehWab_ps_unmatched_LOGIT_FS1.csv')

logitUndajusted1 <- trimming.funct(logitUndajusted1)
logitUndajusted2 <- trimming.funct(logitUndajusted2)
logitUndajusted3 <- trimming.funct(logitUndajusted3)
logitUndajusted4 <- trimming.funct(logitUndajusted4)

# nearest neighbor(feature selected)

m_out_logit_fs_nn1 <- matchit(formula = forumla1, data = logitUndajusted1, method = "nearest",ratio=1,distance = logitUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE)
m_out_logit_fs_nn2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "nearest",ratio=1, distance = logitUndajusted2$propensity_score,caliper = caliper2, replace =  TRUE)
m_out_logit_fs_nn3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "nearest",ratio=1, distance = logitUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_logit_fs_nn4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "nearest",ratio=1, distance = logitUndajusted4$propensity_score,caliper = caliper4, replace =  TRUE)

# genetic matching (feature selected)
m_out_logit_fs_gen1 <- matchit(formula = forumla1, data = logitUndajusted1,method = "genetic",distance = logitUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_logit_fs_gen2 <- matchit(formula = forumla1, data = logitUndajusted2, method = "genetic", distance = logitUndajusted2$propensity_score,caliper = caliper2, replace =  TRUE,pop.size = 50)
m_out_logit_fs_gen3 <- matchit(formula = forumla2, data = logitUndajusted3, method = "genetic", distance = logitUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_logit_fs_gen4 <- matchit(formula = forumla2, data = logitUndajusted4, method = "genetic", distance = logitUndajusted4$propensity_score,caliper = caliper4, replace =  TRUE,pop.size = 50)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# CART

# Get comparison selected data

cartUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswCps_lalonde_ps_unmatched_CART.csv')
cartUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswPsid_lalonde_ps_unmatched_CART.csv')
cartUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswCps_dehWab_ps_unmatched_CART.csv')
cartUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/cart/unmatched/nswPsid_dehWab_ps_unmatched_CART.csv')

cartUndajusted1 <- trimming.funct(cartUndajusted1)
cartUndajusted2 <- trimming.funct(cartUndajusted2)
cartUndajusted3 <- trimming.funct(cartUndajusted3)
cartUndajusted4 <- trimming.funct(cartUndajusted4)

caliper1 = sd(cartUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(cartUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(cartUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(cartUndajusted4$propensity_score, na.rm = FALSE)*0.25

# nearest neighbor 

m_out_cart_nn1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "nearest",ratio=1,distance = cartUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE)
m_out_cart_nn2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "nearest",ratio=1, distance = cartUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE)
m_out_cart_nn3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "nearest",ratio=1, distance = cartUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_cart_nn4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "nearest",ratio=1, distance = cartUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE)

# genetic matching

m_out_cart_gen1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "genetic",distance = cartUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_cart_gen2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "genetic", distance = cartUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_cart_gen3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "genetic", distance = cartUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_cart_gen4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "genetic", distance = cartUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

# Get featured selected data

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

# nearest neighbor(feature selected)

m_out_cart_fs_nn1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "genetic",distance = cartUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_cart_fs_nn1 <- matchit(formula = forumla1, data = cartUndajusted2, method = "genetic", distance = cartUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_cart_fs_nn1 <- matchit(formula = forumla2, data = cartUndajusted3, method = "genetic", distance = cartUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_cart_fs_nn1 <- matchit(formula = forumla2, data = cartUndajusted4, method = "genetic", distance = cartUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

# genetic matching (feature selected)

m_out_cart_fs_gen1 <- matchit(formula = forumla1, data = cartUndajusted1, method = "genetic",distance = cartUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_cart_fs_gen2 <- matchit(formula = forumla1, data = cartUndajusted2, method = "genetic", distance = cartUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_cart_fs_gen3 <- matchit(formula = forumla2, data = cartUndajusted3, method = "genetic", distance = cartUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_cart_fs_gen4 <- matchit(formula = forumla2, data = cartUndajusted4, method = "genetic", distance = cartUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Get comparison  data

# Forest

forestUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswCps_lalonde_ps_unmatched_FOREST.csv')
forestUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswPsid_lalonde_ps_unmatched_FOREST.csv')
forestUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswCps_dehWab_ps_unmatched_FOREST.csv')
forestUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/forest/unmatched/nswPsid_dehWab_ps_unmatched_FOREST.csv')

forestUndajusted1 <- trimming.funct(forestUndajusted1)
forestUndajusted2 <- trimming.funct(forestUndajusted2)
forestUndajusted3 <- trimming.funct(forestUndajusted3)
forestUndajusted4 <- trimming.funct(forestUndajusted4)

caliper1 = sd(forestUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(forestUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(forestUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(forestUndajusted4$propensity_score, na.rm = FALSE)*0.25

# nearest neighbor 

m_out_forest_nn1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "nearest",ratio=1,distance = forestUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE)
m_out_forest_nn2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "nearest",ratio=1, distance = forestUndajusted2$propensity_score,caliper =caliper2,replace =  TRUE)
m_out_forest_nn3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "nearest",ratio=1, distance = forestUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_forest_nn4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "nearest",ratio=1, distance = forestUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE)

# genetic matching 
m_out_forest_gen1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "genetic",distance = forestUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_forest_gen2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "genetic", distance = forestUndajusted2$propensity_score,caliper =caliper2,replace =  TRUE,pop.size = 50)
m_out_forest_gen3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "genetic", distance = forestUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_forest_gen4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "genetic", distance = forestUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

# Get featured selected data

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

#nearest neighbor   

m_out_forest1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "nearest",ratio=5,distance = forestUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE)
m_out_forest2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "nearest",ratio=5, distance = forestUndajusted2$propensity_score,caliper =caliper2,replace =  TRUE)
m_out_forest3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "nearest",ratio=5, distance = forestUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_forest4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "nearest",ratio=5, distance = forestUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE)


#Genetic matching 

m_out_forest_gen1 <- matchit(formula = forumla1, data = forestUndajusted1, method = "genetic",distance = forestUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_forest_gen2 <- matchit(formula = forumla1, data = forestUndajusted2, method = "genetic", distance = forestUndajusted2$propensity_score,caliper =caliper2,replace =  TRUE,pop.size = 50)
m_out_forest_gen3 <- matchit(formula = forumla2, data = forestUndajusted3, method = "genetic", distance = forestUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_forest_gen4 <- matchit(formula = forumla2, data = forestUndajusted4, method = "genetic", distance = forestUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Boost 
boostUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswCps_lalonde_ps_unmatched_BOOST.csv')
boostUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswPsid_lalonde_ps_unmatched_BOOST.csv')
boostUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswCps_dehWab_ps_unmatched_BOOST.csv')
boostUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/boost/unmatched/nswPsid_dehWab_ps_unmatched_BOOST.csv')
# trimming 
boostUndajusted1 <- trimming.funct(boostUndajusted1)
boostUndajusted2 <- trimming.funct(boostUndajusted2)
boostUndajusted3 <- trimming.funct(boostUndajusted3)
boostUndajusted4 <- trimming.funct(boostUndajusted4)
#calipers
caliper1 = sd(boostUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(boostUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(boostUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(boostUndajusted4$propensity_score, na.rm = FALSE)*0.25
#nearest neighbor matching
m_out_boost_nn1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "nearest",ratio=1, distance = boostUndajusted1$propensity_score,caliper =caliper1,replace =  TRUE)
m_out_boost_nn2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "nearest",ratio=1, distance = boostUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE)
m_out_boost_nn3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "nearest",ratio=1, distance = boostUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_boost_nn4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "nearest",ratio=1, distance = boostUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE)
#genetic matching 
m_out_boost_gen1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "genetic", distance = boostUndajusted1$propensity_score,caliper =caliper1,replace =  TRUE,pop.size = 50)
m_out_boost_gen2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "genetic", distance = boostUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_boost_gen3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "genetic", distance = boostUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_boost_gen4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "genetic", distance = boostUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)


# Get featured selected data
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

# Nearest neighbor
m_out_boost_fs_nn1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "nearest",ratio=5, distance = boostUndajusted1$propensity_score,caliper =caliper1,replace =  TRUE)
m_out_boost_fs_nn2 <- matchit(formula = forumla1, data = boostUndajusted2, method = "nearest",ratio=5, distance = boostUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE)
m_out_boost_fs_nn3 <- matchit(formula = forumla2, data = boostUndajusted3, method = "nearest",ratio=5, distance = boostUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_boost_fs_nn4 <- matchit(formula = forumla2, data = boostUndajusted4, method = "nearest",ratio=5, distance = boostUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE)

# Genetic matching 
m_out_boost_fs_gen1 <- matchit(formula = forumla1, data = boostUndajusted1, method = "genetic", distance = boostUndajusted1$propensity_score,caliper =caliper1,replace =  TRUE,pop.size = 50)
m_out_boost_fs_gen1 <- matchit(formula = forumla1, data = boostUndajusted2, method = "genetic", distance = boostUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_boost_fs_gen1 <- matchit(formula = forumla2, data = boostUndajusted3, method = "genetic", distance = boostUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_boost_fs_gen1 <- matchit(formula = forumla2, data = boostUndajusted4, method = "genetic", distance = boostUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
#ANN 
annUndajusted1<- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswCps_lalonde_ps_unmatched_ANN.csv')
annUndajusted2<-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswPsid_lalonde_ps_unmatched_ANN.csv')
annUndajusted3 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswCps_dehWab_ps_unmatched_ANN.csv')
annUndajusted4 <-read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/quasi data/ann/unmatched/nswPsid_dehWab_ps_unmatched_ANN.csv')
#trimming
annUndajusted1 <- trimming.funct(annUndajusted1)
annUndajusted2 <- trimming.funct(annUndajusted2)
annUndajusted3 <- trimming.funct(annUndajusted3)
annUndajusted4 <- trimming.funct(annUndajusted4)
#calipers
caliper1 = sd(annUndajusted1$propensity_score, na.rm = FALSE)*0.25
caliper2 = sd(annUndajusted2$propensity_score, na.rm = FALSE)*0.25
caliper3 = sd(annUndajusted3$propensity_score, na.rm = FALSE)*0.25
caliper4 = sd(annUndajusted4$propensity_score, na.rm = FALSE)*0.25
#nearest neighbor
m_out_ann_nn1 <- matchit(formula = forumla1, data = annUndajusted1, method = "nearest",ratio=1, distance = annUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE)
m_out_ann_nn2 <- matchit(formula = forumla1, data = annUndajusted2, method = "nearest",ratio=1, distance = annUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE)
m_out_ann_nn3 <- matchit(formula = forumla2, data = annUndajusted3, method = "nearest",ratio=1, distance = annUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE)
m_out_ann_nn4 <- matchit(formula = forumla2, data = annUndajusted4, method = "nearest",ratio=1, distance = annUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE)
#genetic matching 
m_out_ann_gen1 <- matchit(formula = forumla1, data = annUndajusted1, method = "genetic", distance = annUndajusted1$propensity_score,caliper = caliper1,replace =  TRUE,pop.size = 50)
m_out_ann_gen2 <- matchit(formula = forumla1, data = annUndajusted2, method = "genetic", distance = annUndajusted2$propensity_score,caliper = caliper2,replace =  TRUE,pop.size = 50)
m_out_ann_gen3 <- matchit(formula = forumla2, data = annUndajusted3, method = "genetic", distance = annUndajusted3$propensity_score,caliper = caliper3,replace =  TRUE,pop.size = 50)
m_out_ann_gen4 <- matchit(formula = forumla2, data = annUndajusted4, method = "genetic", distance = annUndajusted4$propensity_score,caliper = caliper4,replace =  TRUE,pop.size = 50)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# Covariate balance plots 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

# LOGIT - balance plots
balance_plot_logit1 <- love.plot(forumla1,data = logitUndajusted1,
                                 weights = list(genMatch = m_out_logit_nn1,
                                                nearest = m_out_logit_gen1),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                   line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_logit1

balance_plot_logit2 <- love.plot(forumla1,data = logitUndajusted2,
                                 weights = list(genMatch = m_out_logit_nn2,
                                                nearest = m_out_logit_gen2),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde's with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_logit2

balance_plot_logit3 <- love.plot(forumla2,data = logitUndajusted3,
                                 weights = list(genMatch = m_out_logit_nn3,
                                                nearest = m_out_logit_gen3),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_logit3



balance_plot_logit4 <- love.plot(forumla2,data = logitUndajusted4,
                                   weights = list(genMatch = m_out_logit_nn4,
                                                  nearest = m_out_logit_gen4),
                                   stat = c("mean.diffs"),
                                   drop.distance = TRUE, 
                                   var.order = "unadjusted",
                                   abs = FALSE,
                                   line =TRUE, 
                                   stars = "raw",
                                   size = 3.5,
                                   shapes = c("circle filled", "circle filled","circle filled"),
                                   thresholds = c(m = .2),
                                   colors = c("#003366","#E31B23","#006627"),
                                   sample.names = c("unadjusted", "nearest neighbour","genetic"))+
    xlab("D-W with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_logit4


logit_lalonde <-ggarrange(balance_plot_logit1,balance_plot_logit2, 
                          balance_plot_logit3,balance_plot_logit4,
                          ncol = 4,nrow = 1, common.legend = TRUE, legend="right")
logit_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/Balance plots/LOGIT.png', plot = last_plot(),
      width = 22, height = 8 , device = "png", dpi=700)

# CART- balance plots
balance_plot_cart1 <- love.plot(forumla1,data = cartUndajusted1,
                                 weights = list(genMatch = m_out_cart_nn1,
                                                nearest = m_out_cart_gen1),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_cart1

balance_plot_cart2 <- love.plot(forumla1,data = cartUndajusted2,
                                 weights = list(genMatch = m_out_cart_nn2,
                                                nearest = m_out_cart_gen2),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde's with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_cart2

balance_plot_cart3 <- love.plot(forumla2,data = cartUndajusted3,
                                 weights = list(genMatch = m_out_cart_nn3,
                                                nearest = m_out_cart_gen3),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_cart3



balance_plot_cart4 <- love.plot(forumla2,data = cartUndajusted4,
                                 weights = list(genMatch = m_out_cart_nn4,
                                                nearest = m_out_cart_gen4),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_cart4


cart_lalonde <-ggarrange(balance_plot_cart1,balance_plot_cart2, 
                          balance_plot_cart3,balance_plot_cart4,
                          ncol = 4,nrow = 1, common.legend = TRUE, legend="right")
cart_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/Balance plots/cart.png', plot = last_plot(),
       width = 22, height = 8 , device = "png", dpi=700)

# boost- balance plots
balance_plot_boost1 <- love.plot(forumla1,data = boostUndajusted1,
                                weights = list(genMatch = m_out_boost_nn1,
                                               nearest = m_out_boost_gen1),
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line =TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled","circle filled"),
                                thresholds = c(m = .2),
                                colors = c("#003366","#E31B23","#006627"),
                                sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_boost1

balance_plot_boost2 <- love.plot(forumla1,data = boostUndajusted2,
                                weights = list(genMatch = m_out_boost_nn2,
                                               nearest = m_out_boost_gen2),
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line =TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled","circle filled"),
                                thresholds = c(m = .2),
                                colors = c("#003366","#E31B23","#006627"),
                                sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde's with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_boost2

balance_plot_boost3 <- love.plot(forumla2,data = boostUndajusted3,
                                weights = list(genMatch = m_out_boost_nn3,
                                               nearest = m_out_boost_gen3),
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line =TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled","circle filled"),
                                thresholds = c(m = .2),
                                colors = c("#003366","#E31B23","#006627"),
                                sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_boost3



balance_plot_boost4 <- love.plot(forumla2,data = boostUndajusted4,
                                weights = list(genMatch = m_out_boost_nn4,
                                               nearest = m_out_boost_gen4),
                                stat = c("mean.diffs"),
                                drop.distance = TRUE, 
                                var.order = "unadjusted",
                                abs = FALSE,
                                line =TRUE, 
                                stars = "raw",
                                size = 3.5,
                                shapes = c("circle filled", "circle filled","circle filled"),
                                thresholds = c(m = .2),
                                colors = c("#003366","#E31B23","#006627"),
                                sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_boost4


boost_lalonde <-ggarrange(balance_plot_boost1,balance_plot_boost2, 
                         balance_plot_boost3,balance_plot_boost4,
                         ncol = 4,nrow = 1, common.legend = TRUE, legend="right")
boost_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/Balance plots/boost.png', plot = last_plot(),
       width = 22, height = 8 , device = "png", dpi=700)

# forest- balance plots
balance_plot_forest1 <- love.plot(forumla1,data = forestUndajusted1,
                                 weights = list(genMatch = m_out_forest_nn1,
                                                nearest = m_out_forest_gen1),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_forest1

balance_plot_forest2 <- love.plot(forumla1,data = forestUndajusted2,
                                 weights = list(genMatch = m_out_forest_nn2,
                                                nearest = m_out_forest_gen2),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde's with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_forest2

balance_plot_forest3 <- love.plot(forumla2,data = forestUndajusted3,
                                 weights = list(genMatch = m_out_forest_nn3,
                                                nearest = m_out_forest_gen3),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_forest3



balance_plot_forest4 <- love.plot(forumla2,data = forestUndajusted4,
                                 weights = list(genMatch = m_out_forest_nn4,
                                                nearest = m_out_forest_gen4),
                                 stat = c("mean.diffs"),
                                 drop.distance = TRUE, 
                                 var.order = "unadjusted",
                                 abs = FALSE,
                                 line =TRUE, 
                                 stars = "raw",
                                 size = 3.5,
                                 shapes = c("circle filled", "circle filled","circle filled"),
                                 thresholds = c(m = .2),
                                 colors = c("#003366","#E31B23","#006627"),
                                 sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_forest4


forest_lalonde <-ggarrange(balance_plot_forest1,balance_plot_forest2, 
                          balance_plot_forest3,balance_plot_forest4,
                          ncol = 4,nrow = 1, common.legend = TRUE, legend="right")
forest_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/Balance plots/forest.png', plot = last_plot(),
       width = 22, height = 8 , device = "png", dpi=700)

# ann- balance plots
balance_plot_ann1 <- love.plot(forumla1,data = annUndajusted1,
                                  weights = list(genMatch = m_out_ann_nn1,
                                                 nearest = m_out_ann_gen1),
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line =TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled","circle filled"),
                                  thresholds = c(m = .2),
                                  colors = c("#003366","#E31B23","#006627"),
                                  sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_ann1

balance_plot_ann2 <- love.plot(forumla1,data = annUndajusted2,
                                  weights = list(genMatch = m_out_ann_nn2,
                                                 nearest = m_out_ann_gen2),
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line =TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled","circle filled"),
                                  thresholds = c(m = .2),
                                  colors = c("#003366","#E31B23","#006627"),
                                  sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("Lalonde's with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_ann2

balance_plot_ann3 <- love.plot(forumla2,data = annUndajusted3,
                                  weights = list(genMatch = m_out_ann_nn3,
                                                 nearest = m_out_ann_gen3),
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line =TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled","circle filled"),
                                  thresholds = c(m = .2),
                                  colors = c("#003366","#E31B23","#006627"),
                                  sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with CPS controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_ann3



balance_plot_ann4 <- love.plot(forumla2,data = annUndajusted4,
                                  weights = list(genMatch = m_out_ann_nn4,
                                                 nearest = m_out_ann_gen4),
                                  stat = c("mean.diffs"),
                                  drop.distance = TRUE, 
                                  var.order = "unadjusted",
                                  abs = FALSE,
                                  line =TRUE, 
                                  stars = "raw",
                                  size = 3.5,
                                  shapes = c("circle filled", "circle filled","circle filled"),
                                  thresholds = c(m = .2),
                                  colors = c("#003366","#E31B23","#006627"),
                                  sample.names = c("unadjusted", "nearest neighbour","genetic"))+
  xlab("D-W with PSID controls")+
  theme(legend.box.background = element_rect(),
        text = element_text(family = "Times New Roman",size = 18),
        legend.position="none",
        panel.border=element_rect(size=1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.minor.y = element_line(colour="white", size=0.5),
        panel.grid.major.y = element_line(colour="white", size=0.5),
        panel.grid.minor.x= element_line(colour="white", size=0.5),
        panel.grid.major.x= element_line(colour="white", size=0.5))+
  ggtitle("")
balance_plot_ann4


ann_lalonde <-ggarrange(balance_plot_ann1,balance_plot_ann2, 
                           balance_plot_ann3,balance_plot_ann4,
                           ncol = 4,nrow = 1, common.legend = TRUE, legend="right")
ann_lalonde
ggsave('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/Plots/Balance plots/ann.png', plot = last_plot(),
       width = 22, height = 8 , device = "png", dpi=700)

