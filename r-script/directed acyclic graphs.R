
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# load required libraries 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----

library("dagitty")
library("ggdag")
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# global set up and information 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
options(scipen=999)
set.seed(1234)

# KEY 
# 1 - nsw treated + cps control  (lalonde's original sample)
# 2 - nsw treated + psid control (lalonde's original sample)
# 3 - nsw treated + cps control  (dehejia & wahba sub-sample)
# 4 - nsw treated + psid control (dehejia & wahba sub-sample)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
 # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
