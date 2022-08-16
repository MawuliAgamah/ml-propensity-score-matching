
library("ggplot2") # visualization
library("MatchIt") # Matching
library("cobalt") #  covariate balance

# Load and organise matched datasets 

# 1 - nsw lalonde + cps contol
# 2 - nsw lalonde + psid contol
# 3 - nsw dehejia & wahba + cps contol
# 3 - nsw dehejia & wahba + psid contol


# LOGIT matched and unmatched ----

# adjusted logit data set's
LogitMatched1 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswCps_lalonde_LOGIT_psMatched.csv')
logitMatched2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswPsid_lalonde_LOGIT_psMatched.csv')
logitmatched3 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswCps_dehWab_LOGIT_psMatched.csv')
logitMatched4 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswPsid_dehWab_LOGIT_psMatched.csv')

LogitMatched1$comparison_group <- "cps"
logitMatched2$comparison_group <- "psid"
logitmatched3$comparison_group <- "cps"
logitMatched4$comparison_group <- "psid"

LogitMatched1$sample <- "lalonde"
logitMatched2$sample <- "lalonde"
logitmatched3$sample <- "dehejia_wahba"
logitMatched4$sample <- "dehejia_wahba"

# un-adjusted logit data set's
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

Call: 
  matchit(formula = catholic ~ race_white + w3income + p5hmage + 
            p5numpla + w3momed_hsb, data = ecls_nomiss, method = "nearest", 
          distance = "logit", replace = TRUE)

# define lists

logitMatched_list <- list(LogitMatched1 = "lalonde_cps",LogitMatched2 = "lalonde_psid",LogitMatched3 = "dehWab_cps",LogitMatched4 = "dehWab_psid") # dataframes with a list 
logitunadjusted_list <- list(logitUndajusted1 = "lalonde_cps",logitUndajusted2 = "lalonde_psid",logitUndajusted3 = "dehWab_cps",logitUndajusted3 = "dehWab_psid") # dataframes with a list 

Map(cbind, logitMatched_list, model = "logit")    # add column to identify the model used for matching
Map(cbind, logitunadjusted_list, model = "logit") # add column to identify the model used for matching


# CART 


CART_Matched1 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswCps_lalonde_CART_psMatched.csv')
CART_Matched2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswPsid_lalonde_CART_psMatched.csv')
CART_Matched3 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswCps_dehWab_CART_psMatched.csv')
CART_Matched4 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswPsid_dehWab_CART_psMatched.csv')

CART_Matched1$comparison_group <- "cps"
CART_Matched2$comparison_group <- "psid"
CART_Matched3$comparison_group <- "cps"
CART_Matched4$comparison_group <- "psid"

CART_Matched1$sample <- "lalonde"
CART_Matched2$sample <- "lalonde"
CART_Matched3$sample <- "dehejia_wahba"
CART_Matched4$sample <- "dehejia_wahba"

cart_list <- list(CART_Matched1 = "lalonde_cps",CART_Matched2 = "lalonde_psid",CART_Matched3 = "dehWab_cps",CART_Matched4 = "dehWab_psid") # dataframes with a list 

Map(cbind, cart_list, model = "cart") # add column to identify the model used for matching

# Random Forest 
FOREST_Matched1 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/nswCps_lalonde_forest_psMatched.csv')
FOREST_Matched2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/nswPsid_lalonde_forest_psMatched.csv')
FOREST_Matched3 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/nswCps_dehWab_forest_psMatched.csv')
FOREST_Matched4 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/forest/nswPsid_dehWab_forest_psMatched.csv')

forest_list <- list(FOREST_Matched1,FOREST_Matched2,FOREST_Matched3,FOREST_Matched4)

# Boosted trees 
BOOST_Matched1 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/nswCps_lalonde_boost_psMatched.csv')
BOOST_Matched2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/nswPsid_lalonde_boost_psMatched.csv')
BOOST_Matched3 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/nswCps_dehWab_boost_psMatched.csv')
BOOST_Matched4 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/boost/nswPsid_dehWab_boost_psMatched.csv')


# ANN 


colnames(LOGIT_Matched1)

# balance plots 

loveplotfunct = treat ~ age + education. + black + hispanic + married + nodegree + re75 + re78 + propensity_score
love.plot(loveplotfunct, data = LOGIT_Matched4,stat = c("mean.diffs", "variance.ratios"), thresholds = c(m = .1))






# BOXPLOT prep ----
# Create dataset of all propensity scores


logit_list[["lalonde_cps"]]

# Create boxplots 
ggplot(, 
       aes( y=propensity_score,x=as.factor(treat),
            fill=as.factor(treat))) + 
  stat_boxplot(geom ='errorbar',width = 0.2) + 
  geom_boxplot()


LOGIT_Matched1$treat
