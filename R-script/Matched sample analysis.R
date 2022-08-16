
# Load ggplot2
library("ggplot2")


# Load matched datasets 

# 1 - nsw lalonde + cps contol
# 2 - nsw lalonde + psid contol
# 3 - nsw dehejia & wahba + cps contol
# 3 - nsw dehejia & wahba + psid contol

# LOGIT 

LOGIT_Matched1 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswCps_lalonde_LOGIT_psMatched.csv')
LOGIT_Matched2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswPsid_lalonde_LOGIT_psMatched.csv')
LOGIT_Matched3 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswCps_dehWab_LOGIT_psMatched.csv')
LOGIT_Matched4 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/logit/nswPsid_dehWab_LOGIT_psMatched.csv')

LOGIT_Matched1$comparison_group <- "cps"
LOGIT_Matched2$comparison_group <- "psid"
LOGIT_Matched3$comparison_group <- "cps"
LOGIT_Matched4$comparison_group <- "psid"

LOGIT_Matched1$sample <- "lalonde"
LOGIT_Matched2$sample <- "lalonde"
LOGIT_Matched3$sample <- "dehejia_wahba"
LOGIT_Matched4$sample <- "dehejia_wahba"

logit_list <- list(LOGIT_Matched1 = "lalonde_cps",LOGIT_Matched2 = "lalonde_psid",LOGIT_Matched3 = "dehWab_cps",LOGIT_Matched4 = "dehWab_psid") # dataframes with a list 

Map(cbind, logit_list, model = "logit") # add column to indentify the model used for matching

# CART 
CART_Matched1 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswCps_lalonde_CART_psMatched.csv')
CART_Matched2 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswPsid_lalonde_CART_psMatched.csv')
CART_Matched3 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswCps_dehWab_CART_psMatched.csv')
CART_Matched4 <- read.csv('/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/quasi data/cart/nswPsid_dehWab_CART_psMatched.csv')

cart_list <- list(CART_Matched1,CART_Matched2,CART_Matched3,CART_Matched4)

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

boost_list <- list(BOOST_Matched1,BOOST_Matched2,BOOST_Matched3,BOOST_Matched4)

# ANN 



# Create dataset of all propensity scores


logit_lalonde_cps <- rbind(LOGIT_Matched1$treat,LOGIT_Matched1$propensity_score)
logit_lalonde_cps <- rbind(LOGIT_Matched2$treat,LOGIT_Matched2$propensity_score)
logit_dehWab_cps <- rbind(LOGIT_Matched3$treat,LOGIT_Matched3$propensity_score)
logit_dehWab_cps <- rbind(LOGIT_Matched4$treat,LOGIT_Matched4$propensity_score)
logit 

# Create boxplots 
ggplot(BOOST_Matched1, 
       aes( y=propensity_score,x=as.factor(treat),
            fill=as.factor(treat))) + 
  stat_boxplot(geom ='errorbar',width = 0.2) + 
  geom_boxplot()


LOGIT_Matched1$treat
