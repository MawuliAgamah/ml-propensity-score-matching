
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
library("pcalg")
library("pcalg")
#install.packages("/Users/mawuliagamah/gitprojects/causal_inference/r packages/pcalg", repos = NULL, type="source",dependencies=TRUE)


remove.packages("pcalg", lib = "/Users/mawuliagamah/gitprojects/causal_inference/r packages/pcalg")
#install.packages("/Users/mawuliagamah/gitprojects/causal_inference/r packages/graph",repos=NULL, type="source")
#install.packages("/Users/mawuliagamah/gitprojects/causal_inference/r packages/BiocGenerics",repos=NULL, type="source")
#install.packages("/Users/mawuliagamah/gitprojects/causal_inference/r packages/RBGL",repos=NULL, type="source")
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
 # Load the dataset's we need 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
# make it work for one dataset and then roll everything into an objet/function

#class(s) <- "DAG"

# Key variables needed for DAGs
# 1 = individual/unit 
# 2 = outcome variable = wage in 1978 
# 3 = intervention/treatment variable  = treat column




baseline_dataset_1 <- read.csv("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/data/nsw_featureTransformed.csv")
covariates <- colnames(baseline_dataset_1)
covariates[1]

dagified <- dagify(x ~ z,
                   y ~ z,
                   exposure = "x",
                   outcome = "y")
tidy_dagitty(dagified)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ----
