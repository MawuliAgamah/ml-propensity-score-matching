
# ======== load libraries ========  
library("rpart.plot")
library("rpart")
library("readxl")
library("gbm")
library("multcomp")
library("party")
library("haven")
library("dplyr")  
# ========  set working directory ======== 

setwd("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets") 
set.seed(290875)
# ========  Load data set ========  

# Lalonde treated and controls 
nswRe74_treat <- read_xlsx("nswre74_treated.xlsx")
nswRe74_control <- read_xls("nswre74_control.xls")

nswRe74_total <- rbind(nswRe74_treat, nswRe74_control) # drop outcome 
nswRe74_total <- within(nswRe74_total, rm(re78))
nswRe74_total$treat = factor(nswRe74_total$treat, levels = c(0, 1)) # turn treatment into factor 

# CPS control group 
cps_control <- read_dta("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/cps_controls.dta")
cps_control <- within(cps_control, rm(data_id)) # drop id 
cps_control <- within(cps_control, rm(re78))    # drop outcome 
cps_control$treat = factor(cps_control$treat, levels = c(0, 1)) # turn treatment into factor  



# ======== CART ========

# Classification tree 
model <- rpart(treat ~. ,data = nswRe74_total,method = "class")
par(xpd = NA) # otherwise on some devices the text is clipped
plot(model)
text(model, digits = 3)


#estimate probabilities

predicted_probabilities <- model %>% predict(cps_control,"prob")
predicted_classes <- model %>% predict(cps_control,"class") 
predicted_probabilities_df <- unname(unlist(predicted_probabilities)) %>% as.data.frame()
predicted_classes_df <- unname(unlist(predicted_classes)) %>% as.data.frame()








