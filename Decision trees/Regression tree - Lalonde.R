
# ======== load libraries ========  
library("rpart.plot")
library("rpart")
library("readxl")
library("gbm")
library("multcomp")
library("party")
library("haven")
library("dplyr")  
library("caTools")
library("randomForest")
# ========  set working directory ======== 

setwd("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets") 
set.seed(290875)
# ========  Load data sets ========  

# Lalonde

nswRe74_treat <- read_xlsx("nswre74_treated.xlsx")
nswRe74_control <- read_xls("nswre74_control.xls")

nswRe74_total <- rbind(nswRe74_treat, nswRe74_control) # drop outcome 
nswRe74_total <- within(nswRe74_total, rm(re78))
nswRe74_total$treat = factor(nswRe74_total$treat, levels = c(0, 1)) # turn treatment into factor 


# Creating training and testing split of Lalonde Dataset 
#sample <- sample.split(n = nrow(nswRe74_total), size = floor(.75*nrow(nswRe74_total)), replace = F)
#lalonde_train <- nswRe74_total[sample, ]
#lalode_test  <- nswRe74_total[-sample, ]


sample = sample.split(nswRe74_total$treat, SplitRatio = .75)
lalonde_train = subset(nswRe74_total, sample == TRUE)
lalode_test  = subset(nswRe74_total, sample == FALSE)

# CPS control group 
cps_control <- read_dta("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/cps_controls.dta")
cps_control <- within(cps_control, rm(data_id)) # drop id 
cps_control <- within(cps_control, rm(re78))    # drop outcome 
cps_control$treat = factor(cps_control$treat, levels = c(0, 1)) # turn treatment into factor  



# ======== CART ========

model <- rpart(treat ~. ,data = lalonde_train,method = "class")
plot(model)
text(model, digits = 3)

#estimate classes 
predicted.classes <- model %>% predict(lalode_test, type = "class")
# See accuracy 
mean(predicted.classes == lalode_test$treat)


#Store predicted probabilites 
predicted_probabilities <- model %>% predict(lalode_test,"prob")
predicted_classes <- model %>% predict(lalode_test,"class") 
predicted_probabilities_df <- unname(unlist(predicted_probabilities)) %>% as.data.frame()
predicted_classes_df <- unname(unlist(predicted_classes)) %>% as.data.frame()

#Accuracy 
mean(predicted.classes == lalode_test$treat)

# combine data frame with propensity scores 
cps_control_with_ps <- cps_control
cps_control_with_ps$PS <- predicted_probabilities_df[2]

# ======== RANDOM Forest ========
classifier_RF = randomForest(x = lalonde_train[-1],
                             y = lalonde_train$treat,
                             ntree = 500)


y_pred = predict(classifier_RF, newdata = lalode_test[-1])
mean(y_pred == lalode_test$treat)


