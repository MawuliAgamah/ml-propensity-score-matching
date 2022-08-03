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
library("rattle")
library("rpart.plot")
library("RColorBrewer")
library("xgboost")
library("DiagrammeR")

# ========  set working directory ======== 

setwd("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets") 
set.seed(290875)
# ========  Load data sets ========  

# Lalonde

nswRe74_treat <- read_xlsx("nswre74_treated.xlsx")
nswRe74_control <- read_xls("nswre74_control.xls")

nswRe74_total <- rbind(nswRe74_treat, nswRe74_control) 
#nswRe74_total <- within(nswRe74_total, rm(re78))

nswRe74_total$treat = factor(nswRe74_total$treat, levels = c(0, 1)) # turn treatment into factor 


# Creating train and test splits for Lalonde data

sample = sample.split(nswRe74_total$treat, SplitRatio = .75)
lalonde_train = subset(nswRe74_total, sample == TRUE)
lalode_test  = subset(nswRe74_total, sample == FALSE)

assignment_vec = lalode_test %>% select(c("treat"))
lalode_test <- select(lalode_test,-c("treat")) # drop assignment col from test

# CPS control group 
cps_control <- read_dta("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/cps_controls.dta")
cps_control <- within(cps_control, rm(data_id)) # drop id 
cps_control <- within(cps_control, rm(re78))    # drop outcome 
cps_control$treat = factor(cps_control$treat, levels = c(0, 1)) # turn treatment into factor  



# ======== CART ========

model <- rpart(treat ~. ,data = lalonde_train,method = "class")
plot(model)
text(model, digits = 3)


fancyRpartPlot(model, caption = NULL)
#estimate classes 
predicted  <- model %>% predict(lalode_test, type = "class")
# See accuracy 
mean(predicted == factor(assignment_vec$treat))


# Predict probabilities for CPS groups 




#Store predicted probability 

predicted_probabilities <- model %>% predict(lalode_test,"prob")
predicted_classes <- model %>% predict(lalode_test,"class") 
predicted_probabilities_df <- unname(unlist(predicted_probabilities)) %>% as.data.frame()
predicted_classes_df <- unname(unlist(predicted_classes)) %>% as.data.frame()

# Merge estimated propensity scores with the CPS dataframe 
cps_control_with_ps <- cps_control
cps_control_with_ps$PS <- predicted_probabilities_df[2]

# ======== XGBOOST ========

train_features <- select(lalonde_train,-c("treat"))
train_outcome <- select(lalonde_train,c("treat"))
train_outcome <- factor(train_outcome$treat)

features_matrix <- train_features %>% as.matrix()
labels_matrix <- train_outcome %>% as.matrix()

bst <- xgboost(data = features_matrix, label = labels_matrix, max.depth = 8, eta = 2, nthread = 2, nrounds = 5, objective = "binary:logistic")
test_matrix <- lalode_test %>% as.matrix()
pred <- predict(bst, test_matrix)
prediction <- as.numeric(pred > 0.5)
mean(prediction == assignment_vec$treat)
xgb.plot.tree(model = bst)

# ======== RANDOM FOREST ========

train_features <- select(lalonde_train,-c("treat"))
train_outcome <- select(lalonde_train,c("treat"))
train_outcome <- factor(train_outcome$treat)

classifier_RF = randomForest(x = train_features,
                             y = train_outcome,
                             ntree = 500)

y_pred = predict(classifier_RF, newdata = lalode_test)

mean(y_pred == factor(assignment_vec$treat))


