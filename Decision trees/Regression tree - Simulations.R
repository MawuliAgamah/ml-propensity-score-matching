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

simulated_rct_data <- read.csv("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/simulated_experimental_dataset.csv")

simulated_rct_data <- within(simulated_rct_data, rm(X))
simulated_rct_data <- within(simulated_rct_data, rm(outcome)) # drop outcome
simulated_rct_data$assignment = factor(simulated_rct_data$assignment, levels = c(0, 1)) # turn treatment into factor 

sample = sample.split(simulated_rct_data$assignment, SplitRatio = .75)
train = subset(simulated_rct_data, sample == TRUE)
test  = subset(simulated_rct_data, sample == FALSE)



# ======== CART ========

model <- rpart(assignment ~. ,data = train,method = "class")
plot(model)
text(model, digits = 3)

#estimate classes 
predicted.classes <- model %>% predict(test, type = "class")
# See accuracy 
mean(predicted.classes == train$assignment)


#Store predicted probabilites 
predicted_probabilities <- model %>% predict(lalode_test,"prob")
predicted_classes <- model %>% predict(lalode_test,"class") 
predicted_probabilities_df <- unname(unlist(predicted_probabilities)) %>% as.data.frame()
predicted_classes_df <- unname(unlist(predicted_classes)) %>% as.data.frame()



# combine data frame with propensity scores 
cps_control_with_ps <- cps_control
cps_control_with_ps$PS <- predicted_probabilities_df[2]

# ======== RANDOM Forest ========
classifier_RF = randomForest(x = train[-1],
                             y = train$assignment,
                             ntree = 100)


y_pred = predict(classifier_RF, newdata = test[-1])
mean(y_pred == train$assignment)
