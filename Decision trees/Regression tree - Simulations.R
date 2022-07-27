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

simulated_rct_data <- read.csv("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets/simulated_experimental_dataset.csv")
simulated_rct_data <- select(simulated_rct_data,-c("X")) # Drop ID col


sample = sample.split(simulated_rct_data$assignment, SplitRatio = .75)
train = subset(simulated_rct_data, sample == TRUE)
test  = subset(simulated_rct_data, sample == FALSE)

assignment_vec = test %>% select(c("assignment"))
test <- select(test,-c("assignment")) # drop assignment col from test

# ======== CART ========

model <- rpart(assignment ~. ,data = train,method = "class")
plot(model)
text(model, digits = 1)

#estimate classes 
predicted <-predict(model,newdata = test, type = "class")

# Accuracy 
predicted <- predicted %>% as.data.frame()
mean(predicted == assignment_vec)


# ======== RANDOM Forest ========

train_features <- select(train,-c("assignment"))
train_outcome <- select(train,c("assignment"))

train_outcome <- factor(train_outcome$assignment)


classifier_RF = randomForest(x = train_features,y = train_outcome,ntree = 500)

y_pred = predict(classifier_RF, newdata = test)

#Show accuracy 
mean(y_pred == factor(assignment_vec$assignment))

     