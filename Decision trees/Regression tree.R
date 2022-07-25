install.packages("rpart.plot")
library("rpart")
library(rpart.plot)
library(readxl)
library(gbm)

setwd("/Users/mawuliagamah/gitprojects/causal_inference/causal_inference/datasets")

# ==== Load data set ==== 

nswRe74_treat <- read_xlsx("nswre74_treated.xlsx")
nswRe74_control <- read_xls("nswre74_control.xls")
nswRe74_total <- rbind(nswRe74_treat, nswRe74_control)

nswRe74_control <- drop(nswRe74_control$re78)

gb <- rpart(treat ~ ., data = nswRe74_total, method = "anova")
rpart.plot(gb,type = 3 , digits = 3 , fallen.leaves = TRUE)


