#######################################
#
# polyaIrisTest
#
# Sample code to demonstrate how to use polyaTrees.R on the Iris dataset. 
# The probability output will be higher for the "case" classification and lower for the "control" classification if the model can accurately separate the data.
#
# Usage:
#   R --no-save < polyaIrisTest.R
#
# Credit:  Brian Connolly / Pestian Lab
# License:  MIT (https://github.com/PestianLab/polyaTrees/blob/master/LICENSE)
# Copyright (c) 2018-present, Pestian Lab / Cincinnati Children's Hospital Medical Center
#
#######################################
library("e1071")
source("polyaTrees.R")

attach(iris)

# setup variables
train_size<-30
test_size<-20
case_label<-"setosa"
control_label<-"versicolor"

# get "case" training and testing data
setosaIndex <- iris$Species==case_label
setosa <- iris[setosaIndex,]
setosa_training <- head(setosa,train_size)
setosa_testing <- subset(tail(setosa,test_size),select=-Species)

# get "control" training and testing data
versaIndex <- iris$Species==control_label
versa <- iris[versaIndex,]
versa_training <- head(versa,train_size)
versa_testing <- subset(tail(versa,test_size),select=-Species)

# merge case and control training data
iris2_training <- rbind(setosa_training, versa_training)

# merge case and control testing data
iris2_testing <- rbind(setosa_testing, versa_testing)

# define training data for tuning svm
x <- subset(iris2_training, select=-Species)
y <- subset(iris2_training, select=Species)
y <- y$Species

# tune svm to get cost and gamma
svm_tuned <- tune(svm, train.x=x, train.y=y, kernel="radial", ranges=list(cost=test_size^(-1:2), gamma=c(.5,1,2)))

cost  <- as.numeric(svm_tuned$best.parameters["cost"])
gamma <- as.numeric(svm_tuned$best.parameters["gamma"])

# create an svm model
svm_model_after_tune <- svm(x,y,kernel="radial",cost=cost,gamma=gamma)

# predict on training data to get decision values for the training data
pred_after_tune <- predict(svm_model_after_tune,x,decision.values=TRUE)
decision_values <- attr(pred_after_tune, "decision.values")

# predict testing data
prediction_testing <- predict(svm_model_after_tune,iris2_testing,decision.values=TRUE)
prediction_row_names <- attr(prediction_testing, "names")
lvls <- c(case_label,control_label)
prediction_labels <- lvls[prediction_testing]

# get test data decision values
test_decision_values = attr(prediction_testing, "decision.values")

# split out training decision values into "case" and "control"
case_decision_values <- head(decision_values,train_size)
control_decision_values <- tail(decision_values,train_size)

i <- 1

# compute polya tree probability for each test decision value
for (d in test_decision_values) {
   cat("=====================================================\n")
   print(iris2_testing[i,])
   cat("Decision Value: ", d, "\n")
   cat("Prediction: ", prediction_labels[i], "\n")
   polya <- polyaTwoSampleComparison(case_decision_values, control_decision_values, d, mapDistributionsToSoftmax=FALSE, verbose=FALSE)
   cat("Polya Probability: ", polya, "\n")
   cat("=====================================================\n")
   i <- i + 1
}
