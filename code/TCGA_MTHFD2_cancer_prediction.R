# libraries
library(tidyverse)
library(gtools)
library(ggpubr)
library(kernlab)
library(caret)
library(randomForest)
library(ROCR)

# 1. Read input ------------------------------------------------------------------------------------------------------------------

# obtain MTHFD2 sample sheet
MTHFD2_sample_sheet <- read_tsv("../files/MTHFD2_expression_normal_vs_tumor.tsv")
str(MTHFD2_sample_sheet) # 7,048 × 11
colnames(MTHFD2_sample_sheet)

MTHFD2_sample_sheet$SampleType <- factor(MTHFD2_sample_sheet$SampleType, levels = c("Solid Tissue Normal", "Primary Tumor"))
table(MTHFD2_sample_sheet$SampleType)

# Separate breast, lung and colon cancer
breast_MTHFD2_data <- MTHFD2_sample_sheet %>%
  filter(ProjectID == "TCGA-BRCA") %>%
  select(SampleID, SampleType, MTHFD2_expression)
table(breast_MTHFD2_data$SampleType)

lung_MTHFD2_data <- MTHFD2_sample_sheet %>%
  filter(ProjectID == "TCGA-LUAD" | ProjectID == "TCGA-LUSC") %>%
  select(SampleID, SampleType, MTHFD2_expression)
table(lung_MTHFD2_data$SampleType)

colon_MTHFD2_data <- MTHFD2_sample_sheet %>%
  filter(ProjectID == "TCGA-COAD") %>%
  select(SampleID, SampleType, MTHFD2_expression)
table(colon_MTHFD2_data$SampleType)

# 2. Machine learning ----------------------------------------------------------------------------------------------------------

# remove ID
breast_MTHFD2_data_noID <- breast_MTHFD2_data[,-1]
lung_MTHFD2_data_noID <- lung_MTHFD2_data[,-1]
colon_MTHFD2_data_noID <- colon_MTHFD2_data[,-1]

# factor
breast_MTHFD2_data_noID$SampleType <- factor(breast_MTHFD2_data_noID$SampleType, levels = c("Primary Tumor", "Solid Tissue Normal"))
lung_MTHFD2_data_noID$SampleType <- factor(lung_MTHFD2_data_noID$SampleType, levels = c("Primary Tumor", "Solid Tissue Normal"))
colon_MTHFD2_data_noID$SampleType <- factor(colon_MTHFD2_data_noID$SampleType, levels = c("Primary Tumor", "Solid Tissue Normal"))

# ANALYSIS FOR BREAST

# separate training and test
set.seed(12345)
random_ID <- order(runif(nrow(breast_MTHFD2_data_noID)))

breast_train <- breast_MTHFD2_data_noID[random_ID[1:round((1-0.3)
                                       *nrow(breast_MTHFD2_data_noID),0)], ]
breast_test <- breast_MTHFD2_data_noID[random_ID[(1+round((1-0.3)
                                       *nrow(breast_MTHFD2_data_noID),0)):
                                nrow(breast_MTHFD2_data_noID)], ]

dim(breast_train)
dim(breast_test)

# model random forest
set.seed(1234567)
rf_breast <- randomForest(SampleType ~ ., data = breast_train)

# predictions
rf_breast_pred <- predict(rf_breast, breast_test)

# evaluate performance
rf_breast_cm <- confusionMatrix(reference = breast_test$SampleType, rf_breast_pred,
                                positive = "Primary Tumor")
rf_breast_cm

# probabilities
rf_breast_pred_prob <- predict(rf_breast, breast_test, type = "prob")

# get ROC curves
pred_breast <- data.frame(actual_type = breast_test$SampleType,
                        predicted_type_rf = rf_breast_pred,
                        prob_rf = rf_breast_pred_prob[,"Primary Tumor"])

predict_breast_prob_rf <- prediction(pred_breast$prob_rf, labels = pred_breast$actual_type,
                                   label.ordering = c("Solid Tissue Normal", "Primary Tumor"))

perf_breast_rf <- performance(predict_breast_prob_rf, measure = "tpr", x.measure = "fpr") 

pdf("../MTHFD2_cancer_prediction/ROC_breast_rf.pdf", width = 4, height =4)
plot(perf_breast_rf, col = "#E27AD8", lwd = 3, main = "Breast Cancer")
abline(a = 0, b = 1, lwd = 2, lty = 2)
dev.off()

# get AUC
perf.auc_breast_rf <- performance(predict_breast_prob_rf, measure = "auc")
unlist(perf.auc_breast_rf@y.values)

# ANALYSIS FOR LUNG

# separate training and test
set.seed(12345)
random_ID <- order(runif(nrow(lung_MTHFD2_data_noID)))

lung_train <- lung_MTHFD2_data_noID[random_ID[1:round((1-0.3)
                                                          *nrow(lung_MTHFD2_data_noID),0)], ]
lung_test <- lung_MTHFD2_data_noID[random_ID[(1+round((1-0.3)
                                                          *nrow(lung_MTHFD2_data_noID),0)):
                                                   nrow(lung_MTHFD2_data_noID)], ]

dim(lung_train)
dim(lung_test)

# model random forest
set.seed(1234567)
rf_lung <- randomForest(SampleType ~ ., data = lung_train)

# predictions
rf_lung_pred <- predict(rf_lung, lung_test)

# evaluate performance
rf_lung_cm <- confusionMatrix(reference = lung_test$SampleType, rf_lung_pred, 
                              positive = "Primary Tumor")
rf_lung_cm

# probabilities
rf_lung_pred_prob <- predict(rf_lung, lung_test, type = "prob")

# get ROC curves
pred_lung <- data.frame(actual_type = lung_test$SampleType,
                         predicted_type_rf = rf_lung_pred,
                         prob_rf = rf_lung_pred_prob[,"Primary Tumor"])

predict_lung_prob_rf <- prediction(pred_lung$prob_rf, labels = pred_lung$actual_type,
                                    label.ordering = c("Solid Tissue Normal", "Primary Tumor"))

perf_lung_rf <- performance(predict_lung_prob_rf, measure = "tpr", x.measure = "fpr") 

pdf("../MTHFD2_cancer_prediction/ROC_lung_rf.pdf", width = 4, height =4)
plot(perf_lung_rf, col = "#669DDB", lwd = 3, main = "Lung Cancer")
abline(a = 0, b = 1, lwd = 2, lty = 2)
dev.off()

# get AUC
perf.auc_lung_rf <- performance(predict_lung_prob_rf, measure = "auc")
unlist(perf.auc_lung_rf@y.values)

# ANALYSIS FOR COLON

# separate training and test
set.seed(12345)
random_ID <- order(runif(nrow(colon_MTHFD2_data_noID)))

colon_train <- colon_MTHFD2_data_noID[random_ID[1:round((1-0.3)
                                                      *nrow(colon_MTHFD2_data_noID),0)], ]
colon_test <- colon_MTHFD2_data_noID[random_ID[(1+round((1-0.3)
                                                      *nrow(colon_MTHFD2_data_noID),0)):
                                               nrow(colon_MTHFD2_data_noID)], ]

dim(colon_train)
dim(colon_test)

# model random forest
set.seed(1234567)
rf_colon <- randomForest(SampleType ~ ., data = colon_train)

# predictions
rf_colon_pred <- predict(rf_colon, colon_test)

# evaluate performance
rf_colon_cm <- confusionMatrix(reference = colon_test$SampleType, rf_colon_pred, 
                               positive = "Primary Tumor")
rf_colon_cm

# probabilities
rf_colon_pred_prob <- predict(rf_colon, colon_test, type = "prob")

# get ROC curves
pred_colon <- data.frame(actual_type = colon_test$SampleType,
                         predicted_type_rf = rf_colon_pred,
                         prob_rf = rf_colon_pred_prob[,"Primary Tumor"])

predict_colon_prob_rf <- prediction(pred_colon$prob_rf, labels = pred_colon$actual_type,
                                    label.ordering = c("Solid Tissue Normal", "Primary Tumor"))

perf_colon_rf <- performance(predict_colon_prob_rf, measure = "tpr", x.measure = "fpr") 


pdf("../MTHFD2_cancer_prediction/ROC_colon_rf.pdf", width = 4, height =4)
plot(perf_colon_rf, col = "#DFA247", lwd = 3, main = "Colon Cancer")
abline(a = 0, b = 1, lwd = 2, lty = 2)
dev.off()

# get AUC
perf.auc_colon_rf <- performance(predict_colon_prob_rf, measure = "auc")
unlist(perf.auc_colon_rf@y.values)
