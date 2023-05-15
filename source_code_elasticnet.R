rm(list=ls())

library(caret)
library(glmnet)
library(Hmisc)
library(AUC)
library('dplyr')
library(MLmetrics)
library(pROC)
library(e1071)
library(sva)
library(PRROC)
library(doMC) 
registerDoMC(cores = 8)

### read in cell proportion data and responder class assignment
orig_data<-read.table('~/cellproportions_inv_baseline.txt',sep='\t',header=T)
colnames(orig_data)[1] = 'Class'

### set up model parameters
ctrl2 <- trainControl(method = "repeatedcv", repeats = 5, number = 10, classProbs=TRUE, sampling='up', 
                      summaryFunction=twoClassSummary, savePredictions = 'all', allowParallel = TRUE)

###########Elastic Nets model with nested cross-validation ############
ncv<-100;# number of CV rounds
set.seed(5627)
trueClassENET_out<-NULL;
predClassENET_prob<-NULL;
predClassENET_binary<-NULL
auc.enetModel_per_sample_out<-NULL
pr_non_response_enetModel_per_sample_out <-NULL
orig_data_zero=orig_data[orig_data[,1] == 'zero',]
orig_data_one=orig_data[orig_data[,1] == 'one',]
nrA<-nrow(orig_data_zero);
nrB<-nrow(orig_data_one);
ntestA<-0.2*nrA
ntestB<-0.2*nrB
for(i in 1:ncv){# for cv rounds
  ##############################################################################
  # set up training and test split 20:80########################################
  ##############################################################################
  testIndexA<-sample(1:nrA,ntestA)
  testIndexB<-sample(1:nrB,ntestB)
  trainingIndexA<-setdiff(1:nrA,testIndexA)
  trainingIndexB<-setdiff(1:nrB,testIndexB)
  trainData<-rbind(orig_data_zero[trainingIndexA,],orig_data_one[trainingIndexB,])
  testData<-rbind(orig_data_zero[testIndexA,],orig_data_one[testIndexB,])
  ##############################################################################
  # preprocess cell proportion data: scale #####################################
  ##############################################################################
  preProc = preProcess(trainData[,-1], method = c("center", "scale"))
  corrected_train_s = predict(preProc,trainData[,-1])
  corrected_test_s = predict(preProc,testData[,-1])
  dataTrain = data.frame(as.factor(trainData[,1]),corrected_train_s)
  colnames(dataTrain)[1] = 'Class'
  dataTest = data.frame(as.factor(testData[,1]),corrected_test_s)
  colnames(dataTest)[1] = 'Class'
  ##############################################################################
  # optimize model parameters towards ROC ######################################
  ##############################################################################
  fit_binary_enet <- train(Class ~ ., data = dataTrain, method = "glmnet", trControl = ctrl2, 
                           tuneGrid = expand.grid(alpha = seq(0, 1, length = 5),lambda = seq(0.0001, 1, length = 100)), metric='ROC',family="binomial")
  ##############################################################################
  # assess performance of model on kept out test data and save parameters ######
  ##############################################################################
  predClassENET_prob<-rbind(predClassENET_prob,predict(object=fit_binary_enet, newdata=dataTest[,-c(1)], type = 'prob'))
  predClassENET_binary<-rbind(predClassENET_prob,as.numeric(predict(object=fit_binary_enet, newdata=dataTest[,-c(1)])))
  trueClassENET_out<-rbind(trueClassENET_out,matrix(dataTest[,c(1)]))
  roc.enetModel_per_sample <- pROC::roc(factor(dataTest[,c(1)]), enet_Pred_prob$zero)
  auc.enetModel_per_sample_out<-rbind(auc.enetModel_per_sample_out,pROC::auc(roc.enetModel_per_sample))
  pr1 = pr.curve(scores.class0=enet_Pred_prob$zero[dataTest[,c(1)] =="zero"],scores.class1=enet_Pred_prob$zero[dataTest[,c(1)] =="one"], curve=F)
  pr_non_response_enetModel_per_sample_out<-rbind(pr_non_response_enetModel_per_sample_out,cbind(pr1$auc.integral,pr1$auc.davis.goadrich))
  save(pr_non_response_enetModel_per_sample_out,ctrl2,predClassENET_prob,predClassENET_binary,trueClassENET_out, roc.enetModel_per_sample, 
       auc.enetModel_per_sample, auc.enetModel_per_sample_out, file=paste('classification_binary_nonresponse_combine_baseline_cellproportions_inv_glmnet_model_cv_',i,'.RData',sep=''))
}

save(pr_non_response_enetModel_per_sample_out,ctrl2,predClassENET_binary,predClassENET_prob,trueClassENET_out, auc.enetModel_per_sample_out,file='classification_binary_nonresponse_combine_baseline_cellproportions_inv_glmnet_final_model.RData')

q('no')     



