library(caret)
library(reshape2)
library(tidyverse)
library(parallel)
library(parallelMap)
args = commandArgs(T)
filename <- args[1]
geneName <- read.csv(filename,header = FALSE)

## process raw data
parallelStartSocket(cpus = detectCores())
for (j in c(1:length(geneName$V1))){
  acc_glm <- c()
  pvalue_glm <- c()
  AUC_glm <- c()
  gene <- geneName$V1[j]
  geneframe <- data.frame(matrix(ncol = 4))
  
  setwd("/exports/eddie/scratch/s1804132/main/combine")
  Data <- read.csv(paste(gene,"csv",sep = "."),stringsAsFactors = FALSE)
  data <- Data %>% select(-c(X,population))
  data <- data/rowSums(data)
  data[is.na(data)] <- 0
  
  setwd("/exports/eddie/scratch/s1804132/main/Data/winModel/win0/lenWin")
  nor_length <- read.csv(paste(gene,"nor.csv",sep = "_"),header = TRUE)
  nor_len <- nor_length$fraction
  nor <- as.data.frame(t(apply(data, 1, function(x) x/nor_len)))
  nor[is.na(nor)] <- 0
  nor <- do.call(data.frame,lapply(nor, function(x) replace(x, is.infinite(x),0)))
  col_num <- apply(nor,2,function(x) length(unique(x)) == 1)
  
  data <- data[,!col_num]
  x <- cor(data)
  data <- data[,-findCorrelation(x,cutoff = 0.9)]
  if (ncol(data) == 0){next}
  
  data$population <- Data$population
  data$population[which(data$population %in% c("CEU","TSI","FIN","GBR"))] <- "EUR"
  data$population <- as.factor(data$population)

  i = 1 ## revise here if repeating test
  inTrainVal <- createDataPartition (y = data$population,
                                     p = 0.8,
                                     list = FALSE)
  training <- data[inTrainVal,]
  testing <- data[-inTrainVal,]
  up_train <- upSample(x = training[,-ncol(training)],
                       y = training$population)
  names(up_train)[names(up_train) == "Class"] <- "population"
  
  # set control about training methods
  fitControl <- trainControl(method = "cv",
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary)
  
  glm_label <- train (population ~.,
                      data=up_train,
                      method = "glm",
                      family = "binomial",
                      metric = "ROC",
                      control=list(maxit=150),
                      preProcess = c("center","scale"),
                      trControl = fitControl)
  glm_predict <- predict(glm_label, newdata = testing)
  glm_matrix <- confusionMatrix(data = glm_predict,
                                testing$population)
  AUC_glm[i] <- round(glm_label$results$ROC,3)
  acc_glm[i] <- round(glm_matrix$overall[1],3)
  pvalue_glm[i] <- round(glm_matrix$overall[6],3)

  geneframe[1,1] <- as.character(gene)
  geneframe[1,2] <- round(acc_glm,3)
  geneframe[1,3] <- round(AUC_glm,3)
  geneframe[1,4] <- round(pvalue_glm,3)

  write.table(geneframe,
              file=paste("~","res_CompareModel_14.csv",sep = "/"),
              append=TRUE,
              col.names = FALSE)

}