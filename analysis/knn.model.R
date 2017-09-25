### KNN
library(class)
library(lubridate)
library(DMwR)
library(dummies)

set.seed(100)

folds.length = 5;
folds = cvFolds(n = length(unique(dataPurged$trial)), K = folds.length, type = "random")

dataQuarters = dummy(dataPurged$quarter)
dataPurged = cbind(dataPurged, dataQuarters)
dataPurged = dataPurged[,which(names(dataPurged)!= "quarter")]

dataPurged <- within(dataPurged, {
  quarter1 <- as.numeric(quarter1)
  quarter2 <- as.numeric(quarter2)
  quarter3 <- as.numeric(quarter3)
  quarter4<-as.numeric(quarter4)
})

dataKnn = knnImputation(dataPurged[,3:ncol(dataPurged)])

k<-1:10

knn.perf = adply(.data = k, .margins = 1, function(k) { 
  
  fold.accuracy <- rep(0,5)
  for(n in 1:folds.length) {
  
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(data.l$trial)[test.ids]
  
  testdata = dataKnn[which(dataPurged$trial %in% test.trials), ] 
  traindata = dataKnn[-which(dataPurged$trial %in% test.trials), ] 
  
  class = dataPurged[-which(dataPurged$trial %in% test.trials), which(names(dataPurged)=="success")] 
  
  prediction <- knn(train = traindata, test = testdata, cl = class, k = k, prob=TRUE)
  fold.accuracy[n] = mean(prediction == dataPurged[which(dataPurged$trial %in% test.trials), which(names(dataPurged)=="success")])
  }
  return(perf = mean(fold.accuracy))
})


ggplot(data = knn.perf, aes(x = X1,y =V1)) + geom_point() + labs(x = "k", y = "Performance")

