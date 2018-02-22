set.seed(62417)
folds.length = 5;
dataPurged %<>% rename(name = trial) 

folds = cvFolds(n = length(unique(dataPurged$name)), K = folds.length, type = "random")

glmerperf=rep(0, folds.length); glmperf=glmerperf;
data.summary.mat = data.frame(matrix(nrow = folds.length, ncol  = 2))
colnames(data.summary.mat) = c( "AIC","p.value")

#Initializing Variables 
build = TRUE
variable.list = c()
variable.per = c()
variable.aic = c()
variable.min.per = c()
variable.max.per = c()
anova.results = rep(0, folds.length)

fixed.part.1 = "success ~ value + "
fixed.part.2 = " + (1 | name)"

dataPurged$success <- relevel(dataPurged$success, ref = "Transient")
dataPurged %<>% select(-group)

while(build == TRUE) {
  # Step 1: Sets the formula and manipulates the data set
  if(length(variable.list) == 0 )  { # This is the first round
    formula = "success~value + (1 | name)"
    dataPurged %>%
      ungroup() %>%
      gather(key = variable, value = value, -success, -name) %>%
      mutate_at("value", as.numeric) -> data.scaled.l
  } else { # There are already significant variables 
    variable.part = paste(variable.list, collapse = "+")
    formula = paste0(fixed.part.1, variable.part, fixed.part.2)
    formula.null = paste0("success ~", paste0(variable.list, collapse = "+"), fixed.part.2)
    
    dataPurged %>%
      ungroup() %>%
      gather(key = variable, value = value, -success, -name,  -one_of(variable.list)) %>%
      mutate_at("value", as.numeric) -> data.scaled.l
  }
  
  model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
    # K-fold validation for each model based on a single term 
    for(n in 1:folds.length) { 
      test.ids = folds$subsets[folds$which==n]
      test.trials = unique(data.scaled.l$name)[test.ids]
      testdata = x[which(x$name %in% test.trials), ]
      traindata = x[-which(x$name %in% test.trials),]
      GLMER <- lme4::glmer(formula, data = traindata, family="binomial", 
                           control = glmerControl(optimizer="bobyqa"),nAGQ=0)
      data.summary.mat[n,] = c(glance(GLMER)[3], tidy(GLMER)[2,5])
      glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
      glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
      glmerperf[n] <- glmer.ROC$auc
      if (length(variable.list) > 0) { # Compare new model against the null 
        model.null = lme4::glmer(formula.null, data = traindata, family="binomial", 
                                 control = glmerControl(optimizer="bobyqa"),nAGQ=0)
        anova.test = anova(model.null, GLMER)
        anova.results[n] = anova.test$`Pr(>Chisq)`[2]
      }
    }
    anova.mean = mean(anova.results)
    performance=mean(glmerperf)
    min.perf = min(glmerperf)
    max.perf = max(glmerperf)
    data.summary = colMeans(as.data.frame(data.summary.mat))
    return(t(c(per = performance, min.perf = min.perf, max.perf = max.perf, data.summary, anova = anova.mean)))# put the se of the performance 
  })
  variable = model.results %>% arrange(AIC) %>% slice(1)  # Select the variable that is the best 
  if(length(variable.list) == 0 ) { # Is this the first? If yes, add this variable
    variable.list = c(variable.list, variable$variable)
    variable.per = c(variable.per, variable$per)
    variable.aic = c(variable.aic, variable$AIC)
    variable.max.per = c(variable.max.per, variable$max.perf)
    variable.min.per = c(variable.min.per, variable$min.perf)
    print(paste0("Adding : ", variable$variable))
  } else if (variable$anova < .05) { # Determine if this new variable beats the null
    variable.list = c(variable.list, variable$variable)
    variable.per = c(variable.per, variable$per)
    variable.max.per = c(variable.max.per, variable$max.perf)
    variable.min.per = c(variable.min.per, variable$min.perf)
    variable.aic = c(variable.aic, variable$AIC)
    
    print(paste0("Adding : ", variable$variable))
  } else { 
    print(paste0("Final Single Term Model is: ", formula.null))
    build = FALSE
  }
}

