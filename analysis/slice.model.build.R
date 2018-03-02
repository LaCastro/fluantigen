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
#dataPurged %<>% select(-group)

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

<<<<<<< Updated upstream
single.term.results = bind_cols(variable = variable.list, per =  variable.per, max.per = variable.max.per, min.per = variable.min.per, aic = variable.aic)
write.csv(single.term.results, '../results/022318.slicediff.csv', row.names = FALSE)
=======

single.term.results = bind_cols(variable = variable.list, per =  variable.per, max.per = variable.max.per, min.per = variable.min.per, aic = variable.aic)
write.csv(single.term.results, '../results/022218.slice.05.csv', row.names = FALSE)
>>>>>>> Stashed changes

roc = NULL
for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(data.scaled.l$name)[test.ids]
  testdata = dataPurged[which(dataPurged$name %in% test.trials), ]
  traindata = dataPurged[-which(dataPurged$name %in% test.trials),]
  GLMER <- lme4::glmer(formula.null, data = traindata, family="binomial", 
                       control = glmerControl(optimizer="bobyqa"),nAGQ=0)
  data.summary.mat[n,] = c(glance(GLMER)[3], tidy(GLMER)[2,5])
  glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
  glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
  
  roc.values = data.frame(cbind(sen = glmer.ROC$sensitivities,
                                spec = glmer.ROC$specificities,
                                thres = glmer.ROC$thresholds))
  roc.values %<>% mutate(fold = n)
  roc = rbind(roc, roc.values)
}
<<<<<<< Updated upstream
write.csv(roc, "../results/022318.slicediff.roc.csv", row.names = FALSE)
=======
write.csv(roc, "../results/022218.slice.roc.05.csv", row.names = FALSE)




############# Distribution of Time Advance

### Looking at the frequency (randomly, when do I get it)
data %>% 
  filter(variable == "frequency") %>%
  group_by(success) %>%
  summarize(min.iqr = quantile(value,probs = .25),
            med.iqr = median(value),
            max.iqr = quantile(value, probs = .75)) -> data.summary


data %>% 
  filter(variable == "frequency") %>%
  ggplot(aes(value,  fill = success)) + geom_density(alpha = .8) + 
  labs(x = "Relative frequency when sampled", y = "Density", fill = "")  +
  scale_fill_manual(values = c("black", "grey")) -> density.frequency.full
  
data %>% 
  filter(variable == "frequency") %>%
  ggplot(aes(value,  fill = success)) + geom_histogram() + 
  labs(x = "Relative frequency when sampled", y = "Counts", fill = "")  +
  scale_fill_manual(values = c("black", "grey")) -> histo.frequency.full

>>>>>>> Stashed changes
