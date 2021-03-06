################### Doing K-cross validation 
set.seed(62417)
folds.length = 5;
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

single.term.results = bind_cols(variable = variable.list, per =  variable.per, max.per = variable.max.per, min.per = variable.min.per, aic = variable.aic)
write.csv(single.term.results, '../results/022018.3point.10.csv', row.names = FALSE)

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
write.csv(roc, "../results/022018.3point.roc.10.csv", row.names = FALSE)


############# Distribution of Time Advance
subset.analyze %>% 
  select(day.success, success, antigentype, name) -> subset.calculate.day 

freq.t.three %>%
  select(success, antigentype, name, day) %>%
  left_join(subset.calculate.day) %>% 
  filter(success == "Est." & day != 1833) %>%
  mutate(lead.time = day.success - day)  -> lead.data

lead.data %>% 
  ggplot(aes(lead.time)) + geom_histogram(bins = 30)
summary(lead.data$lead.time)

write.csv(lead.data, file= "../results/lead.time/022018.3point.lead.10.csv", row.names = FALSE)


############## STEP TWO: Interactions 
#base.formula = paste0(paste(variable.list, collapse = "+"), "+(1|name)")
# Generates possible combinations 
#combn.df = data.frame(combn(variable.list, m=2))
#comb.list = unname(apply(combn.df, MARGIN = 2, function(x) paste(x, collapse = "*")))

dataPurged$success <- relevel(dataPurged$success, ref = "Transient")
dataPurged %>%
  ungroup() %>%
  dplyr::select(success, name, one_of(variable.list)) -> data.scaled.l

# Initialize Variables 
build = TRUE
anova.results = rep(0, folds.length)
term.placement = length(variable.list) + 2
int.list = c()
int.aic = c()
int.per = c()

while(build == TRUE) {
  contenders = data.frame(per = double(), AIC = double(), p.value = double(), anova = double(), combo = character())
  if(length(int.list) > 0) { 
  comb.list = comb.list[-which(comb.list %in% int.list)] # Update which ones to try 
  }
  
  formula.null = paste0("success~", paste(int.list, collapse = "+"), "+", base.formula)
  
  for(combo in comb.list) {    
    test.formula = paste0("success~", combo,"+",  paste(int.list, collapse = "+"), "+", base.formula)
    for(n in 1:folds.length) { 
      test.ids = folds$subsets[folds$which==n]
      test.trials = unique(data.scaled.l$name)[test.ids]
      testdata = data.scaled.l[which(data.scaled.l$name %in% test.trials), ]
      traindata = data.scaled.l[-which(data.scaled.l$name %in% test.trials),]
      
      model.null <- lme4::glmer(formula.null,data = traindata, family="binomial", 
                                control = glmerControl(optimizer="bobyqa"),nAGQ=0)
      model.int <- lme4::glmer(test.formula, data = traindata, family="binomial", 
                               control = glmerControl(optimizer="bobyqa"),nAGQ=0)
      model.test = anova(model.null, model.int)
      anova.results[n] = model.test$`Pr(>Chisq)`[2]
      data.summary.mat[n,] = c(glance(model.int)[3], tidy(model.int)[term.placement,5])
      glmer.probs <- predict(model.int, newdata=testdata, type="response", allow.new.levels=TRUE)
      glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
      glmerperf[n] <- glmer.ROC$auc
    } 
    if(mean(anova.results) < .05) { # If the variable is significant-keep track  
      info = data.frame(t(c(per = mean(glmerperf), colMeans(as.data.frame(data.summary.mat)), anova = mean(anova.results))))
      info$combo = combo
      contenders = rbind(contenders, info)
    } 
  }
  
  if(nrow(contenders) == 0) { # If none of the variables were 
    print(paste0("Final Model is:", formula.null))
    coefficient.estimates = tidy(model.null)
    build = FALSE
  } else {
    variable = contenders %>% arrange(AIC) %>% slice(1)
    int.list = c(int.list, variable$combo)
    int.per = c(int.per, variable$per)
    int.aic = c(int.aic, variable$AIC)
    print(paste0("Selecting - ", variable$combo))
  }
}
interaction.term.results = bind_cols(variable = int.list, per =  int.per, aic = int.aic)
full.model = rbind(single.term.results , interaction.term.results)
full.model.results = full_join(full.model, coefficient.estimates, by = c("variable" = "term"))

write.csv(full.model.results, '../results/transient.03poplevel.csv')




################################ Looking at ROC 
calculate_SenSpec <- function(model.probs, testdata) {
  performance.data = data.frame(cbind(glmer.probs, truth=testdata$success))
  performance.data %>%
    mutate(predicted = ifelse(glmer.probs < .5, 1,2)) -> performance.data
  data.table = table(Truth = performance.data$truth, Prediction = performance.data$predicted)
  tpr =  sensitivity(data.table, positive = "2")
  tnr =  specificity(data.table, negative = "1")
  return(c(sen = tpr, spec = tnr))
}
calculate_SenSpec(glmer.probs, testdata)


#####################
redo.pop = data.frame(cbind(sen = glmer.ROC$sensitivities,
                              spec = glmer.ROC$specificities,
                              thres = glmer.ROC$thresholds))
write.csv(redo.pop, "~/Dropbox/current_fluantigen/model_csvs/realworld3.roc.csv")

redo.pop %>%
  mutate(false.pos = 1-spec) %>%
  ggplot(aes(x = false.pos, y = sen)) + geom_line(size = 2) +
  labs(x = "False Positive Rate", y = "True Positive Rate") -> new.perf
save_plot(new.perf, filename = "exploratory.figs/new.perf.pdf")


######## New performance
single.term.results %>%
  mutate(id = row_number()) %>%
  ggplot(aes(x = id, y = per)) + geom_line(size =2 ) + 
  labs(x = "# of Model Terms", y = "Average AUC") +
  scale_y_continuous(breaks = seq(from = .55, to =  .95, by = .05), limits = c(.55, .95)) +
  scale_x_continuous(breaks = seq(from = 1, to = 16, by = 2)) -> new.auc

save_plot(new.auc, filename = "exploratory.figs/new.auc.pdf")


trans.pop3 %>%
  ggplot(aes(x  = (1-spec), y = sen)) + geom_line()

roc.data = data.frame(cbind(thres = glmer.ROC$thresholds,
                            sen = glmer.ROC$sensitivities,
                            spec = glmer.ROC$specificities))

### Need to get estimates 
plot_pred_type_distribution <- function(df, threshold) {
  v <- rep(NA, nrow(df))
  v <- ifelse(df$pred >= threshold & df$survived == 1, "TP", v)
  v <- ifelse(df$pred >= threshold & df$survived == 0, "FP", v)
  v <- ifelse(df$pred < threshold & df$survived == 1, "FN", v)
  v <- ifelse(df$pred < threshold & df$survived == 0, "TN", v)
  
  df$pred_type <- v
  
  ggplot(data=df, aes(x=survived, y=pred)) + 
    geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
    geom_jitter(aes(color=pred_type), alpha=0.6) +
    geom_hline(yintercept=threshold, color="red", alpha=0.6) +
    scale_color_discrete(name = "type") +
    labs(title=sprintf("Threshold at %.2f", threshold))
}
set_test_data <- function(data, fold) {
  test.ids = fold$row.index
  testnames = unique(data$name)[test.ids]
  testdata = data[which(data$name %in% testnames),]
  return(testdata)
}
set_train_data <- function(data, fold) {
  test.ids = fold$row.index
  testnames = unique(data$name)[test.ids]
  traindata = data[-which(data$name %in% testnames),]
  return(traindata)
}
