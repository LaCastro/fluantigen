freq.1.subset %>% select(ratio.varSigma, infected, diversity, antigenicDiversity,
                         success,antigentype, name) -> dataPurged
freq.2.subset %>% select(ratio.varR, individual.varSigma,
                         success, antigentype,name) -> dataPurged
freq.3.subset %>% select(ratio.meanR, varR, ratio.mutation,
                         success, antigentype, name) -> dataPurged


growth.1.subset %>% select(ratio.meanR, ratio.mutation,
                           success, antigentype, name, totalI) -> dataPurged
growth.2.subset %>% select(ratio.meanR, varR, entropy, ratio.varR, day,
                           success, antigentype, name, totalI) -> dataPurged
accelerate %>% select(ratio.mutation, day, success, antigentype, name) -> dataPurged

dataPurged  <- within(dataPurged, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
  #freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})

factor.variables = which(sapply(dataPurged,is.factor)==TRUE)
colnames(dataPurged)[-factor.variables] = paste0("accel.", colnames(dataPurged)[-factor.variables])

freq.1 = dataPurged
freq.2 = dataPurged
freq.3 = dataPurged

diff.1 = dataPurged
diff.2 = dataPurged
accel = dataPurged

freq.1 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.1.l
freq.3 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.3.l
freq.2 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.2.l
diff.2 %>% gather(key = variable, value = value, -success, -name, -antigentype) -> diff.2.l
diff.1 %>% gather(key = variable, value = value, -success, -name, -antigentype)  -> diff.1.l
accel %>% gather(key = variable, value = value, -success, -name, -antigentype) -> accel.l

the.whole.data = bind_rows(freq.1.l, freq.2.l, freq.3.l, diff.2.l, diff.1.l, accel.l) %>%
  spread(key = variable, value = value) 
combined.data = the.whole.data


head(combined.data)
##########################################################################################
set.seed(62417)
folds.length = 5;
folds = cvFolds(n = length(unique(combined.data$name)), K = folds.length, type = "random")

glmerperf=rep(0, folds.length); glmperf=glmerperf;
data.summary.mat = data.frame(matrix(nrow = folds.length, ncol  = 2))
colnames(data.summary.mat) = c( "AIC","p.value")

#Initializing Variables 
build = TRUE
variable.list = c()
variable.per = c()
variable.aic = c()
anova.results = rep(0, folds.length)

combined.data$success <- relevel(combined.data$success, ref = "Transient")
combined.data = combined.data %>%  mutate_if(is.numeric, scale)


formula = "success ~ accel.day + accel.ratio.mutation + diff.1.ratio.meanR + diff.1.ratio.mutation + 
diff.2.day + diff.2.entropy + diff.2.ratio.meanR + diff.2.ratio.varR + diff.2.varR + freq.1.antigenicDiversity + 
freq.1.diversity + freq.1.infected + freq.1.ratio.varSigma + freq.2.individual.varSigma + freq.2.ratio.varR +
freq.3.ratio.meanR + freq.3.ratio.mutation + freq.3.varR + (1 | name)"
  
  
# K-fold validation for each model based on a single term 
for(n in 1:folds.length) { 
      test.ids = folds$subsets[folds$which==n]
      test.trials = unique(combined.data$name)[test.ids]
      testdata = combined.data[which(combined.data$name %in% test.trials), ]
      traindata = combined.data[-which(combined.data$name %in% test.trials),]
      GLMER <- lme4::glmer(formula, data = traindata, family="binomial", 
                           control = glmerControl(optimizer="bobyqa"),nAGQ=0)
      data.summary.mat[n,] = c(glance(GLMER)[3], tidy(GLMER)[2,5])
      glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
      #glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
      #glmerperf[n] <- glmer.ROC$auc
    
      ## How to check 
      forensic = data.frame(cbind(probs = glmer.probs, actual = testdata$success))
      forensic = cbind(forensic, name = testdata$name, antigentype=testdata$antigentype)

      forensic %>%
        mutate(characterization = ifelse(probs < .5 & actual == "2", "wrong", "right")) -> forensic 

      forensic.test = left_join(forensic, testdata)
      forensic.test = rbind(forensic.test, forensic.test)
}
      

    
forensic.test %>%
  select(name, antigentype, characterization) %>%
  left_join(combined.data) %>% 
# ggplot(aes( x = characterization, y = diff.2.totalI, color = characterization)) + geom_boxplot() 
  gather(key = variable, value = value, -name, -antigentype, -success,-characterization) %>%
  ggplot(aes(x = characterization, y = value, color = characterization)) +
  geom_boxplot(alpha = .8) + 
  scale_fill_manual(values = c("black", "grey"), labels = c("True", "False")) + 
  facet_wrap(~variable, scales = "free") + 
  scale_color_manual(values = c("black", "grey"), labels = c("True", "False")) +
  theme(strip.text = element_text(size = 8)) +
  labs(y = "Density", x = "Value", fill = "Positive") + guides(color = FALSE) 

-> forensic.unscaled

save_plot(forensic.unscaled, filename = "exploratory.figs/forensic.unscaled.pdf",
          base_height = 8, base_aspect_ratio = 2)
#######################

forensic.test %>%
  filter(characterization == "wrong")

forensic.test %>% 
  gather(key = variable, value = value, -probs, -actual, -name, -antigentype, -success,-characterization) %>%
  ggplot(aes(value, color = characterization, fill = characterization)) +
  geom_density(alpha = .8) + 
  scale_fill_manual(values = c("black", "grey"), labels = c("True", "False")) + 
  facet_wrap(~variable, scales = "free") + 
  scale_color_manual(values = c("black", "grey")) +
  theme(strip.text = element_text(size = 8)) +
  labs(y = "Density", x = "Value", fill = "Positive") + guides(color = FALSE) -> forensic.plot
save_plot(filename = "exploratory.figs/forensic.plot.pdf", forensic.plot, base_height = 8, base_aspect_ratio = 2)

library(lattice)
forensic.test$characterization = as.factor(forensic.test$characterization)





forensic.test %>%
  select(-probs, -actual, -name, -antigentype,-success) %>%
  gather(key = second.variable, value = value, -characterization, -freq.3.varR) %>%
  mutate_at("second.variable", as.factor) %>%
  ggplot(aes(x =  freq.3.varR, y = value, color = characterization)) + 
  geom_point(alpha = .5) + facet_wrap(~second.variable, scales = "free") +
  scale_color_manual(values = c("orange", "black"), labels = c("True", "False")) + 
  labs(color = "Positive") ->  freq.3.varR.plot



qqplot(forensic.test$diff.1.totalI[forensic.test$characterization == "right"],
       forensic.test$diff.1.totalI[forensic.test$characterization == "wrong"])
abline(a = 0, b = 1, lty = 3)

ks.test(forensic.test$diff.2.totalI[forensic.test$characterization == "right"],
        forensic.test$diff.2.totalI[forensic.test$characterization == "wrong"])


shapiro.test(forensic.test$diff.1.totalI[forensic.test$characterization == "wrong"])

wilcox.test(diff.1.totalI ~ characterization, data = forensic.test)
t.test(diff.2.totalI ~ characterization, data = forensic.test)

m <- table(forensic.test$characterization, forensic.test$diff.2.totalI)
chisq.test(m)

