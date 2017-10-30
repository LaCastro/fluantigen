excluded = c("postAntigen", "oriAntigen",  "day.1", "date", "totalN", "totalR",
             "totalCases", "dominant.type", "totalI", "max.I", "min.I")

data = freq.two[, -which(predictorNames%in%excluded)]


data <- within(data, {
  success <- factor(success, levels=c("yes", "no"), labels = c(1,0))
  .id <- factor(.id)
  day <- as.numeric(as.character(day))
    mutLoad <- as.numeric(as.character(mutLoad))
    distance <- as.numeric(as.character(distance))
  simDay <- as.numeric(as.character(simDay))
})

data %>%
  mutate(time.of.year = simDay-floor(x = day/365)) %>%
  mutate(quarter = ifelse(time.of.year < .25, 1,
                          ifelse(time.of.year < .5, 2,
                                 ifelse(time.of.year < .75, 3,4)))) -> dataRaw

dataRaw$quarter = as.factor(dataRaw$quarter)

# Take our any data that's before 3 years
dataRaw %>%
  dplyr::filter(day > (3*365)) -> dataRaw

dataScaled <- dataRaw

factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
dataScaled$netau[!is.finite(dataScaled$netau)] <- NA

dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)



variable.list = c("distance", "mutLoad", "normalize.I", "meanLoad")
p = vector(mode = "list", length = length(variable.list))

for (v in 1:length(variable.list)) {
  variable = variable.list[v]
  dataScaled %>%
    select(success, eval(variable)) -> dataScaled.select
  colnames(dataScaled.select)[2] ="scaled"
  dataScaled.select$scaled = as.numeric(dataScaled.select$scaled)
  
  dataRaw %>%
    select(eval(variable)) -> dataRaw.select
  colnames(dataRaw.select)[1] ="raw"
  
  dataSelect = cbind(dataScaled.select, dataRaw.select)
  unsucessful = which(dataSelect$success == 0)
  sample.no = sample(unsucessful, size = 100, replace = FALSE)
  
  
  plot <- ggplot(dataSelect, aes(x = raw, y = scaled, color = success)) + 
    geom_jitter(height = .15, alpha = .5) +
  labs(title = eval(variable)) + theme(axis.text.x = element_text(size = 8)) + facet_wrap(~success)
  p[[v]] <- plot
}

p

plot.freq.emerge = cowplot::plot_grid(plotlist = p)
cowplot::save_plot(plot.freq.emerge, filename = "exploratory.figs/data.reg.emerge.pdf", base_aspect_ratio = 1.8)

p[[3]]
p[[2]]




########################
library(sjPlot)
library(sjmisc)

set_theme(theme = "forest",
          geom.label.size = 3,
          axis.textsize = .9,
          axis.title.size = .9)

### Random Effects 
sjp.glmer(model.glm, y.offset = .4, sort.est = "(Intercept)")
sjp.glmer(GLMER, y.offset = .4, sort.est = "(Intercept)")


### Fixed Effects 
sjp.glmer(GLMER, type = "fe")
sjp.glmer(model.glm, type = "fe")

sjp.glmer(GLMER, type = "eff", show.ci = TRUE)


#########

freq.explore = rbind(data.frame(freq = ".01", freq.one),
                     data.frame(freq = ".02", freq.two),
                     data.frame(freq = ".05", freq.five))


freq.explore %>%
  mutate(ratio.mutation = individual.meanMut/meanLoad,
         ratio.meanR = individual.meanR/meanR,
         ratio.varR = individual.varR/varR,
         ratio.meanBeta = individual.meanBeta/meanBeta,
         ratio.varBeta = individual.varBeta/varBeta,
         ratio.meanSigma = individual.meanSigma/meanSigma,
         ratio.varSigma = individual.varSigma/varSigma) -> freq.explore

freq.explore %>%
  filter(day != 1833) %>%
  filter(freq == ".01") -> freq.first.check

# calculating the ratios 
individual = c("individual.meanMut", "individual.varMut", "individual.meanR", 
               "individual.varR", "individual.meanBeta", "individual.varBeta", 
               "individual.meanSigma","individual.varSigma")
popAntigen = c("antigenicTypes", "dominant.freq", "diversity", "tmrca", "antigenicDiversity")
viralFitness = c("meanLoad", "meanR", "varR", "meanBeta", "varBeta", "meanSigma","varSigma", "covBetaSigma")
popDynamics = c("netau", "serialInterval", "totalI")
ratios = c("ratio.mutation", 'ratio.meanR', 'ratio.varR', 'ratio.meanBeta', 'ratio.varBeta', 'ratio.meanSigma', 'ratio.varSigma')

freq.explore %>%
  gather(key = variable, value = value, -freq, - antigentype, -success, -name,-day) %>%
  filter(variable %in% ratios) %>%
  filter(day != 1833) %>%
  ggplot(aes(x = freq, y = value, color = success)) +
  facet_wrap(~variable, scales = "free") + 
  geom_boxplot() +
  scale_color_manual(values = c("orange", "purple")) + 
  labs(x = "Frequency Point", color  = "Antigen Fate") -> ratios.plot



save_plot(ratios.plot, filename = "exploratory.figs/ratios.pdf", base_height = 8, base_aspect_ratio = 1.5)

freq.explore %>%
  ggplot(aes(x = individual.varSigma , y = varSigma, color = success)) +
  facet_wrap(~freq) + geom_point() + 
  scale_color_manual(values  = c("orange", "purple")) + geom_smooth()

freq.explore %>%
  group_by(antigentype, name) %>%
  filter(day != 1833) %>%
  gather(key = variable, value = value, -freq, -antigentype,-success,-name) %>%
  arrange(antigentype) %>%
  group_by(name, antigentype, variable) %>%
  summarise(first.growth = first_growth(value),
            second.growth = second_growth(value),
            success = success[1]) %>%
  ungroup() %>%
  mutate_at("name", as.factor) -> diff.df


diff.df %>%
  select(-second.growth) %>%
  spread(key = variable, first.growth) -> diff.first.growth

individual = c("individual.meanMut", "individual.varMut", "individual.meanR", 
               "individual.varR", "individual.meanBeta", "individual.varBeta", 
               "individual.meanSigma","individual.varSigma")
popAntigen = c("antigenicTypes", "dominant.freq", "diversity", "tmrca", "antigenicDiversity")
viralFitness = c("meanLoad", "meanR", "varR", "meanBeta", "varBeta", "meanSigma","varSigma", "covBetaSigma")
popDynamics = c("netau", "serialInterval", "totalI")


diff.df %>%
  filter(variable %in% ratios) %>%
  gather(phase, value, first.growth, second.growth) %>%
  ggplot(aes(phase, value, color = success)) + geom_boxplot() + facet_wrap(~variable, scales = "free") + 
  scale_color_manual(values = c("orange", "purple")) -> ratios.diff.plot

save_plot(ratios.diff.plot, filename = "exploratory.figs/ratios.diff.plot.pdf",
          base_height = 8, base_aspect_ratio = 1.8)




diff.df %>%

  
first_growth = function(x) {x[2]-x[1]}
second_growth = function(x){x[3]-x[2]}