trials = c("tropics_100", "tropics_20")

freq.02 %>%
  filter(.id %in% trials) -> freq.2.sub

freq.05 %>%
  filter(.id %in% trials) -> freq.5.sub

freq.10 %>%
  filter(.id %in% trials) -> freq.10.sub

freq.15 %>%
  filter(.id %in% trials ) -> freq.15.sub

freq.sub = rbind(data.frame(frequency = "2",  freq.2.sub),
                 data.frame(frequency = "5", freq.5.sub),
                 data.frame(frequency = "10", freq.10.sub), 
                 data.frame(frequency = "15", freq.15.sub))

antigen.data %>%
  filter(.id %in% trials) -> antigen.sub


###### Now get the antigen data so can plot 
success.types = return_success_types(antigen.data)
antigen.freq.success.df = filter_frequencies_success(success.types = success.types, 
                                                     antigen.frequencies = antigen.freq.sub)
success.types %>%
  filter(.id %in% trials) -> success.sub

antigen.frequencies %>%
  filter(.id %in% trials) -> antigen.freq.sub


antigen.freq.success.df %>%
  filter(day < 7300) -> antigen.freq.sucess.df

antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = n_distinct(antigentype)) -> num.transitions

antigen.freq.success.df = ddply(antigen.freq.success.df,.variables = ".id", function(antigen) {
  trial = antigen$.id[1]
  num.distinct.cluster = num.transitions[which(num.transitions$.id == trial), "num.transitions"]
  unique.antigens = unique(antigen$antigentype)
  cluster.pair = data.frame(cbind(unique.antigens, seq(1:num.distinct.cluster$num.transitions)))
  left_join(x = antigen, y = cluster.pair, by = c("antigentype" = "unique.antigens")) -> antigen
  colnames(antigen)[6] = "cluster.number"
  return(antigen)
})

antigen.freq.success.df %>%
  select(.id, antigentype, cluster.number) -> cluster.key

freq.sub$.id = as.factor(freq.sub$.id)
freq.sub$postAntigen = as.factor(freq.sub$postAntigen)

cluster.key 
  
indices = unique(cluster.key[,c('.id', 'antigentype')])
str(cluster.key)
cluster.key %>%
  distinct(.id, antigentype, .keep_all = TRUE) -> cluster.key

cluster.key$.id = as.factor(cluster.key$.id)
cluster.key$antigentype=as.factor(cluster.key$antigentype)

freq.sub %>%
  inner_join(cluster.key, by = c("postAntigen" = "antigentype", ".id" = ".id")) -> freq.sub



max.color = max(num.transitions$num.transitions)
myColors = set_my_colors(max.color)

## Need to first go in and fill all the missing values and then combine
ant.freq.success.l = ddply(.data = antigen.freq.success.df, .variables = ".id", function(sim) {
  sim %>%
    distinct(day, cluster.number, .keep_all = TRUE)%>%
    dplyr::select(-antigentype) -> part1
  part1 %>%
    spread(key = cluster.number, value = frequency, fill = 0) -> step2 
  step2 %>%
    gather(key = cluster.number, value = frequency, -1, -2, - infected) -> antigen.freq.long
  return(antigen.freq.long)
})

ant.freq.success.l$cluster.number = as.factor(ant.freq.success.l$cluster.number)


freq.sub %>%
  filter(cluster.number == 4) %>%
  select(frequency, postAntigen, day, .id) -> cluster.four
cluster.four$day=as.numeric(cluster.four$day)

ant.freq.success.l %>%
  mutate(year = day/365) %>%
  #mutate(prevalence = frequency*infected) %>%
  mutate(prevalence = frequency) %>%
  filter(prevalence > 0) %>%
  filter(.id %in% trials ) %>%
  ggplot(aes(x = year, y = prevalence, fill = cluster.number)) +
  geom_area(color = "black", aes(color = antigentype, fill = cluster.number)) +
  #geom_line(aes(x = year, y = infected), color = "black") + 
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") +
  coord_cartesian(xlim = c(0,10)) + 
  #scale_x_continuous(breaks = seq(1:20)) +
  guides(col = FALSE) + guides(fill = FALSE) +
  geom_vline(data = cluster.four, aes(xintercept = day/365, linetype = frequency), size = 1.5)
#-> prev.plot

ant.freq.success.l %>%
  mutate(year = day/365) %>%
  mutate(prevalence = frequency*infected) %>%
  #mutate(prevalence = frequency) %>%
  filter(prevalence > 0) %>%
  filter(.id %in% trials ) %>%
  ggplot(aes(x = year, y = prevalence, fill = cluster.number)) +
  geom_area(color = "black", aes(color = antigentype, fill = cluster.number)) +
  #geom_line(aes(x = year, y = infected), color = "black") + 
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") +
  coord_cartesian(xlim = c(0,10)) + 
  #scale_x_continuous(breaks = seq(1:20)) +
  guides(col = FALSE) + guides(fill = FALSE) +
  geom_vline(data = cluster.four, aes(xintercept = day/365, linetype = frequency), size = 1.5)

prev.plot + ggplot(cluster.two, aes(day)) + geom_vline()
