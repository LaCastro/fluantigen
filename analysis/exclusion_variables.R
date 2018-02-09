## Exclusion 


# What I believe was taken out at this point 
# vif.excluded = c("freq.2.diversity", "freq.2.antigenicDiversity", "freq.2.entropy",
#                 "freq.3.diversity", "freq.3.antigenicDiversity", "freq.2.dominant.freq",
#                 "freq.2.infected", "freq.1.entropy", "freq2.meanR", "freq.2.meanR",
#                 "freq.2.totalS", "freq.3.dominant.freq", "freq.2.tmrca")


###### 1/22 ################################################### 
# Snapshot 
vif.excluded = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR", "individual.meanSigma",
                 "ratio.meanSigma", "meanSigma", "covBetaSigma", "meanBeta", "individual.varR", "individual.varBeta", "ratio.varR")

# Growth 
vif.excluded = c("individual.meanMut", "individual.varBeta", "individual.varR", "ratio.meanR", "meanR", "ratio.meanBeta",
                 "individual.meanR", "ratio.meanSigma", "varSigma", "individual.varMut")

# Combined 
vif.excluded = c("freq.2.ratio.meanBeta", "freq.1.ratio.meanBeta", "freq.2.meanR", "diff.1.meanSigma",
                 "freq.1.meanR", "diff.1.individual.meanR", "diff.1.individual.meanSigma", "diff.2.individual.meanSigma",
                 "freq.2.varSigma", "freq.2.ratio.meanR", "freq.3.varSigma", "diff.1.ratio.varR", "accel.ratio.mutation",
                 "accel.ratio.varR", "diff.1.individual.meanBeta", "accel.ratio.varBeta", "accel.individual.meanBeta",
                 "freq.1.ratio.mutation", "freq.2.ratio.mutation", "accel.meanLoad", "freq.3.ratio.mutation")



##### 1/24 ######################################################

## Exclusion based on collinearity after everything has been combined and up until VIF < 5 
vif.excluded = c("freq.1.individual.meanMut", "freq.1.individual.meanBeta", "freq.2.individual.meanMut", "freq.3.individual.meanMut",
                 "diff.1.individual.meanMut", "diff.1.individual.meanBeta", "diff.2.individual.meanMut", "accel.individual.meanMut",
                 "diff.1.individual.varMut", "diff.1.individual.varBeta", "diff.1.individual.varR", "accel.individual.varMut",
                 "accel.individual.varBeta", "accel.individual.varR", "freq.2.individual.varMut", "freq.3.individual.varMut", 
                 "freq.1.individual.varBeta", "freq.1.individual.varR", "diff.2.individual.varMut", "diff.1.ratio.meanR",
                 "accel.ratio.meanR", "freq.2.individual.meanBeta", "freq.3.individual.meanBeta", "diff.1.individual.meanR",
                 "diff.1.ratio.meanSigma", "accel.ratio.meanSigma", "freq.2.ratio.meanSigma", "freq.3.ratio.meanSigma",
                 "diff.1.meanR", "diff.2.individual.meanR", "freq.2.individual.meanR", "freq.1.individual.meanR", "accel.ratio.meanBeta",
                 "accel.meanR", "freq.2.meanR", "freq.3.meanR", "freq.1.individual.meanSigma", "freq.2.individual.meanSigma", 
                 "freq.3.individual.meanSigma", "diff.2.ratio.meanR", "freq.2.ratio.meanR", "accel.individual.meanBeta", "accel.meanSigma",
                 "freq.3.individual.meanR", "diff.2.ratio.meanBeta", "freq.1.ratio.meanBeta", "freq.3.ratio.meanBeta", "accel.individual.meanSigma",
                 "freq.1.ratio.meanR", "diff.2.meanR", "diff.2.ratio.meanSigma", "diff.1.varSigma", "freq.1.meanSigma", "freq.2.meanSigma",
                 "freq.3.meanSigma", "freq.1.ratio.meanSigma", "accel.varSigma", "freq.2.varSigma", "freq.3.varSigma", "accel.individual.meanR",
                 "accel.varR", "freq.1.covBetaSigma", "freq.2.covBetaSigma", "freq.3.covBetaSigma", "diff.2.varSigma", "freq.1.meanBeta",
                 "freq.2.meanBeta", "freq.3.meanBeta", "accel.ratio.varBeta", "freq.2.individual.varBeta", "freq.2.individual.varR",
                 "freq.3.individual.varBeta", "freq.3.individual.varR", "freq.1.individual.varMut", "diff.1.ratio.varR", "diff.2.individual.varBeta",
                 "diff.2.individual.varR", "accel.ratio.varR", "accel.ratio.mutation", "freq.2.ratio.varR", "freq.3.ratio.varR",
                 "freq.1.ratio.varR", 'diff.1.meanLoad', "freq.1.ratio.mutation", "freq.2.ratio.mutation", "freq.2.ratio.meanBeta",
                 "diff.2.ratio.mutation", "diff.2.ratio.varR", "diff.1.ratio.varSigma", "diff.1.ratio.mutation", "freq.1.entropy",
                 "freq.2.entropy", "freq.3.entropy", "accel.diversity", "accel.antigenicDiversity", "freq.1.varSigma", "accel.ratio.varSigma",
                 "diff.1.entropy", "diff.1.dominant.freq", "diff.1.antigenicTypes")

## Exclusion based on collinearity after everything has been combined and up until VIF < 5 and
# incorporating the metric of I/N

################### 1/29 -- Realized was pulling in wrong data 
vif.excluded = c("freq.2.individual.meanMut", "freq.2.individual.meanBeta", "freq.1.individual.meanMut", "freq.3.individual.meanMut",
                 "freq.1.individual.meanBeta", "gp.1.individual.meanBeta", "gp.1.individual.meanMut", "gp.2.individual.meanBeta",
                 "accel.individual.meanMut", "freq.2.individual.varMut", "freq.2.individual.varBeta", "freq.2.individual.varR", "freq.3.individual.varMut",
                 "gp.1.individual.varMut", "accel.individual.varMut", "freq.2.meanBeta", "freq.1.individual.varBeta", "freq.1.individual.varR",
                 "gp.1.individual.varBeta", "gp.1.individual.varR","gp.2.individual.varBeta", "gp.2.individual.varR", "freq.2.ratio.meanBeta",
                 "freq.2.meanSigma", "freq.3.meanBeta", "freq.2.meanR", "freq.2.individual.meanR", "gp.2.individual.varMut", "freq.2.ratio.meanSigma",
                 "freq.2.individual.meanSigma", "freq.2.ratio.meanR", "freq.3.ratio.meanBeta", "freq.3.individual.meanBeta", "freq.3.meanSigma",
                 "freq.1.individual.meanR", "gp.1.individual.meanR","gp.2.individual.meanR", "freq.3.ratio.meanSigma", "gp.1.ratio.meanSigma",
                 "accel.ratio.meanSigma", "freq.3.meanR", "gp.1.meanR", "accel.ratio.meanBeta", "freq.1.individual.meanSigma", "accel.meanR",
                 "accel.individual.meanBeta", "freq.3.individual.meanSigma","accel.meanSigma", "accel.ratio.meanR", "freq.3.ratio.meanR", "gp.2.ratio.meanR",
                 "freq.2.covBetaSigma", "freq.2.varSigma", "freq.3.individual.meanR", "gp.2.ratio.meanBeta", "accel.individual.meanSigma", "freq.2.individual.varSigma",
                 "accel.individual.meanR", "freq.1.ratio.meanSigma", "gp.2.meanR", "freq.2.varR", "freq.1.covBetaSigma", "freq.1.meanSigma", 
                 "gp.2.ratio.meanSigma", "freq.2.meanLoad", "freq.1.varSigma", "gp.1.varSigma", "gp.2.varSigma", "gp.1.ratio.meanR", "freq.3.covBetaSigma",
                 "freq.3.individual.varSigma", "accel.varSigma","freq.2.diversity", "freq.2.antigenicDiversity", "freq.1.meanBeta", "freq.2.ratio.mutation",
                 "freq.3.meanLoad", "freq.2.ratio.varBeta", "freq.2.entropy", "freq.2.ratio.varR", "freq.3.diversity", "freq.3.antigenicDiversity",
                 "freq.2.antigenicTypes", "freq.2.dominant.freq", "freq.2.varBeta","freq.2.prop.I", "freq.2.totalI", "freq.2.totalS", "freq.1.individual.varMut",
                 "freq.1.ratio.mutation","accel.individual.varBeta", "accel.individual.varR", "freq.3.individual.varBeta", "freq.3.individual.varR",
                 "freq.3.entropy", "accel.ratio.mutation", "freq.1.ratio.varBeta", "gp.1.ratio.varBeta", "gp.2.ratio.varBeta", "freq.3.antigenicTypes",
                 "freq.1.antigenicTypes", "freq.3.dominant.freq", "gp.1.antigenicDiversity", "freq.1.varBeta", "freq.3.prop.I", "freq.3.totalI", "freq.3.totalS",
                 "gp.1.prop.I", "gp.1.totalI", "gp.1.totalS", "gp.1.meanLoad", "gp.1.diversity", "freq.2.ratio.varSigma", "freq.3.ratio.varR","gp.1.ratio.varR",
                 "accel.ratio.varR", "freq.1.ratio.meanBeta", "gp.1.antigenicTypes", "freq.2.tmrca", "gp.2.ratio.mutation", "freq.3.varSigma", "gp.1.ratio.mutation",
                 "gp.1.entropy", "gp.1.dominant.freq", "freq.1.entropy", "gp.2.ratio.varR", "accel.antigenicDiversity", "gp.2.antigenicTypes", "gp.1.meanBeta")


################################### 
#### Real World 2/7/18

vif.excluded = c("freq2.diversity", "freq2.antigenicDiversity", "freq2.entropy", "freq2.antigenicTypes",
                 "freq3.diversity", "freq3.antigenicDiversity", "freq2.infected", "freq2.dominant.freq",
                 "freq3.entropy", "freq2.meanR", "freq3.antigenicTypes", "freq1.antigenicTypes", "freq2.totalI",
                 "freq2.totalS", "freq3.dominant.freq", "freq2.tmrca", "freq3.infected", "gp1.infected",
                 "gp1.antigenicTypes", "accel.antigenicDiversity", "gp1.diversity", "gp1.dominant.freq",
                 "freq1.dominant.freq", "gp1.entropy", "freq1.meanR", "gp1.meanR", "gp2.antigenicTypes")

####################################################################
# Freq 1
eliminated.variables = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR",
                         "individual.meanSigma", "ratio.meanSigma",  "meanSigma",   "covBetaSigma",  "meanBeta",
                         "individual.varR", "individual.varBeta", "ratio.varR", "antigenicTypes",  "ratio.meanBeta", "entropy")             
# Freq 2 


#########################################################################
## Poster Figures 


