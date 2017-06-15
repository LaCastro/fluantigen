rm(list=ls())

if(grepl('laurencastro', Sys.info()['login'])) {
  setwd('~/Documents/projects/fluantigen/analysis/')
  savepath <- "../data/"
}

if(grepl('lacastro', Sys.info()['login'])) {
  setwd('/home1/03131/lacastro/')
}



## First give number of simulations you want to run
N = seq(from = 1, to = 12, by = 1)
type.sim = "tropics"

sink('../launcher/run_sims.txt') # creates a text file in the launcher

for(n in N) {
    startCmd <- "java -jar trackAntigen.jar"
    dirName <- paste0(" ", type.sim, "_" , n)
    parmFile <- " parameters_tropics.yml"
    outputDir <- " data"
    grbCmd <- " -XX:+UseSerialGC"
    memoryCmd <- " -Xmx5G "
    mainCmd <- " Mutantigen"
    
    full_cmd <- paste0(startCmd, dirName, parmFile, outputDir, grbCmd, memoryCmd, mainCmd)
    #full_cmd <- paste0(startCmd, dirName, parmFile, outputDir,
    #                   grbCmd,mainCmd)
    
    cat(full_cmd) # puts it in the text file
    cat('\n') # starts a new line
  }
sink() # shows that this is done 

