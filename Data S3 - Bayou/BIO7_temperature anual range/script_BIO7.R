# https://github.com/uyedaj/EQG2016/wiki
library(bayou)
library(ape)

data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)


mydata <- data.frame(data_tot$original.name[is.na(data_tot$bio7_Mean)==FALSE],data_tot$bio7_Mean[is.na(data_tot$bio7_Mean)==FALSE])
mydata2 <- data.frame(data_tot$original.name[is.na(data_tot$bio7_Mean)==FALSE], data_tot$squaredSEM_bioclim7[is.na(data_tot$bio7_Mean)==FALSE])
mydata[,1]  

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, F)
tree <- drop.tip(tree, setdiff(tree$tip.label, data_tot$original.name[is.na(data_tot$bio7_Mean)==FALSE]))

#Ntip(tree)

#bio7
dat.bio7 <- mydata[,2][is.na(mydata[,1])==FALSE]
names(dat.bio7) <- as.character(mydata[,1][is.na(mydata[,1])==FALSE])
se.bio7 <- mydata2[,2][is.na(mydata2[,1])==FALSE]
names(se.bio7) <- as.character(mydata2[,1][is.na(mydata2[,1])==FALSE])


priorOU.bio7 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                         param=list(dalpha=list(scale=0.1),
                                    dsig2=list(scale=0.1),
                                    dk=list(lambda=trunc(length(tree$edge.length)*2.5/100), kmax=trunc(length(tree$edge.length)*5/100)), #35 #71
                                    dsb=list(bmax=1, prob=1), 
                                    dtheta=list(mean=mean(dat.bio7), 
                                                sd=1.5*sd(dat.bio7))),
                         plot.prior = TRUE)

#Run 1
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_1 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_1$run(1000000) # Run the MCMC
chainOU.bio7_1 <- mcmcOU.bio7_1$load()
summary.bio7_1 <- summary(chainOU.bio7_1)


#Run 2
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_2 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_2$run(1000000) # Run the MCMC
chainOU.bio7_2 <- mcmcOU.bio7_2$load()
summary.bio7_2 <- summary(chainOU.bio7_2)

#Run 3
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_3 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_3$run(1000000) # Run the MCMC
chainOU.bio7_3 <- mcmcOU.bio7_3$load()
summary.bio7_3 <- summary(chainOU.bio7_3)

#Run 4
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_4 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_4$run(1000000) # Run the MCMC
chainOU.bio7_4 <- mcmcOU.bio7_4$load()
summary.bio7_4 <- summary(chainOU.bio7_4)

#Run 5
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_5 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                              new.dir=TRUE, outname="./modelOU_r05", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_5$run(1000000) # Run the MCMC
chainOU.bio7_5 <- mcmcOU.bio7_5$load()
summary.bio7_5 <- summary(chainOU.bio7_5)

#Run 6
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_6 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_6$run(1000000) # Run the MCMC
chainOU.bio7_6 <- mcmcOU.bio7_6$load()
summary.bio7_6 <- summary(chainOU.bio7_6)


#Run 7
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_7 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                            new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_7$run(1000000) # Run the MCMC
chainOU.bio7_7 <- mcmcOU.bio7_7$load()
summary.bio7_7 <- summary(chainOU.bio7_7)

#Run 8
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_8 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_8$run(1000000) # Run the MCMC
chainOU.bio7_8 <- mcmcOU.bio7_8$load()
summary.bio7_8 <- summary(chainOU.bio7_8)

#Run 9
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_9 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_9$run(1000000) # Run the MCMC
chainOU.bio7_9 <- mcmcOU.bio7_9$load()
summary.bio7_9 <- summary(chainOU.bio7_9)

#Run 10
startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio7(startpars.bio7)

mcmcOU.bio7_10 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, startpar=startpars.bio7,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio7_10$run(1000000) # Run the MCMC
chainOU.bio7_10 <- mcmcOU.bio7_10$load()
summary.bio7_10 <- summary(chainOU.bio7_10)


list_chains <- list(chainOU.bio7_1, chainOU.bio7_2, chainOU.bio7_3, chainOU.bio7_4, chainOU.bio7_5,
                    chainOU.bio7_6, chainOU.bio7_7, chainOU.bio7_8, chainOU.bio7_9, chainOU.bio7_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.bio7 <- summary(combine )

par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.bio7, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.bio7, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.bio7, burnin = 0, combine, pp.cutoff = 0.7)
