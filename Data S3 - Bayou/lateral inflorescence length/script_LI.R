
library(bayou)

data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)


mydata <- data.frame(data_tot$original.name[data_tot$LI.length.mp >= 0],data_tot$LI.length.mp[data_tot$LI.length.mp >= 0])
mydata2 <- data.frame(data_tot$original.name[data_tot$LI.length.mp >= 0], data_tot$LI.length.1.4[data_tot$LI.length.mp >= 0])
mydata[,1]  

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, T)
tree <- drop.tip(tree, setdiff(tree$tip.label, data_tot$original.name[data_tot$LI.length.mp >= 0]))

Ntip(tree)

#LI
dat.LI <- mydata[,2][is.na(mydata[,1])==FALSE]
names(dat.LI) <- as.character(mydata[,1][is.na(mydata[,1])==FALSE])
se.LI <- mydata2[,2][is.na(mydata2[,1])==FALSE]
names(se.LI) <- as.character(mydata2[,1][is.na(mydata2[,1])==FALSE])


priorOU.LI <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dlnorm"), 
                         param=list(dalpha=list(scale=0.1),
                                    dsig2=list(scale=0.1),
                                    dk=list(lambda=trunc(length(tree$edge.length)*5/100), kmax=trunc(length(tree$edge.length)*10/100)), #66 #132
                                    dsb=list(bmax=1, prob=1), 
                                    dtheta=list(meanlog=log(mean(dat.LI)), 
                                                sdlog=1)),
                         plot.prior = TRUE)

#Run 1
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_1 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_1$run(1000000) # Run the MCMC
chainOU.LI_1 <- mcmcOU.LI_1$load()
summary.LI_1 <- summary(chainOU.LI_1)


#Run 2
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_2 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_2$run(1000000) # Run the MCMC
chainOU.LI_2 <- mcmcOU.LI_2$load()
summary.LI_2 <- summary(chainOU.LI_2)

#Run 3
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_3 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_3$run(1000000) # Run the MCMC
chainOU.LI_3 <- mcmcOU.LI_3$load()
summary.LI_3 <- summary(chainOU.LI_3)

#Run 4
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_4 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_4$run(1000000) # Run the MCMC
chainOU.LI_4 <- mcmcOU.LI_4$load()
summary.LI_4 <- summary(chainOU.LI_4)

#Run 5
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_5 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                              new.dir=TRUE, outname="./modelOU_r05", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_5$run(1000000) # Run the MCMC
chainOU.LI_5 <- mcmcOU.LI_5$load()
summary.LI_5 <- summary(chainOU.LI_5)

#Run 6
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_6 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_6$run(1000000) # Run the MCMC
chainOU.LI_6 <- mcmcOU.LI_6$load()
summary.LI_6 <- summary(chainOU.LI_6)


#Run 7
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_7 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                            new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_7$run(1000000) # Run the MCMC
chainOU.LI_7 <- mcmcOU.LI_7$load()
summary.LI_7 <- summary(chainOU.LI_7)

#Run 8
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_8 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_8$run(1000000) # Run the MCMC
chainOU.LI_8 <- mcmcOU.LI_8$load()
summary.LI_8 <- summary(chainOU.LI_8)

#Run 9
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_9 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_9$run(1000000) # Run the MCMC
chainOU.LI_9 <- mcmcOU.LI_9$load()
summary.LI_9 <- summary(chainOU.LI_9)

#Run 10
startpars.LI <- priorSim(priorOU.LI, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.LI(startpars.LI)

mcmcOU.LI_10 <- bayou.makeMCMC(tree, dat.LI, SE=se.LI, prior=priorOU.LI, startpar=startpars.LI,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.LI_10$run(1000000) # Run the MCMC
chainOU.LI_10 <- mcmcOU.LI_10$load()
summary.LI_10 <- summary(chainOU.LI_10)


list_chains <- list(chainOU.LI_1, chainOU.LI_2, chainOU.LI_3, chainOU.LI_4, chainOU.LI_5,
                    chainOU.LI_6, chainOU.LI_7, chainOU.LI_8, chainOU.LI_9, chainOU.LI_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.LI <- summary(combine )


par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.LI, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.LI, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.LI, burnin = 0, combine, pp.cutoff = 0.7)
