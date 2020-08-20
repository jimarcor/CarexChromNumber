
library(bayou)

data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)


mydata <- data.frame(data_tot$original.name[is.na(data_tot$bio12_Mean)==FALSE],data_tot$bio12_Mean[is.na(data_tot$bio12_Mean)==FALSE])
mydata2 <- data.frame(data_tot$original.name[is.na(data_tot$bio12_Mean)==FALSE], data_tot$squaredSEM_bioclim12[is.na(data_tot$bio12_Mean)==FALSE])
mydata[,1]  

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, T)
tree <- drop.tip(tree, setdiff(tree$tip.label, data_tot$original.name[is.na(data_tot$bio12_Mean)==FALSE]))

#Ntip(tree)

#bio12
dat.bio12 <- mydata[,2][is.na(mydata[,1])==FALSE]
names(dat.bio12) <- as.character(mydata[,1][is.na(mydata[,1])==FALSE])
se.bio12 <- mydata2[,2][is.na(mydata2[,1])==FALSE]
names(se.bio12) <- as.character(mydata2[,1][is.na(mydata2[,1])==FALSE])


priorOU.bio12 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                         param=list(dalpha=list(scale=0.1),
                                    dsig2=list(scale=0.1),
                                    dk=list(lambda=trunc(length(tree$edge.length)*2.5/100), kmax=trunc(length(tree$edge.length)*5/100)), #35 #71
                                    dsb=list(bmax=1, prob=1), 
                                    dtheta=list(mean=mean(dat.bio12), 
                                                sd=1.5*sd(dat.bio12))),
                         plot.prior = TRUE)

#Run 1
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_1 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_1$run(1000000) # Run the MCMC
chainOU.bio12_1 <- mcmcOU.bio12_1$load()
summary.bio12_1 <- summary(chainOU.bio12_1)


#Run 2
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_2 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_2$run(1000000) # Run the MCMC
chainOU.bio12_2 <- mcmcOU.bio12_2$load()
summary.bio12_2 <- summary(chainOU.bio12_2)

#Run 3
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_3 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_3$run(1000000) # Run the MCMC
chainOU.bio12_3 <- mcmcOU.bio12_3$load()
summary.bio12_3 <- summary(chainOU.bio12_3)

#Run 4
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_4 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_4$run(1000000) # Run the MCMC
chainOU.bio12_4 <- mcmcOU.bio12_4$load()
summary.bio12_4 <- summary(chainOU.bio12_4)

#Run 5
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_5 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                              new.dir=TRUE, outname="./modelOU_r05", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_5$run(1000000) # Run the MCMC
chainOU.bio12_5 <- mcmcOU.bio12_5$load()
summary.bio12_5 <- summary(chainOU.bio12_5)

#Run 6
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_6 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_6$run(1000000) # Run the MCMC
chainOU.bio12_6 <- mcmcOU.bio12_6$load()
summary.bio12_6 <- summary(chainOU.bio12_6)


#Run 7
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_7 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                            new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_7$run(1000000) # Run the MCMC
chainOU.bio12_7 <- mcmcOU.bio12_7$load()
summary.bio12_7 <- summary(chainOU.bio12_7)

#Run 8
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_8 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_8$run(1000000) # Run the MCMC
chainOU.bio12_8 <- mcmcOU.bio12_8$load()
summary.bio12_8 <- summary(chainOU.bio12_8)

#Run 9
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_9 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_9$run(1000000) # Run the MCMC
chainOU.bio12_9 <- mcmcOU.bio12_9$load()
summary.bio12_9 <- summary(chainOU.bio12_9)

#Run 10
startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio12(startpars.bio12)

mcmcOU.bio12_10 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, startpar=startpars.bio12,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio12_10$run(1000000) # Run the MCMC
chainOU.bio12_10 <- mcmcOU.bio12_10$load()
summary.bio12_10 <- summary(chainOU.bio12_10)


list_chains <- list(chainOU.bio12_1, chainOU.bio12_2, chainOU.bio12_3, chainOU.bio12_4, chainOU.bio12_5,
                    chainOU.bio12_6, chainOU.bio12_7, chainOU.bio12_8, chainOU.bio12_9, chainOU.bio12_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.bio12 <- summary(combine )


par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.bio12, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.bio12, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.bio12, burnin = 0, combine, pp.cutoff = 0.7)

