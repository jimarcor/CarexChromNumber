
library(bayou)

data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)


mydata <- data.frame(data_tot$original.name[is.na(data_tot$bio1_Mean)==FALSE],data_tot$bio1_Mean[is.na(data_tot$bio1_Mean)==FALSE])
mydata2 <- data.frame(data_tot$original.name[is.na(data_tot$bio1_Mean)==FALSE], data_tot$squaredSEM_bioclim1[is.na(data_tot$bio1_Mean)==FALSE])
mydata[,1]  

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, T)
tree <- drop.tip(tree, setdiff(tree$tip.label, data_tot$original.name[is.na(data_tot$bio1_Mean)==FALSE]))

#Ntip(tree)

#bio1
dat.bio1 <- mydata[,2][is.na(mydata[,1])==FALSE]
names(dat.bio1) <- as.character(mydata[,1][is.na(mydata[,1])==FALSE])
se.bio1 <- mydata2[,2][is.na(mydata2[,1])==FALSE]
names(se.bio1) <- as.character(mydata2[,1][is.na(mydata2[,1])==FALSE])


priorOU.bio1 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                         param=list(dalpha=list(scale=0.1),
                                    dsig2=list(scale=0.1),
                                    dk=list(lambda=trunc(length(tree$edge.length)*2.5/100), kmax=trunc(length(tree$edge.length)*5/100)), #35 #71
                                    dsb=list(bmax=1, prob=1), 
                                    dtheta=list(mean=mean(dat.bio1), 
                                                sd=1.5*sd(dat.bio1))),
                         plot.prior = TRUE)

#Run 1
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_1 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_1$run(1000000) # Run the MCMC
chainOU.bio1_1 <- mcmcOU.bio1_1$load()
summary.bio1_1 <- summary(chainOU.bio1_1)


#Run 2
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_2 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_2$run(1000000) # Run the MCMC
chainOU.bio1_2 <- mcmcOU.bio1_2$load()
summary.bio1_2 <- summary(chainOU.bio1_2)

#Run 3
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_3 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_3$run(1000000) # Run the MCMC
chainOU.bio1_3 <- mcmcOU.bio1_3$load()
summary.bio1_3 <- summary(chainOU.bio1_3)

#Run 4
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_4 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_4$run(1000000) # Run the MCMC
chainOU.bio1_4 <- mcmcOU.bio1_4$load()
summary.bio1_4 <- summary(chainOU.bio1_4)

#Run 5
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_5 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                              new.dir=TRUE, outname="./modelOU_r05", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_5$run(1000000) # Run the MCMC
chainOU.bio1_5 <- mcmcOU.bio1_5$load()
summary.bio1_5 <- summary(chainOU.bio1_5)

#Run 6
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_6 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_6$run(1000000) # Run the MCMC
chainOU.bio1_6 <- mcmcOU.bio1_6$load()
summary.bio1_6 <- summary(chainOU.bio1_6)


#Run 7
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_7 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                            new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_7$run(1000000) # Run the MCMC
chainOU.bio1_7 <- mcmcOU.bio1_7$load()
summary.bio1_7 <- summary(chainOU.bio1_7)

#Run 8
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_8 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_8$run(1000000) # Run the MCMC
chainOU.bio1_8 <- mcmcOU.bio1_8$load()
summary.bio1_8 <- summary(chainOU.bio1_8)

#Run 9
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_9 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_9$run(1000000) # Run the MCMC
chainOU.bio1_9 <- mcmcOU.bio1_9$load()
summary.bio1_9 <- summary(chainOU.bio1_9)

#Run 10
startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio1(startpars.bio1)

mcmcOU.bio1_10 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, startpar=startpars.bio1,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio1_10$run(1000000) # Run the MCMC
chainOU.bio1_10 <- mcmcOU.bio1_10$load()
summary.bio1_10 <- summary(chainOU.bio1_10)


list_chains <- list(chainOU.bio1_1, chainOU.bio1_2, chainOU.bio1_3, chainOU.bio1_4, chainOU.bio1_5,
                    chainOU.bio1_6, chainOU.bio1_7, chainOU.bio1_8, chainOU.bio1_9, chainOU.bio1_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.bio1 <- summary(combine )


par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.bio1, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.bio1, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.bio1, burnin = 0, combine, pp.cutoff = 0.7)
