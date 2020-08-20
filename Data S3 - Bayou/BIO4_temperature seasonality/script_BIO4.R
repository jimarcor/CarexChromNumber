
library(bayou)

data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)


mydata <- data.frame(data_tot$original.name[is.na(data_tot$bio4_Mean)==FALSE],data_tot$bio4_Mean[is.na(data_tot$bio4_Mean)==FALSE])
mydata2 <- data.frame(data_tot$original.name[is.na(data_tot$bio4_Mean)==FALSE], data_tot$squaredSEM_bioclim4[is.na(data_tot$bio4_Mean)==FALSE])
mydata[,1]  

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, T)
tree <- drop.tip(tree, setdiff(tree$tip.label, data_tot$original.name[is.na(data_tot$bio4_Mean)==FALSE]))

#Ntip(tree)

#bio4
dat.bio4 <- mydata[,2][is.na(mydata[,1])==FALSE]
names(dat.bio4) <- as.character(mydata[,1][is.na(mydata[,1])==FALSE])
se.bio4 <- mydata2[,2][is.na(mydata2[,1])==FALSE]
names(se.bio4) <- as.character(mydata2[,1][is.na(mydata2[,1])==FALSE])


priorOU.bio4 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                         param=list(dalpha=list(scale=0.1),
                                    dsig2=list(scale=0.1),
                                    dk=list(lambda=trunc(length(tree$edge.length)*2.5/100), kmax=trunc(length(tree$edge.length)*5/100)), #35 #71
                                    dsb=list(bmax=1, prob=1), 
                                    dtheta=list(mean=mean(dat.bio4), 
                                                sd=1.5*sd(dat.bio4))),
                         plot.prior = TRUE)

#Run 1
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_1 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_1$run(1000000) # Run the MCMC
chainOU.bio4_1 <- mcmcOU.bio4_1$load()
summary.bio4_1 <- summary(chainOU.bio4_1)


#Run 2
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_2 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_2$run(1000000) # Run the MCMC
chainOU.bio4_2 <- mcmcOU.bio4_2$load()
summary.bio4_2 <- summary(chainOU.bio4_2)

#Run 3
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_3 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_3$run(1000000) # Run the MCMC
chainOU.bio4_3 <- mcmcOU.bio4_3$load()
summary.bio4_3 <- summary(chainOU.bio4_3)

#Run 4
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_4 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_4$run(1000000) # Run the MCMC
chainOU.bio4_4 <- mcmcOU.bio4_4$load()
summary.bio4_4 <- summary(chainOU.bio4_4)

#Run 5
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_5 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                              new.dir=TRUE, outname="./modelOU_r05", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_5$run(1000000) # Run the MCMC
chainOU.bio4_5 <- mcmcOU.bio4_5$load()
summary.bio4_5 <- summary(chainOU.bio4_5)

#Run 6
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_6 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_6$run(1000000) # Run the MCMC
chainOU.bio4_6 <- mcmcOU.bio4_6$load()
summary.bio4_6 <- summary(chainOU.bio4_6)


#Run 7
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_7 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                            new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_7$run(1000000) # Run the MCMC
chainOU.bio4_7 <- mcmcOU.bio4_7$load()
summary.bio4_7 <- summary(chainOU.bio4_7)

#Run 8
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_8 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_8$run(1000000) # Run the MCMC
chainOU.bio4_8 <- mcmcOU.bio4_8$load()
summary.bio4_8 <- summary(chainOU.bio4_8)

#Run 9
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_9 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_9$run(1000000) # Run the MCMC
chainOU.bio4_9 <- mcmcOU.bio4_9$load()
summary.bio4_9 <- summary(chainOU.bio4_9)

#Run 10
startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.bio4(startpars.bio4)

mcmcOU.bio4_10 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, startpar=startpars.bio4,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.bio4_10$run(1000000) # Run the MCMC
chainOU.bio4_10 <- mcmcOU.bio4_10$load()
summary.bio4_10 <- summary(chainOU.bio4_10)


list_chains <- list(chainOU.bio4_1, chainOU.bio4_2, chainOU.bio4_3, chainOU.bio4_4, chainOU.bio4_5,
                    chainOU.bio4_6, chainOU.bio4_7, chainOU.bio4_8, chainOU.bio4_9, chainOU.bio4_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.bio4 <- summary(combine )


par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.bio4, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.bio4, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.bio4, burnin = 0, combine, pp.cutoff = 0.7)
