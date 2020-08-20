
library(bayou)

data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)


mydata <- data.frame(data_tot$original.name[data_tot$Culm.length.mp > 0],data_tot$Culm.length.mp[data_tot$Culm.length.mp > 0])
mydata2 <- data.frame(data_tot$original.name[data_tot$Culm.length.mp > 0], data_tot$Culm.length.1.4[data_tot$Culm.length.mp > 0])
mydata[,1]  

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, T)
tree <- drop.tip(tree, setdiff(tree$tip.label, data_tot$original.name[data_tot$Culm.length.mp > 0]))

Ntip(tree)

#culm
dat.culm <- mydata[,2][is.na(mydata[,1])==FALSE]
names(dat.culm) <- as.character(mydata[,1][is.na(mydata[,1])==FALSE])
se.culm <- mydata2[,2][is.na(mydata2[,1])==FALSE]
names(se.culm) <- as.character(mydata2[,1][is.na(mydata2[,1])==FALSE])


priorOU.culm <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dlnorm"), 
                         param=list(dalpha=list(scale=0.1),
                                    dsig2=list(scale=0.1),
                                    dk=list(lambda=trunc(length(tree$edge.length)*2.5/100), kmax=trunc(length(tree$edge.length)*5/100)), #35 #71
                                    dsb=list(bmax=1, prob=1), 
                                    dtheta=list(meanlog=log(mean(dat.culm)), 
                                                sdlog=0.5)),
                         plot.prior = TRUE)

#Run 1
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_1 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_1$run(1000000) # Run the MCMC
chainOU.culm_1 <- mcmcOU.culm_1$load()
summary.culm_1 <- summary(chainOU.culm_1)


#Run 2
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_2 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_2$run(1000000) # Run the MCMC
chainOU.culm_2 <- mcmcOU.culm_2$load()
summary.culm_2 <- summary(chainOU.culm_2)


#Run 3
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_3 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_3$run(1000000) # Run the MCMC
chainOU.culm_3 <- mcmcOU.culm_3$load()
summary.culm_3 <- summary(chainOU.culm_3)


#Run 4
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_4 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_4$run(1000000) # Run the MCMC
chainOU.culm_4 <- mcmcOU.culm_4$load()
summary.culm_4 <- summary(chainOU.culm_4)


#Run 5
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_5 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r005", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_5$run(1000000) # Run the MCMC
chainOU.culm_5 <- mcmcOU.culm_5$load()
summary.culm_5 <- summary(chainOU.culm_5)


#Run 6
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_6 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_6$run(1000000) # Run the MCMC
chainOU.culm_6 <- mcmcOU.culm_6$load()
summary.culm_6 <- summary(chainOU.culm_6)


#Run 7
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_7 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_7$run(1000000) # Run the MCMC
chainOU.culm_7 <- mcmcOU.culm_7$load()
summary.culm_7 <- summary(chainOU.culm_7)


#Run 8
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_8 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_8$run(1000000) # Run the MCMC
chainOU.culm_8 <- mcmcOU.culm_8$load()
summary.culm_8 <- summary(chainOU.culm_8)


#Run 9
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_9 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_9$run(1000000) # Run the MCMC
chainOU.culm_9 <- mcmcOU.culm_9$load()
summary.culm_9 <- summary(chainOU.culm_9)


#Run 10
startpars.culm <- priorSim(priorOU.culm, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.culm(startpars.culm)

mcmcOU.culm_10 <- bayou.makeMCMC(tree, dat.culm, SE=se.culm, prior=priorOU.culm, startpar=startpars.culm,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.culm_10$run(1000000) # Run the MCMC
chainOU.culm_10 <- mcmcOU.culm_10$load()
summary.culm_10 <- summary(chainOU.culm_10)


list_chains <- list(chainOU.culm_1, chainOU.culm_2, chainOU.culm_3, chainOU.culm_4, chainOU.culm_5,
					chainOU.culm_6, chainOU.culm_7, chainOU.culm_8, chainOU.culm_9, chainOU.culm_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.culm <- summary(combine )


par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.culm, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.culm, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.culm, burnin = 0, combine, pp.cutoff = 0.7)
