
library(bayou)

tree <- read.tree("Carex2n.tree")
tree <- ladderize(tree, T)
data_tot <- read.table("data_tot.csv", sep=";", dec=".", header=T)
mydata <- data.frame(data_tot[,"original.name"],data_tot[,"mean"])
mydata2 <- data.frame(data_tot[,"original.name"], data_tot[,"weighted_squaredSEM"])
mydata[,1]  
Ntip(tree)

#2n
dat.2n <- mydata[,2]
names(dat.2n) <- as.character(mydata[,1])
se.2n <- mydata2[,2]
names(se.2n) <- as.character(mydata2[,1])

priorOU.2n <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1),
                                       dsig2=list(scale=0.1),
                                       dk=list(lambda=107, kmax=trunc(length(tree$edge.length)*10/100)), #107 sections represented
                                       dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.2n), 
                                                   sd=1.5*sd(dat.2n))),
                            plot.prior = TRUE)



#Run 1
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_1 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                            new.dir=TRUE, outname="./modelOU_r001", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_1$run(2500000) # Run the MCMC
chainOU.2n_1 <- mcmcOU.2n_1$load()
summary.2n_1 <- summary(chainOU.2n_1)

save(chainOU.2n_1, file = "mcmcOU.2n.2halfMrun_1.Rdata")

#Run 2
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_2 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                            new.dir=TRUE, outname="./modelOU_r002", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_2$run(2500000) # Run the MCMC
chainOU.2n_2 <- mcmcOU.2n_2$load()
summary.2n_2 <- summary(chainOU.2n_2)

save(chainOU.2n_2, file = "mcmcOU.2n.2halfMrun_2.Rdata")

#Run 3
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_3 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                              new.dir=TRUE, outname="./modelOU_r003", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_3$run(2500000) # Run the MCMC
chainOU.2n_3 <- mcmcOU.2n_3$load()
summary.2n_3 <- summary(chainOU.2n_3)

save(chainOU.2n_3, file = "mcmcOU.2n.2halfMrun_3.Rdata")

#Run 4
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_4 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                              new.dir=TRUE, outname="./modelOU_r004", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_4$run(2500000) # Run the MCMC
chainOU.2n_4 <- mcmcOU.2n_4$load()
summary.2n_4 <- summary(chainOU.2n_4)

save(chainOU.2n_4, file = "mcmcOU.2n.2halfMrun_4.Rdata")

#Run 5
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_5 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                              new.dir=TRUE, outname="./modelOU_r05", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_5$run(2500000) # Run the MCMC
chainOU.2n_5 <- mcmcOU.2n_5$load()
summary.2n_5 <- summary(chainOU.2n_5)

save(chainOU.2n_5, file = "mcmcOU.2n.2halfMrun_5.Rdata")

#Run 6
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_6 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                            new.dir=TRUE, outname="./modelOU_r006", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_6$run(2500000) # Run the MCMC
chainOU.2n_6 <- mcmcOU.2n_6$load()
summary.2n_6 <- summary(chainOU.2n_6)

save(chainOU.2n_6, file = "mcmcOU.2n.2halfMrun_6.Rdata")

#Run 7
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_7 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                            new.dir=TRUE, outname="./modelOU_r007", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_7$run(2500000) # Run the MCMC
chainOU.2n_7 <- mcmcOU.2n_7$load()
summary.2n_7 <- summary(chainOU.2n_7)

save(chainOU.2n_7, file = "mcmcOU.2n.2halfMrun_7.Rdata")

#Run 8
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_8 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                              new.dir=TRUE, outname="./modelOU_r008", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_8$run(2500000) # Run the MCMC
chainOU.2n_8 <- mcmcOU.2n_8$load()
summary.2n_8 <- summary(chainOU.2n_8)

save(chainOU.2n_8, file = "mcmcOU.2n.2halfMrun_8.Rdata")

#Run 9
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_9 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                              new.dir=TRUE, outname="./modelOU_r009", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_9$run(2500000) # Run the MCMC
chainOU.2n_9 <- mcmcOU.2n_9$load()
summary.2n_9 <- summary(chainOU.2n_9)

save(chainOU.2n_9, file = "mcmcOU.2n.2halfMrun_9.Rdata")

#Run 10
startpars.2n <- priorSim(priorOU.2n, tree, plot=T, cex=0.1)$pars[[1]]
priorOU.2n(startpars.2n)

mcmcOU.2n_10 <- bayou.makeMCMC(tree, dat.2n, SE=se.2n, prior=priorOU.2n, startpar=startpars.2n,
                              new.dir=TRUE, outname="./modelOU_r010", plot.freq=NULL) # Set up the MCMC

mcmcOU.2n_10$run(2500000) # Run the MCMC
chainOU.2n_10 <- mcmcOU.2n_10$load()
summary.2n_10 <- summary(chainOU.2n_10)

save(chainOU.2n_10, file = "mcmcOU.2n.2halfMrun_10.Rdata")

list_chains <- list(chainOU.2n_1, chainOU.2n_2, chainOU.2n_3, chainOU.2n_4, chainOU.2n_5,
                    chainOU.2n_6, chainOU.2n_7, chainOU.2n_8, chainOU.2n_9, chainOU.2n_10)
combine <- combine.chains(list_chains, burnin.prop = 0.3 )

summary.2n <- summary(combine )



par(mfrow=c(3, 5) )
plot(combine, auto.layout=FALSE)

par(mfrow=c(1, 1) )
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.3, cex = 0.1) 
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.5, cex = 0.1)
plotSimmap.mcmc(combine, burnin = 0, pp.cutoff = 0.7, cex = 0.1)
plotBranchHeatMap(tree, chain = combine, "theta", burnin = 0, pal = terrain.colors, cex = 0.1)
phenogram.density(tree, dat.2n, burnin = 0, combine, pp.cutoff = 0.3)
phenogram.density(tree, dat.2n, burnin = 0, combine, pp.cutoff = 0.5)
phenogram.density(tree, dat.2n, burnin = 0, combine, pp.cutoff = 0.7)

