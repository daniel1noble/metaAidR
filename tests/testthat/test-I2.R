context("Tests I2 function")

test_that("Check that output is correct", {
	library(phytools)
	library(metafor)
	library(MCMCglmm)
	treeSim <- rtree(n = 30)
	
	#Simulate some meta-analytic data
	set.seed(69)
	n = 1000
	spN = 100
	stdyN =50 

	       V <- rgamma(n, shape = 1, rate = 10)
	 muEs <- rnorm(n, 3, sqrt(V))
	    spp <- rep(rnorm(spN, 0, 0.5), each = n / spN)
	    stdy <- rep(rnorm(stdyN, 0, 4), each = n / stdyN)
	         e <- rnorm(n, 0, 1)

	es <- muEs + spp + stdy + e
	data <- data.frame(es, spp = as.factor(spp), stdy = as.factor(stdy), V)

	# Run in metafor
	metaFor <- rma.mv(es, V, random = list(~1|spp, ~1|stdy), data = data)

	#Run in MCMCglmm
	metaMCMC <- MCMCglmm(es~1, random = ~spp + stdy, mev = V, data = data)

	het <- I2(metaMCMC, v = data$V, re.list=list(phylo = "animal", spp = "spp", stdy = "stdy"))
})