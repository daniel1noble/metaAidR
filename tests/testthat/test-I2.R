context("Tests I2 function")

test_that("Check that output is correct", {
	library(metafor)
	library(MCMCglmm)
	
	#Simulate some meta-analytic data
		set.seed(123)
		n = 1000
		spN = 100
		stdyN =50 

	       V <- rgamma(n, shape = 1, rate = 10)
	 muEs <- rnorm(n, 3, sqrt(V))
	    spp <- rep(rnorm(spN, 0, 0.5), each = n / spN)
	    stdy <- rep(rnorm(stdyN, 0, 4), each = n / stdyN)
	         e <- rnorm(n, 0, 1)

	es <- muEs + spp + stdy + e
	data <- data.frame(es, spp = as.factor(spp), stdy = as.factor(stdy), V, obs = 1:length(es))

	# True I2
	wi <- 1/V  #weight
	Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))
	
	I2spp <- (0.5^2) / (4^2 + 0.5^2 + Vw + 1)
	I2st <- (4^2) / (4^2 + 0.5^2 + Vw + 1)
	I2T <- (4^2 + 0.5^2 + 1) / (4^2 + 0.5^2 + Vw + 1)

	#This model estimates a residual variance
	metaFor <- metafor::rma.mv(es, V, random = list(~1|spp, ~1|stdy, ~1|obs), struct="UN",data = data)

	#Run in MCMCglmm
	metaMCMC <- MCMCglmm::MCMCglmm(es~1, random = ~spp + stdy, mev = V, data = data)

	MetaF <- I2(metaFor, v = data$V, re.list=list(spp = "spp", stdy = "stdy"))
	MCMC <- I2(metaMCMC, v = data$V, re.list=list(spp = "spp", stdy = "stdy"))

	expect_equal(as.numeric(MetaF[rownames(MetaF) == "spp",][1]) ,  I2spp, tolerance = 0.005, info = "faildMetaspp")
	expect_equal(as.numeric(MetaF[rownames(MetaF) == "stdy",][1]) ,  I2st, tolerance = 0.05, info = "faildMetastudy")
	expect_equal(as.numeric(MetaF[rownames(MetaF) == "total",][1]) ,  I2T, tolerance = 0.05, info = "faildMetatot")

	expect_equal(as.numeric(MCMC[rownames(MCMC) == "spp",][1]) ,  I2spp, tolerance = 0.005, info = "faildMCspp")
	expect_equal(as.numeric(MCMC[rownames(MCMC) == "stdy",][1]) ,  I2st, tolerance = 0.05, info = "faildMCstudy")
	expect_equal(as.numeric(MCMC[rownames(MCMC) == "total",][1]) ,  I2T, tolerance = 0.05, info = "faildMCTot")
})