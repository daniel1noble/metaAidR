context("Tests I2 function")

test_that("Check that output is correct", {
	library(metafor)
	library(MCMCglmm)
	library(lme4)
	
	#Simulate some meta-analytic data
	set.seed(123)
	n = 1000
	spN = 100
	stdyN =50 
	rgN = 10
	SDsp = 0.5
	SDsty = 4
	SDrg = 3
	SDe = 1

	V <- rgamma(n, shape = 1, rate = 10)
	muEs <- rnorm(n, 3, sqrt(V))
	spp <- rep(rnorm(spN, 0, SDsp), each = n / spN)
	stdy <- rep(rnorm(stdyN, 0, SDsty), each = n / stdyN)
	rg <- rep(rnorm(rgN, 0, SDrg), each = n / rgN)
	e <- rnorm(n, 0, SDe)

	#Simulated data 1
	es <- muEs + spp + stdy + e
	data <- data.frame(es, spp = as.factor(spp), stdy = as.factor(stdy), V, obs = 1:length(es))

	#Simulated data 2
	es2 <- muEs + spp + stdy + rg + e
	data2 <- data.frame(es, spp = as.factor(spp), stdy = as.factor(stdy), rg = as.factor(rg), V, obs = 1:length(es))

	# True I2
	wi <- 1/V  #weight
	Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))
	
	I2spp <- (SDsp^2) / (SDsty^2 + SDsp^2 + Vw + SDe)
	
	I2st <- (SDsty^2) / (SDsty^2 + SDsp^2 + Vw + SDe)
	I2st2 <- (SDsty^2) / (SDsty^2 + SDrg^2 + SDsp^2 + Vw + SDe)

	I2T <- (SDsty^2 + SDsp^2 + SDe) / (SDsty^2 + SDsp^2 + Vw + SDe)

	# Run with metafor
	 metaFor <- metafor::rma.mv(es, V, random = list(~1|spp, ~1|stdy, ~1|obs), data = data)
	metaFor2 <- metafor::rma.mv(es, V, random = list(~1|spp, ~1|stdy, ~1|rg, ~1|obs), data = data2)
	
	#Run in MCMCglmm
	 metaMCMC <- MCMCglmm::MCMCglmm(es~1, random = ~spp + stdy, mev = V, data = data)
	metaMCMC2 <- MCMCglmm::MCMCglmm(es~1, random = ~spp + stdy + rg, mev = V, data = data2)

	#Test should fail
	metaFor3 <- metafor::rma.mv(es, V, random = list(~1|spp, ~1|stdy), data = data)

	# Wrong model!
	lmER <- lme4::lmer(es ~1 + (1|spp) + (1|stdy), weights = V, data = data)

	# Testing the I2 function
	 MetaF <- I2(metaFor, v = data$V)
	MetaF2 <- I2(metaFor2, v = data$V)

	 MCMC <- I2(metaMCMC, v = data$V, phylo = FALSE)
	MCMC2 <- I2(metaMCMC2, v = data$V, phylo = FALSE)

	#Tests
	expect_equal(as.numeric(MetaF[rownames(MetaF) == "spp",][1]) ,  
		 		 I2spp, tolerance = 0.01, info = "faildMetaspp")
	expect_equal(as.numeric(MetaF[rownames(MetaF) == "stdy",][1]) ,  
					I2st, tolerance = 0.05, info = "faildMetastudy")
	expect_equal(as.numeric(MetaF[rownames(MetaF) == "total",][1]) ,  
					I2T, tolerance = 0.05, info = "faildMetatot")
	
	expect_equal(dim(MetaF2),  c(4,3), info = "dimMetaforFail")
	expect_error(I2(metaFor3, v = data$V), info = "failedMetafor_no_obs")
	expect_error(I2(lmER, v = data$V), info = "failedWrongModel")

	expect_equal(as.numeric(MCMC[rownames(MCMC) == "spp",][1]) ,  
				I2spp, tolerance = 0.005, info = "faildMCspp")
	expect_equal(as.numeric(MCMC[rownames(MCMC) == "stdy",][1]) ,  
				I2st, tolerance = 0.05, info = "faildMCstudy")
	expect_equal(as.numeric(MCMC[rownames(MCMC) == "total",][1]) ,  
				I2T, tolerance = 0.05, info = "faildMCTot")
	expect_equal(dim(MCMC2),  c(4,3), info = "dimMCMCFail")
})