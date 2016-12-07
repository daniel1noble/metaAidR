
#' @title I^2 Function 
#' @description Function for calculating commonly reported  I^2 measures from MCMCglmm model objects.  See Nakagawa and Santos (2012) for detailed explanation on the various types of I^2
#' @param model The MCMCglmm model object
#' @param v The vector of sampling variance for each effect size. 
#' @param sims The number of simulations used for calculating distribution for metafor objects.
#' @param re.list A list of variance estimates one wishes to derive I^2 estimates for. These could include phylogenetic heritability (H^2), species ($I^{2}_{sp}$) and study ($I^{2}_{stdy}$) or often called ($tau^{2}$). At the moment only three types are provided: phylogenetic (phylo), species (spp) and study (stdy). The names of your specific variance component do not matter, but these names should be specified in the re.list argument in the respective argument.
#' @return A data.frame containing the relevant I^2 measures along with the 95 percent credible intervals for each element listed in the re.list argument 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @references Nakagawa, S. and Santos, E.S.A. (2012) Methodological issues and advances in biological meta-analysis. Evolutionary Ecology, 26:1253-1274.
#' @export

I2 <- function(model, v, sims = 1500, re.list = list(phylo = "animal", spp = "species", stdy = "study")){
	
	if(class(model) != "MCMCglmm" & class(model) != "metafor"){
		stop("The model object is not of class 'MCMCglmm' or 'metafor'")
	}

	if(class(model) == "MCMCglmm" ){
		#Calculate measurement error
		  wi <- 1/v  #weight
		Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))
		
		# Get posterior distribution 
		post <- model$VCV

		#Extract names of re.list.
		phylo <- post[,grep(re.list$phylo, colnames(post))]
		   spp  <- post[,grep(re.list$spp, colnames(post))]
		   stdy <- post[,grep(re.list$stdy, colnames(post))]
		        r <- post[,grep("units", colnames(post))]

		# Calculate the proportion of heterogeneity explained by variance components
		I2.phylo <- (phylo) / (phylo + spp + stdy + r)
		I2.stdy <- (stdy) / (phylo + spp + stdy + r + Vw)
		I2.spp <- (spp) / (phylo + spp + stdy + r + Vw)
		I2.t <- (phylo + spp + stdy + r) / (phylo + spp + stdy + r + Vw)

		tmpMatrix <- Matrix::cBind(I2.phylo, I2.stdy, I2.spp, I2.t)
		        mode <- MCMCglmm::posterior.mode(coda::as.mcmc(tmpMatrix))
		             CI <- coda::HPDinterval(coda::as.mcmc(tmpMatrix))

		p.mcmc = ifelse(round(CI[,"lower"], digits =3) != 0, "*", "")
		
	  return(data.frame(post.mode = round(mode, digits = 3), CI = round(CI, digits = 3), p.mcmc = p.mcmc))
  	}

  	if(class(model) == "metafor"){
  		#Monte Carlo Simulations
  		 wi <- 1/v  #weight
		Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))

		# General form of I2 for meta-regression; From http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

		W <- diag(1/v)
		X <- model.matrix(model)
		P <- W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W # Not clear what P is estimating....if residual variance its way off from what is should be, 1. I believe that this is the amount of variance explained by the fixed effects. We could estimate this by including an observational random effect and estimating residual variance, but user would need to do this during model fitting. Turns out that this appears to be both sampling variance and the fixed effect variance!

		# From metafor extract the important statistics
  		sigma2 <- matrix(model$sigma2, nrow = 1, ncol = length(model$sigma2))
  		colnames(sigma2) <- model$s.names

  		#For each variance estimate simulate data
  		Sims <- as.data.frame(sapply(sigma2, function(x) simulate(x, sims = 1000)))
		colnames(Sims) <- colnames(sigma2) 
		Vt <- rowSums(cBind(Sims, Vw))
		
		I2_MC <- Sims / Vt
		I_CI <- ldply(lapply(I2_MC, function(x) quantile(x, c(0.025, 0.975))))

  	}

}


#' @title Parametric simulation 
#' @description Function for calculating I2 estimates using parametric simulations of model estimates taken from metafor. Note that the effectiveness of these simulations depends on the accuracy of model variance estimates.
#' @param estimate The estimate (i.e. variance) from a metafor model
#' @param sims The number of simulations 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
  simulate <- function(estimate, sims){
  		set.seed(07)
  		tmp <- data.frame(num = rep(1:sims, each = 1000), y = rnorm(1000*sims, 0, sqrt(estimate)))
  		Var <- plyr::ddply(tmp, .(num), summarise, var = var(y))
  		return(as.numeric(Var$var))
  	}