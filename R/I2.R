
#' @title I^2 Function 
#' @description Function for calculating commonly reported  I^2 measures from MCMCglmm model objects.  See Nakagawa and Santos (2012) for detailed explanation on the various types of I^2
#' @param model The MCMCglmm model object
#' @param v The vector of sampling variance for each effect size. 
#' @param re.list A list of variance estimates one wishes to derive I^2 estimates for. These could include phylogenetic heritability (H^2), species ($I^{2}_{sp}$) and study ($I^{2}_{stdy}$) or often called ($tau^{2}$). At the moment only three types are provided: phylogenetic (phylo), species (spp) and study (stdy). The names of your specific variance component do not matter, but these names should be specified in the re.list argument in the respective argument.
#' @return A data.frame containing the relevant I^2 measures along with the 95 percent credible intervals for each element listed in the re.list argument 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @references Nakagawa, S. and Santos, E.S.A. (2012) Methodological issues and advances in biological meta-analysis. Evolutionary Ecology, 26:1253-1274.
#' @export

I2 <- function(model, v, re.list = list(phylo = "animal", spp = "species", stdy = "study")){
	
	if(class(model) != "MCMCglmm"){
	stop("The model object is not of class 'MCMCglmm'")
	}

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

	tmpMatrix <- cbind(I2.phylo, I2.stdy, I2.spp, I2.t)
	        mode <- posterior.mode(as.mcmc(tmpMatrix))
	             CI <- HPDinterval(as.mcmc(tmpMatrix))

	p.mcmc = ifelse(round(CI[,"lower"], digits =3) != 0, "*", "")
	
  return(data.frame(post.mode = round(mode, digits = 3), CI = round(CI, digits = 3), p.mcmc = p.mcmc))
}