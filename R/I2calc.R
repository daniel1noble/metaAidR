
#' @title $I^2$ Function 
#' @description Function for calculating \eqn{I^2} from MCMCglmm model objects according to Nakagawa and Santos (2012)
#' @param model The MCMCglmm model object
#' @param v The vector of sampling variance for each effect size. 
#' @param group The specific variance components one wishes to extract. Note that mev might be included in the model or as part of ginverse and these will be excluded if group = FALSE
#' @param re.list A list of variance estimates one wishes to derive $I^2$ estimates for. These could include phylogenetic heritability (H^{2}), species ($I^2_{sp}$) and study (\eqn{I^2_{stdy}} or often called ($\tau^{2}$). 
#' @return A data.frame containing the relevant I^2 measures along with the 95$%$ credible intervals for each element listed in the re.list argument 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @references Nakagawa, S. and Santos, E.S.A. (2012) Methodological issues and advances in biological meta-analysis. Evolutionary Ecology, 26:1253-1274.
#' @export

I2calc <- function(model, v, group = FALSE, re.list = list(phylo = "phylo", spp = "spp", stdy = "Paper_No")){
	
	if(class(model) != "MCMCglmm"){
	stop("The model object is not of class 'MCMCglmm'")
	}

	#calculate measurement error
	wi <- 1/v  #weight
	Vm <- Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))
	
	# Get posterior distribution 
	if(group == FALSE){
	post <- model$VCV[,!colnames(model$VCV) %in% c("esID", "sqrt(mev):sqrt(mev).meta")]
	}else{
	post <- model$VCV[, colnames(model$VCV)[grep(group, colnames(model$VCV))]]
	}

	#Extract names of re.list. 
	phylo <- post[,grep(re.list$phylo, colnames(post))]
	   spp  <- post[,grep(re.list$spp, colnames(post))]
	   stdy <- post[,grep(re.list$stdy, colnames(post))]
	        r <- post[,grep("units", colnames(post))]
 

	# Calculate the proportion of heterogeneity explained by "phylogeny"
	I2.phylo <- phylo/(phylo+spp+stdy+r)

	I2.stdy<-(stdy)/(phylo+spp+stdy+r+Vw)

	I2.spp<-(spp)/(phylo+spp+stdy+r+Vw)

	I2.t <- (phylo+spp+stdy+r)/(phylo+spp+stdy+r+Vw)

	tmpMatrix <- cbind(I2.phylo, I2.stdy, I2.spp, I2.t)
	mode <- posterior.mode(as.mcmc(tmpMatrix))
	CI <- HPDinterval(as.mcmc(tmpMatrix))

	p.mcmc = ifelse(round(CI[,"lower"], digits =3) != 0, "*", "")
	
  return(data.frame(post.mode = round(mode, digits = 3), CI = round(CI, digits = 3), p.mcmc = p.mcmc))
}