
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
	
	if(class(model) != "MCMCglmm" & class(model) != "rma.mv" & class(model) != "rma"){
		stop("The model object is not of class 'MCMCglmm' or 'metafor'")
		}

		wi <- 1/v  #weight
		Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))

	if(class(model) == "MCMCglmm" ){
		# Get posterior distribution 
		post <- model$VCV[,-match(c("sqrt(mev):sqrt(mev).meta"), colnames(model$VCV))]

		#Calculate total variance
		         VT <- rowSums(Matrix::cBind(post, Vw))
		          Vt <- rowSums(post)  # remove Vw
		
		# For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
		  I2_re <- post / VT

		  if("phylo" %in% names(re.list)){
		  	I2_phylo <- post[,grep(re.list$phylo), colnames(sigma2)] / Vt
		   }else{
		  	I2_phylo <- FALSE
		  }

		 I2_total  <- Vt / VT
		
		if(I2_phylo == FALSE){
			tmpMatrix <- Matrix::cBind(I2_re, total = I2_total)
		} else{
		            	tmpMatrix <- Matrix::cBind(I2_re, I2_phylo, total = I2_total)
		}
		   
		   mode <- MCMCglmm::posterior.mode(coda::as.mcmc(tmpMatrix))
		   CI <- coda::HPDinterval(coda::as.mcmc(tmpMatrix))
		    colnames(CI) <- c("2.5% CI", "97.5% CI")

		    Est_Table <- Matrix::cBind(Est. = mode, CI)
	return(round_df(Est_Table, digits = 4))
  	}

  	if(class(model) ==  "rma.mv" | class(model) ==  "rma"){
  		#Monte Carlo Simulations
		# From metafor extract the important statistics
  		sigma2 <- matrix(model$sigma2, nrow = 1, ncol = length(model$sigma2))
  		colnames(sigma2) <- model$s.names

  		#For each variance estimate use Monte Carlo simulation of data
  		Sims <- data.frame(sapply(sigma2, function(x) simMonteCarlo(x, sims = 1000)))
		colnames(Sims) <- colnames(sigma2) 
		
		#Calculate total variance
		         VT <- rowSums(Matrix::cBind(Sims, Vw))
		          Vt <- rowSums(Sims)  # remove Vw
		
		# For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
		  I2_re       <- Sims / VT

		  if("phylo" %in% names(re.list)){
		  	I2_phylo <- Sims[,grep(re.list$phylo), colnames(sigma2)] / Vt
		   }else{
		  	I2_phylo <- FALSE
		  }

		  I2_total   <- Vt / VT

		  if(I2_phylo == FALSE){
			tmpMatrix <- Matrix::cBind(I2_re[,-match("obs", colnames(I2_re))], total = I2_total)
		} else{
		            	tmpMatrix <- Matrix::cBind(I2_re[,-match("obs", colnames(I2_re))], phylo = I2_phylo, total = I2_total)
		}

		      I_CI <- plyr::ldply(lapply(tmpMatrix, function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)))
		I2_table <- round_df(Matrix::cBind(Est = colMeans(tmpMatrix), I_CI[,-1]), digits = 4)

	return(I2_table)
  	}

}


#' @title Parametric simulation 
#' @description Function for calculating I2 estimates using parametric simulations of model estimates taken from metafor. Note that the effectiveness of these simulations depends on the accuracy of model variance estimates.
#' @param estimate The estimate (i.e. variance) from a metafor model
#' @param sims The number of simulations 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
  simMonteCarlo <- function(estimate, sims){
  		set.seed(07)
  		tmp <- data.frame(num = rep(1:sims, each = 5000), y = stats::rnorm(5000*sims, 0, sqrt(estimate)))
  		Var <- dplyr::summarise(dplyr::group_by(tmp, num), var = stats::var(y))
  		return(as.numeric(Var$var))
  	}


 # Function for  rounding a data frame
round_df <- function(x, digits) {
    numeric_columns <- sapply(x, class) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}
