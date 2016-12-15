
#' @title I^2 Function 
#' @description Function for calculating commonly reported  I^2 measures from MCMCglmm model objects.  See Nakagawa and Santos (2012) for detailed explanation on the various types of I^2
#' @param model The MCMCglmm or metafor model object. Note that if using a metafor model object an observation level random effect must be used to calculate the residual variance. This should be input as "~1|obs" in the random effect list. 
#' @param v The vector of sampling variance for each effect size. 
#' @param sims The number of simulations used for calculating confidence intervals on I^2 estimates for metafor objects.
#' @param phylo A character string with the name of the phylogenetic random effect. Defaults to FALSE meaning that no phylogenetic heritability is calculated. 
#' @return A data.frame containing the relevant I^2 measures along with the 95 percent confidence / credible intervals.
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @references Nakagawa, S. and Santos, E.S.A. (2012) Methodological issues and advances in biological meta-analysis. Evolutionary Ecology, 26:1253-1274.
#' @export

I2 <- function(model, v, sims = 1500, phylo = FALSE){
	
	if(class(model) != "MCMCglmm" && class(model) != "rma.mv" && class(model) != "rma"){
		stop("The model object is not of class 'MCMCglmm' or 'metafor'")
		}

		wi <- 1/v  #weight
		Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))

	if("MCMCglmm" %in% class(model)){
		# Get posterior distribution 
		post <- model$VCV[,-match(c("sqrt(mev):sqrt(mev).meta"), colnames(model$VCV))]

		#Calculate total variance
		         VT <- rowSums(Matrix::cBind(post, Vw))
		          Vt <- rowSums(post)  # remove Vw
		
		# For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
		  I2_re <- post / VT
		  I2_total  <- Vt / VT

		if(phylo == FALSE){
			tmpMatrix <- Matrix::cBind(I2_re, total = I2_total)
		}else{
		  	I2_phylo <- post[,match(phylo, colnames(sigma2))] / Vt
		  	tmpMatrix <- Matrix::cBind(I2_re, I2_phylo, total = I2_total)
		  	}

		   mode <- MCMCglmm::posterior.mode(coda::as.mcmc(tmpMatrix))
		   CI <- coda::HPDinterval(coda::as.mcmc(tmpMatrix))
		    colnames(CI) <- c("2.5% CI", "97.5% CI")

		    I2_Table <- as.data.frame(Matrix::cBind(I2_Est. = mode[-match("units", names(mode))], CI[-match("units", rownames(CI)),]))
		    class(I2_Table) <- c("metaAidR", "data.frame")

	return(round_df(I2_Table, digits = 4))
  	}

  	if("rma.mv" %in% class(model) | "rma" %in% class(model)){
  		#Monte Carlo Simulations
		# From metafor extract the important statistics
  		sigma2 <- matrix(model$sigma2, nrow = 1, ncol = length(model$sigma2))
  		colnames(sigma2) <- model$s.names
  		sigmaN <- model$s.nlevels

  		if("obs" %in% colnames(sigma2) == FALSE){
  			stop("The metafor object does not contain a residual variance estimate. Please include an observation-level random effect (~1|obs) when fitting model")
  		}

  		#For each variance estimate use Monte Carlo simulation of data
  		Sims <- data.frame(mapply(function(x,y) simMonteCarlo(x, y, sims = sims), sigma2, sigmaN))
		colnames(Sims) <- colnames(sigma2) 
		
		#Calculate total variance
		VT <- rowSums(Matrix::cBind(Sims, Vw))
		Vt <- rowSums(Sims)  # remove Vw
		
		# For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
		 I2_re       <- Sims / VT
		 I2_total   <- Vt / VT

		  if(phylo == FALSE){
		  	tmpMatrix <- Matrix::cBind(I2_re[,-match("obs", colnames(I2_re))], total = I2_total)
		   }else{
		  	I2_phylo <- Sims[, match(phylo, colnames(sigma2))] / Vt
		  	 tmpMatrix <- Matrix::cBind(I2_re[,-match("obs", colnames(I2_re))], phylo = I2_phylo, total = I2_total)
		  }

		CI <- lapply(tmpMatrix, function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE))
		I_CI <- as.data.frame(do.call(rbind, CI))
		colnames(I_CI) <- c("2.5% CI", "97.5% CI")
		I2_table <- Matrix::cBind(I2_Est. = colMeans(tmpMatrix), I_CI )
		
		class(I2_table) <- c("metaAidR", "data.frame")

	return(round_df(I2_table, digits = 4))
  	}

}


#' @title Parametric simulation 
#' @description Function for calculating I2 estimates using parametric simulations of model estimates taken from metafor. Note that the effectiveness of these simulations depends on the accuracy of model variance estimates.
#' @param estimate The estimate (i.e. variance) from a metafor model
#' @param sims The number of simulations 
#' @param n The sample size used in estimating the variance 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
  simMonteCarlo <- function(estimate, n, sims){
  		set.seed(07)
  		tmp <- data.frame(num = rep(1:sims, each = n), y = stats::rnorm(n*sims, 0, sqrt(estimate)))
  		Var <- dplyr::summarise(dplyr::group_by(tmp, num), var = stats::var(y))
  		return(as.numeric(Var$var))
  	}


 # Function for  rounding a data frame
round_df <- function(x, digits) {
    numeric_columns <- sapply(x, class) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}
