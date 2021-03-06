
#' @title Folded normal distribution
#' @description Applying the folded normal distribution using the posterior parameter estimates of effect sizes that are normally distributed. This function will allow one to understand the overall magnitude of effect regardless of effect size direction. 
#' @param mu The posterior distribution of mean estimates from an MCMCglmm object. Alternatively, you can give it the mean and sd of a normal distribution and it will provide the mean (only used if type = "raw").
#' @param sd The standard deviation of the posterior distribution of variance-covariance matrix from an MCMCglmm object or just the sd of a normal distribution. 
#' @param type Indicates whether the posterior mean (i.e. "mean") or mode (i.e. "mode") should be returned. Alternatively, if one knows the mean and the sd of you can use type = "raw" to get the mean and variance of the folded normal distribution.
#' @return Posterior mean or mode along with the 95 percent credible intervals (i.e. the highest posterior density interval - HPDinterval). Alternatively if type = "raw" it will return the mean and variance for the folded normal distribution.
#' @references Morrisey, M.B. 2016. Meta-analysis of magnitudes, differences and variation in evolutionary parameters. Journal of Evolutionary Biology, 29:1882-1904
#' @references Morrisey, M.B. 2016. Rejoinder: Further considerations for meta-analysis of transformed quantities such as absolute values. Journal of Evolutionary Biology, 29:1922-1931
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' set.seed(60)
#' x <- rnorm(1000, 1, 1)
#' mean(x)
#' sd(x)
#' # Calculate the mean of the folded norm from raw x
#' foldnorm(mu = mean(x), sd = sd(x), type = "raw")
#' #Should be close to: 
#' mean(abs(x))
#' var(abs(x))
#' @export

foldnorm <-  function(mu, sd, type = c("mean", "mode", "raw")){
	
     postfnorm <- stats::dnorm(mu, 0, sd)*2*sd^2 + mu*(2*stats::pnorm(mu, 0, sd) -1)

	if(type == "raw"){
		    est <- postfnorm
		    var.fnorm <- mu^2 + sd^2 - (sd*sqrt(2/pi)*exp((-1*mu^2)/(2*sd^2)) + mu*(1-2*stats::pnorm(-1*mu/sd, 0, 1)))^2
		    est <- data.frame(Mean=est, Variance = var.fnorm)
	}

	if(type == "mean"){
		mean <- mean(postfnorm)
		    CI <- coda::HPDinterval(postfnorm)
		    est <- Matrix::cBind(mean = mean, CI)
	}

	if(type == "mode"){
		mode <- MCMCglmm::posterior.mode(postfnorm)
		     CI <- coda::HPDinterval(postfnorm)
		est <- Matrix::cBind(mode = mode, CI)
	}
     return(est)
}