
#' Applying the folded normal distribution using the posterior parameter estimates of effect sizes that are normally distributed. This function will allow one to understand the overall magnitude of effect regardless of effect size direction. 
#' @param mu The posterior distribution of mean estimates from an MCMCglmm object.
#' @param sd The standard deviation of the posterior distribution of variance-covariance matrix from an MCMCglmm object.
#' @param type Indicates whether the posterior mean (i.e. "mean") or mode (i.e. "mode") should be returned.
#' @return Posterior mean or mode along with the 95% credible intervals (i.e. the highest posterior density interval - HPDinterval) 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @examples
#' dir <- paste0(system.file("extdata",  package="ibutt"), "/")
#' data <- ibuttProc(dir, date = "2015-11-06", time = c("10:00", "11:00"), list = TRUE)
#' @export

foldnorm <-  function(mu, sd, type = c("mean", "mode")){
	
	postfnorm <- dnorm(mu, 0, sd)*2*sd^2 + mu*(2*pnorm(mu, 0, sd) -1)

	if(type == "mean"){
		mean <- mean(postfnorm)
		    CI <- coda::HPDinterval(postfnorm)
		    est <- cBind(mean = mean, CI)
	}

	if(type == "mode"){
		mode <- MCMCglmm::posterior.mode(postfnorm)
		     CI <- coda::HPDinterval(postfnorm)
		est <- cBind(mode = mode, CI)
	}
     return(est)
}