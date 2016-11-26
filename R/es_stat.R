#' @title es_stat: Function for calculating effect size statistics
#' @description Function for calculating common effect size statistics.  Effect size statistics include both the effect size itself and its sampling error. This function estimates both commonly used effect sizes (Hedges' d & g, Zr, lnOR, lnRR) along with less commonly used (lnHR) or newly developed effect sizes (i.e. for variance: lnCVR, lnVR, lnSD). For two-group comparisons, directionally is dependent on the group used for m1 or p1. In all cases, m1 - m2 or p1 - p2 is contrasted.
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param p1 Proportion in group 1. Used with type = "lnOR" only.
#' @param p2 Proportion in group 2. Used with type = "lnOR" only.
#' @param r Correlation coefficient. Used with type = "r" only .
#' @param nr Sample size used for estimating the correlation coefficient. Used with type "r" only.
#' @param type The type specifies the specific effect statistics one wishes to calculate. Types include: "d" = Hedges' d, "g" = bias corrected Hedges' d, "Zr" = Fisher's z-transformed correlation coefficient,"lnOR" = log odds ratio, "lnHR" = log hazards ratio, "lnRR" = log response ratio, "lnCVR" = log coefficient of variation ratio, "lnVR" = log variance ratio. 
#' @return Function returns the effect size and its sampling variance in a matrix (two column, n rows). The arguments can be vectors. 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

es_stat <- function(m1, m2, sd1, sd2, n1, n2, p1, p2, r, nr, type = c("d", "g", "Zr","lnOR", "lnRR", "lnCVR", "lnVR", "lnHR")){

	if(type == "g"){
		g     <- hedge(m1, m2, sd1, sd2, n1, n2, type = "g")
		v_g <- var_hedge(g, n1, n2, type = "g")
		return(cbind(g, v_g))
	}

	if(type == "d"){
		d     <- hedge(m1, m2, sd1, sd2, n1, n2, type = "d")
		v_d <- var_hedge(d, n1, n2, type = "d")	
		return(cbind(d, v_d))
	}

	if(type == "Zr"){
		return(Zr_es(r, nr))
	}

	if(type == "lnOR"){
		return(lOR_es(p1, p2, n1, n2))
	}

	if(type == "lnRR"){
		return(lnRR_es(m1, m2, sd1, sd2, n1, n2))
	}

}
	
#' @title hedge
#' @description Function for calculating Hedges' d or biased corrected Hedges' g. 
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param type  Whether Hedges' d (type = "d") and biased corrected Hedges' d or Hedges' g (type = "g") is the be used. 
#' @references Borenstein et al. 2009. Effect sizes based on means. In: Introduction to Meta-Analysis (eds. Borenstein, M., Hedges, L.V., Higgins, J.P.T. and Rothstein, H.R.). pg 21-32.
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

hedge <- function(m1, m2,  sd1, sd2, n1, n2, type = c("d", "g")){
	if(type == "g"){
	J <- 1 - (3 / ((4 * (n1 + n2 - 2)) -1))
	sd_pool<- (((n1 - 1) * (sd1^2)) + ((n2 - 1) * (sd2^2)))/(n1 + n2 - 2)
	h_g <- ((m1 - m2) / sqrt(sd_pool))*J
	}

	if(type == "d"){
	sd_pool<- (((n1 - 1) * (sd1^2)) + ((n2 - 1) * (sd2^2)))/(n1 + n2 - 2)
	h_g <- ((m1 - m2) / sqrt(sd_pool))
	}
return(h_g)
}

#' @title var_hedge
#' @description Function for calculating the sampling variance for Hedges' d or biased corrected Hedges' g. 
#' @param hedge Calculation of Hedges g or d
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param type  Whether Hedges' d (type = "d") and biased corrected Hedges' d or Hedges' g (type = "g") is the be used. 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
var_hedge <- function(hedge, n1, n2, type = c("d", "g")){
	if(type == "g"){
	J <- 1 - (3 / ((4 * (n1 + n2 - 2)) -1))
	vg <- ((n1 + n2)/(n1 * n2)) + ((hedge^2)/(2 * (n1 + n2 - 2)))
	v <- vg * (J^2)
	}

	if(type == "d"){
	v <- ((n1 + n2)/(n1 * n2)) + ((hedge^2)/(2 * (n1 + n2 - 2)))
	}
return(v)
}	

#' @title Zr_es
#' @description Function for calculating Fisher's z-transformed correlation. 
#' @param r Correlation coefficient
#' @param n Sample size
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

Zr_es <- function(r, n){
	z    <- 0.5*(log((1+r)/(1-r)))
	Vz <- 1 / (n -3)
	return(cbind(z, Vz))
}

#' @title lnOR
#' @description Function for calculating log odds ratio effect size statistics
#' @param p1 Proportion in group 1
#' @param p2 Proportion in group 2
#' @param n1 Sample size in group 1
#' @param n2 Sample size in group 2
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

lOR_es <- function(p1, p2, n1, n2){
	lnOR    <- log(p1 / (1-p1)) - log(p2 / (1-p2))
	v_lnOR   <- (1 / (p1*n1)) + (1/ ((1-p1)*n1)) + (1 / (p2*n2)) + (1/ ((1-p2)*n2))
return(cbind(lnOR, v_lnOR))
}

#' @title lnRR
#' @description Function for calculating log response ratio
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

lnRR_es <- function(m1, m2,  sd1, sd2, n1, n2){
     sd_pool<- (((n1 - 1) * (sd1^2)) + ((n2 - 1) * (sd2^2)))/(n1 + n2 - 2)
	lnRR <- log(m1) - log(m2)
	VlnRR <- (sd_pool^2)*((1/(n1*(m1^2)) + 1/(n2*(m2^2))))
return(cbind(lnRR,VlnRR))
}

#' @title lnCVR
#' @description Function for calculating log coefficient of variation ratio
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @references Nakagawa et al. 2015. 
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @export

lnCVR_es<-function(m1, m2, sd1, sd2, n1, n2){
	lnCVR <-(log(sd1) - log(m1) + 1 / (2*(n1 - 1))) - (log(sd2) - log(m2) + 1 / (2*(n2 - 1)))
	return(lnCVR)	
}


#' @title Sampling variance for lnCVR
#' @description Function for calculating sampling variance of log coefficient of variation ratio under different assumptions about mean-variance relationships. 
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param Equal.E.C.Corr Logical indicating whether the the correlation between mean and variance is equal in both group groups.
#' @param repeated.control Logical indicating whether there are repeated control groups used in the data set.
#' @param Control.IDs A vector indicating which rows share a common control
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @references Nakagawa et al. 2015. 
#' @export
var_lnCVR<-function(m1, m2, sd1, sd2, n1, n2, Equal.E.C.Corr=TRUE, repeated.control = FALSE, Control.IDs){

	if(repeated.control == TRUE){
		mean.control.for.cor<-m1[match(unique(Control.IDs), Control.IDs)]
		sd.control.for.cor<-sd1[match(unique(Control.IDs), Control.IDs)]
	}else{
		mean.control.for.cor<-m1
		sd.control.for.cor<-sd1	
		}
	
	if(Equal.E.C.Corr==TRUE){
		mvcorr<-stats::cor.test(log(c(mean.control.for.cor, m2)), log(c(sd.control.for.cor, sd2)))$estimate
			
		S2<- sd1^2 / (n1 * (m1^2)) + 1 / (2 * (n1 - 1)) - 2 * mvcorr * sqrt((sd1^2 / (n1 * (m1^2))) * (1 / (2 * (n1 - 1)))) + sd2^2 / (n2 * (m2^2)) + 1 / (2 * (n2 - 1)) - 2 * mvcorr * sqrt((sd2^2 / (n2 * (m2^2))) * (1 / (2 * (n2 - 1))))
		}else{
		Cmvcorr<-stats::cor.test(log(mean.control.for.cor), log(sd.control.for.cor))$estimate
		Emvcorr<-stats::cor.test(log(m2), (sd2))$estimate
	
		S2<- sd1^2 / (n1 * (m1^2)) + 1 / (2 * (n1 - 1)) - 2 * Cmvcorr * sqrt((sd1^2 / (n1 * (m1^2))) * (1 / (2 * (n1 - 1)))) + sd2^2 / (n2 * (m2^2)) + 1 / (2 * (n2 - 1)) - 2 * Emvcorr * sqrt((sd2^2 / (n2 * (m2^2))) * (1 / (2 * (n2 - 1))))	
	}
return(S2)	
}

#' @title lnSD_es
#' @description Function for calculating log standard deviation effect size statistics
#' @param sd Standard deviation of group
#' @param n Sample size of group
#' @references Nakagawa et al. 2015. 
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @export
lnSD_es<-function(sd, n){
	lnSD <- log(sd) + (1 / (2 * (n - 1)))
	v_lnSD <- (1 / (2 * (n - 1)))	
	return(cbind(lnSD, v_lnSD))
}

#' @title lnVR_es
#' @description Function for calculating log variation ratio effect size statistics.
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @references Nakagawa et al. 2015. 

lnVR_es<-function(sd1, sd2, n1, n2){
	lnVR <- log(sd1 / sd2) + (1 / (2 * (n1 - 1))) - (1 / (2 * (n2 - 1)))
	v_lnVR <- (1 / (2 * (n1 - 1))) + (1 / (2 * (n2 - 1)))
	return(cbind(lnVR, v_lnVR))
}
