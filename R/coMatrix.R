#' @title Covariance and correlation matrix function
#' @description Function for generating various covariance and correlation matrices used to control for sources of non-independence in meta-analyses. Currently only shared-control and general forms of within-study covariances (e.g. shared traits - for sensitivity analysis) can be conducted. More complex matrices are not yet implemented. 
#' @param data The data object that a correlation or covariance matrix needs to be generated for. 
#' @param V Character strong indicating the column name of the known sampling error variance.
#' @param var1 The vector describing the first variable each effect size is derived from. For example, this could be study ID.
#' @param var2 An optional second variable that will further split the data. For example, by trait type. If var2 is not provided it is assumed only var1 is of interest.
#' @param cor The known or hypothesized correlation between effect sizes resulting from within-study non-independence. cor values could be derived from raw data or from known correlations taken from the literature. Alternatively, multiple cor values can be used to derive several matrices each of which can be used in a sensitivity analysis.
#' @param type Whether a general within study blocked matrix is needed ("ws" - note this could be based on more than one variable), or a shared control covariance matrix is needed ("sc").
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
cov_matrix <- function(data, V, var1, var2 = FALSE, cor, type = c("sc", "ws")){
	if(type == "ws"){
	if(var2 == FALSE){
		spltP <- split(data, list(data[, var1]))
	} else{
		spltP <- split(data, list(data[, var1], data[, var2]))
	}
		covMat <- lapply(spltP, function(x) covMatrix(x, es_var = V, cor = cor))
            		ws_cov <- Matrix::bdiag(covMat)
            		rownames(ws_cov) <- colnames(ws_cov) <- rownames(ws_cov)
		return(ws_cov)
	}
}

#' @title Generate a covariance matrix
#' @description This is a sub-function used by cov_matrix which generates a covariance matrix for a sub-set of a data frame containing the effect size sampling variance and a known correlation matrix
#' @param data The data object that a correlation or covariance matrix needs to be generated for. 
#' @param es_var Character string indicating the column name of the known sampling error variance.
#' @param cor The known or hypothesized correlation between effect sizes resulting from within-study non-independence. cor values could be derived from raw data or from known correlations taken from the literature. Alternatively, multiple cor values can be used to derive several matrices each of which can be used in a sensitivity analysis.
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
covMatrix <- function(data, es_var, cor){
	 tmp <- expand.grid(sqrt(data[, es_var]), sqrt(data[, es_var]))
	 ## NEED TO CHANGE below. This won't work in all cases. We need 1 only in cases where we know for sure the row is the same. 
	  tmp$cor <- ifelse(tmp$Var1 == tmp$Var2, 1, 0.5) 
	 tmp$cov <- tmp$cor * tmp$Var1 * tmp$Var2
	  corMat <- matrix(tmp$cov , nrow = nrow(data), ncol = nrow(data))
  return(corMat)
}