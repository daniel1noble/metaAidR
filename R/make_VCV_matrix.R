
#' @title Covariance and correlation matrix function basing on shared level ID
#' @description Function for generating simple covariance and correlation matrices 
#' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#' @param V Name of the variable (as a string – e.g, "V1") containing effect size variances variances
#' @param cluster Name of the variable (as a string – e.g, "V1") indicating which effects belong to the same cluster. Same value of 'cluster' are assumed to be nonindependent (correlated).
#' @param obs Name of the variable (as a string – e.g, "V1") containing individual IDs for each value in the V (Vector of variances). If this parameter is missing, label will be labelled with consecutive integers starting from 1.
#' @param rho Known or assumed correlation value among effect sizes sharing same 'cluster' value. Default value is 0.5.
#' @param type Optional logical parameter indicating whether a full variance-covariance matrix (default or "vcv") is needed or a correlation matrix ("cor") for the non-independent blocks of variance values.
#' @export

make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  type <- match.arg(type)
  if (missing(data)) {
    stop("Must specify dataframe via 'data' argument.")
  }
  if (missing(V)) {
    stop("Must specify name of the variance variable via 'V' argument.")
  }
  if (missing(cluster)) {
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  }
  if (missing(obs)) {
    obs <- 1:length(V)   
  }
  if (missing(type)) {
    type <- "vcv" 
  }
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  # find start and end coordinates for the subsets
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(utils::combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  return(new_matrix)
}


#Value: Labelled full variance-covariance or correlation matrice of the size and labels matching initial dataframe will be returned 

## TO DO: Implement a shared control equations for SMD and lnRR, maybe similar to eho = 0.5, but we have exact calculations so probably worth implementing