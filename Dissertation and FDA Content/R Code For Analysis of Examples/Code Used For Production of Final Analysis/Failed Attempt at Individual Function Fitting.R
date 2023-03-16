fit_individ_func <- function(time_points,obs_matrix,basis,AR1_Weighting=FALSE,penalty=2){
  
  GCV_func <- function(log_lambda,basis,observations,time_points,penalty,weight=FALSE){
    
    if (is.matrix(weight) == FALSE){
      
      lambda <- 10^log_lambda
      
      fd_par_obj <- fdPar(basis,penalty,lambda)
      
      smoothbasisobj <- smooth.basis(time_points,observations,fd_par_obj)
      
      return(sum(smoothbasisobj$gcv))
      
    } else{
      
      lambda <- 10^log_lambda
      
      fd_par_obj <- fdPar(basis,penalty,lambda)
      
      smoothbasisobj <- smooth.basis(time_points,observations,fd_par_obj,wtvec = weight)
      
      return(sum(smoothbasisobj$gcv))
      
    }
  }
  
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  rep_dims <-  colnames(obs_matrix)
  
  if (AR1_Weighting == FALSE){
    for (i in 1:ncol(obs_matrix)){
      
      dim_i <- rep_dims[i]
      
      optimised_function <- optimise(GCV_func,lower=-10,upper=10,basis=basis,observations=obs_matrix[,i],time_points=time_points,penalty=penalty)
      
      minimum_log_lambda <- optimised_function$minimum
      
      minimum_lambda <- 10^minimum_log_lambda
      
      fd_par_obj <- fdPar(basis,2,minimum_lambda)
      
      fd_obj <- smooth.basis(time_points,obs_matrix[,i],fd_par_obj)$fd
      
      coefs_i <- fd_obj$coefs
      
      colnames(coefs_i) <- dim_i
      
      if (i==1){
        coefs <- coefs_i
      } else{
        coefs <- cbind(coefs, coefs_i)
      }
    } }
  
  
  if (AR1_Weighting == TRUE){
    for (i in 1:ncol(obs_matrix)){
      
      dim_i <- rep_dims[i]
      
      cor_coef <- acf(obs_matrix[,i],plot = FALSE)$acf[2]
      
      Weight_matrix <- solve(ar1_cor(length(obs_matrix[,i]),rho=cor_coef))
      
      optimised_function <- optimise(GCV_func,lower=-10,upper=10,basis=basis,observations=obs_matrix[,i],time_points=time_points,penalty=penalty,weight=Weight_matrix)
      
      minimum_log_lambda <- optimised_function$minimum
      
      minimum_lambda <- 10^minimum_log_lambda
      
      fd_par_obj <- fdPar(basis,2,minimum_lambda)
      
      fd_obj <- smooth.basis(time_points,obs_matrix[,i],fd_par_obj)$fd
      
      coefs_i <- fd_obj$coefs
      
      colnames(coefs_i) <- dim_i
      
      if (i==1){
        coefs <- coefs_i
      } else{
        coefs <- cbind(coefs, coefs_i)
      }
    }}
  
  
  
  merged_fd <- fd(coef = coefs, basisobj =basis)
  merged_fd$fdnames$reps <- rep_dims
  
  return(merged_fd)
  
}