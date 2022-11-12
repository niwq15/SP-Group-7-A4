## The 'newt' function has inputs, including theta, func, grad, hess


newt <- function(theta,
                 func,
                 grad,
                 hess=NULL,
                 ...,
                 tol=1e-8,
                 fscale=1,
                 maxit=100,
                 max.half=20,
                 eps=1e-6) {
  ## use the 'func' function and 'theta' to get the objective function
  f <- func(theta,...) 
  
  ## use the 'grad' function and 'theta' to get the gradient vector
  gradval <- grad(theta,...)  
  
  ## check whether the objective and derivatives are not finite
  ## ? f, gradval is single value or a vector
  if (is.finite(f)==FALSE) {
    break
    return("The objective function is not finite at the initial theta value.") 
  }
  
  if (is.finite(gradval)==FALSE) {
    break
    return("The derivatives are not finite at the initial theta value.") 
  }
  
  ## If the objective or derivatives are finite, do the following steps.
  ## If the Hessian matrix is not supplied, we use an approximation by finite
  ## differencing of the gradient vector
  if (hess=NULL) {
    len <- length(gradval) ## the length of the gradient vector
    Hfd <- matrix(0, len,len) ## finite difference Hessian
    for (i in 1:length(theta)) {## loop over parameters
      th1 <- theta; th1[i] <- th1[i]+ eps ## increase th1[i] by eps
      grad1 <- grad(th1,...) ## compute resulting derivative values
      Hfd[i,] <- (grad1 - gradval)/eps ## approximate -dl/dth[i]
    }
    
    ##make sure the matrix is symmetric
    Hfd <- (t(Hfd)+Hfd)/2
    
  }## now we have Hfd, an approximate Hessian matrix 
   
}
