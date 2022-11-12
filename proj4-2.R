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
    ## the length of the gradient vector
    len <- length(gradval) 
    
    ## finite difference Hessian
    Hfd <- matrix(0, len,len)
    
    for (i in 1:length(theta)) {## loop over parameters
      ## increase th1[i] by eps
      th1 <- theta; th1[i] <- th1[i]+ eps 
      
      ## compute resulting derivative values
      grad1 <- grad(th1,...)
      
      ## approximate -dl/dth[i]
      Hfd[i,] <- (grad1 - gradval)/eps 
    }
    
    ##make sure the matrix is symmetric
    Hfd <- (t(Hfd)+Hfd)/2
    
  }## now we have Hfd, an approximate Hessian matrix
  
}



### Hessian inverse test
# create hessian
hess <- array(c(3,0,1,0,3,0,1,0,3),dim=c(3,3))

# save eigenvalues
eigensH <- eigen(hess)$values

## check if hessian is positive definite
# for positive definite all eigenvalues greater than 0
if (any(abs(eigensH) <= 0)){
  
  print("Hessian not positive definite")
  # change this print to a return in the function
  
} else {
  
  # assign cholesky upper triangle
  chess <- chol(hess)
  
  # assign the hessian inverse (from notes)
  Hi <- backsolve(chess,forwardsolve(t(chess),diag(c(rep(1,length(eigensH))))))
  
  print(Hi)
  # remove this print in the function
}