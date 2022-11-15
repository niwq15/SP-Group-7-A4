## The 'newt' function has inputs, including theta, func, grad, hess


newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  f <- func(theta,...) ## use the 'func' function and 'theta' to get the objective function
  gradval <- grad(theta,...) ## use the 'grad' function and 'theta' to get the gradient vector 
  ## check whether the objective and derivatives are not finite
  ## ? f, gradval is single value or a vector
  if (is.finite(f)==FALSE) {
    break
    return("The objective function is not finite at the initial theta value.") 
  }
  if (any(is.finite(gradval)==FALSE)) {## gradval is a vector
    break
    return("The derivatives are not finite at the initial theta value.") 
  }
  ## If the objective or derivatives are finite, do the following steps 
  ## If the Hessian matrix is not supplied, we use an approximation by finite differencing of the gradient vector
  if (hess=NULL) {
    len <- length(gradval) ## the length of the gradient vector 
    Hfd <- matrix(0, len,len) ## finite difference Hessian
    for (i in 1:length(theta)) {## loop over parameters
      th1 <- theta; th1[i] <- th1[i]+ eps ## increase th1[i] by eps
      grad1 <- grad(th1,...) ## compute resulting derivative values
      Hfd[i,] <- (grad1 - gradval)/eps ## approximate -dl/dth[i]
    }
    hess <- (t(Hfd)+Hfd)/2 ##make sure the matrix is symmetric and store into 'hess' for later use
  }## now we have Hfd, an approximate Hessian matrix
  
  ## Check if the (approximated) Hessian matrix is positive definite, if so, use Cholesky decomposition to calculate the inverse
  ## The hessian matrix is positive definite if it has only positive eigenvalues
  eigensH <- eigen(hess)$values ## caculate the eigenvalues of 'hess'
  if (any(eigensH <= 0)){## if there exists some non-positive eigenvalues
    
    print("Warning: the Hessian is not positive definite.")
  } else {## if all the eigenvalues are positive
    
    chess <- chol(hess) ## solve with cholesky
    Hi <- backsolve(chess,forwardsolve(t(chess),diag(c(rep(1,length(eigensH)))))) ## calculate the hessian inverse
    ## print(Hi) ##result can be checked with solve(hess)
  }

  ## After each try, if a better theta value is found, use 'theta' to store the value
  n_iter <- 0; n_half <- 0
  
  while (n_iter < maxit) {
    theta2 <- theta - Hi%*%gradval ## use the Newton's formula
    f2 <- func(theta2,...)
    
    while (f2 > f) {
      n_half <- n_half + 1 ## counter goes up
      theta2 <-  theta - ( 0.5 ^(n_half) )*( Hi %*% gradval )
      f2 <- func(theta2)
      
      if (n_half >= max.half) {
        stop("The step fails to improve the objective")
        ## should return the value corresponding to initial theta value
      }
    }
    
    while (f2 < f) {
      n_iter <- n_iter + 1 ## count the number of iterations
      bound <- tol*abs(f2)*fscale
      if (any(abs(gradval)) < bound) {
       theta <- theta2 ## use 'theta' to store the optimal theta value
       f <- func(theta,...)
       iter <- n_iter
       g <- grad(theta,...) 
       #Hi <- 
       return(f, theta, iter, g, Hi)
      }
    }
    }
    
}
