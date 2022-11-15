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
  
  if (any(is.finite(gradval)==FALSE)) {
    break
    return("The derivatives are not finite at the initial theta value.") 
  }
  
  ## If the objective or derivatives are finite, do the following steps.
  ## If the Hessian matrix is not supplied, we use an approximation by finite
  ## differencing of the gradient vector
  if (is.null(hess)) {
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
  }
    
  ## now we have Hfd, an approximate Hessian matrix
  
  ## empty vector to store the points
  xval <- c() 
  
  ## use the 'theta' values as initial x0
  xval <- theta
  len <- length(xval)
  
  ## loop over 2 to max iterations
  for (i in 2:maxit) {
    gradval <- grad(xval)
    
    ## cholesky decomposition to get hessian inverse
    chess <- chol(hess(xval))
    Hi <- backsolve(chess,forwardsolve(t(chess),diag(rep(1,len))))
    
    ## use formula
    xvalcheck <- xval - Hi %*% gradval
    
    ## if the step fails to reduce the objective
    while (norm(xvalcheck - xval,"F") > 0) {
      
      xvalcheck <- xval - (1/2) * Hi %*% gradval
    }
    
    ## If we have an answer close enough to last time then stop
    if ( norm(xval - xvalcheck, "F") < tol) {
      
      ## get final values
      theta <- xvalcheck
      f <- func(theta, ...)
      g <- grad(theta, ...)
      chess <- chol(hess(theta))
      Hi <- backsolve(chess,forwardsolve(chess,diag(rep(1,len))))
      iter <- i
      
      return(f, theta, iter, g, Hi)
    }
  }
}



rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


### Hessian inverse test
# create hessian
hess <- array(c(3,0,1,0,3,0,1,0,3),dim=c(3,3))

# save eigenvalues
eigensH <- eigen(hess)$values

## check if hessian is positive definite
# for positive definite all eigenvalues greater than 0
if (any(eigensH <= 0)){ ## if matrix is not positive definite
  
  print("Hessian not positive definite")
  # change this print to a return in the function
  
} else { ## matrix is positive definite
  
  # assign cholesky upper triangle
  chess <- chol(hess)
  
  # assign the hessian inverse (from notes)
  Hi <- backsolve(chess,forwardsolve(t(chess),diag(c(rep(1,length(eigensH))))))
  
  print(Hi)
  # remove this print in the function, result can be checked with 
  # print(solve(hess))
}